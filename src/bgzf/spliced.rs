//! Line reading that keeps the compressed blocks a line came from.
//!
//! When a VCF line is megabytes long it spans dozens of BGZF blocks, and a
//! block whose uncompressed content contains no newline lies entirely inside
//! one line. If the merge copies that stretch of the line into its output
//! verbatim, the *compressed* block is equally valid output — BGZF blocks are
//! self-contained gzip members with their own CRCs, and records do not need
//! to align to block boundaries. This reader inflates a file exactly as a
//! normal reader would, but tags each line with the raw compressed blocks
//! that sit wholly inside it, so `merge` can splice them into its output and
//! skip re-deflating those bytes.
//!
//! One thread per file reads, inflates and splits lines, feeding a bounded
//! channel — the same latency-hiding shape as the pooled BGZF reader, minus
//! the fan-out (a merge has one such reader per input, which is parallelism
//! enough).

use std::io;
use std::sync::mpsc::{sync_channel, Receiver};
use std::thread::JoinHandle;

use crate::io::open_raw_blocks;

/// A compressed block lying wholly inside one line.
pub struct InteriorBlock {
    /// Offset within the line's text where this block's content begins.
    pub line_off: usize,
    /// Uncompressed length of the block.
    pub ulen: usize,
    /// The complete BGZF member, exactly as read from the input.
    pub raw: Vec<u8>,
}

impl InteriorBlock {
    #[inline]
    pub fn end(&self) -> usize {
        self.line_off + self.ulen
    }
}

/// One line plus the compressed blocks interior to it. `text` excludes the
/// trailing newline; interior blocks never contain one by definition.
pub struct LineUnit {
    pub text: Vec<u8>,
    pub blocks: Vec<InteriorBlock>,
}

pub struct SplicedLineReader {
    rx: Receiver<Result<LineUnit, String>>,
    cur: LineUnit,
    handle: Option<JoinHandle<()>>,
}

impl SplicedLineReader {
    pub fn open(path: &str) -> io::Result<SplicedLineReader> {
        let mut blocks = open_raw_blocks(path)?;
        let path_owned = path.to_string();
        // Depth bounds memory at a few lines per file; wide-cohort lines are
        // megabytes each, so this is the number that keeps ten inputs from
        // holding a hundred lines of text between them.
        let (tx, rx) = sync_channel::<Result<LineUnit, String>>(4);
        let handle = std::thread::spawn(move || {
            let mut d = libdeflater::Decompressor::new();
            let mut tmp: Vec<u8> = Vec::with_capacity(1 << 16);
            let mut cur_text: Vec<u8> = Vec::new();
            let mut cur_blocks: Vec<InteriorBlock> = Vec::new();
            loop {
                let blk = match blocks.next_block() {
                    Ok(Some(b)) => b,
                    Ok(None) => break,
                    Err(e) => {
                        let _ = tx.send(Err(format!("{path_owned}: {e}")));
                        return;
                    }
                };
                if blk.isize == 0 {
                    continue; // EOF marker (or an empty block)
                }
                tmp.clear();
                if let Err(e) = blk.inflate_into(&mut d, &mut tmp) {
                    let _ = tx.send(Err(format!("{path_owned}: {e}")));
                    return;
                }
                if memchr::memchr(b'\n', &tmp).is_none() {
                    // The whole block belongs to the current line: this is the
                    // spliceable case, so keep the compressed bytes.
                    cur_blocks.push(InteriorBlock {
                        line_off: cur_text.len(),
                        ulen: tmp.len(),
                        raw: blk.bytes,
                    });
                    cur_text.extend_from_slice(&tmp);
                    continue;
                }
                // Boundary block: complete lines end (and the next begins)
                // inside it, so its bytes travel as text only.
                let mut start = 0usize;
                while let Some(p) = memchr::memchr(b'\n', &tmp[start..]) {
                    let nl = start + p;
                    cur_text.extend_from_slice(&tmp[start..nl]);
                    let unit = LineUnit {
                        text: std::mem::take(&mut cur_text),
                        blocks: std::mem::take(&mut cur_blocks),
                    };
                    if tx.send(Ok(unit)).is_err() {
                        return; // consumer hung up
                    }
                    start = nl + 1;
                }
                cur_text.extend_from_slice(&tmp[start..]);
            }
            if !cur_text.is_empty() {
                let unit = LineUnit { text: cur_text, blocks: cur_blocks };
                let _ = tx.send(Ok(unit));
            }
        });
        Ok(SplicedLineReader {
            rx,
            cur: LineUnit { text: Vec::new(), blocks: Vec::new() },
            handle: Some(handle),
        })
    }

    /// Move to the next line; false at end of input.
    pub fn advance(&mut self) -> Result<bool, String> {
        match self.rx.recv() {
            Ok(Ok(unit)) => {
                self.cur = unit;
                Ok(true)
            }
            Ok(Err(e)) => Err(e),
            Err(_) => Ok(false),
        }
    }

    #[inline]
    pub fn line(&self) -> &[u8] {
        &self.cur.text
    }

    /// Take ownership of the current line's text (the reader's copy becomes
    /// empty). Offsets in `take_blocks` refer to this text.
    pub fn take_line(&mut self) -> Vec<u8> {
        std::mem::take(&mut self.cur.text)
    }

    pub fn take_blocks(&mut self) -> Vec<InteriorBlock> {
        std::mem::take(&mut self.cur.blocks)
    }
}

impl Drop for SplicedLineReader {
    fn drop(&mut self) {
        // Unblock the producer if we stop early, then reap it.
        while self.rx.try_recv().is_ok() {}
        drop(std::mem::replace(&mut self.rx, {
            let (_t, r) = sync_channel(1);
            r
        }));
        if let Some(h) = self.handle.take() {
            let _ = h.join();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    /// Round-trip a synthetic file with lines both longer and shorter than a
    /// block, checking text reassembly and that interior blocks are exactly
    /// the newline-free ones and map back onto the text.
    #[test]
    fn lines_reassemble_and_interior_blocks_map_back() {
        let dir = std::env::temp_dir().join("vcfr_spliced_test");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("t.vcf.gz");
        let mut lines: Vec<Vec<u8>> = Vec::new();
        lines.push(b"##header line".to_vec());
        for i in 0..6u64 {
            let mut l = format!("chr1\t{}\tlong\t", i + 1).into_bytes();
            // ~200KB lines span several 64KB blocks.
            while l.len() < 200_000 + (i as usize * 7919) {
                l.extend_from_slice(format!("\tF{:x}", i * 31 + l.len() as u64 % 97).as_bytes());
            }
            lines.push(l);
            lines.push(format!("chr1\t{}\tshort", i + 100).into_bytes());
        }
        {
            let f = std::fs::File::create(&path).unwrap();
            let mut w = crate::bgzf::BgzfWriter::new(f, 0, 1, false);
            for l in &lines {
                w.write_all(l).unwrap();
                w.write_all(b"\n").unwrap();
            }
            w.finish().unwrap();
        }

        let mut r = SplicedLineReader::open(path.to_str().unwrap()).unwrap();
        let mut d = libdeflater::Decompressor::new();
        let mut got = 0usize;
        let mut spliceable = 0usize;
        while r.advance().unwrap() {
            let blocks = r.take_blocks();
            let text = r.take_line();
            assert_eq!(text, lines[got], "line {got} text mismatch");
            for b in &blocks {
                // Each interior block's content must equal the text slice it
                // claims to cover, and must contain no newline.
                let rb = crate::bgzf::RawBlockReader::new(std::io::Cursor::new(&b.raw[..]))
                    .next_block()
                    .unwrap()
                    .unwrap();
                let mut back = Vec::new();
                rb.inflate_into(&mut d, &mut back).unwrap();
                assert_eq!(back.len(), b.ulen);
                assert_eq!(&back[..], &text[b.line_off..b.end()]);
                assert!(!back.contains(&b'\n'));
                spliceable += 1;
            }
            got += 1;
        }
        assert_eq!(got, lines.len());
        assert!(spliceable >= 12, "long lines should carry interior blocks, saw {spliceable}");
        let _ = std::fs::remove_file(&path);
    }
}
