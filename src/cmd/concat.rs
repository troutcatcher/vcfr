//! `vcfr concat` — concatenate VCFs that hold the same samples.

use std::io::Write;

use crate::io::{detect_compression, open_raw_blocks, open_reader, Compression};
use crate::vcf::header::{union_meta, Header};

use super::{ignore_broken_pipe, resolve_threads, split_threads, OutputOpts};

#[derive(clap::Args, Debug)]
pub struct ConcatArgs {
    /// VCFs to concatenate, in output order
    #[arg(value_name = "FILE")]
    pub inputs: Vec<String>,

    /// Copy compressed blocks straight through without recompressing.
    /// Requires BGZF input and output and byte-identical headers.
    #[arg(short = 'n', long = "naive")]
    pub naive: bool,

    /// Read the file list from FILE, one path per line
    #[arg(short = 'F', long = "file-list", value_name = "FILE")]
    pub file_list: Option<String>,

    #[arg(long = "threads", default_value_t = 0, value_name = "N")]
    pub threads: usize,

    #[command(flatten)]
    pub out: OutputOpts,
}

pub fn run(a: &ConcatArgs) -> Result<(), String> {
    let inputs = collect_inputs(a)?;
    let threads = resolve_threads(a.threads);
    if a.naive {
        return naive(&inputs, a, threads);
    }

    let (rthreads, wthreads) = split_threads(threads, a.out.bgzf()?);
    let mut w = a.out.writer(wthreads)?;

    // Read every header first so an incompatibility is reported before we
    // start emitting output.
    let mut headers = Vec::with_capacity(inputs.len());
    for p in &inputs {
        let mut r = open_reader(p, 1).map_err(|e| e.to_string())?;
        headers.push(Header::read(&mut r).map_err(|e| format!("{p}: {e}"))?);
    }
    for (p, h) in inputs.iter().zip(&headers).skip(1) {
        if h.samples != headers[0].samples {
            return Err(format!(
                "{p}: sample columns differ from {} — concat requires identical samples in identical order (did you mean `vcfr merge`?)",
                inputs[0]
            ));
        }
    }
    let merged = union_meta(&headers);
    let mut buf = Vec::with_capacity(1 << 16);
    merged.write(&mut buf, &headers[0].samples, false);
    w.write_all(&buf).map_err(|e| e.to_string())?;

    let res = (|| -> std::io::Result<()> {
        for p in &inputs {
            let mut r = open_reader(p, rthreads)?;
            // Skip this file's own header.
            while r.advance()? {
                if !r.line().starts_with(b"#") {
                    w.write_all(r.line())?;
                    w.write_all(b"\n")?;
                    break;
                }
            }
            while r.advance()? {
                let line = r.line();
                if line.is_empty() {
                    continue;
                }
                w.write_all(line)?;
                w.write_all(b"\n")?;
            }
        }
        Ok(())
    })();
    ignore_broken_pipe(res)?;
    ignore_broken_pipe(w.finish())
}

/// Block-level concatenation.
///
/// Every input is BGZF, so once the header has been skipped the remaining
/// blocks are already valid output: they are copied verbatim, never inflated.
/// Only the one block straddling the end of the header is re-encoded.
fn naive(inputs: &[String], a: &ConcatArgs, threads: usize) -> Result<(), String> {
    if !a.out.bgzf()? {
        return Err("--naive requires BGZF output (-O z)".into());
    }
    for p in inputs {
        if detect_compression(p).map_err(|e| e.to_string())? != Compression::Bgzf {
            return Err(format!("{p}: --naive requires BGZF-compressed input"));
        }
    }

    let mut reference: Option<Vec<u8>> = None;
    let mut w = a.out.writer(threads)?;

    for p in inputs.iter() {
        let mut blocks = open_raw_blocks(p).map_err(|e| e.to_string())?;
        let mut d = libdeflater::Decompressor::new();
        let mut header = Vec::new();
        let mut tail: Vec<u8> = Vec::new();
        let mut found = false;

        // Inflate blocks until the header ends, keeping whatever record bytes
        // shared that block.
        while let Some(b) = blocks.next_block().map_err(|e| format!("{p}: {e}"))? {
            let base = header.len();
            b.inflate_into(&mut d, &mut header).map_err(|e| format!("{p}: {e}"))?;
            if let Some(end) = header_end(&header, base) {
                tail = header.split_off(end);
                found = true;
                break;
            }
            if header.len() > 64 << 20 {
                return Err(format!("{p}: no #CHROM line in the first 64 MiB"));
            }
        }
        if !found {
            return Err(format!("{p}: missing #CHROM header line"));
        }
        match &reference {
            None => {
                reference = Some(header.clone());
                // The first file's header is the output header.
                w.write_all(&header).map_err(|e| e.to_string())?;
            }
            Some(r) => {
                if *r != header {
                    return Err(format!(
                        "{p}: header differs from {}; --naive needs byte-identical headers",
                        inputs[0]
                    ));
                }
            }
        }
        // Re-encode the partial block, then pass the rest through untouched.
        if !tail.is_empty() {
            w.write_all(&tail).map_err(|e| e.to_string())?;
        }
        while let Some(b) = blocks.next_block().map_err(|e| format!("{p}: {e}"))? {
            if b.isize == 0 {
                continue; // BGZF EOF marker
            }
            w.write_raw_block(b.bytes).map_err(|e| e.to_string())?;
        }
    }
    ignore_broken_pipe(w.finish())
}

/// Offset just past the newline ending the `#CHROM` line, if it is present in
/// `buf` at or after `from`.
fn header_end(buf: &[u8], from: usize) -> Option<usize> {
    // The #CHROM line may have started in an earlier block, so rescan from the
    // beginning of the last unterminated line.
    let start = buf[..from].iter().rposition(|&b| b == b'\n').map_or(0, |p| p + 1);
    let mut i = start;
    while i < buf.len() {
        let nl = memchr::memchr(b'\n', &buf[i..])? + i;
        if buf[i..].starts_with(b"#CHROM") {
            return Some(nl + 1);
        }
        if !buf[i..].starts_with(b"#") {
            return Some(i);
        }
        i = nl + 1;
    }
    None
}

fn collect_inputs(a: &ConcatArgs) -> Result<Vec<String>, String> {
    let mut v = a.inputs.clone();
    if let Some(f) = &a.file_list {
        let text = std::fs::read_to_string(f).map_err(|e| format!("{f}: {e}"))?;
        v.extend(
            text.lines()
                .map(|l| l.trim().to_string())
                .filter(|l| !l.is_empty() && !l.starts_with('#')),
        );
    }
    if v.is_empty() {
        return Err("no input files".into());
    }
    Ok(v)
}
