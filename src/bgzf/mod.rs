//! BGZF (blocked gzip) support with multi-threaded (de)compression.
//!
//! BGZF is gzip with an extra `BC` subfield that records the size of each
//! self-contained deflate block. Because every block stands alone, blocks can be
//! inflated and deflated on independent threads and, for `concat --naive`, copied
//! verbatim without touching the deflate stream at all.

pub mod reader;
pub mod writer;

use std::io::{self, Read};

#[allow(unused_imports)]
pub use reader::BgzfReader;
pub use writer::BgzfWriter;

/// Fixed part of a BGZF block header (before the extra field).
pub const GZIP_HEADER_LEN: usize = 12;
/// Largest legal BGZF block, header and footer included.
pub const MAX_BLOCK_SIZE: usize = 65536;
/// Uncompressed payload per block we emit. Matches htslib and leaves room for
/// the worst-case deflate expansion inside a 64 KiB block.
pub const UNCOMPRESSED_BLOCK_SIZE: usize = 0xff00;

/// The 28-byte empty block every BGZF file ends with.
pub const EOF_BLOCK: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

fn invalid(msg: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg.to_string())
}

/// One BGZF block, still compressed.
pub struct RawBlock {
    /// Every byte of the block, exactly as it appeared in the file.
    pub bytes: Vec<u8>,
    /// Offset of the deflate stream within `bytes`.
    pub payload_start: usize,
    pub payload_end: usize,
    /// Size of the block once inflated.
    pub isize: u32,
}

impl RawBlock {
    #[inline]
    pub fn payload(&self) -> &[u8] {
        &self.bytes[self.payload_start..self.payload_end]
    }

    pub fn inflate_into(&self, d: &mut libdeflater::Decompressor, out: &mut Vec<u8>) -> io::Result<()> {
        if self.isize == 0 {
            return Ok(());
        }
        let base = out.len();
        out.resize(base + self.isize as usize, 0);
        let n = d
            .deflate_decompress(self.payload(), &mut out[base..])
            .map_err(|e| invalid(&format!("corrupt BGZF block: {e}")))?;
        if n != self.isize as usize {
            return Err(invalid("BGZF block size mismatch"));
        }
        Ok(())
    }
}

/// Reads BGZF blocks off a byte stream without decompressing them.
pub struct RawBlockReader<R: Read> {
    inner: R,
}

impl<R: Read> RawBlockReader<R> {
    pub fn new(inner: R) -> Self {
        RawBlockReader { inner }
    }

    /// Reads exactly `n` bytes, returning `false` on a clean end of stream.
    fn read_exact_eof(&mut self, buf: &mut [u8]) -> io::Result<bool> {
        let mut n = 0;
        while n < buf.len() {
            match self.inner.read(&mut buf[n..])? {
                0 => {
                    if n == 0 {
                        return Ok(false);
                    }
                    return Err(invalid("truncated BGZF block"));
                }
                k => n += k,
            }
        }
        Ok(true)
    }

    pub fn next_block(&mut self) -> io::Result<Option<RawBlock>> {
        let mut bytes = vec![0u8; GZIP_HEADER_LEN];
        if !self.read_exact_eof(&mut bytes)? {
            return Ok(None);
        }
        if bytes[0] != 0x1f || bytes[1] != 0x8b || bytes[2] != 8 {
            return Err(invalid("not a gzip stream"));
        }
        if bytes[3] & 0x04 == 0 {
            return Err(invalid("gzip stream is not BGZF (no extra field)"));
        }
        let xlen = u16::from_le_bytes([bytes[10], bytes[11]]) as usize;
        bytes.resize(GZIP_HEADER_LEN + xlen, 0);
        if !self.read_exact_eof(&mut bytes[GZIP_HEADER_LEN..])? {
            return Err(invalid("truncated BGZF extra field"));
        }
        let bsize = find_bsize(&bytes[GZIP_HEADER_LEN..]).ok_or_else(|| invalid("BGZF BC subfield missing"))?;
        let total = bsize as usize + 1;
        if total < GZIP_HEADER_LEN + xlen + 8 {
            return Err(invalid("nonsensical BGZF block size"));
        }
        let head = bytes.len();
        bytes.resize(total, 0);
        if !self.read_exact_eof(&mut bytes[head..])? {
            return Err(invalid("truncated BGZF block body"));
        }
        let isize = u32::from_le_bytes([
            bytes[total - 4],
            bytes[total - 3],
            bytes[total - 2],
            bytes[total - 1],
        ]);
        Ok(Some(RawBlock {
            payload_start: head,
            payload_end: total - 8,
            bytes,
            isize,
        }))
    }
}

/// Scans the gzip extra field for the `BC` subfield holding `BSIZE - 1`.
fn find_bsize(extra: &[u8]) -> Option<u16> {
    let mut i = 0;
    while i + 4 <= extra.len() {
        let si1 = extra[i];
        let si2 = extra[i + 1];
        let slen = u16::from_le_bytes([extra[i + 2], extra[i + 3]]) as usize;
        if i + 4 + slen > extra.len() {
            return None;
        }
        if si1 == b'B' && si2 == b'C' && slen == 2 {
            return Some(u16::from_le_bytes([extra[i + 4], extra[i + 5]]));
        }
        i += 4 + slen;
    }
    None
}

/// True if `head` (at least 18 bytes of a file) looks like a BGZF block header.
pub fn is_bgzf_header(head: &[u8]) -> bool {
    head.len() >= 18
        && head[0] == 0x1f
        && head[1] == 0x8b
        && head[2] == 8
        && head[3] & 0x04 != 0
        && {
            let xlen = u16::from_le_bytes([head[10], head[11]]) as usize;
            head.len() >= 12 + xlen && find_bsize(&head[12..12 + xlen]).is_some()
        }
}

pub fn is_gzip_header(head: &[u8]) -> bool {
    head.len() >= 2 && head[0] == 0x1f && head[1] == 0x8b
}

/// Deflate `data` into a complete BGZF block using vcfr's own preset-code
/// encoder (`src/deflate`). CRC comes from crc32fast rather than libdeflate.
pub fn deflate_block_preset(
    codes: &crate::deflate::CodeSet,
    m: &mut crate::deflate::Matcher,
    data: &[u8],
    out: &mut Vec<u8>,
) {
    debug_assert!(data.len() <= UNCOMPRESSED_BLOCK_SIZE);
    out.clear();
    out.extend_from_slice(&[
        0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0x00, b'B', b'C', 0x02, 0x00, 0, 0,
    ]);
    crate::deflate::compress_into(codes, m, data, out);
    let crc = crc32fast::hash(data);
    out.extend_from_slice(&crc.to_le_bytes());
    out.extend_from_slice(&(data.len() as u32).to_le_bytes());
    let total = out.len();
    debug_assert!(total <= MAX_BLOCK_SIZE);
    out[16..18].copy_from_slice(&((total - 1) as u16).to_le_bytes());
}

/// Build a complete BGZF block with the high-effort per-block-Huffman encoder.
pub fn deflate_block_high(m: &mut crate::deflate::HighMatcher, data: &[u8], out: &mut Vec<u8>) {
    debug_assert!(data.len() <= UNCOMPRESSED_BLOCK_SIZE);
    out.clear();
    out.extend_from_slice(&[
        0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0x00, b'B', b'C', 0x02, 0x00, 0, 0,
    ]);
    crate::deflate::compress_high_into(m, data, out);
    let crc = crc32fast::hash(data);
    out.extend_from_slice(&crc.to_le_bytes());
    out.extend_from_slice(&(data.len() as u32).to_le_bytes());
    let total = out.len();
    debug_assert!(total <= MAX_BLOCK_SIZE);
    out[16..18].copy_from_slice(&((total - 1) as u16).to_le_bytes());
}

/// Deflate `data` into a complete BGZF block.
pub fn deflate_block(c: &mut libdeflater::Compressor, data: &[u8], out: &mut Vec<u8>) {
    debug_assert!(data.len() <= UNCOMPRESSED_BLOCK_SIZE);
    let header_len = 18;
    out.clear();
    out.resize(header_len, 0);
    out[..header_len].copy_from_slice(&[
        0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0x00, b'B', b'C', 0x02, 0x00, 0, 0,
    ]);
    let bound = c.deflate_compress_bound(data.len());
    out.resize(header_len + bound, 0);
    let n = match c.deflate_compress(data, &mut out[header_len..]) {
        Ok(n) => n,
        // Incompressible input: fall back to stored deflate blocks, which is
        // what `deflate_compress_bound` accounts for.
        Err(e) => panic!("deflate failed: {e}"),
    };
    out.truncate(header_len + n);
    let crc = libdeflater::crc32(data);
    out.extend_from_slice(&crc.to_le_bytes());
    out.extend_from_slice(&(data.len() as u32).to_le_bytes());
    let total = out.len();
    debug_assert!(total <= MAX_BLOCK_SIZE);
    let bsize = (total - 1) as u16;
    out[16..18].copy_from_slice(&bsize.to_le_bytes());
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lines::BufSource;
    use std::io::Write;

    fn roundtrip(data: &[u8], threads: usize, level: u32) -> Vec<u8> {
        let mut buf = Vec::new();
        {
            let mut w = writer::BgzfWriter::new(&mut buf, threads, level, false);
            w.write_all(data).unwrap();
            w.finish().unwrap();
        }
        let mut r = reader::BgzfReader::new(std::io::Cursor::new(buf.clone()), threads);
        let mut got = Vec::new();
        while let Some(chunk) = r.next_chunk().unwrap() {
            got.extend_from_slice(&chunk);
        }
        assert_eq!(got, data, "round trip mismatch");
        buf
    }

    #[test]
    fn preset_engine_roundtrips_and_interoperates() {
        // Level 0 = vcfr's own encoder. The stream must come back through the
        // BGZF reader (which inflates with libdeflate) byte-for-byte.
        let mut data = Vec::new();
        let mut x = 7u64;
        while data.len() < 5 * UNCOMPRESSED_BLOCK_SIZE + 1234 {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            data.extend_from_slice(
                format!("chr2\t{}\trs{}\tC\tT\t99\tPASS\tDP=42\tGT\t0/1\n", x % 999999, x % 777)
                    .as_bytes(),
            );
        }
        let mut buf = Vec::new();
        {
            let mut w = writer::BgzfWriter::new(&mut buf, 2, 0, false);
            w.write_all(&data).unwrap();
            w.finish().unwrap();
        }
        assert!(is_bgzf_header(&buf[..64]));
        assert!(buf.ends_with(&EOF_BLOCK));
        let mut r = reader::BgzfReader::new(std::io::Cursor::new(buf), 2);
        let mut got = Vec::new();
        while let Some(chunk) = r.next_chunk().unwrap() {
            got.extend_from_slice(&chunk);
        }
        assert_eq!(got, data);
    }

    #[test]
    fn roundtrips_across_block_boundaries() {
        for len in [0usize, 1, 1024, UNCOMPRESSED_BLOCK_SIZE - 1, UNCOMPRESSED_BLOCK_SIZE, UNCOMPRESSED_BLOCK_SIZE + 1, 5 * UNCOMPRESSED_BLOCK_SIZE + 7] {
            let data: Vec<u8> = (0..len).map(|i| (i % 251) as u8).collect();
            for threads in [1usize, 4] {
                roundtrip(&data, threads, 6);
            }
        }
    }

    #[test]
    fn output_is_a_valid_bgzf_stream() {
        let data: Vec<u8> = (0..200_000).map(|i| (i % 97) as u8).collect();
        let encoded = roundtrip(&data, 3, 1);
        assert!(is_bgzf_header(&encoded[..64]));
        assert!(encoded.ends_with(&EOF_BLOCK), "must end with the EOF marker");
        // Every block must declare a size that walks us exactly to the end.
        let mut off = 0usize;
        let mut blocks = 0;
        while off < encoded.len() {
            let xlen = u16::from_le_bytes([encoded[off + 10], encoded[off + 11]]) as usize;
            let bsize = find_bsize(&encoded[off + 12..off + 12 + xlen]).unwrap() as usize + 1;
            assert!(bsize <= MAX_BLOCK_SIZE);
            off += bsize;
            blocks += 1;
        }
        assert_eq!(off, encoded.len());
        assert!(blocks > 3, "expected several blocks, got {blocks}");
    }

    #[test]
    fn incompressible_data_still_fits_one_block() {
        // Pseudo-random bytes deflate to slightly more than their input size.
        let mut x = 0x2545F4914F6CDD1Du64;
        let data: Vec<u8> = (0..UNCOMPRESSED_BLOCK_SIZE)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                x as u8
            })
            .collect();
        roundtrip(&data, 1, 12);
    }
}
