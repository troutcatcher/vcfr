//! Zero-copy line splitting over a chunked byte source.
//!
//! `LineReader` hands out `&[u8]` slices that point straight into the buffer
//! produced by the underlying source; a line is only ever copied when it happens
//! to straddle a chunk boundary.

use std::io::{self, Read};

/// A source of successive byte chunks (a decompressed BGZF stream, a plain file, …).
pub trait BufSource {
    fn next_chunk(&mut self) -> io::Result<Option<Vec<u8>>>;
    /// Hand a spent buffer back so it can be refilled instead of reallocated.
    fn recycle(&mut self, _buf: Vec<u8>) {}
}

const PLAIN_CHUNK: usize = 512 * 1024;

/// `BufSource` over any `Read` (uncompressed input).
pub struct PlainSource<R: Read> {
    inner: R,
    spare: Option<Vec<u8>>,
}

impl<R: Read> PlainSource<R> {
    pub fn new(inner: R) -> Self {
        PlainSource { inner, spare: None }
    }
}

impl<R: Read> BufSource for PlainSource<R> {
    fn next_chunk(&mut self) -> io::Result<Option<Vec<u8>>> {
        let mut buf = self.spare.take().unwrap_or_else(|| vec![0u8; PLAIN_CHUNK]);
        buf.resize(PLAIN_CHUNK, 0);
        let mut n = 0;
        while n < buf.len() {
            match self.inner.read(&mut buf[n..])? {
                0 => break,
                k => n += k,
            }
        }
        if n == 0 {
            self.spare = Some(buf);
            return Ok(None);
        }
        buf.truncate(n);
        Ok(Some(buf))
    }

    fn recycle(&mut self, buf: Vec<u8>) {
        if buf.capacity() >= PLAIN_CHUNK {
            self.spare = Some(buf);
        }
    }
}

pub struct LineReader<S: BufSource> {
    src: S,
    buf: Vec<u8>,
    pos: usize,
    /// Holds a line that spans two chunks.
    carry: Vec<u8>,
    use_carry: bool,
    start: usize,
    end: usize,
    eof: bool,
}

impl<S: BufSource> LineReader<S> {
    pub fn new(src: S) -> Self {
        LineReader {
            src,
            buf: Vec::new(),
            pos: 0,
            carry: Vec::new(),
            use_carry: false,
            start: 0,
            end: 0,
            eof: false,
        }
    }

    /// The line produced by the last successful `advance`, without its newline.
    #[inline]
    pub fn line(&self) -> &[u8] {
        if self.use_carry {
            &self.carry
        } else {
            &self.buf[self.start..self.end]
        }
    }

    /// Read the next line into internal state. Returns `false` at end of input.
    pub fn advance(&mut self) -> io::Result<bool> {
        self.use_carry = false;
        self.carry.clear();
        loop {
            if self.pos < self.buf.len() {
                match memchr::memchr(b'\n', &self.buf[self.pos..]) {
                    Some(rel) => {
                        let s = self.pos;
                        let nl = self.pos + rel;
                        self.pos = nl + 1;
                        if self.carry.is_empty() {
                            let mut e = nl;
                            if e > s && self.buf[e - 1] == b'\r' {
                                e -= 1;
                            }
                            self.start = s;
                            self.end = e;
                        } else {
                            self.carry.extend_from_slice(&self.buf[s..nl]);
                            if self.carry.last() == Some(&b'\r') {
                                self.carry.pop();
                            }
                            self.use_carry = true;
                        }
                        return Ok(true);
                    }
                    None => {
                        self.carry.extend_from_slice(&self.buf[self.pos..]);
                        self.pos = self.buf.len();
                    }
                }
            }
            if self.eof {
                // Trailing data with no final newline.
                if !self.carry.is_empty() {
                    if self.carry.last() == Some(&b'\r') {
                        self.carry.pop();
                    }
                    self.use_carry = true;
                    return Ok(true);
                }
                return Ok(false);
            }
            let old = std::mem::take(&mut self.buf);
            if !old.is_empty() {
                self.src.recycle(old);
            }
            match self.src.next_chunk()? {
                Some(b) => {
                    self.buf = b;
                    self.pos = 0;
                }
                None => self.eof = true,
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Hands out the chunks it was built with, so line boundaries can be made
    /// to fall anywhere.
    struct Chunks(Vec<Vec<u8>>);
    impl BufSource for Chunks {
        fn next_chunk(&mut self) -> io::Result<Option<Vec<u8>>> {
            if self.0.is_empty() {
                Ok(None)
            } else {
                Ok(Some(self.0.remove(0)))
            }
        }
    }

    fn lines_of(chunks: &[&str]) -> Vec<String> {
        let mut r = LineReader::new(Chunks(chunks.iter().map(|c| c.as_bytes().to_vec()).collect()));
        let mut out = Vec::new();
        while r.advance().unwrap() {
            out.push(String::from_utf8(r.line().to_vec()).unwrap());
        }
        out
    }

    #[test]
    fn splits_lines_in_one_chunk() {
        assert_eq!(lines_of(&["a\nbb\nccc\n"]), ["a", "bb", "ccc"]);
    }

    #[test]
    fn joins_lines_across_chunks() {
        assert_eq!(lines_of(&["ab", "cd\ne", "f\n"]), ["abcd", "ef"]);
        assert_eq!(lines_of(&["abc", "", "def\n"]), ["abcdef"]);
    }

    #[test]
    fn handles_missing_final_newline() {
        assert_eq!(lines_of(&["a\nb"]), ["a", "b"]);
        assert_eq!(lines_of(&["a\nb", "c"]), ["a", "bc"]);
    }

    #[test]
    fn strips_carriage_returns_even_across_chunks() {
        assert_eq!(lines_of(&["a\r\nb\r\n"]), ["a", "b"]);
        assert_eq!(lines_of(&["ab\r", "\ncd\r\n"]), ["ab", "cd"]);
        assert_eq!(lines_of(&["ab\r"]), ["ab"]);
    }

    #[test]
    fn empty_input_yields_nothing() {
        assert!(lines_of(&[]).is_empty());
        assert_eq!(lines_of(&["\n\n"]), ["", ""]);
    }
}
