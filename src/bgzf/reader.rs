//! Multi-threaded BGZF decompression.

use std::io::{self, Read};

use super::{RawBlock, RawBlockReader};
use crate::lines::BufSource;
use crate::pool::OrderedPool;

/// How many blocks are inflated by a single job. Bigger batches amortise the
/// channel traffic; 32 blocks is ~2 MiB of output per job.
const BLOCKS_PER_JOB: usize = 32;

pub struct BgzfReader<R: Read + Send + 'static> {
    blocks: Option<RawBlockReader<R>>,
    pool: OrderedPool<io::Result<Vec<u8>>>,
    exhausted: bool,
}

impl<R: Read + Send + 'static> BgzfReader<R> {
    pub fn new(inner: R, threads: usize) -> Self {
        BgzfReader {
            blocks: Some(RawBlockReader::new(inner)),
            pool: OrderedPool::new(threads, threads * 3),
            exhausted: false,
        }
    }

    /// Queue jobs until the pool is full or the input runs out.
    fn pump(&mut self) -> io::Result<()> {
        while !self.exhausted && self.pool.outstanding() < self.pool.capacity() {
            let src = self.blocks.as_mut().unwrap();
            let mut batch: Vec<RawBlock> = Vec::with_capacity(BLOCKS_PER_JOB);
            let mut total = 0usize;
            for _ in 0..BLOCKS_PER_JOB {
                match src.next_block()? {
                    Some(b) => {
                        total += b.isize as usize;
                        batch.push(b);
                    }
                    None => {
                        self.exhausted = true;
                        break;
                    }
                }
            }
            if batch.is_empty() {
                break;
            }
            self.pool.submit(move || {
                let mut d = libdeflater::Decompressor::new();
                let mut out = Vec::with_capacity(total);
                for b in &batch {
                    b.inflate_into(&mut d, &mut out)?;
                }
                Ok(out)
            });
        }
        Ok(())
    }
}

impl<R: Read + Send + 'static> BufSource for BgzfReader<R> {
    fn next_chunk(&mut self) -> io::Result<Option<Vec<u8>>> {
        loop {
            self.pump()?;
            match self.pool.next() {
                // A batch of pure EOF blocks decompresses to nothing; skip it
                // rather than signalling end of input to the caller.
                Some(Ok(v)) if v.is_empty() => continue,
                Some(Ok(v)) => return Ok(Some(v)),
                Some(Err(e)) => return Err(e),
                None => return Ok(None),
            }
        }
    }
}

/// Single-threaded BGZF decompression, used where the extra threads would not
/// pay for themselves (reading just a header, or many files at once in `merge`).
pub struct SerialBgzfReader<R: Read> {
    blocks: RawBlockReader<R>,
    d: libdeflater::Decompressor,
}

impl<R: Read> SerialBgzfReader<R> {
    pub fn new(inner: R) -> Self {
        SerialBgzfReader {
            blocks: RawBlockReader::new(inner),
            d: libdeflater::Decompressor::new(),
        }
    }
}

impl<R: Read> BufSource for SerialBgzfReader<R> {
    fn next_chunk(&mut self) -> io::Result<Option<Vec<u8>>> {
        let mut out = Vec::new();
        loop {
            let mut any = false;
            for _ in 0..BLOCKS_PER_JOB {
                match self.blocks.next_block()? {
                    Some(b) => {
                        any = true;
                        b.inflate_into(&mut self.d, &mut out)?;
                    }
                    None => break,
                }
            }
            if !out.is_empty() {
                return Ok(Some(out));
            }
            // A run of empty blocks (a mid-stream BGZF EOF marker) is not the
            // end of the input; keep going until a block actually yields bytes.
            if !any {
                return Ok(None);
            }
        }
    }
}
