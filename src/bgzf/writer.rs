//! Multi-threaded BGZF compression.

use std::io::{self, Write};

use super::{deflate_block, EOF_BLOCK, UNCOMPRESSED_BLOCK_SIZE};
use crate::pool::OrderedPool;

pub struct BgzfWriter<W: Write> {
    inner: W,
    pool: OrderedPool<Vec<u8>>,
    buf: Vec<u8>,
    level: i32,
    finished: bool,
}

impl<W: Write> BgzfWriter<W> {
    pub fn new(inner: W, threads: usize, level: u32) -> Self {
        BgzfWriter {
            inner,
            pool: OrderedPool::new(threads, threads * 3),
            buf: Vec::with_capacity(UNCOMPRESSED_BLOCK_SIZE),
            level: level.clamp(1, 12) as i32,
            finished: false,
        }
    }

    /// Drain finished blocks, optionally waiting until the pool has room.
    fn drain(&mut self, until_empty: bool) -> io::Result<()> {
        while self.pool.outstanding() >= self.pool.capacity()
            || (until_empty && self.pool.outstanding() > 0)
        {
            match self.pool.next() {
                Some(b) => self.inner.write_all(&b)?,
                None => break,
            }
        }
        Ok(())
    }

    fn flush_block(&mut self) -> io::Result<()> {
        if self.buf.is_empty() {
            return Ok(());
        }
        let data = std::mem::replace(&mut self.buf, Vec::with_capacity(UNCOMPRESSED_BLOCK_SIZE));
        let level = self.level;
        self.pool.submit(move || {
            let mut c = libdeflater::Compressor::new(
                libdeflater::CompressionLvl::new(level).expect("valid compression level"),
            );
            let mut out = Vec::with_capacity(data.len() / 2 + 64);
            deflate_block(&mut c, &data, &mut out);
            out
        });
        self.drain(false)
    }

    /// Append an already-formed BGZF block verbatim. Used by `concat --naive`.
    pub fn write_raw_block(&mut self, block: Vec<u8>) -> io::Result<()> {
        self.flush_block()?;
        self.drain(true)?;
        self.inner.write_all(&block)
    }

    /// Finish the stream: flush pending data and append the BGZF EOF block.
    pub fn finish(&mut self) -> io::Result<()> {
        if self.finished {
            return Ok(());
        }
        self.finished = true;
        self.flush_block()?;
        self.drain(true)?;
        self.inner.write_all(&EOF_BLOCK)?;
        self.inner.flush()
    }
}

impl<W: Write> Write for BgzfWriter<W> {
    fn write(&mut self, mut data: &[u8]) -> io::Result<usize> {
        let n = data.len();
        while !data.is_empty() {
            let room = UNCOMPRESSED_BLOCK_SIZE - self.buf.len();
            let take = room.min(data.len());
            self.buf.extend_from_slice(&data[..take]);
            data = &data[take..];
            if self.buf.len() == UNCOMPRESSED_BLOCK_SIZE {
                self.flush_block()?;
            }
        }
        Ok(n)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.flush_block()?;
        self.drain(true)?;
        self.inner.flush()
    }
}

impl<W: Write> Drop for BgzfWriter<W> {
    fn drop(&mut self) {
        let _ = self.finish();
    }
}
