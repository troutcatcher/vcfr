//! Multi-threaded BGZF compression.

use std::cell::RefCell;
use std::io::{self, Write};
use std::sync::Arc;

use super::{deflate_block, deflate_block_high, deflate_block_preset, EOF_BLOCK, UNCOMPRESSED_BLOCK_SIZE};
use crate::deflate::{CodeSet, HighMatcher, Matcher};
use crate::pool::OrderedPool;

thread_local! {
    static MATCHER: RefCell<Matcher> = RefCell::new(Matcher::default());
    static HIGH_MATCHER: RefCell<HighMatcher> = RefCell::new(HighMatcher::default());
}

/// Which encoder turns a block into deflate bits.
enum Engine {
    /// libdeflate at a given level (1-12).
    Lib(i32),
    /// vcfr's preset-code encoder (level 0). Codes are trained on early
    /// blocks and refreshed occasionally; every job holds an Arc to the set
    /// it was submitted with, so retraining never affects blocks in flight.
    Preset { codes: Option<Arc<CodeSet>>, blocks: u64 },
    /// vcfr's high-effort encoder (`--codec rust`, levels 1-6): per-block
    /// optimal Huffman over a lazy hash-chain matcher.
    High { chain: usize, nice: usize, lazy: usize },
}

pub struct BgzfWriter<W: Write> {
    inner: W,
    /// `None` when compressing inline on the calling thread.
    pool: Option<OrderedPool<Vec<u8>>>,
    inline: Option<libdeflater::Compressor>,
    engine: Engine,
    scratch: Vec<u8>,
    buf: Vec<u8>,
    finished: bool,
}

impl<W: Write> BgzfWriter<W> {
    /// `workers` compression threads; `0` deflates inline on the calling
    /// thread, so a caller asking for one thread really gets one thread.
    /// Level 0 selects vcfr's own preset-code encoder; 1-12 are libdeflate
    /// unless `rust_codec` routes them to the built-in high-effort encoder.
    pub fn new(inner: W, workers: usize, level: u32, rust_codec: bool) -> Self {
        let engine = if level == 0 {
            Engine::Preset { codes: None, blocks: 0 }
        } else if rust_codec {
            let (chain, nice, lazy) = crate::deflate::effort_for_level(level);
            Engine::High { chain, nice, lazy }
        } else {
            Engine::Lib(level.clamp(1, 12) as i32)
        };
        BgzfWriter {
            inner,
            // Depth is blocks in flight between the feeding thread and the
            // deflate workers. At 3 per worker (~0.75 MiB on 4 workers) any
            // hiccup on the feeding side starves them: measured 2.6 of 4 cores
            // at level 1 and 3.7 at level 6. At 16 per worker (~4 MiB) both
            // reach ~3.85; 48 buys nothing more.
            pool: (workers >= 1).then(|| OrderedPool::new(workers, workers * 16)),
            inline: (workers == 0 && level != 0 && !rust_codec).then(|| {
                libdeflater::Compressor::new(
                    libdeflater::CompressionLvl::new(level.clamp(1, 12) as i32)
                        .expect("valid compression level"),
                )
            }),
            engine,
            scratch: Vec::with_capacity(UNCOMPRESSED_BLOCK_SIZE),
            buf: Vec::with_capacity(UNCOMPRESSED_BLOCK_SIZE),
            finished: false,
        }
    }

    /// Drain finished blocks, optionally waiting until the pool has room.
    fn drain(&mut self, until_empty: bool) -> io::Result<()> {
        let pool = match &mut self.pool {
            Some(p) => p,
            None => return Ok(()),
        };
        while pool.outstanding() >= pool.capacity() || (until_empty && pool.outstanding() > 0) {
            match pool.next() {
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
        let level = match &mut self.engine {
            Engine::Lib(l) => *l,
            Engine::Preset { codes, blocks } => {
                // Train on the first block, retrain once clear of the VCF
                // header, then occasionally in case the data drifts. Jobs in
                // flight keep the Arc they were submitted with.
                if codes.is_none() || *blocks == 8 || *blocks % 65536 == 0 {
                    let set = MATCHER
                        .with(|m| CodeSet::train(&self.buf, &mut m.borrow_mut()));
                    *codes = Some(Arc::new(set));
                }
                *blocks += 1;
                let codes = Arc::clone(codes.as_ref().unwrap());
                if self.pool.is_none() {
                    MATCHER.with(|m| {
                        deflate_block_preset(&codes, &mut m.borrow_mut(), &self.buf, &mut self.scratch)
                    });
                    self.inner.write_all(&self.scratch)?;
                    self.buf.clear();
                    return Ok(());
                }
                let data =
                    std::mem::replace(&mut self.buf, Vec::with_capacity(UNCOMPRESSED_BLOCK_SIZE));
                self.pool.as_mut().unwrap().submit(move || {
                    MATCHER.with(|m| {
                        let mut out = Vec::with_capacity(data.len() / 2 + 64);
                        deflate_block_preset(&codes, &mut m.borrow_mut(), &data, &mut out);
                        out
                    })
                });
                return self.drain(false);
            }
            Engine::High { chain, nice, lazy } => {
                let (chain, nice, lazy) = (*chain, *nice, *lazy);
                if self.pool.is_none() {
                    HIGH_MATCHER.with(|m| {
                        let mut m = m.borrow_mut();
                        m.set_effort(chain, nice, lazy);
                        deflate_block_high(&mut m, &self.buf, &mut self.scratch);
                    });
                    self.inner.write_all(&self.scratch)?;
                    self.buf.clear();
                    return Ok(());
                }
                let data =
                    std::mem::replace(&mut self.buf, Vec::with_capacity(UNCOMPRESSED_BLOCK_SIZE));
                self.pool.as_mut().unwrap().submit(move || {
                    HIGH_MATCHER.with(|m| {
                        let mut m = m.borrow_mut();
                        m.set_effort(chain, nice, lazy);
                        let mut out = Vec::with_capacity(data.len() / 2 + 64);
                        deflate_block_high(&mut m, &data, &mut out);
                        out
                    })
                });
                return self.drain(false);
            }
        };
        if let Some(c) = &mut self.inline {
            deflate_block(c, &self.buf, &mut self.scratch);
            self.inner.write_all(&self.scratch)?;
            self.buf.clear();
            return Ok(());
        }
        let data = std::mem::replace(&mut self.buf, Vec::with_capacity(UNCOMPRESSED_BLOCK_SIZE));
        self.pool.as_mut().expect("pool or inline compressor").submit(move || {
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

    /// Splice an already-formed BGZF block into the stream *through* the
    /// ordered pool: pending buffered text is flushed as its own (possibly
    /// short) block, then the raw block travels as an identity job, keeping
    /// its place among in-flight deflate jobs without stalling them the way
    /// `write_raw_block`'s full drain does. Used by merge's block splicing,
    /// which does this hundreds of times per output line.
    pub fn splice_block(&mut self, block: Vec<u8>) -> io::Result<()> {
        self.flush_block()?;
        match &mut self.pool {
            Some(pool) => {
                pool.submit(move || block);
                self.drain(false)
            }
            None => self.inner.write_all(&block),
        }
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
