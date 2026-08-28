pub mod concat;
pub mod merge;
pub mod view;

use std::io;

use crate::io::{path_implies_bgzf, Writer};

/// Output options shared by every subcommand.
#[derive(clap::Args, Debug, Clone)]
pub struct OutputOpts {
    /// Write to FILE instead of standard output
    #[arg(short = 'o', long = "output", value_name = "FILE")]
    pub output: Option<String>,

    /// Output type: v = uncompressed VCF, z = BGZF-compressed VCF
    #[arg(short = 'O', long = "output-type", value_name = "v|z")]
    pub output_type: Option<char>,

    /// BGZF compression level: 1-12 (libdeflate), or 0 for vcfr's built-in
    /// encoder — the fastest option, with a ratio between levels 1 and 2
    #[arg(short = 'l', long = "compression-level", default_value_t = 6, value_name = "N")]
    pub level: u32,

    /// Codec for BGZF output: "lib" (libdeflate, default) or "rust" — vcfr's
    /// pure-Rust encoders throughout: level 0 as always, levels 1-6 mapped to
    /// the high-effort per-block-Huffman encoder
    #[arg(long = "codec", value_name = "lib|rust", default_value = "lib")]
    pub codec: String,
}

impl OutputOpts {
    pub fn bgzf(&self) -> Result<bool, String> {
        match self.output_type {
            Some('z') | Some('b') => Ok(true),
            Some('v') | Some('u') => Ok(false),
            Some(c) => Err(format!("unknown output type '{c}' (expected v or z)")),
            None => Ok(path_implies_bgzf(self.output.as_deref())),
        }
    }

    pub fn writer(&self, threads: usize) -> Result<Writer, String> {
        let bgzf = self.bgzf()?;
        let rust_codec = match self.codec.as_str() {
            "lib" => false,
            "rust" => true,
            c => return Err(format!("unknown codec '{c}' (expected lib or rust)")),
        };
        if rust_codec && self.level > 6 {
            return Err("--codec rust supports levels 0-6".to_string());
        }
        Writer::create(self.output.as_deref(), bgzf, threads, self.level, rust_codec)
            .map_err(|e| e.to_string())
    }
}

/// Split a thread budget between the reader and the writer.
///
/// Both numbers count *worker* threads, so `--threads 1` really is one thread:
/// zero workers means the main thread does that half of the job itself.
///
/// Deflating costs several times more than inflating -- on this codebase's
/// benchmark, compressing a stream costs roughly six times what decompressing
/// it does -- so with BGZF output the writer is what saturates the machine and
/// gets the whole budget. Inflation workers are latency-hiding rather than
/// throughput-consuming: they spend most of a run blocked waiting for the
/// consumer, so charging them against the writer's share just starves it.
/// With plain output the writer needs no workers at all.
pub fn split_threads(total: usize, bgzf_out: bool) -> (usize, usize) {
    if total <= 1 {
        return (0, 0);
    }
    if !bgzf_out {
        return (total, 0);
    }
    ((total / 3).max(1), total)
}

/// Thread budget for a k-way merge.
///
/// Reader workers are latency-hiding rather than throughput-consuming: they
/// spend most of the run blocked, and their only job is to keep inflation off
/// the main thread, which is the merge's serial bottleneck. So every input gets
/// at least one, and they are not charged against the writer's share.
pub fn merge_threads(total: usize, inputs: usize, bgzf_out: bool) -> (usize, usize) {
    if total <= 1 {
        return (0, 0);
    }
    let per_file = ((total / 2) / inputs.max(1)).max(1);
    (per_file, if bgzf_out { total } else { 0 })
}

/// Resolve `--threads 0` to the machine's parallelism.
pub fn resolve_threads(t: usize) -> usize {
    if t == 0 {
        std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1)
    } else {
        t
    }
}

/// Turn a broken-pipe error (`vcfr view … | head`) into a clean exit.
pub fn ignore_broken_pipe(r: io::Result<()>) -> Result<(), String> {
    match r {
        Ok(()) => Ok(()),
        Err(e) if e.kind() == io::ErrorKind::BrokenPipe => Ok(()),
        Err(e) => Err(e.to_string()),
    }
}

/// A `write_all` call site whose surrounding function returns `Result<...,
/// String>` rather than an `io::Result` closure `ignore_broken_pipe` can wrap
/// wholesale — merge's write points sit inside a closure that also carries
/// stringified parse errors from unrelated sources, so the closure itself
/// can't be `io::Result`. Checking the error kind here, before it is
/// stringified and that information is lost, keeps a downstream `| head`
/// a silent, successful exit instead of a spurious "Broken pipe" failure.
pub fn write_or_broken_pipe(w: &mut crate::io::Writer, buf: &[u8]) -> Result<(), String> {
    use std::io::Write;
    match w.write_all(buf) {
        Ok(()) => Ok(()),
        Err(e) if e.kind() == io::ErrorKind::BrokenPipe => Ok(()),
        Err(e) => Err(e.to_string()),
    }
}

/// Emit the byte range `[a, b)` of `line`, splicing any interior compressed
/// block lying wholly inside the range through the BGZF writer untouched.
/// Text between splices accumulates in `buf`, which is flushed (becoming its
/// own, possibly short, deflate block) before each splice so stream order is
/// preserved. With no spliceable block in range this degenerates to
/// `buf.extend_from_slice(&line[a..b])`. Shared by `merge` and `view -s`.
pub(crate) fn splice_range(
    line: &[u8],
    blocks: &mut [crate::bgzf::spliced::InteriorBlock],
    a: usize,
    b: usize,
    buf: &mut Vec<u8>,
    w: &mut crate::io::Writer,
) -> io::Result<()> {
    use std::io::Write;
    let mut cur = a;
    for k in blocks.iter_mut() {
        if k.line_off < a || k.end() > b || k.raw.is_empty() {
            continue;
        }
        debug_assert!(k.line_off >= cur);
        buf.extend_from_slice(&line[cur..k.line_off]);
        if !buf.is_empty() {
            w.write_all(buf)?;
            buf.clear();
        }
        w.splice_block(std::mem::take(&mut k.raw))?;
        cur = k.end();
    }
    buf.extend_from_slice(&line[cur..b]);
    Ok(())
}

/// Same treatment for `Writer::write_raw_block`, used by `concat --naive`'s
/// block-copy loop.
pub fn write_raw_block_or_broken_pipe(w: &mut crate::io::Writer, block: Vec<u8>) -> Result<(), String> {
    match w.write_raw_block(block) {
        Ok(()) => Ok(()),
        Err(e) if e.kind() == io::ErrorKind::BrokenPipe => Ok(()),
        Err(e) => Err(e.to_string()),
    }
}
