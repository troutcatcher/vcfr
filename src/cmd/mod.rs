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

    /// BGZF compression level (1-12)
    #[arg(short = 'l', long = "compression-level", default_value_t = 6, value_name = "N")]
    pub level: u32,
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
        Writer::create(self.output.as_deref(), bgzf, threads, self.level).map_err(|e| e.to_string())
    }
}

/// Split a thread budget between the reader and the writer.
///
/// Both numbers count *worker* threads, so `--threads 1` really is one thread:
/// zero workers means the main thread does that half of the job itself.
///
/// Deflating costs several times more than inflating, so when the output is
/// BGZF most of the budget goes to the writer; with plain output the writer
/// needs none at all.
pub fn split_threads(total: usize, bgzf_out: bool) -> (usize, usize) {
    if total <= 1 {
        return (0, 0);
    }
    if !bgzf_out {
        return (total, 0);
    }
    let read = (total / 3).max(1);
    (read, (total - read).max(1))
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
