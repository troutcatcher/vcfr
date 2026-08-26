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
