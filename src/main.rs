//! vcfr — a fast VCF toolkit.

mod bgzf;
mod cmd;
mod io;
mod lines;
mod pool;
mod vcf;

use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(
    name = "vcfr",
    version,
    about = "Fast VCF merging, concatenation and subsetting",
    long_about = "vcfr reads plain, gzip and BGZF VCFs, decompresses and recompresses \
BGZF across all cores, and only parses the record fields a command actually needs."
)]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Subset samples and variants, and convert between VCF and BGZF
    View(cmd::view::ViewArgs),
    /// Concatenate files that hold the same samples
    Concat(cmd::concat::ConcatArgs),
    /// Merge files that hold different samples
    Merge(cmd::merge::MergeArgs),
}

fn main() {
    let cli = Cli::parse();
    let r = match &cli.cmd {
        Cmd::View(a) => cmd::view::run(a),
        Cmd::Concat(a) => cmd::concat::run(a),
        Cmd::Merge(a) => cmd::merge::run(a),
    };
    if let Err(e) = r {
        eprintln!("vcfr: {e}");
        std::process::exit(1);
    }
}
