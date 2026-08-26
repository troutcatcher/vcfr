//! Deterministic VCF generator used by the test and benchmark scripts.
//!
//!   cargo run --release --example gen -- --sites 100000 --samples 500 > big.vcf
//!
//! Sites are derived from `--site-seed` alone, so two runs that differ only in
//! `--gt-seed`/`--sample-prefix` describe the same variants for different
//! cohorts — exactly the input `merge` expects.

use std::io::{BufWriter, Write};

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        // splitmix64
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
}

fn arg(name: &str, default: &str) -> String {
    let args: Vec<String> = std::env::args().collect();
    for (i, a) in args.iter().enumerate() {
        if a == name {
            return args.get(i + 1).cloned().unwrap_or_else(|| default.to_string());
        }
    }
    default.to_string()
}

const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

fn main() {
    let sites: u64 = arg("--sites", "100000").parse().unwrap();
    let samples: usize = arg("--samples", "100").parse().unwrap();
    let contigs: usize = arg("--contigs", "4").parse().unwrap();
    let site_seed: u64 = arg("--site-seed", "1").parse().unwrap();
    let gt_seed: u64 = arg("--gt-seed", "42").parse().unwrap();
    let prefix = arg("--sample-prefix", "S");
    // Fraction (in percent) of sites to omit, so merged inputs overlap partially.
    let drop_pct: u64 = arg("--drop-pct", "0").parse().unwrap();
    let drop_seed: u64 = arg("--drop-seed", "7").parse().unwrap();

    let mut out = BufWriter::with_capacity(1 << 20, std::io::stdout().lock());
    let contig_len = 250_000_000u64;

    writeln!(out, "##fileformat=VCFv4.2").unwrap();
    writeln!(out, "##FILTER=<ID=PASS,Description=\"All filters passed\">").unwrap();
    writeln!(out, "##FILTER=<ID=LowQual,Description=\"Low quality\">").unwrap();
    for c in 1..=contigs {
        writeln!(out, "##contig=<ID=chr{c},length={contig_len}>").unwrap();
    }
    writeln!(out, "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">").unwrap();
    writeln!(out, "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">").unwrap();
    writeln!(out, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total depth\">").unwrap();
    writeln!(out, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">").unwrap();
    writeln!(out, "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">").unwrap();
    writeln!(out, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">").unwrap();
    writeln!(out, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">").unwrap();
    write!(out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT").unwrap();
    for i in 0..samples {
        write!(out, "\t{prefix}{i:05}").unwrap();
    }
    writeln!(out).unwrap();

    let mut gt_rng = Rng(gt_seed);
    let mut drop_rng = Rng(drop_seed);
    let per_contig = sites / contigs as u64;
    let mut line: Vec<u8> = Vec::with_capacity(64 + samples * 24);

    for c in 0..contigs {
        let mut pos = 0u64;
        let mut srng = Rng(site_seed.wrapping_mul(1_000_003).wrapping_add(c as u64));
        for _ in 0..per_contig {
            pos += 1 + srng.below(400);
            let r = srng.next();
            let ref_i = (r & 3) as usize;
            let refb = BASES[ref_i];
            let kind = r >> 8 & 0xff;
            let (reference, alt): (Vec<u8>, Vec<u8>) = if kind < 205 {
                // SNP, occasionally triallelic. ALT alleles are picked from the
                // bases the site does not already use so they are always distinct.
                let a1_i = (ref_i + 1 + ((r >> 16) % 3) as usize) % 4;
                let a1 = BASES[a1_i];
                if kind < 12 {
                    let rest: Vec<usize> = (0..4).filter(|i| *i != ref_i && *i != a1_i).collect();
                    let a2 = BASES[rest[((r >> 18) & 1) as usize]];
                    (vec![refb], vec![a1, b',', a2])
                } else {
                    (vec![refb], vec![a1])
                }
            } else if kind < 230 {
                // insertion
                let ins = BASES[((r >> 20) & 3) as usize];
                (vec![refb], vec![refb, ins, ins])
            } else {
                // deletion
                let d = BASES[((r >> 22) & 3) as usize];
                (vec![refb, d], vec![refb])
            };
            let n_alt = alt.iter().filter(|&&b| b == b',').count() + 1;

            let keep = drop_pct == 0 || drop_rng.below(100) >= drop_pct;
            if !keep {
                // Still consume genotype randomness so cohorts stay comparable.
                for _ in 0..samples {
                    gt_rng.next();
                }
                continue;
            }

            line.clear();
            let mut ac = vec![0u32; n_alt];
            let mut an = 0u32;
            let mut gts = String::with_capacity(samples * 20);
            for _ in 0..samples {
                let g = gt_rng.next();
                let miss = (g & 0x3f) == 0;
                if miss {
                    gts.push_str("\t./.:.:.:.");
                    continue;
                }
                let a = ((g >> 8) % 10 == 0) as usize * (1 + ((g >> 12) as usize % n_alt));
                let b = ((g >> 20) % 6 == 0) as usize * (1 + ((g >> 24) as usize % n_alt));
                an += 2;
                if a > 0 {
                    ac[a - 1] += 1;
                }
                if b > 0 {
                    ac[b - 1] += 1;
                }
                let dp = 10 + (g >> 32) % 40;
                let gq = 20 + (g >> 40) % 79;
                let ad_ref = dp / 2;
                gts.push('\t');
                gts.push_str(&a.to_string());
                gts.push('/');
                gts.push_str(&b.to_string());
                gts.push_str(":");
                gts.push_str(&ad_ref.to_string());
                for _ in 0..n_alt {
                    gts.push(',');
                    gts.push_str(&(dp - ad_ref).to_string());
                }
                gts.push(':');
                gts.push_str(&dp.to_string());
                gts.push(':');
                gts.push_str(&gq.to_string());
            }

            let qual = 20 + (r >> 40) % 900;
            let filter = if (r >> 50) % 20 == 0 { "LowQual" } else { "PASS" };
            let acs: Vec<String> = ac.iter().map(|v| v.to_string()).collect();
            write!(
                out,
                "chr{}\t{}\trs{}\t{}\t{}\t{}\t{}\tAC={};AN={};DP={}\tGT:AD:DP:GQ{}\n",
                c + 1,
                pos,
                (c as u64 + 1) * 100_000_000 + pos,
                String::from_utf8_lossy(&reference),
                String::from_utf8_lossy(&alt),
                qual,
                filter,
                acs.join(","),
                an,
                an * 15,
                gts
            )
            .unwrap();
        }
    }
    out.flush().unwrap();
}
