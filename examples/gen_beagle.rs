//! Deterministic generator for Beagle-5-style imputation output.
//!
//! Beagle VCFs differ from generic call sets in ways that exercise different
//! merge paths: every genotype is phased (`0|1`, never missing), FORMAT is
//! GT:DS:GP with a Number=G field, INFO carries DR2/AF and the IMP flag,
//! QUAL is ".", and sites are biallelic SNPs on a shared reference panel —
//! so cohorts imputed against the same panel have an identical site universe.
//!
//!   cargo run --release --example gen_beagle -- --sites 150000 --samples 300 \
//!       --sample-prefix A --gt-seed 11 > cohortA.vcf
//!
//! `--gt-only` drops DS/GP and the DR2/AF/IMP INFO fields, leaving FORMAT=GT
//! alone — the shape of a phasing-only pipeline (Beagle run without genotype
//! likelihoods, or output already stripped of dosage/probability fields).

use std::io::{BufWriter, Write};

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
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
    args.iter()
        .position(|a| a == name)
        .and_then(|i| args.get(i + 1).cloned())
        .unwrap_or_else(|| default.to_string())
}

const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

fn main() {
    let sites: u64 = arg("--sites", "150000").parse().unwrap();
    let samples: usize = arg("--samples", "300").parse().unwrap();
    let contigs: usize = arg("--contigs", "2").parse().unwrap();
    let site_seed: u64 = arg("--site-seed", "1").parse().unwrap();
    let gt_seed: u64 = arg("--gt-seed", "42").parse().unwrap();
    let prefix = arg("--sample-prefix", "S");
    // Nonzero: omit a fraction of sites, for cohorts not sharing every marker.
    let drop_pct: u64 = arg("--drop-pct", "0").parse().unwrap();
    let drop_seed: u64 = arg("--drop-seed", "7").parse().unwrap();
    let gt_only = std::env::args().any(|a| a == "--gt-only");

    let mut out = BufWriter::with_capacity(1 << 20, std::io::stdout().lock());
    writeln!(out, "##fileformat=VCFv4.2").unwrap();
    writeln!(out, "##filedate=20260827").unwrap();
    writeln!(out, "##source=\"beagle.27May26.abc.jar\"").unwrap();
    for c in 1..=contigs {
        writeln!(out, "##contig=<ID=chr{c}>").unwrap();
    }
    if !gt_only {
        writeln!(out, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Estimated ALT Allele Frequencies\">").unwrap();
        writeln!(out, "##INFO=<ID=DR2,Number=A,Type=Float,Description=\"Dosage R-Squared: estimated squared correlation between estimated REF dose [P(RA) + 2*P(RR)] and true REF dose\">").unwrap();
        writeln!(out, "##INFO=<ID=IMP,Number=0,Type=Flag,Description=\"Imputed marker\">").unwrap();
    }
    writeln!(out, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">").unwrap();
    if !gt_only {
        writeln!(out, "##FORMAT=<ID=DS,Number=A,Type=Float,Description=\"estimated ALT dose [P(RA) + 2*P(AA)]\">").unwrap();
        writeln!(out, "##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Estimated Genotype Probability\">").unwrap();
    }
    write!(out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT").unwrap();
    for i in 0..samples {
        write!(out, "\t{prefix}{i:05}").unwrap();
    }
    writeln!(out).unwrap();

    let mut gt_rng = Rng(gt_seed);
    let mut drop_rng = Rng(drop_seed);
    let per_contig = sites / contigs as u64;

    for c in 0..contigs {
        let mut pos = 0u64;
        let mut srng = Rng(site_seed.wrapping_mul(1_000_003).wrapping_add(c as u64));
        for _ in 0..per_contig {
            pos += 1 + srng.below(1500);
            let r = srng.next();
            let ref_i = (r & 3) as usize;
            let alt_i = (ref_i + 1 + ((r >> 8) % 3) as usize) % 4;
            // Site allele frequency shapes the cohort genotypes.
            let af_pm = 1 + srng.below(500); // per-mille, 0.001..0.5
            let imputed = srng.below(100) < 85;
            let dr2 = if imputed { 30 + srng.below(70) } else { 95 + srng.below(5) };

            // Every srng draw happens before the keep check, or cohorts with
            // different drop seeds would desync and describe different sites.
            let has_id = srng.below(3) == 0;
            let keep = drop_pct == 0 || drop_rng.below(100) >= drop_pct;
            if !keep {
                for _ in 0..samples {
                    gt_rng.next();
                }
                continue;
            }

            let id = if has_id { format!("rs{}", (c as u64 + 1) * 10_000_000 + pos) } else { ".".to_string() };
            if gt_only {
                write!(
                    out,
                    "chr{}\t{}\t{}\t{}\t{}\t.\tPASS\t.\tGT",
                    c + 1, pos, id, BASES[ref_i] as char, BASES[alt_i] as char,
                )
                .unwrap();
            } else {
                write!(
                    out,
                    "chr{}\t{}\t{}\t{}\t{}\t.\tPASS\tDR2={}.{:02};AF={}.{:04}{}\tGT:DS:GP",
                    c + 1,
                    pos,
                    id,
                    BASES[ref_i] as char,
                    BASES[alt_i] as char,
                    dr2 / 100,
                    dr2 % 100,
                    af_pm / 1000,
                    af_pm % 1000 * 10,
                    if imputed { ";IMP" } else { "" },
                )
                .unwrap();
            }
            let _ = dr2;

            for _ in 0..samples {
                let g = gt_rng.next();
                // Two haplotypes drawn at the site AF; everything phased.
                let a = (g % 1000 < af_pm) as u32;
                let b = ((g >> 20) % 1000 < af_pm) as u32;
                if gt_only {
                    write!(out, "\t{a}|{b}").unwrap();
                    continue;
                }
                let dose = a + b;
                // GP concentrated on the called genotype, Beagle-style rounding.
                let noise = (g >> 40) % 6; // 0.00..0.05
                let p_call = 100 - noise;
                let (p0, p1, p2) = match dose {
                    0 => (p_call, noise, 0),
                    1 => (noise / 2, p_call, noise - noise / 2),
                    _ => (0, noise, p_call),
                };
                // DS = P(het) + 2*P(hom-alt), i.e. the integer dose nudged by
                // the probability mass the call did not get.
                let ds_centi = match dose {
                    0 => noise,
                    1 => 100 - (noise as i64 - 2 * (noise as i64 / 2)).unsigned_abs() % 100 + noise / 2 * 2 - noise / 2 * 2,
                    _ => 200 - noise,
                };
                let ds_centi = if dose == 1 { 100 + noise / 2 - (noise - noise / 2) } else { ds_centi };
                write!(
                    out,
                    "\t{a}|{b}:{}.{:02}:{}.{:02},{}.{:02},{}.{:02}",
                    ds_centi / 100,
                    ds_centi % 100,
                    p0 / 100, p0 % 100,
                    p1 / 100, p1 % 100,
                    p2 / 100, p2 % 100,
                )
                .unwrap();
            }
            writeln!(out).unwrap();
        }
    }
    out.flush().unwrap();
}
