//! `vcfr view` — subset samples and variants.

use std::io::Write;

use crate::io::open_reader;
use crate::vcf::header::{select_samples, Header};
use crate::vcf::record::*;
use crate::vcf::RegionSet;

use super::{ignore_broken_pipe, resolve_threads, split_threads, OutputOpts};

#[derive(clap::Args, Debug)]
pub struct ViewArgs {
    /// Input VCF (plain, gzip or BGZF); "-" reads standard input
    #[arg(value_name = "FILE")]
    pub input: String,

    /// Comma-separated samples to keep; prefix with ^ to exclude instead
    #[arg(short = 's', long = "samples", value_name = "LIST")]
    pub samples: Option<String>,

    /// File of sample names, one per line; prefix the path with ^ to exclude
    #[arg(short = 'S', long = "samples-file", value_name = "FILE")]
    pub samples_file: Option<String>,

    /// Regions to keep, e.g. chr1,chr2:1000-2000. vcfr always streams, so this
    /// behaves like `bcftools view -t` and needs no index.
    #[arg(short = 'r', long = "regions", visible_alias = "targets", short_alias = 't', value_name = "LIST")]
    pub regions: Option<String>,

    /// Regions file: CHROM[<TAB>BEG[<TAB>END]], 1-based inclusive
    #[arg(short = 'R', long = "regions-file", visible_alias = "targets-file", short_alias = 'T', value_name = "FILE")]
    pub regions_file: Option<String>,

    /// Keep only these variant types: snps,indels,mnps,other,ref
    #[arg(short = 'v', long = "types", value_name = "LIST")]
    pub types: Option<String>,

    /// Drop these variant types
    #[arg(short = 'V', long = "exclude-types", value_name = "LIST")]
    pub exclude_types: Option<String>,

    /// Keep only records whose FILTER is one of these (use PASS or .)
    #[arg(short = 'f', long = "apply-filters", value_name = "LIST")]
    pub apply_filters: Option<String>,

    /// Keep only records with these IDs (comma-separated)
    #[arg(long = "ids", value_name = "LIST")]
    pub ids: Option<String>,

    /// File of IDs to keep, one per line
    #[arg(long = "ids-file", value_name = "FILE")]
    pub ids_file: Option<String>,

    /// Minimum number of alleles (REF included)
    #[arg(short = 'm', long = "min-alleles", value_name = "N")]
    pub min_alleles: Option<usize>,

    /// Maximum number of alleles (REF included)
    #[arg(short = 'M', long = "max-alleles", value_name = "N")]
    pub max_alleles: Option<usize>,

    /// Minimum QUAL (long form only: bcftools spells -q as minimum allele frequency)
    #[arg(long = "min-qual", value_name = "F")]
    pub min_qual: Option<f64>,

    /// Maximum QUAL
    #[arg(long = "max-qual", value_name = "F")]
    pub max_qual: Option<f64>,

    /// Minimum non-reference allele count over the kept samples
    #[arg(long = "min-ac", value_name = "N")]
    pub min_ac: Option<u64>,

    /// Drop the FORMAT and sample columns
    #[arg(short = 'G', long = "drop-genotypes")]
    pub drop_genotypes: bool,

    /// Print only the header
    #[arg(short = 'h', long = "header-only")]
    pub header_only: bool,

    /// Suppress the header
    #[arg(short = 'H', long = "no-header")]
    pub no_header: bool,

    /// Do not recompute INFO/AC and INFO/AN after subsetting samples
    #[arg(short = 'I', long = "no-update")]
    pub no_update: bool,

    #[arg(long = "threads", default_value_t = 0, value_name = "N")]
    pub threads: usize,

    #[command(flatten)]
    pub out: OutputOpts,
}

pub fn run(a: &ViewArgs) -> Result<(), String> {
    let threads = resolve_threads(a.threads);
    let (rthreads, wthreads) = split_threads(threads, a.out.bgzf()?);

    let mut rdr = open_reader(&a.input, rthreads).map_err(|e| e.to_string())?;
    let hdr = Header::read(&mut rdr).map_err(|e| e.to_string())?;

    let sample_spec = match (&a.samples, &a.samples_file) {
        (Some(_), Some(_)) => return Err("use either --samples or --samples-file".into()),
        (Some(s), None) => Some(s.clone()),
        (None, Some(f)) => Some(read_sample_file(f)?),
        (None, None) => None,
    };
    let selected: Option<Vec<usize>> = match &sample_spec {
        Some(s) => Some(select_samples(&hdr.samples, s)?),
        None => None,
    };
    let out_samples: Vec<String> = match &selected {
        Some(sel) => sel.iter().map(|&i| hdr.samples[i].clone()).collect(),
        None => hdr.samples.clone(),
    };

    let mut regions = match (&a.regions, &a.regions_file) {
        (Some(_), Some(_)) => return Err("use either --regions or --regions-file".into()),
        (Some(s), None) => Some(RegionSet::parse(s)?),
        (None, Some(f)) => Some(RegionSet::from_file(f)?),
        (None, None) => None,
    };
    if let Some(rs) = regions.as_mut() {
        rs.set_contig_order(|c| hdr.contig_rank(c.as_bytes()));
    }

    let want_types = a.types.as_deref().map(parse_types).transpose()?;
    let drop_types = a.exclude_types.as_deref().map(parse_types).transpose()?;
    let filters: Option<Vec<Vec<u8>>> = a
        .apply_filters
        .as_deref()
        .map(|s| s.split(',').map(|f| f.as_bytes().to_vec()).collect());
    let ids: Option<std::collections::HashSet<Vec<u8>>> = match (&a.ids, &a.ids_file) {
        (Some(s), None) => Some(s.split(',').map(|i| i.as_bytes().to_vec()).collect()),
        (None, Some(f)) => Some(
            std::fs::read_to_string(f)
                .map_err(|e| format!("{f}: {e}"))?
                .lines()
                .map(|l| l.trim().as_bytes().to_vec())
                .filter(|l| !l.is_empty())
                .collect(),
        ),
        (Some(_), Some(_)) => return Err("use either --ids or --ids-file".into()),
        (None, None) => None,
    };

    let mut w = a.out.writer(wthreads)?;
    let mut buf: Vec<u8> = Vec::with_capacity(1 << 16);

    if !a.no_header {
        hdr.write(&mut buf, &out_samples, a.drop_genotypes);
        w.write_all(&buf).map_err(|e| e.to_string())?;
        buf.clear();
    }
    if a.header_only {
        return finish(w);
    }

    // AC/AN only need recomputing when the sample set actually changed.
    let update_info = !a.no_update
        && !a.drop_genotypes
        && selected.as_ref().map_or(false, |s| s.len() != hdr.samples.len());
    let needs_samples = selected.is_some() || a.drop_genotypes || update_info || a.min_ac.is_some();

    let mut rec = Record::new();
    let mut gts: Vec<i32> = Vec::with_capacity(4);
    let mut info_buf: Vec<u8> = Vec::with_capacity(256);
    let all_samples: Vec<usize> = (0..hdr.samples.len()).collect();

    let res = (|| -> std::io::Result<()> {
        while rdr.advance()? {
            let line = rdr.line();
            if line.is_empty() || line[0] == b'#' {
                continue;
            }
            if needs_samples {
                rec.split(line);
            } else {
                rec.split_prefix(line, COL_INFO + 1);
            }
            if rec.ncols() < 8 {
                continue;
            }

            let chrom = rec.get(line, COL_CHROM);
            let pos = parse_u64(rec.get(line, COL_POS)).unwrap_or(0);

            if let Some(rs) = &regions {
                let reference = rec.get(line, COL_REF);
                let end = pos + reference.len().max(1) as u64 - 1;
                if rs.past_all(hdr.contig_rank(chrom), pos) {
                    break;
                }
                if !rs.contains(chrom, pos, end) {
                    continue;
                }
            }
            if let Some(set) = &ids {
                let id = rec.get(line, COL_ID);
                if !id.split(|&b| b == b';').any(|i| set.contains(i)) {
                    continue;
                }
            }
            if want_types.is_some() || drop_types.is_some() {
                let t = site_types(rec.get(line, COL_REF), rec.get(line, COL_ALT));
                if let Some(w) = want_types {
                    if t & w == 0 {
                        continue;
                    }
                }
                if let Some(d) = drop_types {
                    if t & d != 0 {
                        continue;
                    }
                }
            }
            if a.min_alleles.is_some() || a.max_alleles.is_some() {
                let n = count_alts(rec.get(line, COL_ALT)) + 1;
                if a.min_alleles.map_or(false, |m| n < m) || a.max_alleles.map_or(false, |m| n > m) {
                    continue;
                }
            }
            if a.min_qual.is_some() || a.max_qual.is_some() {
                let q = rec.get(line, COL_QUAL);
                let qv = std::str::from_utf8(q).ok().and_then(|s| s.parse::<f64>().ok());
                match qv {
                    Some(v) => {
                        if a.min_qual.map_or(false, |m| v < m) || a.max_qual.map_or(false, |m| v > m)
                        {
                            continue;
                        }
                    }
                    None => {
                        if a.min_qual.is_some() {
                            continue;
                        }
                    }
                }
            }
            if let Some(fs) = &filters {
                let f = rec.get(line, COL_FILTER);
                if !f.split(|&b| b == b';').any(|v| fs.iter().any(|x| x == v)) {
                    continue;
                }
            }

            // Fast path: nothing about the record changes, so pass the bytes on.
            if !needs_samples {
                w.write_all(line)?;
                w.write_all(b"\n")?;
                continue;
            }

            let sel: &[usize] = selected.as_deref().unwrap_or(&all_samples);
            let counts = if update_info || a.min_ac.is_some() {
                let n_alt = count_alts(rec.get(line, COL_ALT));
                let gi = if rec.ncols() > COL_FORMAT {
                    gt_index(rec.get(line, COL_FORMAT)).unwrap_or(usize::MAX)
                } else {
                    usize::MAX
                };
                if gi == usize::MAX {
                    None
                } else {
                    Some(count_alleles(line, &rec, sel, n_alt, gi, &mut gts))
                }
            } else {
                None
            };
            if let (Some(min_ac), Some(c)) = (a.min_ac, counts.as_ref()) {
                if c.ac.iter().copied().max().unwrap_or(0) < min_ac {
                    continue;
                }
            }

            buf.clear();
            if a.drop_genotypes {
                buf.extend_from_slice(rec.span(line, COL_CHROM, COL_INFO));
            } else if let Some(c) = counts.filter(|_| update_info) {
                buf.extend_from_slice(rec.span(line, COL_CHROM, COL_FILTER));
                buf.push(b'\t');
                info_buf.clear();
                rewrite_info(
                    rec.get(line, COL_INFO),
                    &[
                        ("AC", Some(join_u64(&c.ac))),
                        ("AN", Some(c.an.to_string())),
                    ],
                    &mut info_buf,
                );
                buf.extend_from_slice(&info_buf);
                append_samples(&mut buf, line, &rec, sel);
            } else {
                buf.extend_from_slice(rec.span(line, COL_CHROM, COL_INFO));
                append_samples(&mut buf, line, &rec, sel);
            }
            buf.push(b'\n');
            w.write_all(&buf)?;
        }
        Ok(())
    })();

    ignore_broken_pipe(res)?;
    finish(w)
}

fn append_samples(buf: &mut Vec<u8>, line: &[u8], rec: &Record, sel: &[usize]) {
    if rec.ncols() <= COL_FORMAT {
        return;
    }
    buf.push(b'\t');
    buf.extend_from_slice(rec.get(line, COL_FORMAT));
    for &s in sel {
        let c = FIRST_SAMPLE + s;
        buf.push(b'\t');
        if c < rec.ncols() {
            buf.extend_from_slice(rec.get(line, c));
        } else {
            buf.push(b'.');
        }
    }
}

fn finish(mut w: crate::io::Writer) -> Result<(), String> {
    ignore_broken_pipe(w.finish())
}

fn read_sample_file(path: &str) -> Result<String, String> {
    let (invert, p) = match path.strip_prefix('^') {
        Some(rest) => (true, rest),
        None => (false, path),
    };
    let text = std::fs::read_to_string(p).map_err(|e| format!("{p}: {e}"))?;
    let names: Vec<&str> = text
        .lines()
        .map(|l| l.split_whitespace().next().unwrap_or(""))
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .collect();
    Ok(format!("{}{}", if invert { "^" } else { "" }, names.join(",")))
}

