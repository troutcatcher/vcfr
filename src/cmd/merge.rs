//! `vcfr merge` — combine files that hold different samples at the same sites.

use std::borrow::Cow;
use std::collections::HashMap;
use std::io::Write;

use crate::io::{open_reader, Reader};
use crate::vcf::header::{union_meta, Header, Number};
use crate::vcf::record::*;

use super::{ignore_broken_pipe, merge_threads, resolve_threads, OutputOpts};

#[derive(clap::Args, Debug)]
pub struct MergeArgs {
    /// VCFs to merge; each must be coordinate-sorted
    #[arg(value_name = "FILE")]
    pub inputs: Vec<String>,

    /// Read the file list from FILE, one path per line
    #[arg(short = 'F', long = "file-list", value_name = "FILE")]
    pub file_list: Option<String>,

    /// Which records at one position may be combined into a multiallelic
    /// record: snps, indels, both, all or none
    #[arg(short = 'm', long = "merge", default_value = "both", value_name = "STRING")]
    pub merge: String,

    /// Rename duplicated sample names instead of failing
    #[arg(long = "force-samples")]
    pub force_samples: bool,

    /// Fill genotypes missing from an input with 0/0 rather than ./.
    #[arg(short = '0', long = "missing-to-ref")]
    pub missing_to_ref: bool,

    /// Do not recompute INFO/AC and INFO/AN
    #[arg(short = 'I', long = "no-update")]
    pub no_update: bool,

    /// How to combine INFO fields present in several files:
    /// KEY:sum|avg|min|max|join|first,... Use "-" to always take the first value.
    #[arg(short = 'i', long = "info-rules", default_value = "DP:sum,DP4:sum", value_name = "RULES")]
    pub info_rules: String,

    #[arg(long = "threads", default_value_t = 0, value_name = "N")]
    pub threads: usize,

    #[command(flatten)]
    pub out: OutputOpts,
}

/// `-m`: which records sharing a position may become one multiallelic record.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum MergeMode {
    /// Only records with identical alleles are combined.
    None,
    /// SNP records combine with each other; everything else needs identical alleles.
    Snps,
    /// Indel records combine with each other.
    Indels,
    /// SNPs combine with SNPs and indels with indels, but never with each other.
    Both,
    /// Everything at a position with the same REF becomes one record.
    All,
}

impl MergeMode {
    fn parse(s: &str) -> Result<MergeMode, String> {
        Ok(match s {
            "none" => MergeMode::None,
            "snps" => MergeMode::Snps,
            "indels" => MergeMode::Indels,
            "both" => MergeMode::Both,
            "all" => MergeMode::All,
            _ => return Err(format!("unknown -m mode '{s}' (snps, indels, both, all or none)")),
        })
    }

    /// Output rank of a group when a position yields several records. SNPs
    /// come before indels, except under `-m indels`, where the class being
    /// merged leads. Matches bcftools' ordering.
    fn rank(&self, c: Class) -> u8 {
        match (self, c) {
            (MergeMode::Indels, Class::Indel) => 0,
            (MergeMode::Indels, Class::Snp) => 1,
            (_, Class::Snp) => 0,
            (_, Class::Indel) => 1,
            (_, Class::Other) => 2,
        }
    }

    /// May a record of class `b` join a group whose class is `a`?
    fn forces(&self, a: Class, b: Class) -> bool {
        match self {
            MergeMode::All => true,
            MergeMode::None => false,
            MergeMode::Snps => a == Class::Snp && b == Class::Snp,
            MergeMode::Indels => a == Class::Indel && b == Class::Indel,
            MergeMode::Both => a == b && a != Class::Other,
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Class {
    Snp,
    Indel,
    Other,
}

fn classify(reference: &[u8], alt: &[u8]) -> Class {
    let t = site_types(reference, alt);
    if t & T_INDEL != 0 {
        Class::Indel
    } else if t & T_SNP != 0 {
        Class::Snp
    } else {
        Class::Other
    }
}

fn alt_list(alt: &[u8]) -> Vec<&[u8]> {
    if alt == b"." || alt.is_empty() {
        Vec::new()
    } else {
        alt.split(|&b| b == b',').collect()
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Rule {
    Sum,
    Avg,
    Min,
    Max,
    Join,
    First,
}

/// `--info-rules`, matching bcftools' defaults so shared INFO fields combine
/// the same way.
#[derive(Default)]
struct InfoRules(HashMap<String, Rule>);

impl InfoRules {
    fn parse(spec: &str) -> Result<InfoRules, String> {
        let mut m = HashMap::new();
        if spec != "-" {
            for item in spec.split(',').filter(|s| !s.is_empty()) {
                let (k, v) = item
                    .split_once(':')
                    .ok_or_else(|| format!("bad --info-rules entry '{item}', expected KEY:METHOD"))?;
                let rule = match v {
                    "sum" => Rule::Sum,
                    "avg" => Rule::Avg,
                    "min" => Rule::Min,
                    "max" => Rule::Max,
                    "join" => Rule::Join,
                    "first" => Rule::First,
                    _ => return Err(format!("unknown --info-rules method '{v}'")),
                };
                m.insert(k.to_string(), rule);
            }
        }
        Ok(InfoRules(m))
    }
    fn get(&self, key: &[u8]) -> Option<Rule> {
        std::str::from_utf8(key).ok().and_then(|k| self.0.get(k)).copied()
    }
}

/// Contig ordering shared by every input: header contigs first, then any
/// contig seen only in the records, in order of first appearance.
struct RankMap {
    ranks: HashMap<Vec<u8>, usize>,
    next: usize,
}

impl RankMap {
    fn new(contigs: &[String]) -> RankMap {
        let mut ranks = HashMap::new();
        for (i, c) in contigs.iter().enumerate() {
            ranks.insert(c.as_bytes().to_vec(), i);
        }
        RankMap { next: contigs.len(), ranks }
    }
    fn rank(&mut self, chrom: &[u8]) -> usize {
        if let Some(&r) = self.ranks.get(chrom) {
            return r;
        }
        let r = self.next;
        self.next += 1;
        self.ranks.insert(chrom.to_vec(), r);
        r
    }
}

struct Stream {
    path: String,
    rdr: Reader,
    line: Vec<u8>,
    rec: Record,
    valid: bool,
    rank: usize,
    pos: u64,
    n_samples: usize,
    /// Index of this file's first sample in the output sample list.
    sample_off: usize,
}

impl Stream {
    #[inline]
    fn get(&self, col: usize) -> &[u8] {
        self.rec.get(&self.line, col)
    }

    fn advance(&mut self, ranks: &mut RankMap) -> Result<(), String> {
        loop {
            let more = self.rdr.advance().map_err(|e| format!("{}: {e}", self.path))?;
            if !more {
                self.valid = false;
                return Ok(());
            }
            let l = self.rdr.line();
            if l.is_empty() || l[0] == b'#' {
                continue;
            }
            self.line.clear();
            self.line.extend_from_slice(l);
            self.rec.split(&self.line);
            if self.rec.ncols() < 8 {
                continue;
            }
            let pos = parse_u64(self.rec.get(&self.line, COL_POS)).ok_or_else(|| {
                format!("{}: bad POS '{}'", self.path, String::from_utf8_lossy(self.rec.get(&self.line, COL_POS)))
            })?;
            let rank = ranks.rank(self.rec.get(&self.line, COL_CHROM));
            if self.valid && (rank, pos) < (self.rank, self.pos) {
                return Err(format!(
                    "{}: input is not coordinate-sorted at {}:{pos}",
                    self.path,
                    String::from_utf8_lossy(self.rec.get(&self.line, COL_CHROM))
                ));
            }
            self.rank = rank;
            self.pos = pos;
            self.valid = true;
            return Ok(());
        }
    }
}

/// Reusable buffers so the per-record path does not allocate.
#[derive(Default)]
struct Ctx {
    alts: Vec<Vec<u8>>,
    /// Per member: old allele index -> new allele index.
    maps: Vec<Vec<i32>>,
    identity: Vec<bool>,
    fmt_keys: Vec<Vec<u8>>,
    /// Per member: union key index -> that member's subfield index.
    key_pos: Vec<Vec<Option<usize>>>,
    gts: Vec<i32>,
    info_buf: Vec<u8>,
    ac: Vec<u64>,
}

pub fn run(a: &MergeArgs) -> Result<(), String> {
    let inputs = collect_inputs(a)?;
    if inputs.len() < 2 {
        return Err("merge needs at least two input files".into());
    }
    let threads = resolve_threads(a.threads);
    // Inflation must not run on the main thread: assembling merged records is
    // the serial bottleneck, so every input gets its own inflation worker.
    let (per_file, wthreads) = merge_threads(threads, inputs.len(), a.out.bgzf()?);

    let mut headers = Vec::with_capacity(inputs.len());
    let mut readers = Vec::with_capacity(inputs.len());
    for p in &inputs {
        let mut r = open_reader(p, per_file).map_err(|e| e.to_string())?;
        headers.push(Header::read(&mut r).map_err(|e| format!("{p}: {e}"))?);
        readers.push(r);
    }

    // Output sample list.
    let mut out_samples: Vec<String> = Vec::new();
    let mut seen: HashMap<String, usize> = HashMap::new();
    let mut offsets = Vec::with_capacity(inputs.len());
    for (fi, h) in headers.iter().enumerate() {
        offsets.push(out_samples.len());
        for s in &h.samples {
            let name = match seen.get(s) {
                None => s.clone(),
                Some(_) if a.force_samples => format!("{}:{s}", fi + 1),
                Some(_) => {
                    return Err(format!(
                        "duplicate sample '{s}' in {}; pass --force-samples to rename it",
                        inputs[fi]
                    ))
                }
            };
            seen.insert(name.clone(), fi);
            out_samples.push(name);
        }
    }

    let mut merged_hdr = union_meta(&headers);
    merged_hdr.ensure_meta(
        "FORMAT",
        "GT",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    );
    if !a.no_update {
        merged_hdr.ensure_meta("INFO", "AC", "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">");
        merged_hdr.ensure_meta("INFO", "AN", "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">");
    }
    merged_hdr.has_format_col = !out_samples.is_empty();

    let mut ranks = RankMap::new(&merged_hdr.contigs);
    let mut streams: Vec<Stream> = Vec::with_capacity(inputs.len());
    for ((p, r), (h, off)) in inputs
        .iter()
        .zip(readers)
        .zip(headers.iter().zip(offsets.iter()))
    {
        streams.push(Stream {
            path: p.clone(),
            rdr: r,
            line: Vec::with_capacity(1 << 12),
            rec: Record::new(),
            valid: false,
            rank: 0,
            pos: 0,
            n_samples: h.samples.len(),
            sample_off: *off,
        });
        let s = streams.last_mut().unwrap();
        s.advance(&mut ranks)?;
    }

    let mut w = a.out.writer(wthreads)?;
    let mut buf = Vec::with_capacity(1 << 16);
    merged_hdr.write(&mut buf, &out_samples, false);
    w.write_all(&buf).map_err(|e| e.to_string())?;

    let rules = InfoRules::parse(&a.info_rules)?;
    let mode = MergeMode::parse(&a.merge)?;
    let mut ctx = Ctx::default();
    let n_out = out_samples.len();
    let mut members: Vec<usize> = Vec::with_capacity(streams.len());

    let res = (|| -> Result<(), String> {
        loop {
            let mut best: Option<(usize, u64)> = None;
            for s in &streams {
                if s.valid {
                    let k = (s.rank, s.pos);
                    best = Some(match best {
                        Some(b) if b <= k => b,
                        _ => k,
                    });
                }
            }
            let key = match best {
                Some(k) => k,
                None => break,
            };
            members.clear();
            for (i, s) in streams.iter().enumerate() {
                if s.valid && (s.rank, s.pos) == key {
                    members.push(i);
                }
            }
            // Records that disagree on REF describe different events. Beyond
            // that, -m decides which ALTs may share one multiallelic record;
            // records whose alleles are already a subset of the group's join
            // regardless, since that creates no new multiallelic.
            let mut done = vec![false; members.len()];
            let mut groups: Vec<(Class, Vec<usize>)> = Vec::new();
            for gi in 0..members.len() {
                if done[gi] {
                    continue;
                }
                done[gi] = true;
                let head = members[gi];
                let reference = streams[head].get(COL_REF).to_vec();
                let mut alts: Vec<Vec<u8>> = alt_list(streams[head].get(COL_ALT))
                    .into_iter()
                    .map(|a| a.to_vec())
                    .collect();
                let mut class = classify(&reference, streams[head].get(COL_ALT));
                let mut g = vec![head];
                for gj in gi + 1..members.len() {
                    if done[gj] {
                        continue;
                    }
                    let s = &streams[members[gj]];
                    if s.get(COL_REF) != &reference[..] {
                        continue;
                    }
                    let their = alt_list(s.get(COL_ALT));
                    let their_class = classify(&reference, s.get(COL_ALT));
                    let novel = their.iter().filter(|a| !alts.iter().any(|x| x == *a)).count();
                    // No new allele, or the mode says this class may combine.
                    let ok = novel == 0
                        || their.len() >= alts.len() + novel
                        || mode.forces(class, their_class);
                    if !ok {
                        continue;
                    }
                    done[gj] = true;
                    for a in their {
                        if !alts.iter().any(|x| x == a) {
                            alts.push(a.to_vec());
                        }
                    }
                    if class == Class::Other {
                        class = their_class;
                    }
                    g.push(members[gj]);
                }
                groups.push((class, g));
            }
            groups.sort_by_key(|(c, _)| mode.rank(*c));
            for (_, g) in &groups {
                buf.clear();
                emit(&streams, g, n_out, a, &merged_hdr, &rules, &mut ctx, &mut buf)?;
                w.write_all(&buf).map_err(|e| e.to_string())?;
            }
            for &i in &members {
                streams[i].advance(&mut ranks)?;
            }
        }
        Ok(())
    })();
    res?;
    ignore_broken_pipe(w.finish())
}

#[allow(clippy::too_many_arguments)]
fn emit(
    streams: &[Stream],
    members: &[usize],
    n_out: usize,
    a: &MergeArgs,
    hdr: &Header,
    rules: &InfoRules,
    ctx: &mut Ctx,
    out: &mut Vec<u8>,
) -> Result<(), String> {
    let first = &streams[members[0]];

    // Almost always every input describes the site identically. When it does,
    // there is nothing to renumber and no FORMAT union to take, so the whole
    // remapping apparatus below can be skipped and each input's sample columns
    // copied across in one memcpy.
    let first_alt = first.get(COL_ALT);
    let first_fmt: &[u8] = if first.rec.ncols() > COL_FORMAT { first.get(COL_FORMAT) } else { b"" };
    let uniform = members.iter().skip(1).all(|&m| {
        let s = &streams[m];
        s.get(COL_ALT) == first_alt
            && (if s.rec.ncols() > COL_FORMAT { s.get(COL_FORMAT) } else { b"" }) == first_fmt
    });

    let n_alt;
    let gt_slot;
    let n_keys;
    ctx.identity.clear();
    if uniform {
        n_alt = count_alts(first_alt);
        gt_slot = gt_index(first_fmt);
        n_keys = if first_fmt.is_empty() { 0 } else { first_fmt.iter().filter(|&&b| b == b':').count() + 1 };
        ctx.identity.resize(members.len(), true);
    } else {
        // --- ALT union and per-member allele renumbering ---------------------
        ctx.alts.clear();
        ctx.maps.clear();
        for &m in members {
            let s = &streams[m];
            let alt = s.get(COL_ALT);
            let mut map = vec![0i32; 1];
            let mut ident = true;
            if alt != b"." && !alt.is_empty() {
                for (old, al) in alt.split(|&b| b == b',').enumerate() {
                    let idx = match ctx.alts.iter().position(|x| x == al) {
                        Some(i) => i,
                        None => {
                            ctx.alts.push(al.to_vec());
                            ctx.alts.len() - 1
                        }
                    };
                    if idx != old {
                        ident = false;
                    }
                    map.push(idx as i32 + 1);
                }
            }
            ctx.maps.push(map);
            ctx.identity.push(ident);
        }
        n_alt = ctx.alts.len();
        // A member only passes through unchanged if it also carries every allele.
        for mi in 0..members.len() {
            ctx.identity[mi] = ctx.identity[mi] && ctx.maps[mi].len() == n_alt + 1;
        }

        // --- FORMAT union ----------------------------------------------------
        ctx.fmt_keys.clear();
        let mut has_gt = false;
        for &m in members {
            let s = &streams[m];
            if s.rec.ncols() <= COL_FORMAT {
                continue;
            }
            for k in s.get(COL_FORMAT).split(|&b| b == b':') {
                if k == b"GT" {
                    has_gt = true;
                    continue;
                }
                if !ctx.fmt_keys.iter().any(|x| x == k) {
                    ctx.fmt_keys.push(k.to_vec());
                }
            }
        }
        if has_gt || ctx.fmt_keys.is_empty() {
            ctx.fmt_keys.insert(0, b"GT".to_vec());
        }
        ctx.key_pos.clear();
        for &m in members {
            let s = &streams[m];
            let own: Vec<&[u8]> = if s.rec.ncols() > COL_FORMAT {
                s.get(COL_FORMAT).split(|&b| b == b':').collect()
            } else {
                Vec::new()
            };
            ctx.key_pos.push(
                ctx.fmt_keys
                    .iter()
                    .map(|k| own.iter().position(|o| *o == &k[..]))
                    .collect(),
            );
        }
        gt_slot = ctx.fmt_keys.iter().position(|k| k == b"GT");
        n_keys = ctx.fmt_keys.len();
    }

    // --- sample columns --------------------------------------------------
    ctx.ac.clear();
    ctx.ac.resize(n_alt, 0);
    let mut an = 0u64;
    let update = !a.no_update;
    let mut ploidy = 0usize;

    // Determine the ploidy used to pad samples that are absent from a file.
    if let Some(slot) = gt_slot {
        'outer: for (mi, &m) in members.iter().enumerate() {
            let s = &streams[m];
            if s.rec.ncols() <= FIRST_SAMPLE {
                continue;
            }
            if let Some(si) = if uniform { Some(slot) } else { ctx.key_pos[mi][slot] } {
                for j in 0..s.n_samples {
                    let col = FIRST_SAMPLE + j;
                    if col >= s.rec.ncols() {
                        break;
                    }
                    let g = subfield(s.get(col), si);
                    let p = g.iter().filter(|&&b| b == b'/' || b == b'|').count() + 1;
                    if p > 0 && g != b"." {
                        ploidy = p;
                        break 'outer;
                    }
                }
            }
        }
    }
    if ploidy == 0 {
        ploidy = 2;
    }

    // Allele counts have to be known before INFO is written, and INFO precedes
    // the genotypes on the line. Counting them in their own pass lets the
    // genotype columns go straight into the output buffer below, instead of
    // being staged in a second buffer and copied again.
    if update {
        for (mi, &m) in members.iter().enumerate() {
            let s = &streams[m];
            let si = match gt_slot.and_then(|slot| if uniform { Some(slot) } else { ctx.key_pos[mi][slot] }) {
                Some(si) => si,
                None => continue,
            };
            for j in 0..s.n_samples {
                let col = FIRST_SAMPLE + j;
                if col >= s.rec.ncols() {
                    break;
                }
                parse_gt(subfield(s.get(col), si), &mut ctx.gts);
                for &al in ctx.gts.iter() {
                    if al < 0 {
                        continue;
                    }
                    an += 1;
                    let mapped = if ctx.identity[mi] {
                        al
                    } else {
                        ctx.maps[mi].get(al as usize).copied().unwrap_or(-1)
                    };
                    if mapped > 0 && (mapped as usize) <= n_alt {
                        ctx.ac[mapped as usize - 1] += 1;
                    }
                }
            }
        }
        if a.missing_to_ref {
            // Samples no input covers are called 0/0, so they contribute alleles.
            let covered: usize = members.iter().map(|&m| streams[m].n_samples).sum();
            an += (n_out.saturating_sub(covered) * ploidy) as u64;
        }
    }

    // --- fixed columns ----------------------------------------------------
    out.extend_from_slice(first.get(COL_CHROM));
    out.push(b'\t');
    out.extend_from_slice(first.get(COL_POS));
    out.push(b'\t');
    write_ids(out, streams, members);
    out.push(b'\t');
    out.extend_from_slice(first.get(COL_REF));
    out.push(b'\t');
    if uniform {
        out.extend_from_slice(first_alt);
    } else if n_alt == 0 {
        out.push(b'.');
    } else {
        for (i, al) in ctx.alts.iter().enumerate() {
            if i > 0 {
                out.push(b',');
            }
            out.extend_from_slice(al);
        }
    }
    out.push(b'\t');
    write_qual(out, streams, members);
    out.push(b'\t');
    write_filter(out, streams, members);
    out.push(b'\t');

    // --- INFO -------------------------------------------------------------
    ctx.info_buf.clear();
    // When AC/AN are recomputed they are appended last, as bcftools does.
    let skip: &[&[u8]] = if update { &[b"AC", b"AN"] } else { &[] };
    merge_info(
        &mut ctx.info_buf,
        streams,
        members,
        &ctx.maps,
        &ctx.identity,
        n_alt,
        hdr,
        rules,
        skip,
    );
    if update {
        if ctx.info_buf == b"." {
            ctx.info_buf.clear();
        } else {
            ctx.info_buf.push(b';');
        }
        ctx.info_buf.extend_from_slice(b"AN=");
        ctx.info_buf.extend_from_slice(an.to_string().as_bytes());
        ctx.info_buf.extend_from_slice(b";AC=");
        ctx.info_buf.extend_from_slice(join_u64(&ctx.ac).as_bytes());
    }
    out.extend_from_slice(&ctx.info_buf);

    // --- FORMAT + samples --------------------------------------------------
    // With no samples anywhere the record is sites-only and has no FORMAT.
    if n_out > 0 {
        out.push(b'\t');
        if uniform {
            out.extend_from_slice(first_fmt);
        } else {
            for (i, k) in ctx.fmt_keys.iter().enumerate() {
                if i > 0 {
                    out.push(b':');
                }
                out.extend_from_slice(k);
            }
        }

        let mut next_out_sample = 0usize;
        for (mi, &m) in members.iter().enumerate() {
            let s = &streams[m];
            while next_out_sample < s.sample_off {
                push_missing(out, n_keys, gt_slot, ploidy, a.missing_to_ref);
                next_out_sample += 1;
            }
            let passthrough = ctx.identity[mi]
                && s.rec.ncols() > COL_FORMAT
                && (uniform || format_matches(s.get(COL_FORMAT), &ctx.fmt_keys));
            if passthrough && s.rec.ncols() == FIRST_SAMPLE + s.n_samples && s.n_samples > 0 {
                // Every sample column is reproduced verbatim, so copy the whole
                // run of them, tab separators included, in one memcpy.
                out.push(b'\t');
                out.extend_from_slice(s.rec.span(&s.line, FIRST_SAMPLE, s.rec.ncols() - 1));
            } else {
                for j in 0..s.n_samples {
                    let col = FIRST_SAMPLE + j;
                    let field: &[u8] = if col < s.rec.ncols() { s.get(col) } else { b"." };
                    out.push(b'\t');
                    if passthrough {
                        out.extend_from_slice(field);
                    } else {
                        write_sample(
                            out,
                            field,
                            &ctx.fmt_keys,
                            &ctx.key_pos[mi],
                            &ctx.maps[mi],
                            ctx.identity[mi],
                            n_alt,
                            hdr,
                            ploidy,
                        );
                    }
                }
            }
            next_out_sample += s.n_samples;
        }
        while next_out_sample < n_out {
            push_missing(out, n_keys, gt_slot, ploidy, a.missing_to_ref);
            next_out_sample += 1;
        }
    }
    out.push(b'\n');
    Ok(())
}

fn format_matches(fmt: &[u8], keys: &[Vec<u8>]) -> bool {
    let mut it = fmt.split(|&b| b == b':');
    for k in keys {
        match it.next() {
            Some(x) if x == &k[..] => {}
            _ => return false,
        }
    }
    it.next().is_none()
}

fn push_missing(out: &mut Vec<u8>, n_keys: usize, gt_slot: Option<usize>, ploidy: usize, to_ref: bool) {
    out.push(b'\t');
    for i in 0..n_keys {
        if i > 0 {
            out.push(b':');
        }
        if Some(i) == gt_slot {
            for p in 0..ploidy {
                if p > 0 {
                    out.push(b'/');
                }
                out.push(if to_ref { b'0' } else { b'.' });
            }
        } else {
            out.push(b'.');
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn write_sample(
    out: &mut Vec<u8>,
    field: &[u8],
    keys: &[Vec<u8>],
    key_pos: &[Option<usize>],
    map: &[i32],
    identity: bool,
    n_alt: usize,
    hdr: &Header,
    ploidy: usize,
) {
    for (ki, key) in keys.iter().enumerate() {
        if ki > 0 {
            out.push(b':');
        }
        let si = match key_pos[ki] {
            Some(si) => si,
            None => {
                if key == b"GT" {
                    for p in 0..ploidy {
                        if p > 0 {
                            out.push(b'/');
                        }
                        out.push(b'.');
                    }
                } else {
                    out.push(b'.');
                }
                continue;
            }
        };
        let v = subfield(field, si);
        if key == b"GT" {
            if identity {
                out.extend_from_slice(v);
            } else {
                remap_gt(v, map, out);
            }
            continue;
        }
        if identity {
            out.extend_from_slice(v);
            continue;
        }
        match hdr.format.get(std::str::from_utf8(key).unwrap_or("")) {
            Some(Number::A) => remap_values(v, map, n_alt, 1, out),
            Some(Number::R) => remap_values(v, map, n_alt + 1, 0, out),
            // Per-genotype values cannot be renumbered without knowing the
            // full genotype ordering; drop them rather than emit garbage.
            Some(Number::G) => out.push(b'.'),
            _ => out.extend_from_slice(v),
        }
    }
}

/// Renumber the allele indices of a GT field, keeping phasing separators.
fn remap_gt(gt: &[u8], map: &[i32], out: &mut Vec<u8>) {
    let mut cur: i64 = -1;
    let mut seen = false;
    let flush = |out: &mut Vec<u8>, cur: i64, seen: bool| {
        if !seen {
            out.push(b'.');
            return;
        }
        match map.get(cur as usize) {
            Some(&n) if n >= 0 => out.extend_from_slice(n.to_string().as_bytes()),
            _ => out.push(b'.'),
        }
    };
    for &b in gt {
        match b {
            b'0'..=b'9' => {
                if !seen {
                    cur = 0;
                    seen = true;
                }
                cur = cur * 10 + (b - b'0') as i64;
            }
            b'/' | b'|' => {
                flush(out, cur, seen);
                out.push(b);
                cur = -1;
                seen = false;
            }
            _ => {
                seen = false;
                cur = -1;
            }
        }
    }
    flush(out, cur, seen);
}

/// Re-index a comma-separated per-allele list. `offset` is 1 for Number=A
/// (values start at the first ALT) and 0 for Number=R.
fn remap_values(v: &[u8], map: &[i32], n_new: usize, offset: usize, out: &mut Vec<u8>) {
    let mut slots: Vec<&[u8]> = vec![b"."; n_new];
    for (i, val) in v.split(|&b| b == b',').enumerate() {
        let old_allele = i + offset;
        if let Some(&new) = map.get(old_allele) {
            if new >= 0 && new as usize >= offset {
                let target = new as usize - offset;
                if target < n_new {
                    slots[target] = val;
                }
            }
        }
    }
    for (i, s) in slots.iter().enumerate() {
        if i > 0 {
            out.push(b',');
        }
        out.extend_from_slice(s);
    }
}

fn write_ids(out: &mut Vec<u8>, streams: &[Stream], members: &[usize]) {
    let mut any = false;
    let mut seen: Vec<&[u8]> = Vec::new();
    for &m in members {
        let id = streams[m].get(COL_ID);
        if id == b"." || id.is_empty() {
            continue;
        }
        for part in id.split(|&b| b == b';') {
            if part.is_empty() || seen.contains(&part) {
                continue;
            }
            seen.push(part);
            if any {
                out.push(b';');
            }
            out.extend_from_slice(part);
            any = true;
        }
    }
    if !any {
        out.push(b'.');
    }
}

fn write_qual(out: &mut Vec<u8>, streams: &[Stream], members: &[usize]) {
    let mut best: Option<(f64, &[u8])> = None;
    for &m in members {
        let q = streams[m].get(COL_QUAL);
        if let Some(v) = std::str::from_utf8(q).ok().and_then(|s| s.parse::<f64>().ok()) {
            if best.map_or(true, |(b, _)| v > b) {
                best = Some((v, q));
            }
        }
    }
    match best {
        Some((_, raw)) => out.extend_from_slice(raw),
        None => out.push(b'.'),
    }
}

fn write_filter(out: &mut Vec<u8>, streams: &[Stream], members: &[usize]) {
    let mut seen: Vec<&[u8]> = Vec::new();
    let mut saw_pass = false;
    for &m in members {
        for f in streams[m].get(COL_FILTER).split(|&b| b == b';') {
            if f == b"PASS" {
                saw_pass = true;
            } else if f != b"." && !f.is_empty() && !seen.contains(&f) {
                seen.push(f);
            }
        }
    }
    if seen.is_empty() {
        out.extend_from_slice(if saw_pass { b"PASS" } else { b"." });
        return;
    }
    for (i, f) in seen.iter().enumerate() {
        if i > 0 {
            out.push(b';');
        }
        out.extend_from_slice(f);
    }
}

/// One INFO key at this site, gathered across the records being merged.
struct InfoEntry<'a> {
    key: &'a [u8],
    /// `None` for a flag.
    first: Option<Cow<'a, [u8]>>,
    /// Values from the remaining records. Only filled for keys whose result
    /// actually depends on more than the first record, which keeps this empty
    /// — and therefore unallocated — for most keys.
    extra: Vec<Cow<'a, [u8]>>,
    rule: Option<Rule>,
    per_allele: bool,
}

/// Combine the INFO columns of every record at this site.
///
/// A key listed in `--info-rules` is folded across all contributing records,
/// and a per-allele key takes each slot from the first record that has one;
/// every other key is simply the first record's value. Per-allele values are
/// renumbered onto the merged ALT list before any of that, so the rules see
/// comparable vectors.
#[allow(clippy::too_many_arguments)]
fn merge_info<'a>(
    out: &mut Vec<u8>,
    streams: &'a [Stream],
    members: &[usize],
    maps: &[Vec<i32>],
    identity: &[bool],
    n_alt: usize,
    hdr: &Header,
    rules: &InfoRules,
    skip: &[&[u8]],
) {
    let mut entries: Vec<InfoEntry<'a>> = Vec::new();
    for (mi, &m) in members.iter().enumerate() {
        let info = streams[m].get(COL_INFO);
        if info == b"." || info.is_empty() {
            continue;
        }
        for field in info.split(|&b| b == b';') {
            if field.is_empty() {
                continue;
            }
            let (key, val) = match memchr::memchr(b'=', field) {
                Some(p) => (&field[..p], Some(&field[p + 1..])),
                None => (field, None),
            };
            if skip.iter().any(|k| *k == key) {
                continue;
            }
            let number = hdr.info.get(std::str::from_utf8(key).unwrap_or("")).copied();
            // Per-genotype values cannot be renumbered; drop them when the
            // record's alleles were renumbered.
            if !identity[mi] && number == Some(Number::G) {
                continue;
            }
            let per_allele = matches!(number, Some(Number::A) | Some(Number::R));
            // Borrow the value unless this record's alleles were renumbered.
            let value: Option<Cow<'a, [u8]>> = val.map(|v| {
                if identity[mi] {
                    Cow::Borrowed(v)
                } else {
                    let mut b = Vec::with_capacity(v.len());
                    match number {
                        Some(Number::A) => remap_values(v, &maps[mi], n_alt, 1, &mut b),
                        Some(Number::R) => remap_values(v, &maps[mi], n_alt + 1, 0, &mut b),
                        _ => b.extend_from_slice(v),
                    }
                    Cow::Owned(b)
                }
            });
            match entries.iter_mut().find(|e| e.key == key) {
                Some(e) => {
                    if e.rule.is_some() || e.per_allele {
                        if let Some(v) = value {
                            e.extra.push(v);
                        }
                    }
                }
                None => entries.push(InfoEntry {
                    key,
                    first: value,
                    extra: Vec::new(),
                    rule: rules.get(key),
                    per_allele,
                }),
            }
        }
    }

    let mut first_out = true;
    let mut all: Vec<&[u8]> = Vec::new();
    for e in &entries {
        if !first_out {
            out.push(b';');
        }
        first_out = false;
        out.extend_from_slice(e.key);
        let v0 = match &e.first {
            None => continue, // flag
            Some(v) => v,
        };
        out.push(b'=');
        if e.extra.is_empty() {
            out.extend_from_slice(v0);
            continue;
        }
        all.clear();
        all.push(v0);
        all.extend(e.extra.iter().map(|c| &c[..]));
        match e.rule {
            Some(r) => fold_into(&all, r, out),
            // Per-allele values describe different alleles in each file, so
            // fill each slot from the first file that has it.
            None => fold_slots_into(&all, out),
        }
    }
    if first_out {
        out.push(b'.');
    }
}

/// Fold several comma-separated value vectors into one, element by element.
fn fold_into(vals: &[&[u8]], rule: Rule, out: &mut Vec<u8>) {
    if rule == Rule::Join {
        for (i, v) in vals.iter().enumerate() {
            if i > 0 {
                out.push(b',');
            }
            out.extend_from_slice(v);
        }
        return;
    }
    // Scalars are the overwhelmingly common case (DP, DP4 aside), and folding
    // them needs no per-value splitting.
    if vals.iter().all(|v| !v.contains(&b',')) {
        fold_one(vals, rule, out);
        return;
    }
    let cols: Vec<Vec<&[u8]>> = vals.iter().map(|v| v.split(|&b| b == b',').collect()).collect();
    let width = cols.iter().map(|c| c.len()).max().unwrap_or(0);
    let mut slot: Vec<&[u8]> = Vec::with_capacity(cols.len());
    for i in 0..width {
        if i > 0 {
            out.push(b',');
        }
        slot.clear();
        slot.extend(cols.iter().filter_map(|c| c.get(i).copied()));
        fold_one(&slot, rule, out);
    }
}

/// Fold one position's worth of values.
fn fold_one(items: &[&[u8]], rule: Rule, out: &mut Vec<u8>) {
    let present: Vec<&[u8]> = items
        .iter()
        .copied()
        .filter(|t| *t != b"." && !t.is_empty())
        .collect();
    if present.is_empty() {
        out.push(b'.');
        return;
    }
    // Integers stay integers; anything else is folded as a float.
    let ints: Option<Vec<i64>> = present
        .iter()
        .map(|t| std::str::from_utf8(t).ok().and_then(|s| s.parse::<i64>().ok()))
        .collect();
    if let Some(v) = ints {
        let r = match rule {
            Rule::Sum => v.iter().sum::<i64>(),
            Rule::Min => *v.iter().min().unwrap(),
            Rule::Max => *v.iter().max().unwrap(),
            Rule::Avg => {
                let avg = v.iter().sum::<i64>() as f64 / v.len() as f64;
                out.extend_from_slice(fmt_f64(avg).as_bytes());
                return;
            }
            _ => v[0],
        };
        out.extend_from_slice(r.to_string().as_bytes());
        return;
    }
    let fs: Vec<f64> = present
        .iter()
        .filter_map(|t| std::str::from_utf8(t).ok().and_then(|s| s.parse::<f64>().ok()))
        .collect();
    if fs.is_empty() {
        out.extend_from_slice(present[0]);
        return;
    }
    let v = match rule {
        Rule::Sum => fs.iter().sum::<f64>(),
        Rule::Min => fs.iter().cloned().fold(f64::INFINITY, f64::min),
        Rule::Max => fs.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
        Rule::Avg => fs.iter().sum::<f64>() / fs.len() as f64,
        _ => fs[0],
    };
    out.extend_from_slice(fmt_f64(v).as_bytes());
}

/// Take each per-allele slot from the first record that supplies one.
fn fold_slots_into(vals: &[&[u8]], out: &mut Vec<u8>) {
    let cols: Vec<Vec<&[u8]>> = vals.iter().map(|v| v.split(|&b| b == b',').collect()).collect();
    let width = cols.iter().map(|c| c.len()).max().unwrap_or(0);
    for i in 0..width {
        if i > 0 {
            out.push(b',');
        }
        let pick = cols
            .iter()
            .filter_map(|c| c.get(i).copied())
            .find(|t| *t != b"." && !t.is_empty());
        out.extend_from_slice(pick.unwrap_or(b"."));
    }
}

fn fmt_f64(v: f64) -> String {
    if v.fract() == 0.0 && v.abs() < 1e15 {
        format!("{}", v as i64)
    } else {
        let s = format!("{v:.6}");
        s.trim_end_matches('0').trim_end_matches('.').to_string()
    }
}

fn collect_inputs(a: &MergeArgs) -> Result<Vec<String>, String> {
    let mut v = a.inputs.clone();
    if let Some(f) = &a.file_list {
        let text = std::fs::read_to_string(f).map_err(|e| format!("{f}: {e}"))?;
        v.extend(
            text.lines()
                .map(|l| l.trim().to_string())
                .filter(|l| !l.is_empty() && !l.starts_with('#')),
        );
    }
    Ok(v)
}
