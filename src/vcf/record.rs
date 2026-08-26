//! Record-level primitives.
//!
//! Records are never fully parsed. A line is split into column offsets and only
//! the fields a command actually needs are decoded, so passing a record through
//! unchanged costs one `memchr` scan and one `write_all`.

use memchr::memchr_iter;

pub const COL_CHROM: usize = 0;
pub const COL_POS: usize = 1;
pub const COL_ID: usize = 2;
pub const COL_REF: usize = 3;
pub const COL_ALT: usize = 4;
pub const COL_QUAL: usize = 5;
pub const COL_FILTER: usize = 6;
pub const COL_INFO: usize = 7;
pub const COL_FORMAT: usize = 8;
pub const FIRST_SAMPLE: usize = 9;

/// Column boundaries within a record line.
///
/// `starts[i]` is the offset of column `i`; a sentinel one past the end lets
/// `end(i)` be computed as `starts[i + 1] - 1` without a branch.
#[derive(Default, Clone)]
pub struct Record {
    starts: Vec<u32>,
}

impl Record {
    pub fn new() -> Record {
        Record { starts: Vec::with_capacity(64) }
    }

    /// Index every column of `line`.
    #[inline]
    pub fn split(&mut self, line: &[u8]) {
        self.starts.clear();
        self.starts.push(0);
        for p in memchr_iter(b'\t', line) {
            self.starts.push(p as u32 + 1);
        }
        self.starts.push(line.len() as u32 + 1);
    }

    /// Index only the first `max` columns; the final entry spans the remainder
    /// of the line. Used when sample columns are irrelevant.
    #[inline]
    pub fn split_prefix(&mut self, line: &[u8], max: usize) {
        self.starts.clear();
        self.starts.push(0);
        for p in memchr_iter(b'\t', line) {
            if self.starts.len() > max {
                break;
            }
            self.starts.push(p as u32 + 1);
        }
        self.starts.push(line.len() as u32 + 1);
    }

    #[inline]
    pub fn ncols(&self) -> usize {
        self.starts.len() - 1
    }

    #[inline]
    pub fn get<'a>(&self, line: &'a [u8], i: usize) -> &'a [u8] {
        let s = self.starts[i] as usize;
        let e = (self.starts[i + 1] - 1) as usize;
        &line[s..e.min(line.len())]
    }

    /// The slice covering columns `i..=j`, separators included.
    #[inline]
    pub fn span<'a>(&self, line: &'a [u8], i: usize, j: usize) -> &'a [u8] {
        let s = self.starts[i] as usize;
        let e = (self.starts[j + 1] - 1) as usize;
        &line[s..e.min(line.len())]
    }

}

#[inline]
pub fn parse_u64(b: &[u8]) -> Option<u64> {
    if b.is_empty() {
        return None;
    }
    let mut v: u64 = 0;
    for &c in b {
        if !c.is_ascii_digit() {
            return None;
        }
        v = v.wrapping_mul(10).wrapping_add((c - b'0') as u64);
    }
    Some(v)
}

pub const T_SNP: u8 = 1;
pub const T_INDEL: u8 = 2;
pub const T_MNP: u8 = 4;
pub const T_OTHER: u8 = 8;
pub const T_REF: u8 = 16;

pub fn parse_types(spec: &str) -> Result<u8, String> {
    let mut m = 0;
    for t in spec.split(',').filter(|s| !s.is_empty()) {
        m |= match t {
            "snps" | "snp" => T_SNP,
            "indels" | "indel" => T_INDEL,
            "mnps" | "mnp" => T_MNP,
            "other" => T_OTHER,
            "ref" => T_REF,
            _ => return Err(format!("unknown variant type '{t}'")),
        };
    }
    Ok(m)
}

fn allele_type(reference: &[u8], alt: &[u8]) -> u8 {
    if alt.is_empty() || alt == b"." {
        return T_REF;
    }
    if alt == b"*" || alt[0] == b'<' || alt.contains(&b'[') || alt.contains(&b']') {
        return T_OTHER;
    }
    if !alt.iter().all(|c| matches!(c, b'A' | b'C' | b'G' | b'T' | b'N' | b'a' | b'c' | b'g' | b't' | b'n')) {
        return T_OTHER;
    }
    if alt == reference {
        return T_REF;
    }
    if reference.len() == alt.len() {
        if reference.len() == 1 {
            T_SNP
        } else {
            // A single substitution inside equal-length alleles is still a SNP.
            let diff = reference.iter().zip(alt).filter(|(a, b)| a != b).count();
            if diff == 1 {
                T_SNP
            } else {
                T_MNP
            }
        }
    } else {
        T_INDEL
    }
}

/// Union of the types of every ALT allele at a site.
pub fn site_types(reference: &[u8], alt_col: &[u8]) -> u8 {
    if alt_col == b"." {
        return T_REF;
    }
    let mut m = 0;
    for a in alt_col.split(|&b| b == b',') {
        m |= allele_type(reference, a);
    }
    m
}

pub fn count_alts(alt_col: &[u8]) -> usize {
    if alt_col == b"." || alt_col.is_empty() {
        0
    } else {
        alt_col.iter().filter(|&&b| b == b',').count() + 1
    }
}

/// Position of `GT` among the colon-separated keys of a FORMAT column.
pub fn gt_index(format_col: &[u8]) -> Option<usize> {
    format_col.split(|&b| b == b':').position(|k| k == b"GT")
}

/// The `n`-th colon-separated subfield of a sample column.
#[inline]
pub fn subfield(sample: &[u8], n: usize) -> &[u8] {
    let mut start = 0;
    let mut idx = 0;
    for (i, &b) in sample.iter().enumerate() {
        if b == b':' {
            if idx == n {
                return &sample[start..i];
            }
            idx += 1;
            start = i + 1;
        }
    }
    if idx == n {
        &sample[start..]
    } else {
        b""
    }
}

/// Decode a GT field into allele indices; `-1` stands for a missing allele.
/// Returns `true` if any separator was `|`.
#[inline]
pub fn parse_gt(gt: &[u8], out: &mut Vec<i32>) -> bool {
    out.clear();
    let mut phased = false;
    let mut cur: i32 = -1;
    let mut seen_digit = false;
    for &b in gt {
        match b {
            b'0'..=b'9' => {
                if !seen_digit {
                    cur = 0;
                    seen_digit = true;
                }
                cur = cur.saturating_mul(10).saturating_add((b - b'0') as i32);
            }
            b'/' | b'|' => {
                phased |= b == b'|';
                out.push(if seen_digit { cur } else { -1 });
                cur = -1;
                seen_digit = false;
            }
            b'.' => {
                seen_digit = false;
                cur = -1;
            }
            _ => {}
        }
    }
    out.push(if seen_digit { cur } else { -1 });
    phased
}

/// Accumulate allele counts straight from a sample column.
///
/// This runs once per sample per site — 125 million times over a three-way
/// merge of 250k sites — so it reads the GT subfield in place: it stops at the
/// first `:` rather than slicing the subfield out first, and it never
/// materialises the allele list.
///
/// `map` renumbers alleles when the record's ALTs were merged into a longer
/// list; `None` means the numbering already matches.
#[inline]
pub fn count_gt(field: &[u8], map: Option<&[i32]>, n_alt: usize, ac: &mut [u64], an: &mut u64) {
    let mut i = 0;
    while i < field.len() {
        let b = field[i];
        if b == b':' {
            return;
        }
        if b.is_ascii_digit() {
            let mut v = (b - b'0') as usize;
            i += 1;
            while i < field.len() {
                let c = field[i];
                if !c.is_ascii_digit() {
                    break;
                }
                v = v * 10 + (c - b'0') as usize;
                i += 1;
            }
            // A called allele counts towards AN whether or not it survives the
            // renumbering.
            *an += 1;
            let a = match map {
                None => v,
                Some(m) => match m.get(v) {
                    Some(&x) if x > 0 => x as usize,
                    _ => continue,
                },
            };
            if a > 0 && a <= n_alt {
                ac[a - 1] += 1;
            }
        } else {
            i += 1;
        }
    }
}

/// Allele counts over the selected sample columns.
pub struct AlleleCounts {
    pub ac: Vec<u64>,
    pub an: u64,
}

pub fn count_alleles(
    line: &[u8],
    rec: &Record,
    samples: &[usize],
    n_alt: usize,
    gt_idx: usize,
    scratch: &mut Vec<i32>,
) -> AlleleCounts {
    let mut ac = vec![0u64; n_alt];
    let mut an = 0u64;
    for &s in samples {
        let col = FIRST_SAMPLE + s;
        if col >= rec.ncols() {
            continue;
        }
        let field = rec.get(line, col);
        let gt = if gt_idx == 0 {
            match memchr::memchr(b':', field) {
                Some(p) => &field[..p],
                None => field,
            }
        } else {
            subfield(field, gt_idx)
        };
        parse_gt(gt, scratch);
        for &a in scratch.iter() {
            if a < 0 {
                continue;
            }
            an += 1;
            if a > 0 && (a as usize) <= n_alt {
                ac[a as usize - 1] += 1;
            }
        }
    }
    AlleleCounts { ac, an }
}

/// Rewrite an INFO column, replacing the listed keys and appending any that
/// were absent. `updates` values of `None` drop the key.
pub fn rewrite_info(info: &[u8], updates: &[(&str, Option<String>)], out: &mut Vec<u8>) {
    let mut done = vec![false; updates.len()];
    let mut first = true;
    if info != b"." {
        for field in info.split(|&b| b == b';') {
            if field.is_empty() {
                continue;
            }
            let key = match memchr::memchr(b'=', field) {
                Some(p) => &field[..p],
                None => field,
            };
            let mut replaced = false;
            for (i, (k, v)) in updates.iter().enumerate() {
                if key == k.as_bytes() {
                    done[i] = true;
                    replaced = true;
                    if let Some(val) = v {
                        if !first {
                            out.push(b';');
                        }
                        out.extend_from_slice(k.as_bytes());
                        out.push(b'=');
                        out.extend_from_slice(val.as_bytes());
                        first = false;
                    }
                    break;
                }
            }
            if !replaced {
                if !first {
                    out.push(b';');
                }
                out.extend_from_slice(field);
                first = false;
            }
        }
    }
    for (i, (k, v)) in updates.iter().enumerate() {
        if done[i] {
            continue;
        }
        if let Some(val) = v {
            if !first {
                out.push(b';');
            }
            out.extend_from_slice(k.as_bytes());
            out.push(b'=');
            out.extend_from_slice(val.as_bytes());
            first = false;
        }
    }
    if first {
        out.push(b'.');
    }
}

pub fn join_u64(vals: &[u64]) -> String {
    let mut s = String::new();
    for (i, v) in vals.iter().enumerate() {
        if i > 0 {
            s.push(',');
        }
        s.push_str(&v.to_string());
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    const LINE: &[u8] = b"chr1\t100\trs1\tA\tG,T\t50\tPASS\tAC=1,2;AN=6\tGT:AD\t0/1:5,6,0\t1|2:1,2,3\t./.:.";

    #[test]
    fn splits_columns_without_copying() {
        let mut r = Record::new();
        r.split(LINE);
        assert_eq!(r.ncols(), 12);
        assert_eq!(r.get(LINE, COL_CHROM), b"chr1");
        assert_eq!(r.get(LINE, COL_POS), b"100");
        assert_eq!(r.get(LINE, COL_ALT), b"G,T");
        assert_eq!(r.get(LINE, COL_INFO), b"AC=1,2;AN=6");
        assert_eq!(r.get(LINE, FIRST_SAMPLE + 2), b"./.:.");
        assert_eq!(r.span(LINE, COL_CHROM, COL_ID), b"chr1\t100\trs1");
    }

    #[test]
    fn prefix_split_stops_early() {
        let mut r = Record::new();
        r.split_prefix(LINE, COL_INFO + 1);
        assert_eq!(r.get(LINE, COL_INFO), b"AC=1,2;AN=6");
        assert!(r.ncols() <= COL_INFO + 2);
    }

    #[test]
    fn classifies_variant_types() {
        assert_eq!(site_types(b"A", b"G"), T_SNP);
        assert_eq!(site_types(b"A", b"AT"), T_INDEL);
        assert_eq!(site_types(b"AT", b"A"), T_INDEL);
        assert_eq!(site_types(b"AT", b"GC"), T_MNP);
        assert_eq!(site_types(b"AT", b"GT"), T_SNP, "one substitution is a SNP");
        assert_eq!(site_types(b"A", b"<DEL>"), T_OTHER);
        assert_eq!(site_types(b"A", b"*"), T_OTHER);
        assert_eq!(site_types(b"A", b"."), T_REF);
        assert_eq!(site_types(b"A", b"G,AT"), T_SNP | T_INDEL);
    }

    #[test]
    fn counts_alts_and_finds_gt() {
        assert_eq!(count_alts(b"G"), 1);
        assert_eq!(count_alts(b"G,T"), 2);
        assert_eq!(count_alts(b"."), 0);
        assert_eq!(gt_index(b"GT:AD"), Some(0));
        assert_eq!(gt_index(b"AD:DP:GT"), Some(2));
        assert_eq!(gt_index(b"AD:DP"), None);
    }

    #[test]
    fn reads_subfields() {
        assert_eq!(subfield(b"0/1:5,6:9", 0), b"0/1");
        assert_eq!(subfield(b"0/1:5,6:9", 1), b"5,6");
        assert_eq!(subfield(b"0/1:5,6:9", 2), b"9");
        assert_eq!(subfield(b"0/1:5,6:9", 3), b"");
    }

    #[test]
    fn parses_genotypes() {
        let mut v = Vec::new();
        assert!(!parse_gt(b"0/1", &mut v));
        assert_eq!(v, [0, 1]);
        assert!(parse_gt(b"1|2", &mut v));
        assert_eq!(v, [1, 2]);
        parse_gt(b"./.", &mut v);
        assert_eq!(v, [-1, -1]);
        parse_gt(b"12/3", &mut v);
        assert_eq!(v, [12, 3]);
        parse_gt(b"1", &mut v);
        assert_eq!(v, [1], "haploid");
        parse_gt(b"0/1/1", &mut v);
        assert_eq!(v, [0, 1, 1], "triploid");
    }

    #[test]
    fn counts_alleles_over_selected_samples() {
        let mut r = Record::new();
        r.split(LINE);
        let mut scratch = Vec::new();
        // Allele 1 is carried by both called samples, allele 2 by one.
        let c = count_alleles(LINE, &r, &[0, 1, 2], 2, 0, &mut scratch);
        assert_eq!(c.ac, [2, 1]);
        assert_eq!(c.an, 4, "the ./. sample contributes nothing");
        let c = count_alleles(LINE, &r, &[0], 2, 0, &mut scratch);
        assert_eq!(c.ac, [1, 0]);
        assert_eq!(c.an, 2);
    }

    #[test]
    fn rewrites_info_in_place_and_appends_new_keys() {
        let mut out = Vec::new();
        rewrite_info(b"AC=1;AN=2;DP=9", &[("AC", Some("5".into())), ("AN", Some("6".into()))], &mut out);
        assert_eq!(out, b"AC=5;AN=6;DP=9");

        out.clear();
        rewrite_info(b"DP=9;FLAG", &[("AC", Some("5".into()))], &mut out);
        assert_eq!(out, b"DP=9;FLAG;AC=5");

        out.clear();
        rewrite_info(b"DP=9", &[("DP", None)], &mut out);
        assert_eq!(out, b".", "an emptied INFO becomes '.'");

        out.clear();
        rewrite_info(b".", &[("AN", Some("2".into()))], &mut out);
        assert_eq!(out, b"AN=2");
    }
}
