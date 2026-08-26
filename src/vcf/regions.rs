//! Region selection for streaming (index-free) subsetting.

use std::collections::HashMap;
use std::fs;

#[derive(Clone, Copy, Debug)]
pub struct Interval {
    pub start: u64,
    pub end: u64,
}

#[derive(Debug, Default)]
pub struct RegionSet {
    by_chrom: HashMap<String, Vec<Interval>>,
    /// Furthest (contig rank, end) any region reaches; enables early exit on a
    /// coordinate-sorted input.
    stop_after: Option<(usize, u64)>,
}

impl RegionSet {
    /// Parse a comma-separated list such as `chr1,chr2:1000-2000,chr3:500`.
    pub fn parse(spec: &str) -> Result<RegionSet, String> {
        let mut rs = RegionSet::default();
        for item in spec.split(',').filter(|s| !s.is_empty()) {
            let (chrom, iv) = parse_one(item)?;
            rs.add(chrom, iv);
        }
        rs.normalise();
        Ok(rs)
    }

    /// Parse a regions file: `CHROM`, `CHROM<TAB>POS`, or `CHROM<TAB>BEG<TAB>END`
    /// with 1-based inclusive coordinates, matching `bcftools -R`.
    pub fn from_file(path: &str) -> Result<RegionSet, String> {
        let text = fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))?;
        let mut rs = RegionSet::default();
        for (n, line) in text.lines().enumerate() {
            let line = line.trim_end();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            let iv = match f.len() {
                1 => Interval { start: 1, end: u64::MAX },
                2 => {
                    let p = num(f[1], path, n)?;
                    Interval { start: p, end: p }
                }
                _ => Interval {
                    start: num(f[1], path, n)?,
                    end: num(f[2], path, n)?,
                },
            };
            rs.add(f[0].to_string(), iv);
        }
        rs.normalise();
        Ok(rs)
    }

    fn add(&mut self, chrom: String, iv: Interval) {
        self.by_chrom.entry(chrom).or_default().push(iv);
    }

    fn normalise(&mut self) {
        for v in self.by_chrom.values_mut() {
            v.sort_by_key(|i| (i.start, i.end));
            let mut merged: Vec<Interval> = Vec::with_capacity(v.len());
            for iv in v.iter() {
                match merged.last_mut() {
                    Some(last) if iv.start <= last.end.saturating_add(1) => {
                        last.end = last.end.max(iv.end)
                    }
                    _ => merged.push(*iv),
                }
            }
            *v = merged;
        }
    }

    /// Precompute the early-exit point using the header's contig order.
    pub fn set_contig_order(&mut self, rank: impl Fn(&str) -> Option<usize>) {
        let mut best: Option<(usize, u64)> = None;
        for (chrom, ivs) in &self.by_chrom {
            let r = match rank(chrom) {
                Some(r) => r,
                // Unknown contig: we cannot reason about ordering, so never
                // stop early.
                None => return,
            };
            let end = ivs.last().map(|i| i.end).unwrap_or(0);
            best = match best {
                Some(b) if b >= (r, end) => Some(b),
                _ => Some((r, end)),
            };
        }
        self.stop_after = best;
    }

    #[inline]
    pub fn contains(&self, chrom: &[u8], pos: u64, end: u64) -> bool {
        let key = match std::str::from_utf8(chrom) {
            Ok(s) => s,
            Err(_) => return false,
        };
        let ivs = match self.by_chrom.get(key) {
            Some(v) => v,
            None => return false,
        };
        // First interval that could overlap [pos, end].
        let idx = ivs.partition_point(|i| i.end < pos);
        matches!(ivs.get(idx), Some(i) if i.start <= end)
    }

    /// True once a coordinate-sorted stream has moved past every region.
    #[inline]
    pub fn past_all(&self, rank: Option<usize>, pos: u64) -> bool {
        match (self.stop_after, rank) {
            (Some((r, e)), Some(cur)) => cur > r || (cur == r && pos > e),
            _ => false,
        }
    }
}

fn num(s: &str, path: &str, line: usize) -> Result<u64, String> {
    s.replace(',', "")
        .parse::<u64>()
        .map_err(|_| format!("{path}:{}: bad coordinate '{s}'", line + 1))
}

fn parse_one(item: &str) -> Result<(String, Interval), String> {
    // Split on the last ':' so contigs containing ':' still work.
    match item.rfind(':') {
        None => Ok((item.to_string(), Interval { start: 1, end: u64::MAX })),
        Some(p) => {
            let (chrom, rest) = item.split_at(p);
            let rest = &rest[1..];
            let parse = |s: &str| -> Result<u64, String> {
                s.replace(',', "")
                    .parse::<u64>()
                    .map_err(|_| format!("bad region '{item}'"))
            };
            let iv = match rest.split_once('-') {
                None => {
                    let p = parse(rest)?;
                    Interval { start: p, end: p }
                }
                Some((a, "")) => Interval { start: parse(a)?, end: u64::MAX },
                Some(("", b)) => Interval { start: 1, end: parse(b)? },
                Some((a, b)) => Interval { start: parse(a)?, end: parse(b)? },
            };
            Ok((chrom.to_string(), iv))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_region_syntax() {
        let rs = RegionSet::parse("chr1,chr2:100-200,chr3:50,chr4:10-").unwrap();
        assert!(rs.contains(b"chr1", 1, 1));
        assert!(rs.contains(b"chr1", u64::MAX / 2, u64::MAX / 2));
        assert!(!rs.contains(b"chr2", 99, 99));
        assert!(rs.contains(b"chr2", 100, 100));
        assert!(rs.contains(b"chr2", 200, 200));
        assert!(!rs.contains(b"chr2", 201, 201));
        assert!(rs.contains(b"chr3", 50, 50));
        assert!(!rs.contains(b"chr3", 51, 51));
        assert!(rs.contains(b"chr4", 1_000_000, 1_000_000));
        assert!(!rs.contains(b"chr9", 1, 1));
    }

    #[test]
    fn a_record_overlapping_the_region_is_kept() {
        let rs = RegionSet::parse("chr1:100-200").unwrap();
        // A deletion starting before the region but reaching into it.
        assert!(rs.contains(b"chr1", 98, 101));
        assert!(!rs.contains(b"chr1", 95, 99));
    }

    #[test]
    fn merges_overlapping_intervals() {
        let rs = RegionSet::parse("chr1:100-200,chr1:150-300,chr1:301-400").unwrap();
        assert_eq!(rs.by_chrom["chr1"].len(), 1);
        assert!(rs.contains(b"chr1", 400, 400));
    }

    #[test]
    fn early_exit_needs_a_known_contig_order() {
        let mut rs = RegionSet::parse("chr1:100-200").unwrap();
        let order = |c: &str| ["chr1", "chr2"].iter().position(|x| *x == c);
        rs.set_contig_order(order);
        assert!(!rs.past_all(Some(0), 200));
        assert!(rs.past_all(Some(0), 201));
        assert!(rs.past_all(Some(1), 1), "a later contig is past everything");

        // An unknown contig must disable the optimisation entirely.
        let mut rs = RegionSet::parse("chrZ:100-200").unwrap();
        rs.set_contig_order(order);
        assert!(!rs.past_all(Some(5), 10_000));
    }
}
