//! VCF header parsing, merging and sample selection.

use std::collections::HashMap;
use std::io;

use crate::lines::{BufSource, LineReader};

/// The `Number=` attribute of an INFO/FORMAT declaration. It decides whether a
/// value has to be re-ordered when ALT alleles are renumbered during a merge.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Number {
    /// One value per ALT allele.
    A,
    /// One value per allele, REF included.
    R,
    /// One value per genotype.
    G,
    Fixed(u32),
    Unknown,
}

impl Number {
    fn parse(s: &str) -> Number {
        match s {
            "A" => Number::A,
            "R" => Number::R,
            "G" => Number::G,
            "." => Number::Unknown,
            n => n.parse::<u32>().map(Number::Fixed).unwrap_or(Number::Unknown),
        }
    }
}

#[derive(Clone, Debug)]
pub struct Header {
    /// `##` lines, verbatim and in file order.
    pub meta: Vec<Vec<u8>>,
    /// Columns of the `#CHROM` line beyond the fixed nine.
    pub samples: Vec<String>,
    /// Contigs in declaration order, for position ordering checks.
    pub contigs: Vec<String>,
    contig_rank: HashMap<String, usize>,
    pub info: HashMap<String, Number>,
    pub format: HashMap<String, Number>,
    /// True when the `#CHROM` line had a FORMAT column.
    pub has_format_col: bool,
}

impl Default for Header {
    fn default() -> Self {
        Header {
            meta: Vec::new(),
            samples: Vec::new(),
            contigs: Vec::new(),
            contig_rank: HashMap::new(),
            info: HashMap::new(),
            format: HashMap::new(),
            has_format_col: false,
        }
    }
}

impl Header {
    /// Consume header lines from `r`, leaving it positioned at the first record.
    pub fn read<S: BufSource>(r: &mut LineReader<S>) -> io::Result<Header> {
        let mut h = Header::default();
        let mut saw_chrom = false;
        while r.advance()? {
            let line = r.line();
            if line.is_empty() {
                continue;
            }
            if line.starts_with(b"##") {
                h.push_meta(line.to_vec());
            } else if line.starts_with(b"#CHROM") {
                let cols: Vec<&[u8]> = line.split(|&b| b == b'\t').collect();
                if cols.len() > 8 {
                    h.has_format_col = cols[8] == b"FORMAT";
                }
                for c in cols.iter().skip(9) {
                    h.samples.push(String::from_utf8_lossy(c).into_owned());
                }
                saw_chrom = true;
                break;
            } else {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "record found before the #CHROM header line",
                ));
            }
        }
        if !saw_chrom {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "missing #CHROM header line",
            ));
        }
        Ok(h)
    }

    pub fn push_meta(&mut self, line: Vec<u8>) {
        {
            let text = String::from_utf8_lossy(&line);
            if let Some(rest) = text.strip_prefix("##contig=<") {
                if let Some(id) = attr(rest, "ID") {
                    if !self.contig_rank.contains_key(&id) {
                        self.contig_rank.insert(id.clone(), self.contigs.len());
                        self.contigs.push(id);
                    }
                }
            } else if let Some(rest) = text.strip_prefix("##INFO=<") {
                if let (Some(id), Some(n)) = (attr(rest, "ID"), attr(rest, "Number")) {
                    self.info.insert(id, Number::parse(&n));
                }
            } else if let Some(rest) = text.strip_prefix("##FORMAT=<") {
                if let (Some(id), Some(n)) = (attr(rest, "ID"), attr(rest, "Number")) {
                    self.format.insert(id, Number::parse(&n));
                }
            }
        }
        self.meta.push(line);
    }

    pub fn contig_rank(&self, chrom: &[u8]) -> Option<usize> {
        self.contig_rank
            .get(std::str::from_utf8(chrom).ok()?)
            .copied()
    }

    pub fn has_meta_id(&self, key: &str, id: &str) -> bool {
        let needle = format!("##{key}=<");
        self.meta.iter().any(|l| {
            let t = String::from_utf8_lossy(l);
            t.strip_prefix(&needle)
                .and_then(|r| attr(r, "ID"))
                .map_or(false, |v| v == id)
        })
    }

    /// Serialise the header, optionally replacing the sample list.
    pub fn write(&self, out: &mut Vec<u8>, samples: &[String], drop_genotypes: bool) {
        let mut wrote_fileformat = false;
        for l in &self.meta {
            if l.starts_with(b"##fileformat=") {
                if wrote_fileformat {
                    continue;
                }
                wrote_fileformat = true;
            }
            out.extend_from_slice(l);
            out.push(b'\n');
        }
        out.extend_from_slice(b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
        if !drop_genotypes && (!samples.is_empty() || self.has_format_col) {
            out.extend_from_slice(b"\tFORMAT");
            for s in samples {
                out.push(b'\t');
                out.extend_from_slice(s.as_bytes());
            }
        }
        out.push(b'\n');
    }

    /// Add a meta line only if a declaration with the same key/ID is absent.
    pub fn ensure_meta(&mut self, key: &str, id: &str, line: &str) {
        if !self.has_meta_id(key, id) {
            self.push_meta(line.as_bytes().to_vec());
        }
    }
}

/// Structured header lines share a `KEY=<ID=...,Number=...>` shape; pull one
/// attribute out of the part after `<`, honouring double-quoted values.
pub fn attr(body: &str, key: &str) -> Option<String> {
    let bytes = body.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        // Read a key.
        let ks = i;
        while i < bytes.len() && bytes[i] != b'=' && bytes[i] != b',' && bytes[i] != b'>' {
            i += 1;
        }
        if i >= bytes.len() || bytes[i] != b'=' {
            return None;
        }
        let k = &body[ks..i];
        i += 1;
        let vs = i;
        let v = if bytes.get(i) == Some(&b'"') {
            i += 1;
            let inner = i;
            while i < bytes.len() && bytes[i] != b'"' {
                if bytes[i] == b'\\' {
                    i += 1;
                }
                i += 1;
            }
            let v = &body[inner..i.min(bytes.len())];
            i += 1;
            v
        } else {
            while i < bytes.len() && bytes[i] != b',' && bytes[i] != b'>' {
                i += 1;
            }
            &body[vs..i]
        };
        if k == key {
            return Some(v.to_string());
        }
        if i < bytes.len() && bytes[i] == b',' {
            i += 1;
        } else {
            return None;
        }
    }
    None
}

/// Deduplication key for a meta line: `(KEY, ID)` for structured lines,
/// otherwise the whole line.
fn meta_key(line: &[u8]) -> String {
    let t = String::from_utf8_lossy(line);
    if let Some(rest) = t.strip_prefix("##") {
        if let Some(eq) = rest.find('=') {
            let (key, val) = rest.split_at(eq);
            if val.starts_with("=<") {
                if let Some(id) = attr(&val[2..], "ID") {
                    return format!("{key}\u{1}{id}");
                }
            }
        }
    }
    t.into_owned()
}

/// Union of several headers' meta lines: first definition of each ID wins,
/// declaration order of the first header is preserved.
pub fn union_meta(headers: &[Header]) -> Header {
    let mut out = Header::default();
    let mut seen = std::collections::HashSet::new();
    let mut fileformat: Option<Vec<u8>> = None;
    for h in headers {
        for l in &h.meta {
            if l.starts_with(b"##fileformat=") {
                if fileformat.is_none() {
                    fileformat = Some(l.clone());
                }
                continue;
            }
            if seen.insert(meta_key(l)) {
                out.push_meta(l.clone());
            }
        }
        out.has_format_col |= h.has_format_col;
    }
    if let Some(f) = fileformat {
        out.meta.insert(0, f);
    } else {
        out.meta.insert(0, b"##fileformat=VCFv4.2".to_vec());
    }
    out
}

/// Resolve a `-s`/`-S` selection into indices into `all`, in output order.
/// A leading `^` inverts the selection.
pub fn select_samples(all: &[String], spec: &str) -> Result<Vec<usize>, String> {
    let (invert, body) = match spec.strip_prefix('^') {
        Some(rest) => (true, rest),
        None => (false, spec),
    };
    let index: HashMap<&str, usize> = all.iter().enumerate().map(|(i, s)| (s.as_str(), i)).collect();
    let wanted: Vec<&str> = body.split(',').filter(|s| !s.is_empty()).collect();
    let mut picked = Vec::with_capacity(wanted.len());
    for name in &wanted {
        match index.get(name) {
            Some(&i) => picked.push(i),
            None => return Err(format!("sample '{name}' is not present in the header")),
        }
    }
    if !invert {
        return Ok(picked);
    }
    let drop: std::collections::HashSet<usize> = picked.into_iter().collect();
    Ok((0..all.len()).filter(|i| !drop.contains(i)).collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_structured_attributes() {
        let b = r#"ID=DP,Number=1,Type=Integer,Description="Total, depth=x">"#;
        assert_eq!(attr(b, "ID").as_deref(), Some("DP"));
        assert_eq!(attr(b, "Number").as_deref(), Some("1"));
        // A comma or '=' inside a quoted description must not end the value.
        assert_eq!(attr(b, "Description").as_deref(), Some("Total, depth=x"));
        assert_eq!(attr(b, "Missing"), None);
    }

    #[test]
    fn parses_number_attribute() {
        assert_eq!(Number::parse("A"), Number::A);
        assert_eq!(Number::parse("R"), Number::R);
        assert_eq!(Number::parse("G"), Number::G);
        assert_eq!(Number::parse("."), Number::Unknown);
        assert_eq!(Number::parse("3"), Number::Fixed(3));
        assert_eq!(Number::parse("junk"), Number::Unknown);
    }

    fn header_of(text: &str) -> Header {
        let mut r = crate::lines::LineReader::new(crate::lines::PlainSource::new(
            std::io::Cursor::new(text.as_bytes().to_vec()),
        ));
        Header::read(&mut r).unwrap()
    }

    const H: &str = "##fileformat=VCFv4.2\n\
##contig=<ID=chr1,length=10>\n\
##contig=<ID=chr2,length=20>\n\
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"ac\">\n\
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"ad\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\n";

    #[test]
    fn reads_samples_contigs_and_field_numbers() {
        let h = header_of(H);
        assert_eq!(h.samples, ["s1", "s2", "s3"]);
        assert_eq!(h.contigs, ["chr1", "chr2"]);
        assert_eq!(h.contig_rank(b"chr2"), Some(1));
        assert_eq!(h.contig_rank(b"chrX"), None);
        assert_eq!(h.info.get("AC"), Some(&Number::A));
        assert_eq!(h.format.get("AD"), Some(&Number::R));
        assert!(h.has_format_col);
    }

    #[test]
    fn selects_samples_by_name_and_by_exclusion() {
        let all: Vec<String> = ["s1", "s2", "s3"].iter().map(|s| s.to_string()).collect();
        assert_eq!(select_samples(&all, "s3,s1").unwrap(), [2, 0]);
        assert_eq!(select_samples(&all, "^s2").unwrap(), [0, 2]);
        assert!(select_samples(&all, "nope").is_err());
    }

    #[test]
    fn unions_headers_keeping_the_first_definition() {
        let a = header_of(H);
        let b = header_of(
            "##fileformat=VCFv4.3\n\
##contig=<ID=chr2,length=999>\n\
##contig=<ID=chr3,length=30>\n\
##INFO=<ID=AC,Number=1,Type=Integer,Description=\"different\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tz1\n",
        );
        let u = union_meta(&[a, b]);
        assert_eq!(u.contigs, ["chr1", "chr2", "chr3"]);
        assert_eq!(u.info.get("AC"), Some(&Number::A), "first definition wins");
        assert_eq!(u.meta.iter().filter(|l| l.starts_with(b"##fileformat")).count(), 1);
        assert!(u.meta[0].starts_with(b"##fileformat=VCFv4.2"));
        assert!(u.meta.iter().filter(|l| l.starts_with(b"##contig=<ID=chr2")).count() == 1);
    }

    #[test]
    fn writes_the_chrom_line_for_the_given_samples() {
        let h = header_of(H);
        let mut out = Vec::new();
        h.write(&mut out, &["s2".into()], false);
        let text = String::from_utf8(out).unwrap();
        assert!(text.ends_with("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts2\n"));
        let mut out = Vec::new();
        h.write(&mut out, &["s2".into()], true);
        let text = String::from_utf8(out).unwrap();
        assert!(text.ends_with("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"));
    }
}
