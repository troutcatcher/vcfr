//! A DEFLATE encoder specialised for BGZF-sized blocks of VCF text.
//!
//! libdeflate must treat every block as unknown data: count symbol
//! frequencies, build Huffman codes, then emit — per 64 KiB block. VCF blocks
//! all look alike, so this encoder builds one canonical Huffman code set from
//! a sample block ("training") and then compresses every block in a single
//! greedy pass with no counting and no code construction. The output is
//! ordinary DEFLATE: one dynamic block per BGZF member, readable by any
//! inflater.
//!
//! This is not a general-purpose libdeflate replacement — it targets the fast
//! end of the curve on data whose statistics are stable across blocks, and it
//! is benchmarked against libdeflate rather than assumed faster.

mod high;
mod huffman;

pub use high::{compress_high_into, effort_for_level, HighMatcher};

const MAX_MATCH: usize = 258;
const MIN_MATCH: usize = 4;
const MAX_DIST: usize = 32768;
const HASH_BITS: u32 = 15;
/// Bytes hashed per probe position (4, 5 or 6). Longer hashes waste fewer
/// probes on the 4-byte collisions VCF text is full of, at the cost of no
/// longer discovering most 4-byte matches.
const HASH_BYTES: u32 = 6;
/// Interior-of-match insertion steps (short match, long match).
const STEP_SHORT: usize = 2;
const STEP_LONG: usize = 4;
const HASH_SIZE: usize = 1 << HASH_BITS;

/// Length code metadata: (symbol, extra_bits, base_len) per RFC 1951.
const LEN_SYM: [(u32, u32, u32); 29] = [
    (257, 0, 3), (258, 0, 4), (259, 0, 5), (260, 0, 6), (261, 0, 7), (262, 0, 8),
    (263, 0, 9), (264, 0, 10), (265, 1, 11), (266, 1, 13), (267, 1, 15), (268, 1, 17),
    (269, 2, 19), (270, 2, 23), (271, 2, 27), (272, 2, 31), (273, 3, 35), (274, 3, 43),
    (275, 3, 51), (276, 3, 59), (277, 4, 67), (278, 4, 83), (279, 4, 99), (280, 4, 115),
    (281, 5, 131), (282, 5, 163), (283, 5, 195), (284, 5, 227), (285, 0, 258),
];

/// Distance code metadata: (extra_bits, base_dist) per RFC 1951.
const DIST_SYM: [(u32, u32); 30] = [
    (0, 1), (0, 2), (0, 3), (0, 4), (1, 5), (1, 7), (2, 9), (2, 13), (3, 17), (3, 25),
    (4, 33), (4, 49), (5, 65), (5, 97), (6, 129), (6, 193), (7, 257), (7, 385),
    (8, 513), (8, 769), (9, 1025), (9, 1537), (10, 2049), (10, 3073), (11, 4097),
    (11, 6145), (12, 8193), (12, 12289), (13, 16385), (13, 24577),
];

/// Emit-ready code tables trained on one sample block and shared (behind an
/// `Arc`) by every compression job for the rest of the stream.
pub struct CodeSet {
    /// Literal byte -> code | nbits<<16 (code pre-reversed for LSB-first).
    lit: [u32; 256],
    /// Match length (3..=258, index len-3) -> fully combined bits and count:
    /// Huffman code plus extra bits in one value.
    len_bits: [u32; 256],
    len_n: [u8; 256],
    /// Distance symbol -> (revcode, code_len); extra bits appended at emit.
    dist_code: [(u32, u32); 30],
    /// d-1 -> distance symbol, zlib-style two-scale lookup.
    dist_lut: [u8; 512],
    eob: (u32, u32),
    header_words: Vec<u64>,
    header_nbits: u32,
}

/// Per-thread match-finder state, reused across blocks.
///
/// Each bucket packs the two most recent positions for its hash into one u64
/// (newest in the low half), so a probe is one load and an insert is one
/// shifting store — depth-2 matching at close to depth-1 cost.
pub struct Matcher {
    table: Vec<u64>,
    /// Virtual offset of the current block; bumped past MAX_DIST between
    /// blocks so stale table entries can never pass the distance check.
    base: u32,
}

impl Default for Matcher {
    fn default() -> Self {
        Matcher { table: vec![0; HASH_SIZE], base: 1 }
    }
}

#[inline(always)]
fn hash4(v: u32) -> usize {
    (v.wrapping_mul(0x9E37_79B1) >> (32 - HASH_BITS)) as usize
}

/// Hash of HASH_BYTES bytes at `i`; callers guarantee i + 8 <= data.len().
#[inline(always)]
fn hash_at(data: &[u8], i: usize) -> usize {
    if HASH_BYTES == 4 {
        hash4(load4(data, i))
    } else {
        let v = load8(data, i) << ((8 - HASH_BYTES as u64) * 8);
        (v.wrapping_mul(0x9E37_79B9_7F4A_7C15) >> (64 - HASH_BITS as u64)) as usize
    }
}

#[inline(always)]
fn load4(data: &[u8], i: usize) -> u32 {
    debug_assert!(i + 4 <= data.len());
    // SAFETY: every caller guarantees i + 4 <= data.len(); the parse loop
    // bounds i by end_safe = len - 12 and match candidates lie behind i.
    unsafe { (data.as_ptr().add(i) as *const u32).read_unaligned() }.to_le()
}

#[inline(always)]
fn load8(data: &[u8], i: usize) -> u64 {
    debug_assert!(i + 8 <= data.len());
    // SAFETY: match_len only calls this while n + 8 <= max <= len - position.
    unsafe { (data.as_ptr().add(i) as *const u64).read_unaligned() }.to_le()
}

/// Longest common prefix of data[a..] and data[b..], capped at MAX_MATCH,
/// compared eight bytes at a time. The caller has already proven the first
/// four bytes equal, so comparison starts there.
#[inline]
fn match_len(data: &[u8], a: usize, b: usize, limit: usize) -> usize {
    let max = limit.min(MAX_MATCH);
    let mut n = 4;
    while b + n + 8 <= data.len() && n + 8 <= max {
        let x = load8(data, a + n);
        let y = load8(data, b + n);
        let d = x ^ y;
        if d != 0 {
            return (n + (d.trailing_zeros() / 8) as usize).min(max);
        }
        n += 8;
    }
    while n < max && data[a + n] == data[b + n] {
        n += 1;
    }
    n
}

impl CodeSet {
    /// Build a code set from the symbol statistics of `sample`, smoothing so
    /// every symbol that could appear in later blocks has a code.
    pub fn train(sample: &[u8], matcher: &mut Matcher) -> CodeSet {
        let mut lit_freq = [1u64; 286]; // +1 smoothing on every symbol
        let mut dist_freq = [1u64; 30];
        parse::<false>(sample, matcher, None, &mut lit_freq, &mut dist_freq, &mut Vec::new());
        lit_freq[256] += 16; // EOB fires once per block
        Self::from_freqs(&lit_freq, &dist_freq)
    }

    /// Build a code set for exactly these frequencies — no smoothing, so
    /// symbols that never occur get no code. Used per block by the
    /// high-effort encoder, where the counts are this block's own.
    pub fn from_freqs(lit_freq: &[u64; 286], dist_freq: &[u64; 30]) -> CodeSet {
        let lit_lens = huffman::build_lengths(lit_freq, 15);
        let dist_lens = huffman::build_lengths(dist_freq, 15);
        let lit_codes = huffman::assign_codes(&lit_lens);
        let dist_codes = huffman::assign_codes(&dist_lens);
        let (header_words, header_nbits) = huffman::build_header(&lit_lens, &dist_lens);

        let mut cs = CodeSet {
            lit: [0; 256],
            len_bits: [0; 256],
            len_n: [0; 256],
            dist_code: [(0, 0); 30],
            dist_lut: [0; 512],
            eob: (lit_codes[256], lit_lens[256]),
            header_words,
            header_nbits,
        };
        for b in 0..256 {
            cs.lit[b] = lit_codes[b] | (lit_lens[b] << 16);
        }
        for len in 3usize..=258 {
            let s = LEN_SYM.iter().rposition(|&(_, _, base)| base as usize <= len).unwrap();
            let (sym, ebits, base) = LEN_SYM[s];
            let code = lit_codes[sym as usize];
            let n = lit_lens[sym as usize];
            cs.len_bits[len - 3] = code | (((len as u32 - base) as u32) << n);
            cs.len_n[len - 3] = (n + ebits) as u8;
        }
        for s in 0..30 {
            cs.dist_code[s] = (dist_codes[s], dist_lens[s]);
        }
        for dm1 in 0..512usize {
            let d = if dm1 < 256 { dm1 + 1 } else { ((dm1 - 256) << 7) + 1 };
            let s = DIST_SYM.iter().rposition(|&(_, base)| base as usize <= d).unwrap();
            cs.dist_lut[dm1] = s as u8;
        }
        cs
    }

    #[inline(always)]
    fn dist_sym(&self, d: usize) -> usize {
        let dm1 = d - 1;
        self.dist_lut[if dm1 < 256 { dm1 } else { 256 + (dm1 >> 7) }] as usize
    }
}

/// One shared parse for training (EMIT=false: count) and compression
/// (EMIT=true: write bits).
fn parse<const EMIT: bool>(
    data: &[u8],
    m: &mut Matcher,
    codes: Option<&CodeSet>,
    lit_freq: &mut [u64; 286],
    dist_freq: &mut [u64; 30],
    out: &mut Vec<u8>,
) -> (u64, u32) {
    // Refresh the virtual base; wrap the tables long before u32 overflow.
    if m.base > u32::MAX - 3 * (MAX_DIST as u32 + 0x1_0000) {
        m.table.iter_mut().for_each(|x| *x = 0);
        m.base = 1;
    }
    let base = m.base as usize;
    m.base += data.len() as u32 + MAX_DIST as u32 + 1;

    let mut bb: u64 = 0;
    let mut nb: u32 = 0;
    // Flush BEFORE adding: pending bits stay <= 31, and every put is <= 32
    // bits, so the OR below can never shift data off the top of the buffer.
    macro_rules! put {
        ($bits:expr, $n:expr) => {
            if EMIT {
                if nb >= 32 {
                    out.extend_from_slice(&(bb as u32).to_le_bytes());
                    bb >>= 32;
                    nb -= 32;
                }
                bb |= ($bits as u64) << nb;
                nb += $n as u32;
            }
        };
    }

    if EMIT {
        let cs = codes.unwrap();
        // BFINAL=1, BTYPE=10 (dynamic), then the cached header.
        put!(0b101u64, 3);
        let mut left = cs.header_nbits;
        for &w in &cs.header_words {
            let mut word = w;
            let mut take = left.min(64);
            while take > 0 {
                let k = take.min(32);
                put!(word & ((1u64 << k) - 1), k);
                word >>= k;
                take -= k;
            }
            left = left.saturating_sub(64);
        }
    }

    let n = data.len();
    let mut i = 0usize;
    let end_safe = n.saturating_sub(MIN_MATCH + 8); // room for 4-byte load + 8-byte compare
    // Grows through unmatched stretches so probing gets sparser (the skipped
    // positions are still emitted as literals); resets on every match.
    let mut miss_run = 0usize;
    while i < end_safe {
        let cur = load4(data, i);
        let h = hash_at(data, i);
        let bucket = m.table[h];
        m.table[h] = (bucket << 32) | (base + i) as u64;

        let mut best_len = 0usize;
        let mut best_dist = 0usize;
        for cand in [(bucket & 0xffff_ffff) as usize, (bucket >> 32) as usize] {
            // cand == 0 fails cand >= base, since base is always >= 1.
            if cand >= base {
                let dist = base + i - cand;
                if dist <= MAX_DIST {
                    let j = cand - base;
                    if load4(data, j) == cur {
                        let l = match_len(data, j, i, n - i);
                        if l > best_len {
                            best_len = l;
                            best_dist = dist;
                            // A long match will not be beaten by the older,
                            // farther candidate often enough to pay for the probe.
                            if l >= 32 {
                                break;
                            }
                        }
                    }
                }
            }
        }

        if best_len >= MIN_MATCH {
            if EMIT {
                let cs = codes.unwrap();
                put!(cs.len_bits[best_len - 3], cs.len_n[best_len - 3]);
                let s = cs.dist_sym(best_dist);
                let (code, clen) = cs.dist_code[s];
                let (ebits, dbase) = DIST_SYM[s];
                put!(code as u64 | (((best_dist as u32 - dbase) as u64) << clen), clen + ebits);
            } else {
                let s = LEN_SYM.iter().rposition(|&(_, _, b)| b as usize <= best_len).unwrap();
                lit_freq[LEN_SYM[s].0 as usize] += 1;
                let ds = DIST_SYM.iter().rposition(|&(_, b)| b as usize <= best_dist).unwrap();
                dist_freq[ds] += 1;
            }
            // Index the interior of the match sparsely: dense enough to keep
            // finding the repeats VCF is full of, cheap enough not to pay a
            // hash per byte inside long runs.
            let stop = (i + best_len).min(end_safe);
            let step = if best_len < 16 { STEP_SHORT } else { STEP_LONG };
            let mut k = i + 1;
            while k < stop {
                let hh = hash_at(data, k);
                m.table[hh] = (m.table[hh] << 32) | (base + k) as u64;
                k += step;
            }
            i += best_len;
            miss_run = 0;
        } else {
            if EMIT {
                let cs = codes.unwrap();
                let e = cs.lit[data[i] as usize];
                put!(e & 0xffff, e >> 16);
            } else {
                lit_freq[data[i] as usize] += 1;
            }
            i += 1;
            miss_run += 1;
            // Deep in a literal run, emit a few literals without probing at all.
            let skip = miss_run >> 5;
            if skip > 0 {
                let stop = (i + skip).min(end_safe);
                while i < stop {
                    if EMIT {
                        let cs = codes.unwrap();
                        let e = cs.lit[data[i] as usize];
                        put!(e & 0xffff, e >> 16);
                    } else {
                        lit_freq[data[i] as usize] += 1;
                    }
                    i += 1;
                }
            }
        }
    }
    while i < n {
        if EMIT {
            let cs = codes.unwrap();
            let e = cs.lit[data[i] as usize];
            put!(e & 0xffff, e >> 16);
        } else {
            lit_freq[data[i] as usize] += 1;
        }
        i += 1;
    }
    if EMIT {
        let cs = codes.unwrap();
        put!(cs.eob.0, cs.eob.1);
        // Final byte-align flush.
        while nb > 0 {
            out.push(bb as u8);
            bb >>= 8;
            nb = nb.saturating_sub(8);
        }
    }
    (bb, nb)
}

/// Compress `data` as one complete raw DEFLATE stream into `out` (appended).
/// Falls back to stored blocks if the compressed form would be larger.
pub fn compress_into(codes: &CodeSet, m: &mut Matcher, data: &[u8], out: &mut Vec<u8>) {
    let start = out.len();
    let mut dummy_l = [0u64; 286];
    let mut dummy_d = [0u64; 30];
    parse::<true>(data, m, Some(codes), &mut dummy_l, &mut dummy_d, out);
    if out.len() - start > data.len() + 5 * (data.len() / 65535 + 1) {
        // Incompressible: rewrite as stored blocks.
        out.truncate(start);
        store_uncompressed(data, out);
    }
}

/// Emit `data` as DEFLATE stored (BTYPE=00) blocks — the incompressible-input
/// fallback shared by both encoders.
fn store_uncompressed(data: &[u8], out: &mut Vec<u8>) {
    let mut off = 0;
    while off < data.len() || data.is_empty() {
        let chunk = (data.len() - off).min(65535);
        let last = off + chunk == data.len();
        out.push(if last { 1 } else { 0 });
        out.extend_from_slice(&(chunk as u16).to_le_bytes());
        out.extend_from_slice(&(!(chunk as u16)).to_le_bytes());
        out.extend_from_slice(&data[off..off + chunk]);
        off += chunk;
        if data.is_empty() {
            break;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn vcfish(n: usize, seed: u64) -> Vec<u8> {
        let mut x = seed;
        let mut v = Vec::with_capacity(n);
        while v.len() < n {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            v.extend_from_slice(
                format!(
                    "chr1\t{}\trs{}\tA\tG\t{}\tPASS\tAC=1;AN=200\tGT:AD:DP\t0/1:{},{}:{}\n",
                    x % 1_000_000, x % 99999, x % 900, x % 30, x % 30, x % 60
                )
                .as_bytes(),
            );
        }
        v.truncate(n);
        v
    }

    fn roundtrip(data: &[u8]) -> Vec<u8> {
        let mut m = Matcher::default();
        let codes = CodeSet::train(&data[..data.len().min(0xff00)], &mut m);
        let mut out = Vec::new();
        compress_into(&codes, &mut m, data, &mut out);
        let mut d = libdeflater::Decompressor::new();
        let mut back = vec![0u8; data.len()];
        let n = d
            .deflate_decompress(&out, &mut back)
            .expect("libdeflate must accept our stream");
        assert_eq!(n, data.len());
        assert_eq!(&back, data, "round trip mismatch");
        out
    }

    #[test]
    fn roundtrips_vcf_like_text() {
        let data = vcfish(0xff00, 42);
        let out = roundtrip(&data);
        assert!(out.len() < data.len() / 3, "should compress well: {}", out.len());
    }

    #[test]
    fn roundtrips_blocks_the_codes_were_not_trained_on() {
        let mut m = Matcher::default();
        let train = vcfish(0xff00, 1);
        let codes = CodeSet::train(&train, &mut m);
        for seed in 2..12u64 {
            let data = vcfish(0xff00, seed * 7919);
            let mut out = Vec::new();
            compress_into(&codes, &mut m, &data, &mut out);
            let mut d = libdeflater::Decompressor::new();
            let mut back = vec![0u8; data.len()];
            let n = d.deflate_decompress(&out, &mut back).unwrap();
            assert_eq!(n, data.len());
            assert_eq!(back, data);
        }
    }

    #[test]
    fn handles_awkward_inputs() {
        for data in [
            vec![],
            b"x".to_vec(),
            b"abcd".to_vec(),
            vec![b'A'; 70000],
            (0..=255u8).cycle().take(3).collect::<Vec<_>>(),
            vcfish(17, 9),
            vcfish(0xff00, 3),
        ] {
            roundtrip(&data);
        }
    }

    #[test]
    fn incompressible_input_falls_back_to_stored() {
        let mut x = 0x2545F4914F6CDD1Du64;
        let data: Vec<u8> = (0..0xff00)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                x as u8
            })
            .collect();
        let out = roundtrip(&data);
        assert!(out.len() <= data.len() + 20);
    }

    #[test]
    fn bytes_absent_from_training_still_encode() {
        let mut m = Matcher::default();
        let codes = CodeSet::train(b"aaaaaaaabbbbbbbbcccccccc", &mut m);
        let data: Vec<u8> = (0..=255u8).cycle().take(1000).collect();
        let mut out = Vec::new();
        compress_into(&codes, &mut m, &data, &mut out);
        let mut d = libdeflater::Decompressor::new();
        let mut back = vec![0u8; data.len()];
        let n = d.deflate_decompress(&out, &mut back).unwrap();
        assert_eq!(n, data.len());
        assert_eq!(back, data);
    }

    #[test]
    fn distances_never_exceed_the_deflate_window() {
        // A block longer than 32768 with a repeat straddling more than the
        // window: the encoder must re-match locally, not reach back too far.
        let unit = vcfish(40000, 5);
        let mut data = unit.clone();
        data.extend_from_slice(&unit);
        roundtrip(&data);
    }
}
