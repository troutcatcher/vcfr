//! High-effort DEFLATE path: per-block optimal Huffman codes on top of a
//! hash-chain lazy matcher.
//!
//! Where the fast path (`compress_into`) spends nothing on statistics —
//! preset codes, single greedy pass — this one pays for ratio three ways:
//! it tokenizes the whole block first and builds exact Huffman codes for
//! that block's own symbol counts, it walks bounded hash chains instead of
//! a depth-2 bucket, and it matches lazily (a match is deferred when the
//! next position hides a longer one). The output is still one ordinary
//! dynamic-Huffman DEFLATE stream per block.


use super::{load4, match_len, CodeSet, DIST_SYM, LEN_SYM, MAX_DIST};

const MIN_MATCH: usize = 4;
const HASH_BITS_H: u32 = 16;
const HASH_SIZE_H: usize = 1 << HASH_BITS_H;
/// Chain probes per position; the knob that trades speed for ratio.
/// (Runtime fields on HighMatcher during tuning; defaults below.)
const MAX_CHAIN: usize = 64;
/// A match this long is taken immediately — no lazy deferral, no more probes.
const NICE_LEN: usize = 130;
/// Matches at least this long skip the lazy look-ahead (zlib's max_lazy idea).
const MAX_LAZY: usize = 32;

#[inline(always)]
fn hash4_h(v: u32) -> usize {
    (v.wrapping_mul(0x9E37_79B1) >> (32 - HASH_BITS_H)) as usize
}

/// Match token: bit 31 set, length-3 in bits 16..24, dist-1 in bits 0..15.
/// Literals are the bare byte value.
const TOK_MATCH: u32 = 1 << 31;

/// Per-thread state for the high-effort encoder, reused across blocks.
pub struct HighMatcher {
    /// hash -> most recent position, stored as base + pos (0 = empty; a
    /// stale entry from an earlier block fails the `>= base` check).
    head: Vec<u32>,
    /// Block-local chain: prev[pos] -> previous base+pos with the same hash.
    prev: Vec<u32>,
    base: u32,
    tokens: Vec<u32>,
    max_chain: usize,
    nice_len: usize,
    max_lazy: usize,
}

impl Default for HighMatcher {
    fn default() -> Self {
        HighMatcher {
            head: vec![0; HASH_SIZE_H],
            prev: vec![0; 0x1_0000],
            base: 1,
            tokens: Vec::with_capacity(0x1_0000),
            max_chain: MAX_CHAIN,
            nice_len: NICE_LEN,
            max_lazy: MAX_LAZY,
        }
    }
}

/// Effort settings for `--codec rust` levels 1-6, measured on VCF text to
/// track libdeflate's speed/ratio frontier from roughly l3 up to l9:
/// (max_chain, nice_len, max_lazy).
pub fn effort_for_level(level: u32) -> (usize, usize, usize) {
    match level.clamp(1, 6) {
        1 => (8, 32, 16),
        2 => (16, 64, 16),
        3 => (24, 130, 24),
        4 => (64, 130, 32),
        5 => (128, 258, 64),
        _ => (256, 258, 128),
    }
}

impl HighMatcher {
    pub fn set_effort(&mut self, chain: usize, nice: usize, lazy: usize) {
        self.max_chain = chain;
        self.nice_len = nice;
        self.max_lazy = lazy;
    }

    /// Longest match at `i`, walking the chain from `prev[i]` — the caller
    /// guarantees `i` itself was the most recent insertion for its hash, so
    /// starting past it skips the zero-distance self-match. Returns
    /// (len, dist), len 0 if nothing beat `min_beat.max(MIN_MATCH - 1)`.
    #[inline]
    fn find(&self, data: &[u8], i: usize, base: usize, min_beat: usize) -> (usize, usize) {
        let n = data.len();
        let cur = load4(data, i);
        let mut best_len = min_beat.max(MIN_MATCH - 1);
        let mut best_dist = 0usize;
        let mut cand = self.prev[i & 0xffff] as usize;
        let abs_i = base + i;
        for _ in 0..self.max_chain {
            if cand < base {
                break;
            }
            let dist = abs_i - cand;
            if dist > MAX_DIST {
                break;
            }
            let j = cand - base;
            // Cheap rejection: a chain entry that can't beat best_len will
            // differ at data[j + best_len] almost always.
            if j + best_len < n
                && i + best_len < n
                && data[j + best_len] == data[i + best_len]
                && load4(data, j) == cur
            {
                let l = match_len(data, j, i, n - i);
                if l > best_len {
                    best_len = l;
                    best_dist = dist;
                    if l >= self.nice_len {
                        break;
                    }
                }
            }
            cand = self.prev[j & 0xffff] as usize;
        }
        if best_dist == 0 {
            (0, 0)
        } else {
            (best_len, best_dist)
        }
    }

    /// Insert every position in `[*inserted, upto)`, exactly once each — a
    /// double insert would make a position its own chain predecessor and the
    /// walk would spin on it. (A batched AVX2 version of this — pshufb window
    /// extraction, vectorized multiply — measured 25% slower than this scalar
    /// loop: the serial head/prev stores dominate, and the common one-position
    /// call paid a non-inlinable target_feature function call for nothing.)
    #[inline(always)]
    fn insert_upto(&mut self, data: &[u8], base: usize, upto: usize, inserted: &mut usize) {
        while *inserted < upto {
            let i = *inserted;
            let h = hash4_h(load4(data, i));
            self.prev[i & 0xffff] = self.head[h];
            self.head[h] = (base + i) as u32;
            *inserted += 1;
        }
    }
}

/// Static length-3-indexed and distance LUTs for frequency counting.
struct SymLut {
    len_sym: [u8; 256],  // len-3 -> LEN_SYM index
    dist_lut: [u8; 512], // same two-scale scheme as CodeSet::dist_sym
}

fn build_sym_lut() -> SymLut {
    let mut lut = SymLut { len_sym: [0; 256], dist_lut: [0; 512] };
    for len in 3usize..=258 {
        lut.len_sym[len - 3] =
            LEN_SYM.iter().rposition(|&(_, _, base)| base as usize <= len).unwrap() as u8;
    }
    for dm1 in 0..512usize {
        let d = if dm1 < 256 { dm1 + 1 } else { ((dm1 - 256) << 7) + 1 };
        lut.dist_lut[dm1] = DIST_SYM.iter().rposition(|&(_, base)| base as usize <= d).unwrap() as u8;
    }
    lut
}

#[inline(always)]
fn dist_sym_of(lut: &SymLut, d: usize) -> usize {
    let dm1 = d - 1;
    lut.dist_lut[if dm1 < 256 { dm1 } else { 256 + (dm1 >> 7) }] as usize
}

/// Compress `data` as one complete raw DEFLATE stream (appended to `out`),
/// with Huffman codes computed for this block's exact token statistics.
/// Falls back to stored blocks when that would be smaller.
pub fn compress_high_into(m: &mut HighMatcher, data: &[u8], out: &mut Vec<u8>) {
    // ---- pass 1: tokenize with lazy chain matching, counting as we go ----
    if m.base > u32::MAX - 3 * (MAX_DIST as u32 + 0x1_0000) {
        m.head.iter_mut().for_each(|x| *x = 0);
        m.base = 1;
    }
    let base = m.base as usize;
    m.base += data.len() as u32 + MAX_DIST as u32 + 1;

    let lut = build_sym_lut();
    let mut lit_freq = [0u64; 286];
    let mut dist_freq = [0u64; 30];
    let mut tokens = std::mem::take(&mut m.tokens);
    tokens.clear();

    let n = data.len();
    let end_safe = n.saturating_sub(MIN_MATCH + 8);
    let mut i = 0usize;
    let mut inserted = 0usize; // positions below this are in the table
    while i < end_safe {
        m.insert_upto(data, base, i + 1, &mut inserted);
        let (mut len, mut dist) = m.find(data, i, base, 0);
        if len < MIN_MATCH {
            tokens.push(data[i] as u32);
            lit_freq[data[i] as usize] += 1;
            i += 1;
            continue;
        }
        // Lazy: while the match is short, see if i+1 starts a longer one.
        while len < m.max_lazy && i + 1 < end_safe {
            m.insert_upto(data, base, i + 2, &mut inserted);
            let (l2, d2) = m.find(data, i + 1, base, len);
            if l2 > len {
                tokens.push(data[i] as u32);
                lit_freq[data[i] as usize] += 1;
                i += 1;
                len = l2;
                dist = d2;
            } else {
                break;
            }
        }
        tokens.push(TOK_MATCH | (((len - 3) as u32) << 16) | (dist - 1) as u32);
        lit_freq[LEN_SYM[lut.len_sym[len - 3] as usize].0 as usize] += 1;
        dist_freq[dist_sym_of(&lut, dist)] += 1;
        // Insert every interior position: full coverage is what lets the
        // chains keep finding the next repeat.
        m.insert_upto(data, base, (i + len).min(end_safe), &mut inserted);
        i += len;
    }
    while i < n {
        tokens.push(data[i] as u32);
        lit_freq[data[i] as usize] += 1;
        i += 1;
    }
    lit_freq[256] += 1; // EOB

    // ---- pass 2: exact codes for these counts, then replay the tokens ----
    let cs = CodeSet::from_freqs(&lit_freq, &dist_freq);
    let start = out.len();
    emit_tokens(&cs, &tokens, out);

    if out.len() - start > data.len() + 5 * (data.len() / 65535 + 1) {
        out.truncate(start);
        super::store_uncompressed(data, out);
    }
    m.tokens = tokens;
}

fn emit_tokens(cs: &CodeSet, tokens: &[u32], out: &mut Vec<u8>) {
    let mut bb: u64 = 0;
    let mut nb: u32 = 0;
    // Same discipline as the fast path: flush BEFORE the OR, so pending bits
    // stay <= 31 and every put of <= 32 bits fits.
    macro_rules! put {
        ($bits:expr, $n:expr) => {
            if nb >= 32 {
                out.extend_from_slice(&(bb as u32).to_le_bytes());
                bb >>= 32;
                nb -= 32;
            }
            bb |= ($bits as u64) << nb;
            nb += $n as u32;
        };
    }

    put!(0b101u64, 3); // BFINAL=1, BTYPE=10
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

    for &t in tokens {
        if t & TOK_MATCH == 0 {
            let e = cs.lit[t as usize];
            put!(e & 0xffff, e >> 16);
        } else {
            let lm3 = ((t >> 16) & 0xff) as usize;
            let dist = (t & 0xffff) as usize + 1;
            put!(cs.len_bits[lm3], cs.len_n[lm3]);
            let s = cs.dist_sym(dist);
            let (code, clen) = cs.dist_code[s];
            let (ebits, dbase) = DIST_SYM[s];
            put!(code as u64 | (((dist as u32 - dbase) as u64) << clen), clen + ebits);
        }
    }
    put!(cs.eob.0, cs.eob.1);
    while nb > 0 {
        out.push(bb as u8);
        bb >>= 8;
        nb = nb.saturating_sub(8);
    }
}






#[cfg(test)]
mod tests {
    use super::*;

    fn roundtrip_high(data: &[u8]) -> Vec<u8> {
        let mut m = HighMatcher::default();
        let mut out = Vec::new();
        compress_high_into(&mut m, data, &mut out);
        let mut d = libdeflater::Decompressor::new();
        let mut back = vec![0u8; data.len().max(1)];
        let n = d.deflate_decompress(&out, &mut back).expect("must inflate");
        assert_eq!(n, data.len());
        assert_eq!(&back[..n], data);
        out
    }

    #[test]
    fn high_roundtrips_awkward_inputs() {
        for data in [
            vec![],
            b"x".to_vec(),
            b"abcd".to_vec(),
            vec![b'A'; 70000],
            (0..=255u8).cycle().take(3).collect::<Vec<_>>(),
            (0..=255u8).cycle().take(1000).collect::<Vec<_>>(),
        ] {
            roundtrip_high(&data);
        }
    }

    #[test]
    fn high_roundtrips_vcf_like_blocks_and_reuses_state() {
        let mut m = HighMatcher::default();
        let mut x = 42u64;
        for _ in 0..6 {
            let mut data = Vec::new();
            while data.len() < 0xff00 {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                data.extend_from_slice(
                    format!(
                        "chr2\t{}\t.\tC\tT\t.\tPASS\tDR2=0.9{};AF=0.0{}\tGT:DS\t0|0:0.0{}\n",
                        x % 5_000_000, x % 10, x % 1000, x % 10
                    )
                    .as_bytes(),
                );
            }
            data.truncate(0xff00);
            let mut out = Vec::new();
            compress_high_into(&mut m, &data, &mut out);
            let mut d = libdeflater::Decompressor::new();
            let mut back = vec![0u8; data.len()];
            let n = d.deflate_decompress(&out, &mut back).unwrap();
            assert_eq!(n, data.len());
            assert_eq!(back, data);
            assert!(out.len() < data.len() / 3);
        }
    }

    #[test]
    fn high_incompressible_falls_back_to_stored() {
        let mut x = 0x2545F4914F6CDD1Du64;
        let data: Vec<u8> = (0..0xff00)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                x as u8
            })
            .collect();
        let out = roundtrip_high(&data);
        assert!(out.len() <= data.len() + 20);
    }
}
