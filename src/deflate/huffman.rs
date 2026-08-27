//! Canonical Huffman code construction and the RFC 1951 dynamic-block header.
//!
//! Everything here runs once per trained code set, not per block, so it is
//! written for clarity and correctness rather than speed.

/// Build length-limited canonical Huffman code lengths for `freqs`.
///
/// Standard two-queue Huffman followed by a Kraft repair pass when any code
/// exceeds `max_len`. The repair is the classic one: cap the offenders, then
/// lengthen the cheapest short codes until the Kraft sum fits, then shorten
/// greedily while it still fits.
pub fn build_lengths(freqs: &[u64], max_len: u32) -> Vec<u32> {
    let n = freqs.len();
    let mut lens = vec![0u32; n];
    let used: Vec<usize> = (0..n).filter(|&i| freqs[i] > 0).collect();
    match used.len() {
        0 => return lens,
        1 => {
            lens[used[0]] = 1;
            return lens;
        }
        _ => {}
    }

    // Huffman via sorted leaf queue + package queue.
    let mut leaves: Vec<(u64, usize)> = used.iter().map(|&i| (freqs[i], i)).collect();
    leaves.sort();
    // Tree nodes: (weight, left, right); leaves are represented by symbol index.
    #[derive(Clone, Copy)]
    enum Node {
        Leaf(usize),
        Internal(usize, usize),
    }
    let mut nodes: Vec<(u64, Node)> = Vec::with_capacity(2 * leaves.len());
    let mut li = 0usize;
    let mut packages: std::collections::VecDeque<usize> = std::collections::VecDeque::new();
    let take = |nodes: &mut Vec<(u64, Node)>,
                    li: &mut usize,
                    packages: &mut std::collections::VecDeque<usize>|
     -> usize {
        let leaf_w = leaves.get(*li).map(|l| l.0);
        let pack_w = packages.front().map(|&p| nodes[p].0);
        match (leaf_w, pack_w) {
            (Some(lw), Some(pw)) if lw <= pw => {
                let (w, s) = leaves[*li];
                *li += 1;
                nodes.push((w, Node::Leaf(s)));
                nodes.len() - 1
            }
            (Some(_), Some(_)) | (None, Some(_)) => packages.pop_front().unwrap(),
            (Some(_), None) => {
                let (w, s) = leaves[*li];
                *li += 1;
                nodes.push((w, Node::Leaf(s)));
                nodes.len() - 1
            }
            (None, None) => unreachable!(),
        }
    };
    let total = leaves.len();
    let mut remaining = total;
    let mut root = 0usize;
    while remaining > 1 {
        let a = take(&mut nodes, &mut li, &mut packages);
        let b = take(&mut nodes, &mut li, &mut packages);
        let w = nodes[a].0 + nodes[b].0;
        nodes.push((w, Node::Internal(a, b)));
        root = nodes.len() - 1;
        packages.push_back(root);
        remaining -= 1;
    }
    // Depth-first depth assignment.
    let mut stack = vec![(root, 0u32)];
    while let Some((idx, d)) = stack.pop() {
        match nodes[idx].1 {
            Node::Leaf(sym) => lens[sym] = d.max(1),
            Node::Internal(a, b) => {
                stack.push((a, d + 1));
                stack.push((b, d + 1));
            }
        }
    }

    if lens.iter().all(|&l| l <= max_len) {
        return lens;
    }

    // Kraft repair: cap, then fix the sum.
    let unit = 1u64 << max_len; // Kraft contributions scaled by 2^max_len
    for l in lens.iter_mut() {
        if *l > max_len {
            *l = max_len;
        }
    }
    let ksum = |lens: &Vec<u32>| -> u64 {
        lens.iter().filter(|&&l| l > 0).map(|&l| unit >> l).sum()
    };
    // Lengthen the least-frequent short codes until the code is feasible.
    let mut order = used.clone();
    order.sort_by_key(|&i| freqs[i]);
    let mut sum = ksum(&lens);
    'outer: while sum > unit {
        for &i in &order {
            if lens[i] < max_len && lens[i] > 0 {
                sum -= unit >> lens[i];
                lens[i] += 1;
                sum += unit >> lens[i];
                if sum <= unit {
                    break 'outer;
                }
            }
        }
    }
    // Tighten: shorten the most frequent codes while the sum allows it.
    let mut order_desc = order.clone();
    order_desc.reverse();
    let mut changed = true;
    while changed {
        changed = false;
        for &i in &order_desc {
            if lens[i] > 1 {
                let delta = (unit >> (lens[i] - 1)) - (unit >> lens[i]);
                if sum + delta <= unit {
                    lens[i] -= 1;
                    sum += delta;
                    changed = true;
                }
            }
        }
    }
    debug_assert!(sum <= unit);
    lens
}

/// Canonical code assignment (RFC 1951 order), returning bit-reversed codes
/// ready for LSB-first emission.
pub fn assign_codes(lens: &[u32]) -> Vec<u32> {
    let max = lens.iter().copied().max().unwrap_or(0);
    let mut bl_count = vec![0u32; (max + 1) as usize];
    for &l in lens {
        if l > 0 {
            bl_count[l as usize] += 1;
        }
    }
    let mut next = vec![0u32; (max + 2) as usize];
    let mut code = 0u32;
    for bits in 1..=max {
        code = (code + bl_count[(bits - 1) as usize]) << 1;
        next[bits as usize] = code;
    }
    lens.iter()
        .map(|&l| {
            if l == 0 {
                return 0;
            }
            let c = next[l as usize];
            next[l as usize] += 1;
            c.reverse_bits() >> (32 - l)
        })
        .collect()
}

/// Order in which code-length-code lengths appear in the header.
const CL_ORDER: [usize; 19] = [16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15];

/// Serialise the dynamic-block header (everything after BFINAL/BTYPE) for the
/// given litlen and distance code lengths. Returns (bits, total_bit_count),
/// packed LSB-first into u64 words for cheap replay at the start of every block.
pub fn build_header(lit_lens: &[u32], dist_lens: &[u32]) -> (Vec<u64>, u32) {
    // RLE-encode the two length arrays as one stream of CL symbols.
    let nlit = lit_lens.iter().rposition(|&l| l > 0).map_or(257, |p| (p + 1).max(257));
    let ndist = dist_lens.iter().rposition(|&l| l > 0).map_or(1, |p| (p + 1).max(1));
    let all: Vec<u32> = lit_lens[..nlit].iter().chain(dist_lens[..ndist].iter()).copied().collect();

    // (symbol, extra_value, extra_bits)
    let mut rle: Vec<(u32, u32, u32)> = Vec::new();
    let mut i = 0;
    while i < all.len() {
        let v = all[i];
        let mut run = 1;
        while i + run < all.len() && all[i + run] == v {
            run += 1;
        }
        if v == 0 {
            let mut left = run;
            while left >= 11 {
                let take = left.min(138);
                rle.push((18, (take - 11) as u32, 7));
                left -= take;
            }
            if left >= 3 {
                rle.push((17, (left - 3) as u32, 3));
                left = 0;
            }
            for _ in 0..left {
                rle.push((0, 0, 0));
            }
        } else {
            rle.push((v, 0, 0));
            let mut left = run - 1;
            while left >= 3 {
                let take = left.min(6);
                rle.push((16, (take - 3) as u32, 2));
                left -= take;
            }
            for _ in 0..left {
                rle.push((v, 0, 0));
            }
        }
        i += run;
    }

    let mut cl_freq = [0u64; 19];
    for &(s, _, _) in &rle {
        cl_freq[s as usize] += 1;
    }
    let cl_lens = build_lengths(&cl_freq, 7);
    let cl_codes = assign_codes(&cl_lens);
    let hclen = (4..19)
        .rposition(|k| cl_lens[CL_ORDER[k]] > 0)
        .map_or(4, |p| p + 4 + 1)
        .max(4);

    let mut words: Vec<u64> = Vec::new();
    let mut bb: u64 = 0;
    let mut nb: u32 = 0;
    let put = |bits: u64, n: u32, words: &mut Vec<u64>, bb: &mut u64, nb: &mut u32| {
        *bb |= bits << *nb;
        *nb += n;
        if *nb >= 64 {
            words.push(*bb);
            *nb -= 64;
            *bb = if n == *nb { 0 } else { bits >> (n - *nb) };
        }
    };
    put((nlit - 257) as u64, 5, &mut words, &mut bb, &mut nb);
    put((ndist - 1) as u64, 5, &mut words, &mut bb, &mut nb);
    put((hclen - 4) as u64, 4, &mut words, &mut bb, &mut nb);
    for k in 0..hclen {
        put(cl_lens[CL_ORDER[k]] as u64, 3, &mut words, &mut bb, &mut nb);
    }
    for &(s, extra, ebits) in &rle {
        put(cl_codes[s as usize] as u64, cl_lens[s as usize], &mut words, &mut bb, &mut nb);
        if ebits > 0 {
            put(extra as u64, ebits, &mut words, &mut bb, &mut nb);
        }
    }
    let total = words.len() as u32 * 64 + nb;
    if nb > 0 {
        words.push(bb);
    }
    (words, total)
}
