//! Honest micro-benchmark: the built-in Rust encoder vs libdeflate, on the
//! actual blocks of a VCF. Reads a plain VCF from stdin, splits it into
//! BGZF-sized chunks, and measures each encoder end to end.
use std::io::Read;
use std::time::Instant;

#[path = "../src/deflate/mod.rs"]
mod deflate;

fn main() {
    let mut data = Vec::new();
    std::io::stdin().read_to_end(&mut data).unwrap();
    let blocks: Vec<&[u8]> = data.chunks(0xff00).collect();
    let mb = data.len() as f64 / 1e6;
    println!("{} blocks, {:.0} MB", blocks.len(), mb);

    for level in [1i32, 2, 3, 6] {
        let t = Instant::now();
        let mut total = 0usize;
        let mut c = libdeflater::Compressor::new(libdeflater::CompressionLvl::new(level).unwrap());
        let bound = c.deflate_compress_bound(0xff00);
        let mut out = vec![0u8; bound];
        for b in &blocks {
            total += c.deflate_compress(b, &mut out).unwrap();
        }
        let dt = t.elapsed().as_secs_f64();
        println!(
            "libdeflate l{level}: {:>7.1} MB/s   ratio {:.3}   {} bytes",
            mb / dt,
            data.len() as f64 / total as f64,
            total
        );
    }

    let mut m = deflate::Matcher::default();
    // Train on a mid-file block: the first blocks are the VCF header.
    let train_block = blocks[blocks.len() / 2];
    let codes = deflate::CodeSet::train(train_block, &mut m);
    let t = Instant::now();
    let mut total = 0usize;
    let mut out = Vec::with_capacity(0xff00 + 1024);
    for b in &blocks {
        out.clear();
        deflate::compress_into(&codes, &mut m, b, &mut out);
        total += out.len();
    }
    let dt = t.elapsed().as_secs_f64();
    println!(
        "vcfr-rs      : {:>7.1} MB/s   ratio {:.3}   {} bytes",
        mb / dt,
        data.len() as f64 / total as f64,
        total
    );

    // Same encode plus the CRC each real BGZF member needs, isolating the
    // checksum cost that the pure-encode loop above leaves out.
    let t = Instant::now();
    let mut total = 0usize;
    let mut crc_acc = 0u32;
    for b in &blocks {
        out.clear();
        deflate::compress_into(&codes, &mut m, b, &mut out);
        crc_acc ^= crc32fast::hash(b);
        total += out.len();
    }
    let dt = t.elapsed().as_secs_f64();
    println!(
        "vcfr-rs + crc32fast: {:>7.1} MB/s   ({} bytes, crc xor {crc_acc:08x})",
        mb / dt,
        total
    );

    let t = Instant::now();
    let mut crc_acc2 = 0u32;
    for b in &blocks {
        crc_acc2 ^= libdeflater::crc32(b);
    }
    let dt = t.elapsed().as_secs_f64();
    println!("libdeflate crc32 alone: {:>7.1} MB/s (xor {crc_acc2:08x})", mb / dt);
    let t = Instant::now();
    let mut crc_acc3 = 0u32;
    for b in &blocks {
        crc_acc3 ^= crc32fast::hash(b);
    }
    let dt = t.elapsed().as_secs_f64();
    println!("crc32fast alone       : {:>7.1} MB/s (xor {crc_acc3:08x})", mb / dt);
}
