# vcfr

A fast VCF toolkit: **merge**, **concat**, and **subset** VCF files.

`vcfr` does the handful of VCF operations that dominate real pipelines, and
does them quickly. Output is verified against `bcftools` — the differential
suite in `tests/differential.sh` runs 58 paired checks and requires the two
tools to agree exactly.

```
vcfr view   [-s samples] [-r regions] [-v types] ...  file.vcf.gz
vcfr concat [--naive] part1.vcf.gz part2.vcf.gz ...
vcfr merge  cohortA.vcf.gz cohortB.vcf.gz ...
```

## Benchmarks

All numbers below come from **one session on one 4-core machine** (bcftools
1.19 / htslib 1.19, both tools `--threads 4`, best of 2 interleaved rounds,
page cache flushed before each timed run). "cores" is CPU-seconds over wall
seconds — the parallelism actually achieved — because a speedup on its own is
ambiguous: a tool can be faster because it is better, or because it was
handed more cores.

### Generic call set

381 MB BGZF file, 250,000 sites × 500 samples (~2.5 GB of VCF text), the
shape of an ordinary genotyped cohort.

| operation | bcftools | vcfr | speedup | cores b → v | output b / v |
| --- | ---: | ---: | ---: | :---: | :---: |
| `view`: decompress to VCF | 20.1s | 0.66s | **30x** | 1.2 → 3.9 | — |
| `view`: recompress to BGZF | 26.4s | 11.5s | 2.3x † | 3.5 → 4.0 | 398M / 430M |
| `view`: subset 50 of 500 samples | 13.7s | 2.5s | **5.5x** | 1.8 → 3.3 | 45M / 48M |
| `concat`: 4 parts → BGZF | 18.9s | 11.0s | 1.7x † | 3.8 → 4.0 | 398M / 430M |
| `concat --naive` (block copy) | 1.6s | 0.3s | I/O-bound ‡ | — | 430M / 430M |
| `merge`: 3 cohorts → VCF | 25.7s | 1.8s | **14x** | 1.1 → 3.2 | — |
| `merge`: 3 cohorts → BGZF | 27.2s | 10.1s | 2.7x † | 3.2 → 4.0 | 346M / 372M |
| `merge`: 3 cohorts → BGZF, `-l 7` | 27.2s | 15.7s | **1.7x** | 3.2 → 4.0 | 346M / 345M |

### Imputed sequence data (Beagle-style)

Three 300-sample cohorts × 150,000 sites of Beagle output shape — phased
`GT:DS:GP`, `DR2`/`AF`/`IMP` INFO, allele frequencies drawn from the
sequence spectrum (log-uniform on [2e-4, 0.5], ~50% of sites under 1% AF,
mostly hom-ref genotypes; `examples/gen_beagle --af-spectrum seq`). ~3.2 GB
of merged text. The gap is wider than on the generic set because rare-variant
text makes vcfr's copy path and the compressor cheaper per byte, while
bcftools still decodes and re-renders every DS/GP float.

| operation | bcftools | vcfr | speedup | cores b → v | output b / v |
| --- | ---: | ---: | ---: | :---: | :---: |
| `merge`: 3 cohorts → VCF | 37.3s | 1.6s | **23x** | 1.1 → 3.0 | — |
| `merge`: 3 cohorts → BGZF | 39.9s | 6.6s | **6.1x** | 2.2 → 3.9 | 145M / 160M |
| `merge`: 3 cohorts → BGZF, `-l 7` | 39.9s | 12.9s | **3.1x** | 2.2 → 4.0 | 145M / 149M |

### Wide cohorts and block splicing

Ten batches of 220,000 samples each (2.2M samples merged; lines are ~30-44 MB
of text spanning up to ~670 BGZF blocks). Here `merge -Oz` and wide `view -s`
skip most deflate work entirely by **block splicing**: an input block whose
content lies wholly inside a verbatim-copied stretch of a line is passed into
the output still compressed (see below). The subset row drops 3 of 220,000
samples from one batch.

| operation | bcftools | vcfr unspliced | vcfr spliced | speedup | output b / unspliced / spliced |
| --- | ---: | ---: | ---: | ---: | :---: |
| `merge`: 10 × 220k samples → BGZF | 33.7s | 7.8s | 6.4s | **5.2x** | 46M / 50M / 46M |
| `view -S ^3`: drop 3 of 220k → BGZF | 2.9s | 0.6s | 0.4s | **7.0x** | — / 5M / 4M |

Splicing here cuts vcfr's own CPU roughly in half (merge: 15.2 → 8.6
CPU-seconds) and — because copied blocks keep the input's compression — closes
the output-size gap with bcftools to ~0.3%. Decompressed output is
byte-identical with splicing on or off (`VCFR_NO_SPLICE=1` disables it).

### † Compression levels: read the BGZF rows carefully

`vcfr` deflates with libdeflate and `bcftools` with zlib, and the two
disagree about what a level means: at the shared default of 6, libdeflate is
faster but weaker, so the † rows compare unequal work — vcfr produced a
larger file. `-l 7` matches bcftools' default output size almost exactly
(the `-l 7` rows above are the fair matched-ratio comparison). The full
curve on the generic recompress case, same session:

| encoder | time | output |
| --- | ---: | ---: |
| `vcfr -l 0` (built-in Rust encoder) | 3.7s | 476M (+19%) |
| `vcfr -l 1` | 3.9s | 498M (+25%) |
| `vcfr -l 6` (default) | 11.5s | 430M (+8%) |
| `vcfr -l 7` | 17.8s | 398M (±0%) |
| `bcftools -Oz` (zlib 6, default) | 26.4s | 398M |
| `vcfr -l 8` | 37.9s | 372M (−7%) |

Rules of thumb: `-l 0` for scratch and intermediate files (7x faster than
bcftools at +19% size), the default for general use, `-l 7` when the file
should be the size bcftools would have made, `-l 8` for archives.

### ‡ `concat --naive` is I/O-bound

Both tools copy the same compressed blocks, so this is a file copy: ~130x
less CPU than recompressing (measured 1.4s CPU against 184s for 1.7 GB of
parts), with wall time set by the disk and page cache, not the code.

### Scaling notes

Measured separately, same machine: a 4x input (1.5 GB BGZF, ~10 GB text)
reproduces every ratio above within noise, and both tools scale linearly in
input size. Thread scaling on the compression-bound path: vcfr reaches 96%
parallel efficiency at 4 threads against bcftools' 71%, so the head-to-head
gap widens with cores (1.7x at one thread, 2.3x at four on the generic
recompress). On the parse-bound `-Ov` path bcftools stays at ~1.1 cores
regardless of `--threads`, while vcfr scales it. Merge memory is bounded by
line bytes, not file size: ~60 MiB on the generic set, ~750 MiB at 2.2M
samples (whole lines in flight), against bcftools' ~1.15 GB there.

Reproduce with `bench/bench.sh --sites 250000 --samples 500 --reps 3`
(fixtures are generated and cached in `bench/data/`). The harness aborts on
non-zero exits (a failed command must not produce a timing), checks free
disk first, and `sync`s between runs — each of those defends against a
failure mode that produced a wrong number at least once during development.

## Why it is fast

**Nothing is parsed until it is needed.** A record is indexed by scanning
for tab positions; fields a command doesn't touch are never looked at. A
record that survives a filter unchanged costs one `memchr` scan and one
`write_all`. Genotypes are decoded only when actually needed (AC/AN
recomputation, allele renumbering).

**BGZF is (de)compressed on every core.** Blocks are self-contained, so
inflation and deflation are farmed out to a worker pool and reassembled in
order (`src/pool.rs`). The thread budget follows the bottleneck: deflate
costs ~6x inflate, so `-O z` gives the writer the whole budget and `-O v`
gives everything to readers. `--threads 1` really is one thread — the codec
runs inline, not on a pool of one.

**Unchanged compressed blocks are never re-encoded.** Three mechanisms, in
increasing granularity: `concat --naive` copies whole files' blocks verbatim
after checking headers match; and on wide cohorts (≥16,384 samples, where
one line spans many blocks), `merge` and `view` **splice** input blocks
through: a block whose uncompressed content contains no newline lies wholly
inside one line, and when that stretch of the line is copied to the output
verbatim — merge's uniform-site fast path, or a kept run of samples under
`view -s` — the compressed block is equally valid output, since BGZF blocks
are self-contained gzip members with their own CRCs and records need not
align to block boundaries (`src/bgzf/spliced.rs`). Only boundary blocks are
re-deflated. Splicing engages per file from the header's sample count;
scattered or reordered `-s` selections have no spliceable runs and fall back
to plain recompression automatically.

**Two built-in DEFLATE encoders** (`src/deflate`). `-l 0` trains one
canonical Huffman code set on an early block and compresses every block in a
single greedy pass — no counting, no code construction, a 6-byte match hash
because VCF text is saturated with 4-byte collisions (`\t0/0`, `:12,`). It
beats libdeflate level 1 on both axes and is the fastest option. `--codec
rust` levels 1-6 are the opposite trade: a lazy hash-chain matcher plus
exact per-block Huffman codes, tracing libdeflate's speed/ratio frontier
from about l4 to l9 point-for-point (on rare-variant data `--codec rust -l 4`
beat `lib -l 7` end-to-end at matched size; on the generic set it doesn't —
data-dependent, measure on yours). Both emit ordinary DEFLATE that htslib
reads and indexes. Match extension uses AVX2 (32-byte compares behind
runtime detection, worth ~7-8%, SWAR fallback elsewhere); an AVX-512 variant
and vectorized hash insertion both measured slower and were dropped. Levels
1-12 default to libdeflate; decompression is libdeflate throughout.

**`merge` keeps its serial thread minimal.** Combining k sorted streams is
inherently serial, so that thread must not also inflate (each input gets a
reader worker) or assemble records (formatting pools, or splicing writes).
Lines are moved, not copied. At sites where every input agrees on REF, ALT
and FORMAT — essentially every site when cohorts were imputed on one panel —
allele renumbering and FORMAT union are skipped and each input's sample
columns are copied in one memcpy (or spliced compressed, above).

## Correctness

```sh
cargo test                 # 37 unit tests: codecs, splicing, VCF primitives
tests/differential.sh      # 58 paired checks against bcftools 1.19
```

The differential suite covers sample subsetting (including AC/AN
recalculation and its header declarations), variant-type and region filters,
allele renumbering across disagreeing ALT sets, `--info-rules`, all five
`-m` merge modes, BGZF interop in every direction (bgzip, htslib, both
built-in encoders, spliced streams), wide-cohort splicing against the
unspliced path byte-for-byte, and clean `| head` behaviour for every
command. Merged imputation output was additionally verified semantically at
2.2M samples: fixed columns and INFO identical after float normalisation,
all GT columns byte-identical, DS/GP numerically equal.

Deliberate differences from bcftools: INFO key order is preserved rather
than reshuffled (key order carries no meaning; the suite sorts before
comparing), floats are passed through as text rather than re-rendered
(`0.0310` stays `0.0310`), and no `##bcftools_*` provenance headers are
stamped. AC/AN follow bcftools' rules: recomputed on merge only when an
input header declares them; declared in the header when `view -s`
recomputes them.

## Commands

### `vcfr view`

Subsets samples and variants, and converts between plain VCF and BGZF.

| flag | meaning |
| --- | --- |
| `-s, --samples LIST` | keep these samples, in this order; `^` inverts |
| `-S, --samples-file FILE` | same, read from a file (`^file` inverts) |
| `-r, --regions LIST` | `chr1`, `chr2:1000-2000`, `chr3:500`, comma-separated |
| `-R, --regions-file FILE` | `CHROM[<TAB>BEG[<TAB>END]]`, 1-based inclusive |
| `-v, --types LIST` | keep `snps`, `indels`, `mnps`, `other`, `ref` |
| `-V, --exclude-types LIST` | drop those types |
| `-f, --apply-filters LIST` | require one of these FILTER values |
| `-m/-M, --min/max-alleles N` | allele-count bounds (REF included) |
| `--min-qual/--max-qual F` | QUAL bounds |
| `--min-ac N` | minimum non-reference allele count over the kept samples |
| `--ids LIST`, `--ids-file FILE` | keep records with these IDs |
| `-G, --drop-genotypes` | emit sites only |
| `-h/-H` | header only / suppress header |
| `-I, --no-update` | leave `INFO/AC` and `INFO/AN` alone when subsetting |

`vcfr` always streams, so `-r` needs no index and behaves like `bcftools
view -t` (it still stops reading once past the last requested region).
`-q`/`-Q` are deliberately not accepted as short forms of `--min-qual`/
`--max-qual`: in bcftools those letters mean allele *frequency* bounds.

### `vcfr concat`

Concatenates files carrying the same samples in the same order (the one-file
-per-chromosome case). `--naive` copies compressed blocks straight through;
it requires BGZF input and output and byte-identical headers, and refuses to
run otherwise. Files split with `bcftools view -r` are not eligible because
bcftools stamps the region into a header line, making every part's header
differ.

### `vcfr merge`

Merges files carrying *different* samples into one multi-sample file: ALT
alleles are unioned per site, genotypes renumbered onto the merged allele
list, `Number=A`/`R` values re-indexed, absent samples filled with `./.`
(or `0/0` under `-0`).

| flag | meaning |
| --- | --- |
| `-m, --merge MODE` | `snps`, `indels`, `both` (default), `all`, `none` |
| `--info-rules RULES` | `KEY:sum\|avg\|min\|max\|join\|first`, default `DP:sum,DP4:sum` |
| `-0, --missing-to-ref` | fill absent genotypes with `0/0` |
| `--force-samples` | rename duplicate sample names rather than failing |
| `-I, --no-update` | do not recompute `INFO/AC` and `INFO/AN` |

Known limits: REF padding is not performed (records at one position with
different REF strings stay separate — where bcftools instead aborts on such
input, vcfr emits them as separate records); `Number=G` values are dropped
from records whose alleles had to be renumbered, since they cannot be
re-indexed without the full genotype ordering.

### Shared output options

| flag | meaning |
| --- | --- |
| `-o, --output FILE` | output path (default stdout) |
| `-O, --output-type v\|z` | plain VCF or BGZF; inferred from a `.gz` path |
| `-l, --compression-level N` | 0-12, default 6; 0 is the built-in fast encoder |
| `--codec lib\|rust` | libdeflate (default) or the pure-Rust encoders for levels 1-6 |
| `--threads N` | worker threads, `0` = all cores (default) |

Input may be plain VCF, BGZF or ordinary gzip, detected from content rather
than extension; `-` reads stdin. Output is plain VCF or BGZF; BCF is not
supported.

## Building

```sh
cargo build --release      # target/release/vcfr
```

`scripts/build-pgo.sh` builds with profile-guided optimisation (~7% on
`-l 0` end-to-end, byte-identical output; needs `rustup component add
llvm-tools`). Dependencies: `libdeflater`, `memchr`, `flate2` (plain-gzip
input only), `crc32fast`, and `clap`.

## Test data

Two deterministic generators feed the tests and benchmarks:

```sh
cargo run --release --example gen -- --sites 100000 --samples 500 > big.vcf
cargo run --release --example gen_beagle -- --sites 150000 --samples 300 \
    --af-spectrum seq --render beagle --sample-prefix A_ --gt-seed 11 > cohortA.vcf
```

`gen.rs` makes generic call sets; sites depend only on `--site-seed`, so
runs differing in `--gt-seed`/`--sample-prefix` describe the same variants
for different cohorts — exactly what `merge` expects, with `--drop-pct` for
partial overlap. `gen_beagle.rs` mimics Beagle imputation output;
`--af-spectrum seq` draws the rare-variant-heavy sequence spectrum and
`--render beagle` reproduces real Beagle byte habits (trimmed floats,
`0|0:0:1,0,0` hom-ref fields, DR2/heterozygosity-scaled genotype
probabilities), calibrated on a snapshot of production imputed bovine
sequence data. One practical finding from that data, measured with these
tools: dropping FORMAT `DS`+`GP` shrinks a `.vcf.gz` by ~74%; dropping only
`GP` while keeping dosages saves ~38%.
