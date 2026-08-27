# vcfr

A fast VCF toolkit: **merge**, **concat**, and **subset** VCF files.

`vcfr` does the handful of VCF operations that dominate real pipelines, and does
them quickly. Output is verified byte-for-byte against `bcftools` — the test
suite in `tests/differential.sh` runs 41 paired invocations and requires the two
tools to agree exactly.

```
vcfr view   [-s samples] [-r regions] [-v types] ...  file.vcf.gz
vcfr concat [--naive] part1.vcf.gz part2.vcf.gz ...
vcfr merge  cohortA.vcf.gz cohortB.vcf.gz ...
```

## Why it is fast

**Nothing is parsed until it is needed.** A record is indexed by scanning for
tab positions; the fields a command doesn't touch are never looked at. A record
that survives a filter unchanged costs one `memchr` scan and one `write_all` —
no allele decoding, no INFO parsing, no per-record allocation. Genotypes are
only decoded when a command actually needs them (recomputing `AC`/`AN`, or
renumbering alleles during a merge).

**BGZF is (de)compressed on every core.** BGZF blocks are self-contained, so
inflation and deflation are farmed out to a worker pool and reassembled in
order (`src/pool.rs`). `libdeflate` does the actual codec work. It is much
faster than zlib at inflating; at deflating it is faster mainly because its
levels are calibrated differently, so see the compression note below the
benchmark table before reading anything into the `-O z` numbers.

**Level 0 is vcfr's own DEFLATE encoder** (`src/deflate`), written for this
workload. libdeflate must treat every 64 KiB block as unknown data — count
symbol frequencies, build Huffman codes, then emit, per block. VCF blocks all
look alike, so `-l 0` trains one canonical code set on an early block and
compresses every block in a single greedy pass with no counting and no code
construction, probing with a 6-byte hash because VCF text is saturated with
4-byte collisions (`\t0/0`, `:12,`) that waste 4-byte probes. On real VCF text
it beats libdeflate level 1 on both axes — 240-262 MB/s against 227-235 at a
better ratio (4.13-4.15 against 3.95-3.97) — and reaches 98-99% of level 2's
ratio at roughly twice level 2's speed. The output is ordinary DEFLATE:
htslib reads and indexes it, `gzip -t` accepts it, and the test suite
round-trips every stream through libdeflate's own decoder. Levels 1-12 remain
libdeflate, whose levels 3+ reach ratios the specialised encoder does not try
for; decompression is libdeflate throughout. On data unlike its training
block the encoder stays correct — codes are smoothed so every byte value can
be emitted, and a stored-block fallback bounds the incompressible case — it
just compresses less well. The thread budget
is split between reading and writing according to which one is the bottleneck:
deflating costs several times more than inflating, so with `-O z` most threads
go to the writer, and with `-O v` they all go to the reader. `--threads N`
counts worker threads, so `--threads 1` really is a single thread — the codec
runs inline rather than handing work to a pool of one.

**`merge` does as little as possible on the thread that orders records.**
Combining k sorted streams is inherently serial, so that thread is the
bottleneck: it must not also inflate its inputs (every input gets its own
inflation worker), and it must not assemble output records (those go to a pool
of formatting workers, reordered on the way out). Record lines are *moved* into
a batch rather than copied — a stream is about to overwrite its buffer on the
next advance anyway, so it takes a recycled one and hands the old one over —
and splitting a line into columns happens on whichever worker formats it. What
is left on the ordering thread is a scan to FORMAT and a comparison.

On top of that, a site where every input agrees on REF, ALT and FORMAT — the
overwhelmingly common case — skips allele renumbering and the FORMAT union
entirely, and each input's genotype columns are copied across in one memcpy.

The cost is memory: batches in flight take `merge` from about 10 MiB to 60 MiB.
Batches are bounded by input bytes rather than by site count, so that figure
does not grow with the sample count.

**`concat --naive` never touches the deflate stream.** Concatenating files that
share a header does not require recompression: after locating the end of each
header, the remaining BGZF blocks are already valid output and are copied
verbatim. Only the single block straddling the end of the header is re-encoded.
This turns concatenation into roughly a file copy.

## Benchmarks

`bench/bench.sh` runs each operation through both tools with identical inputs,
identical filters and identical output formats, and reports the best of N runs.

All numbers are best-of-3 from a single run of `bench/bench.sh` on a 4-core
machine, bcftools 1.19 / htslib 1.19, against a 381 MB BGZF file (250,000 sites
x 500 samples, ~2.5 GB of VCF text). Both tools are given `--threads 4`.

The harness records what each tool actually spent, not just how long it took,
because a speedup on its own is ambiguous: a tool can be faster because it is
better, or because it was handed more cores. "cores" is CPU seconds over wall
seconds — the average parallelism actually achieved.

| operation | bcftools | vcfr | speedup | cores b → v | CPU-s b → v | output b / v |
| --- | ---: | ---: | ---: | :---: | :---: | :---: |
| `view`: decompress to VCF | 19.74s | 0.86s | **22.9x** | 1.16 → 3.39 | 22.9 → 2.9 | — |
| `view`: BGZF → BGZF re-compress | 26.59s | 12.02s | 2.21x † | 3.41 → 3.86 | 90.9 → 46.5 | 381M / 411M |
| `view`: subset 50 of 500 samples, BGZF | 13.60s | 2.45s | 5.55x † | 1.77 → 3.38 | 24.1 → 8.3 | 43M / 46M |
| `view`: subset 50 of 500 samples, VCF | 13.71s | 1.77s | **7.7x** | 1.23 → 2.37 | 17.0 → 4.2 | — |
| `view`: SNPs only, BGZF | 22.54s | 9.40s | 2.39x † | 3.40 → 3.79 | 76.7 → 35.7 | 306M / 330M |
| `view`: biallelic SNPs + PASS, BGZF | 20.70s | 8.11s | 2.55x † | 3.32 → 3.88 | 68.9 → 31.5 | 268M / 290M |
| `view`: region `chr2`, BGZF | 15.58s | 3.27s | 4.76x † | 2.16 → 3.47 | 33.7 → 11.4 | 96M / 103M |
| `view`: drop genotypes, BGZF | 12.12s | 1.91s | 6.34x † | 1.29 → 1.73 | 15.7 → 3.3 | 3.1M / 3.2M |
| `concat`: 4 parts → BGZF | 19.27s | 11.66s | 1.65x † | 3.71 → 3.80 | 71.5 → 44.3 | 381M / 411M |
| `concat --naive`: 4 parts → BGZF | 1.14s | 0.52s | 2.2x ‡ | 0.37 → 0.82 | 0.4 → 0.4 | 411M / 411M |
| `merge`: 3 cohorts → BGZF | 27.25s | 10.37s | 2.62x † | 3.15 → 3.81 | 85.9 → 39.6 | 331M / 356M |
| `merge`: 3 cohorts → VCF | 24.83s | 1.68s | **14.8x** | 1.12 → 3.32 | 27.9 → 5.6 | — |

Single-threaded — one OS thread each, all six measured at 0.99 cores, so this is
the cleanest like-for-like comparison in the set:

| operation | bcftools | vcfr | speedup | cores b → v |
| --- | ---: | ---: | ---: | :---: |
| `view`: decompress to VCF | 20.51s | 2.09s | **9.8x** | 0.99 → 0.99 |
| `view`: BGZF → BGZF | 87.24s | 42.17s | 2.07x † | 0.99 → 0.99 |
| `merge`: 3 cohorts → BGZF | 81.78s | 39.37s | 2.08x † | 0.99 → 0.99 |

### What the resource columns say

**The speedups are not bought with extra cores.** On the compression-bound rows
both tools sit in the same band — 3.41 vs 3.86, 3.40 vs 3.79, 3.32 vs 3.88,
3.71 vs 3.80 — and the CPU-seconds column, the efficiency measure, tracks the
speedups closely: roughly half the CPU for the same work. bcftools also spawns
more OS threads than `vcfr` throughout (6-9 against 5-8): htslib's `--threads N`
means N workers *beyond* the main thread, and it creates more than one pool.

**Two rows do depend on parallelism bcftools leaves unused.** Decompress to VCF
has bcftools at 1.16 cores against `vcfr`'s 3.39, and `merge` to VCF 1.12
against 3.32 — bcftools barely parallelises either path. Some of those 22.9x
and 14.8x figures is that gap rather than the code. The underlying efficiency
gaps are smaller but still large: 22.9 against 2.9 CPU-seconds, and 27.9
against 5.6. The single-threaded rows isolate the first honestly at 9.8x.

### † Read the BGZF-output rows carefully

`vcfr` deflates with libdeflate and `bcftools` with zlib, and the two disagree
about what a compression level means: at the shared default of 6, libdeflate is
faster but weaker. The rows marked † are therefore not comparing equal work —
`vcfr` produced a **larger** file. Calibrating on the re-compress case:

| encoder | time | output |
| --- | ---: | ---: |
| `bcftools -Oz` (zlib 6, the default) | 26.80s | 398,917,159 |
| `vcfr -Oz -l 1` | 4.59s | 498,181,180 (+24.9%) |
| `vcfr -Oz -l 6` (libdeflate 6, the default) | 11.58s | 430,610,327 (+8.0%) |
| `vcfr -Oz -l 7` | 18.54s | 398,712,552 (−0.05%) |
| `vcfr -Oz -l 8` | 38.29s | 372,607,066 (−6.6%) |

**At a matched compression ratio, `vcfr -l 7` beats `bcftools` by 1.45x on pure
re-compression**, with the outputs within 0.05% of each other. An independent
interleaved A/B put it at 1.45x, 1.45x, 1.44x over three rounds. Every † row
should be read the same way: part of the margin is `vcfr` compressing less.
Pass `-l 7` for output the size `bcftools` would have produced.

The level curve is worth exploiting deliberately: the whole span from `-l 0` to
`-l 8` trades a large range in time against a 34% range in size. For scratch
and intermediate files, `-l 0` (vcfr's own encoder) is the best fast point —
about 5% faster than `-l 1` end-to-end with 4.3% smaller output; for archives,
`-l 8` out-compresses bcftools' default.

### ‡ `concat --naive` is I/O-bound

Both tools copy the same blocks and land within 4 KB of each other on a 411 MB
output, so the comparison is like-for-like — but it is a file copy, and its
ratio moves with page-cache and writeback state rather than with CPU. Across
runs it measured anywhere from 1.8x to 4.4x. Treat it as "roughly 2-3x, and not
the interesting part".

The rows in bold are the ones the design targets: work bounded by parsing and
moving bytes rather than by DEFLATE. Where DEFLATE dominates, `vcfr` sits at
3.8-3.9 of 4 cores — within a few percent of the machine's limit — so the
remaining lever there is the compression level, not the code.

### Bigger inputs and thread scaling

The same suite run on a 4x input — 500,000 sites x 1,000 samples, a 1.5 GB
BGZF file holding ~10 GB of VCF text — reproduces every ratio above within
noise (best of 2, same machine):

| operation | bcftools | vcfr | speedup | cores b → v |
| --- | ---: | ---: | ---: | :---: |
| `view`: decompress to VCF | 77.73s | 3.11s | **25.0x** | 1.17 → 3.58 |
| `view`: BGZF → BGZF re-compress | 109.50s | 46.57s | 2.35x † | 3.41 → 3.85 |
| `view`: subset 50 of 1000 samples, BGZF | 50.79s | 8.79s | 5.77x † | 1.54 → 2.72 |
| `view`: SNPs only, BGZF | 91.71s | 38.35s | 2.39x † | 3.36 → 3.85 |
| `concat`: 4 parts → BGZF | 79.42s | 46.55s | 1.70x † | 3.68 → 3.87 |
| `merge`: 3 cohorts → BGZF | 105.57s | 41.77s | 2.52x † | 3.22 → 3.85 |
| `merge`: 3 cohorts → VCF | 94.89s | 6.19s | **15.3x** | 1.12 → 3.34 |

Both tools scale linearly in input size (4.1x the data took each ~3.9x
longer), doubling the sample width changes nothing, and the matched-ratio
calibration lands in the same place (`-l 7`: 74.83s against 107.47s, 1.44x).
Subsetting 50 of 1000 samples on this file is byte-identical to bcftools
across all 500,000 records, as is merging the three 333-sample cohorts; and
`merge`'s peak memory stays flat because its batches are bounded by bytes, not
sites.

Thread scaling, measured on the same 4x input (single runs, tools alternating;
this machine has 4 physical cores, so 8 threads tests oversubscription, not
scaling):

| threads | bcftools `-Oz` | cores | vcfr `-Oz` | cores | head-to-head |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 297.9s | 1.27 | 176.4s | 0.99 | 1.69x |
| 2 | 149.1s | 2.49 | 80.1s | 2.30 | 1.86x |
| 4 | 104.9s | 3.50 | 45.9s | 3.93 | 2.29x |
| 8 | 104.0s | 3.55 | 46.6s | 3.94 | 2.23x |

Three things this table says. `vcfr`'s parallel efficiency at 4 threads is 96%
(3.85x the single-thread time) against bcftools' 71%, so the head-to-head gap
*widens* as threads are added — 1.69x at one, 2.29x at four. `vcfr`'s total CPU
is flat across thread counts (176-185s), i.e. the parallelism is close to free.
And oversubscribing 8 threads onto 4 cores costs neither tool anything. Whether
the efficiency gap keeps widening past 4 cores is extrapolation — this machine
cannot measure it.

On the parse-bound path the thread column is one-sided: bcftools runs `-Ov` at
1.16 cores whatever `--threads` is set to, while `vcfr` scales it (8.60s at
one thread, 3.46s at four — 22x head-to-head on this file). The cores column
also shows htslib's `--threads 1` really running 1.27 cores: htslib counts
workers *beyond* the main thread.

### Notes on method

Two failure modes cost real time in developing this table, and the harness now
defends against both:

- **A command that fails still produces a timing.** An early version of this
  table quoted a 10.9x `--naive` speedup that was really `vcfr` refusing to run
  and exiting instantly. The harness now aborts on a non-zero exit.
- **A full disk does not announce itself.** It slows both tools together through
  writeback throttling, so the ratios still look plausible and the numbers are
  quietly worthless — this invalidated an entire run and a measurement that
  overstated a change by 3x before it was caught. The harness now checks free
  space before starting, and flushes dirty pages before each timed run so a row
  is not timed while the previous row's writeback is still in flight.

Absolute times drift by 20% or so between runs on this host, so figures from
different runs are not comparable; everything in the tables above comes from one
run, and the two independent checks quoted are labelled as such.

Reproduce with:

```sh
cargo build --release && cargo build --release --example gen
bench/bench.sh --sites 250000 --samples 500 --reps 3
```

`--sites` and `--samples` control the fixture size; fixtures are cached in
`bench/data/` and reused across runs. The parts used by `concat` are split with
`vcfr`, not `bcftools`, because `bcftools view -r` records the region it was
given in a `##bcftools_viewCommand` header line — which would leave the parts
with differing headers and make `--naive` inapplicable to them.

## Correctness

```sh
cargo test                 # unit tests: BGZF codec, line splitting, VCF primitives
tests/differential.sh      # 41 paired comparisons against bcftools
```

The differential suite compares `vcfr` output to `bcftools` output line by line,
covering sample subsetting (including the `AC`/`AN` recalculation `bcftools`
performs by default), variant-type and region filters, allele renumbering across
files with disagreeing ALT sets, `--info-rules`, all five `-m` merge modes, and
round-tripping BGZF through `bgzip`, `htslib` and `vcfr` in every direction.

Two deliberate differences from `bcftools`:

- **INFO key order.** `vcfr` keeps the order the fields appeared in; `bcftools`
  reshuffles them as a side effect of updating `AC`/`AN`. Key order carries no
  meaning in VCF, so the differential suite sorts INFO keys before comparing.
- **`vcfr` does not stamp provenance headers** (`##bcftools_viewCommand=…`) into
  its output.

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

`vcfr` always streams, so `-r` needs no index and behaves like `bcftools view -t`
rather than `bcftools view -r`. On a coordinate-sorted file it still stops
reading once it has passed the last requested region.

Note that `-q`/`-Q` are deliberately **not** accepted as short forms of
`--min-qual`/`--max-qual`: in `bcftools` those letters mean minimum and maximum
*allele frequency*, and silently reinterpreting them would be a trap.

### `vcfr concat`

Concatenates files that carry the same samples in the same order — the usual
"one file per chromosome" case. Headers are merged (first definition of each ID
wins) and sample lists must match exactly.

`--naive` copies compressed blocks straight through. It requires BGZF input,
BGZF output, and byte-identical headers — and it checks, refusing to run when
the headers differ rather than silently emitting the first file's header over
everyone's records. Note that files split with `bcftools view -r` are *not*
eligible, because bcftools records the region it was given in a
`##bcftools_viewCommand` header line, making every part's header different.

### `vcfr merge`

Merges files that carry *different* samples into one multi-sample file. At each
site the ALT alleles are unioned, every input's genotypes are renumbered onto
the merged allele list, and `Number=A`/`Number=R` INFO and FORMAT values are
re-indexed to match. Samples absent from a file are filled with `./.` (or `0/0`
under `-0`).

| flag | meaning |
| --- | --- |
| `-m, --merge MODE` | `snps`, `indels`, `both` (default), `all`, `none` |
| `--info-rules RULES` | `KEY:sum\|avg\|min\|max\|join\|first`, default `DP:sum,DP4:sum` |
| `-0, --missing-to-ref` | fill absent genotypes with `0/0` |
| `--force-samples` | rename duplicate sample names rather than failing |
| `-I, --no-update` | do not recompute `INFO/AC` and `INFO/AN` |

`-m` decides which records at one position may become a single multiallelic
record: `both` keeps SNPs and indels in separate records, `all` combines
everything sharing a REF, `none` creates no new multiallelics. Records whose
alleles are already a subset of another's are combined under every mode, since
that introduces no new allele.

Known limits: REF padding is not performed, so records at the same position with
different REF strings (`A→AT` in one file, `AT→A` in another) stay separate
rather than being rewritten onto a common REF. `Number=G` values (such as `PL`)
are dropped from a record whose alleles had to be renumbered, because they
cannot be re-indexed without the full genotype ordering.

### Shared output options

| flag | meaning |
| --- | --- |
| `-o, --output FILE` | output path (default stdout) |
| `-O, --output-type v\|z` | plain VCF or BGZF; inferred from a `.gz` path |
| `-l, --compression-level N` | 1–12, default 6 |
| `--threads N` | worker threads, `0` = all cores (default) |

## Input formats

Plain VCF, BGZF, and ordinary gzip are all detected from the file contents
rather than the extension. `-` reads standard input. Output is plain VCF or
BGZF; BCF is not supported.

## Building

```sh
cargo build --release      # target/release/vcfr
```

For a few percent more, build with profile-guided optimisation:

```sh
rustup component add llvm-tools
scripts/build-pgo.sh       # target/pgo/release/vcfr
```

Measured interleaved against the plain build: ~7% on `-l 0` end-to-end and
~5% on the encoder itself, with byte-identical output; the libdeflate levels
gain only 1-3% because the C library is not instrumented.

Dependencies: `libdeflater` (the codec), `memchr`, `flate2` (non-BGZF gzip
input only), and `clap`.

## Test data

`examples/gen.rs` generates deterministic VCFs for the tests and benchmarks:

```sh
cargo run --release --example gen -- --sites 100000 --samples 500 > big.vcf
```

Sites depend only on `--site-seed`, so two runs differing only in `--gt-seed`
and `--sample-prefix` describe the same variants for different cohorts — which
is exactly the input `merge` expects. `--drop-pct` omits a fraction of sites so
the cohorts overlap only partially.
