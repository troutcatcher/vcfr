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
benchmark table before reading anything into the `-O z` numbers. The thread budget
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
| `view`: decompress to VCF | 19.40s | 0.88s | **22.1x** | 1.16 → 3.38 | 22.6 → 3.0 | — |
| `view`: BGZF → BGZF re-compress | 26.22s | 11.53s | 2.27x † | 3.43 → 3.73 | 90.1 → 43.1 | 381M / 411M |
| `view`: subset 50 of 500 samples, BGZF | 13.47s | 2.46s | 5.48x † | 1.75 → 3.37 | 23.6 → 8.3 | 43M / 46M |
| `view`: subset 50 of 500 samples, VCF | 13.27s | 1.68s | **7.9x** | 1.22 → 2.46 | 16.3 → 4.1 | — |
| `view`: SNPs only, BGZF | 22.24s | 9.27s | 2.39x † | 3.38 → 3.76 | 75.2 → 34.9 | 306M / 330M |
| `view`: biallelic SNPs + PASS, BGZF | 20.66s | 8.16s | 2.53x † | 3.33 → 3.85 | 69.0 → 31.5 | 268M / 290M |
| `view`: region `chr2`, BGZF | 15.50s | 3.30s | 4.69x † | 2.17 → 3.42 | 33.7 → 11.3 | 96M / 103M |
| `view`: drop genotypes, BGZF | 12.20s | 2.03s | 5.99x † | 1.30 → 1.69 | 15.9 → 3.4 | 3.1M / 3.2M |
| `concat`: 4 parts → BGZF | 19.34s | 11.67s | 1.65x † | 3.73 → 3.71 | 72.2 → 43.4 | 381M / 411M |
| `concat --naive`: 4 parts → BGZF | 1.15s | 0.50s | 2.3x ‡ | 0.35 → 0.61 | 0.4 → 0.3 | 411M / 411M |
| `merge`: 3 cohorts → BGZF | 27.26s | 10.98s | 2.48x † | 3.17 → 3.58 | 86.6 → 39.3 | 331M / 356M |
| `merge`: 3 cohorts → VCF | 25.16s | 1.78s | **14.2x** | 1.12 → 3.38 | 28.3 → 6.0 | — |

Single-threaded — one OS thread each, all six measured at 0.99 cores, so this is
the cleanest like-for-like comparison in the set:

| operation | bcftools | vcfr | speedup | cores b → v |
| --- | ---: | ---: | ---: | :---: |
| `view`: decompress to VCF | 20.70s | 2.20s | **9.4x** | 0.99 → 0.99 |
| `view`: BGZF → BGZF | 88.56s | 42.25s | 2.09x † | 0.99 → 0.99 |
| `merge`: 3 cohorts → BGZF | 83.05s | 39.23s | 2.11x † | 0.99 → 0.99 |

### What the resource columns say

**The speedups are not bought with extra cores.** On the compression-bound rows
both tools sit in the same band — 3.43 vs 3.73, 3.38 vs 3.76, 3.33 vs 3.85 —
and on `concat` they are level (3.73 vs 3.71) with `vcfr` still 1.65x faster.
bcftools also spawns more OS threads than `vcfr` throughout (6-9 against 5-8):
htslib's `--threads N` means N workers *beyond* the main thread, and it creates
more than one pool. The CPU-seconds column is the efficiency measure, and on
those rows it tracks the speedups closely — roughly half the CPU for the same
work.

**Two rows do depend on parallelism bcftools leaves unused.** Decompress to VCF
has bcftools at 1.16 cores against `vcfr`'s 3.38, and `merge` to VCF 1.12
against 3.38 — bcftools barely parallelises either path. Some of those 22.1x
and 14.2x figures is that gap rather than the code. The underlying efficiency
gaps are smaller but still large: 22.6 against 3.0 CPU-seconds, and 28.3 against
6.0. The single-threaded rows isolate the first honestly at 9.4x.

### † Read the BGZF-output rows carefully

`vcfr` deflates with libdeflate and `bcftools` with zlib, and the two disagree
about what a compression level means: at the shared default of 6, libdeflate is
faster but weaker. The rows marked † are therefore not comparing equal work —
`vcfr` produced a **larger** file. Calibrating on the re-compress case:

| encoder | time | output |
| --- | ---: | ---: |
| `bcftools -Oz` (zlib 6, the default) | 27.07s | 398,917,162 |
| `vcfr -Oz -l 6` (libdeflate 6, the default) | 12.75s | 430,610,327 (+8.0%) |
| `vcfr -Oz -l 7` | 18.95s | 398,712,552 (−0.05%) |
| `vcfr -Oz -l 8` | 38.31s | 372,607,066 (−6.6%) |

**At a matched compression ratio, `vcfr -l 7` beats `bcftools` by 1.43x on pure
re-compression**, with the outputs within 0.05% of each other. An independent
interleaved A/B put it at 1.45x, 1.45x, 1.44x over three rounds. Every † row
should be read the same way: part of the margin is `vcfr` compressing less.
Pass `-l 7` for output the size `bcftools` would have produced.

### ‡ `concat --naive` is I/O-bound

Both tools copy the same blocks and land within 4 KB of each other on a 411 MB
output, so the comparison is like-for-like — but it is a file copy, and its
ratio moves with page-cache and writeback state rather than with CPU. Across
runs it measured anywhere from 1.8x to 4.4x (`vcfr` 0.35-0.61s against
bcftools' 1.05-2.03s). Treat it as "roughly 2-3x, and not the interesting part".

The rows in bold are the ones the design targets: work bounded by parsing and
moving bytes rather than by DEFLATE. Where DEFLATE dominates, `vcfr` is close to
the machine's limit — `view -Oz` runs at 3.73 of 4 cores — and no further work
on the record path will move it.

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
