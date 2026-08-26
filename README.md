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

All numbers are best-of-3 on a 4-core machine, bcftools 1.19 / htslib 1.19,
against a 381 MB BGZF file (250,000 sites x 500 samples, ~2.5 GB of VCF text).
Both tools are given `--threads 4`.

`bench/bench.sh` records what each tool actually spent, not just how long it
took, because a speedup on its own is ambiguous: a tool can be faster because
it is better, or because it was handed more cores. "cores" is CPU seconds over
wall seconds — the average parallelism actually achieved.

| operation | bcftools | vcfr | speedup | cores b → v | CPU-s b → v | output b / v |
| --- | ---: | ---: | ---: | :---: | :---: | :---: |
| `view`: decompress to VCF | 13.84s | 0.68s | **20.2x** | 1.18 → 3.47 | 16.4 → 2.4 | — |
| `view`: BGZF → BGZF re-compress | 23.11s | 11.32s | 2.04x † | 3.28 → 3.12 | 76.0 → 35.4 | 381M / 411M |
| `view`: subset 50 of 500 samples, BGZF | 10.27s | 2.25s | 4.56x † | 1.87 → 3.07 | 19.2 → 6.9 | 43M / 46M |
| `view`: subset 50 of 500 samples, VCF | 11.59s | 1.43s | **8.1x** | 1.11 → 2.34 | 13.0 → 3.4 | — |
| `view`: SNPs only, BGZF | 20.49s | 9.57s | 2.14x † | 3.08 → 3.09 | 63.2 → 29.6 | 306M / 330M |
| `view`: biallelic SNPs + PASS, BGZF | 18.09s | 8.39s | 2.15x † | 3.12 → 3.11 | 56.6 → 26.1 | 268M / 290M |
| `view`: region `chr2`, BGZF | 12.44s | 3.32s | 3.74x † | 2.17 → 2.81 | 27.0 → 9.4 | 96M / 103M |
| `view`: drop genotypes, BGZF | 8.77s | 1.78s | 4.93x † | 1.30 → 1.69 | 11.5 → 3.0 | 3.1M / 3.2M |
| `concat`: 4 parts → BGZF | 17.31s | 12.62s | 1.37x † | 3.64 → 2.91 | 63.1 → 36.8 | 381M / 411M |
| `concat --naive`: 4 parts → BGZF | 1.02s | 0.23s | **4.4x** | 0.33 → 0.95 | 0.3 → 0.2 | 411M / 411M |
| `merge`: 3 cohorts → BGZF | 21.89s | 9.45s | 2.31x † | 3.19 → 3.47 | 70.0 → 32.8 | 331M / 356M |
| `merge`: 3 cohorts → VCF | 18.72s | 1.35s | **13.9x** | 1.08 → 3.16 | 20.4 → 4.3 | — |

Single-threaded — one OS thread each, all six measured at 0.98-0.99 cores, so
this is the cleanest like-for-like comparison in the set:

| operation | bcftools | vcfr | speedup | cores b → v |
| --- | ---: | ---: | ---: | :---: |
| `view`: decompress to VCF | 15.22s | 1.93s | **7.9x** | 0.99 → 0.99 |
| `view`: BGZF → BGZF | 75.80s | 35.12s | 2.15x † | 0.98 → 0.99 |
| `merge`: 3 cohorts → BGZF | 67.87s | 32.52s | 2.08x † | 0.99 → 0.99 |

### What the resource columns say

**The speedups are not bought with extra cores.** On every compression-bound
row both tools land within a few percent of each other on parallelism — 3.28 vs
3.12, 3.08 vs 3.09, 3.12 vs 3.11, 3.19 vs 3.47 — and on `concat` bcftools
achieves *more* (3.64 vs 2.91) and is still slower. bcftools also spawns more
OS threads than vcfr throughout (6-9 against 4-6): htslib's `--threads N` means
N workers *beyond* the main thread, and it creates more than one pool. The
CPU-seconds column is the efficiency measure, and on those rows it tracks the
speedups closely — roughly half the CPU for the same work.

**Two rows do depend on parallelism bcftools leaves on the table.** Decompress
to VCF has bcftools at 1.18 cores against vcfr's 3.47, and `merge` to VCF has
1.08 against 3.16 — bcftools barely parallelises either path. Some of those
20.2x and 13.9x figures is that gap rather than the code. The underlying
efficiency gaps are smaller but still large: 16.4 against 2.4 CPU-seconds, and
20.4 against 4.3. The single-threaded rows isolate the first honestly at 7.9x.

`concat --naive` is the mirror image: vcfr uses more of the machine (0.95 cores
against 0.33) on a job that is mostly moving bytes.

### † Read the BGZF-output rows carefully

`vcfr` deflates with libdeflate and `bcftools` with zlib, and the two disagree
about what a compression level means: at the shared default of 6, libdeflate is
faster but weaker. The rows marked † are therefore not comparing equal work —
`vcfr` produced a **larger** file. Calibrating on the re-compress case:

| encoder | time | output |
| --- | ---: | ---: |
| `bcftools -Oz` (zlib 6, the default) | 23.59s | 398,917,161 |
| `vcfr -Oz -l 6` (libdeflate 6, the default) | 12.76s | 430,610,327 (+8.0%) |
| `vcfr -Oz -l 7` | 20.88s | 398,712,552 (−0.05%) |
| `vcfr -Oz -l 8` | 47.28s | 372,607,066 (−6.6%) |

**At a matched compression ratio, `vcfr -l 7` beats `bcftools` by 1.13x on pure
re-compression**, with the outputs within 0.05% of each other. An independent
interleaved A/B run put the same figure at 1.12x. Every † row should be read the
same way: part of the margin is `vcfr` compressing less. Pass `-l 7` for output
the size `bcftools` would have produced, and treat those rows as roughly parity
on the deflate itself plus whatever the read side saves.

The two `concat --naive` figures are directly comparable — both tools copy the
same blocks and land within 4 KB of each other on a 411 MB output — so that
4.4x is like-for-like.

The rows in bold are the ones the design targets: work bounded by parsing and
moving bytes rather than by DEFLATE. Where DEFLATE dominates, `vcfr` is already
at the machine's limit — `merge -Oz` runs at 3.47 of 4 cores, within 6% of the
CPU-bound floor — and no amount of further work on the record path will move it.

### Notes on method

`bench/bench.sh` aborts rather than reporting a time for a command that exited
non-zero: an early version of this table quoted a 10.9x `--naive` speedup that
was really `vcfr` refusing to run and exiting instantly. It also flushes dirty
pages before each timed run, because every row writes hundreds of megabytes and
without that a row is timed while the previous row's writeback is still in
flight — which showed up as the same `bcftools -Oz` invocation taking 23.1s in
one place and 28.2s in another.

The one row above not taken from a single run of the script is the BGZF sample
subset: that row came out at 4.29s in the recorded run against 2.25-2.58s in
seven other measurements, so it was re-measured with the same code and the
consistent figure used.

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
