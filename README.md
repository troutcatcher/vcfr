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

**`merge` keeps inflation off the thread that assembles records.** Combining
k sorted streams is inherently serial, so that thread is the bottleneck and
must not also be inflating its inputs: every input gets its own inflation
worker. Reader threads spend most of a run blocked, so they are not charged
against the writer's share. Above that, a site where every input agrees on REF,
ALT and FORMAT — the overwhelmingly common case — skips allele renumbering and
the FORMAT union entirely, and each input's genotype columns are copied across
in a single memcpy.

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
| `view`: decompress to VCF | 13.71s | 0.70s | **19.6x** | 1.18 → 3.42 | 16.3 → 2.4 | — |
| `view`: BGZF → BGZF re-compress | 24.26s | 11.67s | 2.07x † | 3.08 → 3.06 | 74.9 → 35.8 | 381M / 411M |
| `view`: subset 50 of 500 samples, BGZF | 10.29s | 2.16s | 4.76x † | 1.83 → 3.04 | 18.9 → 6.6 | 43M / 46M |
| `view`: subset 50 of 500 samples, VCF | 9.97s | 1.47s | **6.8x** | 1.23 → 2.36 | 12.3 → 3.5 | — |
| `view`: SNPs only, BGZF | 19.48s | 9.33s | 2.08x † | 3.19 → 3.12 | 62.3 → 29.2 | 306M / 330M |
| `view`: biallelic SNPs + PASS, BGZF | 18.95s | 8.49s | 2.23x † | 3.02 → 3.14 | 57.4 → 26.7 | 268M / 290M |
| `view`: region `chr2`, BGZF | 12.18s | 3.28s | 3.70x † | 2.23 → 2.82 | 27.2 → 9.3 | 96M / 103M |
| `view`: drop genotypes, BGZF | 8.77s | 1.69s | 5.19x † | 1.31 → 1.73 | 11.6 → 2.9 | 3.1M / 3.2M |
| `concat`: 4 parts → BGZF | 17.23s | 11.39s | 1.51x † | 3.64 → 3.13 | 62.8 → 35.7 | 381M / 411M |
| `concat --naive`: 4 parts → BGZF | 0.80s | 0.27s | **2.9x** | 0.41 → 0.98 | 0.3 → 0.3 | 411M / 411M |
| `merge`: 3 cohorts → BGZF | 21.21s | 9.98s | 2.12x † | 3.26 → 3.37 | 69.2 → 33.7 | 331M / 356M |
| `merge`: 3 cohorts → VCF | 18.77s | 3.53s | **5.3x** | 1.11 → 1.45 | 20.9 → 5.1 | — |

Single-threaded — one OS thread each, both measured at 0.98-0.99 cores, so this
is the cleanest like-for-like comparison in the set:

| operation | bcftools | vcfr | speedup | cores b → v |
| --- | ---: | ---: | ---: | :---: |
| `view`: decompress to VCF | 14.57s | 1.96s | **7.4x** | 0.99 → 0.99 |
| `view`: BGZF → BGZF | 73.69s | 35.43s | 2.08x † | 0.99 → 0.98 |
| `merge`: 3 cohorts → BGZF | 67.11s | 33.44s | 2.00x † | 0.99 → 0.99 |

### What the resource columns say

**The speedups are not bought with extra cores.** On every compression-bound
row both tools land within a few percent of each other on parallelism — 3.08 vs
3.06, 3.19 vs 3.12, 3.02 vs 3.14, 3.26 vs 3.37 — and on `concat` bcftools
actually achieves *more* (3.64 vs 3.13) and is still slower. bcftools also
spawns more OS threads than vcfr throughout (6-9 against 4-5): htslib's
`--threads N` means N workers *beyond* the main thread, and it creates more
than one pool. The CPU-seconds column is the efficiency measure, and it tracks
the speedups closely on those rows: roughly half the CPU for the same work.

**One row does depend on parallelism bcftools leaves on the table.** Decompress
to VCF has bcftools at 1.18 cores against vcfr's 3.42 — it barely parallelises
that path. Some of that 19.6x is therefore the gap in parallelism rather than
in the code. The underlying efficiency gap is smaller but still large: 16.3
against 2.4 CPU-seconds. The single-threaded row isolates it honestly at 7.4x.

### † Read the BGZF-output rows carefully

`vcfr` deflates with libdeflate and `bcftools` with zlib, and the two disagree
about what a compression level means: at the shared default of 6, libdeflate is
faster but weaker. The rows marked † are therefore not comparing equal work —
`vcfr` produced a **larger** file. Calibrating on the re-compress case:

| encoder | time | output |
| --- | ---: | ---: |
| `bcftools -Oz` (zlib 6, the default) | 23.08s | 398,917,134 |
| `vcfr -Oz -l 6` (libdeflate 6, the default) | 11.54s | 430,610,327 (+8.0%) |
| `vcfr -Oz -l 7` | 20.61s | 398,712,552 (−0.05%) |
| `vcfr -Oz -l 8` | 45.78s | 372,607,066 (−6.6%) |

**At a matched compression ratio, `vcfr -l 7` beats `bcftools` by 1.12x on pure
re-compression** (interleaved A/B, best of 4, outputs within 0.05% of each
other). Every † row should be read the same way: part of the margin is `vcfr`
compressing less. Pass `-l 7` for output the size `bcftools` would have
produced, and treat those rows as roughly parity on the deflate itself plus
whatever the read side saves.

The two `concat --naive` figures are directly comparable — both tools copy the
same blocks and land within 4 KB of each other on a 411 MB output — so that
2.9x is like-for-like.

The rows in bold are the ones the design targets: work bounded by parsing and
moving bytes rather than by DEFLATE. Decompressing to VCF is fast because
`vcfr` never builds a record object; subsetting samples to VCF is fast for the
same reason; `merge` to VCF is fast because only the genotype fields that
actually change get rebuilt; and `concat --naive` because neither tool
compresses anything, leaving only how fast the bytes move.

`bench/bench.sh` aborts rather than reporting a time for a command that exited
non-zero — an early version of this table quoted a 10.9x `--naive` speedup that
was really `vcfr` refusing to run and exiting instantly.

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
