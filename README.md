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
go to the writer, and with `-O v` they all go to the reader.

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
`--threads 4` is passed to both tools.

| operation | bcftools | vcfr | speedup |
| --- | ---: | ---: | ---: |
| `view`: decompress to VCF | 13.78s | 0.80s | **17.2x** |
| `view`: BGZF -> BGZF re-compress | 24.06s | 12.02s | 2.00x † |
| `view`: subset 50 of 500 samples, BGZF out | 10.11s | 3.39s | 2.98x † |
| `view`: subset 50 of 500 samples, VCF out | 9.93s | 1.39s | **7.1x** |
| `view`: SNPs only, BGZF out | 18.87s | 9.37s | 2.01x † |
| `view`: biallelic SNPs + PASS, BGZF out | 17.59s | 8.35s | 2.10x † |
| `view`: region `chr2`, BGZF out | 11.94s | 3.31s | 3.60x † |
| `view`: drop genotypes, BGZF out | 8.63s | 2.71s | 3.18x † |
| `concat`: 4 parts -> BGZF | 16.71s | 12.61s | 1.32x † |
| `concat --naive`: 4 parts -> BGZF | 0.94s | 0.31s | **3.0x** |
| `merge`: 3 cohorts -> BGZF | 20.58s | 10.40s | 1.97x † |
| `merge`: 3 cohorts -> VCF | 18.89s | 5.50s | **3.4x** |

Single-threaded, so the codec is compared without any parallelism:

| operation | bcftools | vcfr | speedup |
| --- | ---: | ---: | ---: |
| `view`: decompress to VCF | 15.14s | 1.91s | 7.9x |
| `view`: BGZF -> BGZF | 72.08s | 33.76s | 2.13x † |
| `merge`: 3 cohorts -> BGZF | 67.98s | 28.40s | 2.39x † |

### † Read the BGZF-output rows carefully

`vcfr` deflates with libdeflate and `bcftools` with zlib, and the two disagree
about what a compression level means: at the shared default of 6, libdeflate is
faster but weaker. So the rows marked † are not comparing equal work — `vcfr`
produced a **larger** file. Calibrating on the re-compress case:

| encoder | time | output |
| --- | ---: | ---: |
| `bcftools -Oz` (zlib 6, the default) | 23.23s | 398,917,137 |
| `vcfr -Oz -l 6` (libdeflate 6, the default) | 11.65s | 430,610,327 (+8.0%) |
| `vcfr -Oz -l 7` | 21.91s | 398,712,552 (-0.1%) |
| `vcfr -Oz -l 8` | 47.02s | 372,607,066 (-6.6%) |

**At a matched compression ratio, `vcfr -l 7` beats `bcftools` by only 1.06x on
pure re-compression.** Every † row should be read the same way: part of the
margin is `vcfr` compressing less. Pick `-l 7` when you want output the size
`bcftools` would have produced, and treat those rows as roughly parity plus
whatever the read side saves.

The two `concat --naive` figures are directly comparable — both tools copy the
same blocks and land within 4 KB of each other on a 431 MB output — so that 3.0x
is a like-for-like result.

The rows in bold are the honest wins, and they are the ones the design targets:
work bounded by parsing and moving bytes rather than by DEFLATE. Decompressing
to VCF is 17x because `vcfr` never builds a record object; subsetting samples to
VCF is 7x for the same reason; `merge` to VCF is 3.4x because only the genotype
fields that actually change get rebuilt; and `concat --naive` is 3x because
neither tool compresses anything, leaving only how fast the bytes move.

`bench/bench.sh` reports the output size of every BGZF row so this comparison is
visible by default, and aborts rather than reporting a time for a command that
exited non-zero — an early version of this table quoted a 10.9x `--naive`
speedup that was really `vcfr` refusing to run and exiting instantly.

Reproduce with:

```sh
cargo build --release && cargo build --release --example gen
bench/bench.sh --sites 250000 --samples 500 --reps 3
```

`--sites` and `--samples` control the fixture size; the fixtures are cached in
`bench/data/` and reused across runs.

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
