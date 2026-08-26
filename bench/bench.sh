#!/usr/bin/env bash
# Benchmark vcfr against bcftools on the same inputs and the same output type.
#
#   bench/bench.sh [--sites N] [--samples M] [--reps R] [--threads T]
#
# Every pair does identical work: same input file, same filters, same output
# format, output to /dev/null unless the operation needs a real file.
set -euo pipefail

cd "$(dirname "$0")/.."
VCFR=${VCFR:-$PWD/target/release/vcfr}
GEN=${GEN:-$PWD/target/release/examples/gen}
DATA=${DATA:-$PWD/bench/data}

SITES=200000 SAMPLES=500 REPS=3 THREADS=$(nproc)
while [[ $# -gt 0 ]]; do
  case $1 in
    --sites) SITES=$2; shift 2;;
    --samples) SAMPLES=$2; shift 2;;
    --reps) REPS=$2; shift 2;;
    --threads) THREADS=$2; shift 2;;
    *) echo "unknown option $1" >&2; exit 1;;
  esac
done

command -v bcftools >/dev/null || { echo "bcftools is required" >&2; exit 1; }
[[ -x $VCFR ]] || { echo "run: cargo build --release && cargo build --release --example gen" >&2; exit 1; }

mkdir -p "$DATA"
OUT=$(mktemp -d); trap 'rm -rf "$OUT"' EXIT

# ---------------------------------------------------------------- fixtures --
BIG=$DATA/big.$SITES.$SAMPLES.vcf.gz
if [[ ! -f $BIG ]]; then
  echo "generating $SITES sites x $SAMPLES samples ..."
  $GEN --sites "$SITES" --samples "$SAMPLES" --contigs 4 > "$DATA/big.vcf"
  bgzip -@ "$THREADS" -c "$DATA/big.vcf" > "$BIG"
  rm -f "$DATA/big.vcf"
  bcftools index -f "$BIG"
fi
# Split with vcfr, not bcftools: bcftools stamps the region it was given into
# a ##bcftools_viewCommand header line, so the parts would not share a header
# and `concat --naive` (which requires identical headers) could not run at all.
for i in 1 2 3 4; do
  P=$DATA/part$i.$SITES.$SAMPLES.vcf.gz
  [[ -f $P ]] || { $VCFR view --threads "$THREADS" -r "chr$i" -O z -o "$P" "$BIG"; bcftools index -f "$P"; }
done
COHORTS=()
for i in 1 2 3; do
  C=$DATA/cohort$i.$SITES.vcf.gz
  if [[ ! -f $C ]]; then
    $GEN --sites "$SITES" --samples $((SAMPLES / 3)) --contigs 4 \
         --sample-prefix "C${i}_" --gt-seed $((i * 17)) --drop-pct 15 --drop-seed $((i * 5)) \
      | bgzip -@ "$THREADS" -c > "$C"
    bcftools index -f "$C"
  fi
  COHORTS+=("$C")
done
PARTS="$DATA/part1.$SITES.$SAMPLES.vcf.gz $DATA/part2.$SITES.$SAMPLES.vcf.gz $DATA/part3.$SITES.$SAMPLES.vcf.gz $DATA/part4.$SITES.$SAMPLES.vcf.gz"
SAMPLE_SUBSET=$(bcftools query -l "$BIG" | head -50 | paste -sd,)

echo
printf 'input: %s (%s)\n' "$(basename "$BIG")" "$(du -h "$BIG" | cut -f1)"
printf 'threads: %s   reps: %s (best of)\n\n' "$THREADS" "$REPS"

# ------------------------------------------------------------------ timing --
best() { # best <command string> -> seconds
  local b=999999 t
  for _ in $(seq "$REPS"); do
    local s=$(date +%s.%N)
    # A benchmark that times a failing command is worse than no benchmark:
    # refuse to report a number for a command that did not succeed.
    if ! eval "$1" >/dev/null 2>"$OUT/cmd.err"; then
      echo >&2
      echo "benchmark aborted: command failed" >&2
      echo "  $1" >&2
      sed 's/^/  /' "$OUT/cmd.err" >&2
      exit 1
    fi
    local e=$(date +%s.%N)
    t=$(echo "$e - $s" | bc)
    b=$(echo "if ($t < $b) $t else $b" | bc)
  done
  printf '%.2f' "$b"
}

hdr() {
  printf '%-40s %10s %10s %9s  %s\n' "$1" "$2" "$3" "$4" "$5"
}
hdr "operation" "bcftools" "vcfr" "speedup" "output size (bcftools / vcfr)"
hdr "----------------------------------------" "----------" "----------" "---------" "----------------------------"

# row <label> <bcftools cmd> <vcfr cmd>
#
# When both commands wrote $OUT/b.vcf.gz and $OUT/v.vcf.gz the sizes are
# reported too: a BGZF speedup means nothing without knowing how hard each
# tool compressed.
row() {
  local bt vt sizes=""
  rm -f "$OUT/b.vcf.gz" "$OUT/v.vcf.gz"
  bt=$(best "$2"); vt=$(best "$3")
  if [[ -s $OUT/b.vcf.gz && -s $OUT/v.vcf.gz ]]; then
    sizes="$(du -h "$OUT/b.vcf.gz" | cut -f1) / $(du -h "$OUT/v.vcf.gz" | cut -f1)"
  fi
  printf '%-40s %9ss %9ss %8sx  %s\n' "$1" "$bt" "$vt" "$(echo "scale=2; $bt / $vt" | bc)" "$sizes"
}

row "view: decompress to VCF" \
    "bcftools view --threads $THREADS -O v -o /dev/null $BIG" \
    "$VCFR view --threads $THREADS -O v -o /dev/null $BIG"

row "view: BGZF -> BGZF re-compress" \
    "bcftools view --threads $THREADS -O z -o $OUT/b.vcf.gz $BIG" \
    "$VCFR view --threads $THREADS -O z -o $OUT/v.vcf.gz $BIG"

row "view: subset 50/$SAMPLES samples" \
    "bcftools view --threads $THREADS -s $SAMPLE_SUBSET -O z -o $OUT/b.vcf.gz $BIG" \
    "$VCFR view --threads $THREADS -s $SAMPLE_SUBSET -O z -o $OUT/v.vcf.gz $BIG"

row "view: subset samples, VCF out" \
    "bcftools view --threads $THREADS -s $SAMPLE_SUBSET -O v -o /dev/null $BIG" \
    "$VCFR view --threads $THREADS -s $SAMPLE_SUBSET -O v -o /dev/null $BIG"

row "view: SNPs only" \
    "bcftools view --threads $THREADS -v snps -O z -o $OUT/b.vcf.gz $BIG" \
    "$VCFR view --threads $THREADS -v snps -O z -o $OUT/v.vcf.gz $BIG"

row "view: biallelic SNPs, PASS" \
    "bcftools view --threads $THREADS -v snps -m2 -M2 -f PASS -O z -o $OUT/b.vcf.gz $BIG" \
    "$VCFR view --threads $THREADS -v snps -m2 -M2 -f PASS -O z -o $OUT/v.vcf.gz $BIG"

row "view: region chr2 (streaming)" \
    "bcftools view --threads $THREADS -t chr2 -O z -o $OUT/b.vcf.gz $BIG" \
    "$VCFR view --threads $THREADS -r chr2 -O z -o $OUT/v.vcf.gz $BIG"

row "view: drop genotypes (sites only)" \
    "bcftools view --threads $THREADS -G -O z -o $OUT/b.vcf.gz $BIG" \
    "$VCFR view --threads $THREADS -G -O z -o $OUT/v.vcf.gz $BIG"

row "concat: 4 parts -> BGZF" \
    "bcftools concat --threads $THREADS -O z -o $OUT/b.vcf.gz $PARTS" \
    "$VCFR concat --threads $THREADS -O z -o $OUT/v.vcf.gz $PARTS"

row "concat --naive: 4 parts -> BGZF" \
    "bcftools concat --threads $THREADS --naive -O z -o $OUT/b.vcf.gz $PARTS" \
    "$VCFR concat --threads $THREADS --naive -O z -o $OUT/v.vcf.gz $PARTS"

row "merge: 3 cohorts -> BGZF" \
    "bcftools merge --threads $THREADS -O z -o $OUT/b.vcf.gz ${COHORTS[*]}" \
    "$VCFR merge --threads $THREADS -O z -o $OUT/v.vcf.gz ${COHORTS[*]}"

row "merge: 3 cohorts -> VCF" \
    "bcftools merge --threads $THREADS -O v -o /dev/null ${COHORTS[*]}" \
    "$VCFR merge --threads $THREADS -O v -o /dev/null ${COHORTS[*]}"

echo
echo "single-threaded:"
printf '%-40s %10s %10s %9s\n' "----------------------------------------" "----------" "----------" "---------"
row "view: decompress to VCF (1 thread)" \
    "bcftools view -O v -o /dev/null $BIG" \
    "$VCFR view --threads 1 -O v -o /dev/null $BIG"
row "view: BGZF -> BGZF (1 thread)" \
    "bcftools view -O z -o $OUT/b.vcf.gz $BIG" \
    "$VCFR view --threads 1 -O z -o $OUT/v.vcf.gz $BIG"
row "merge: 3 cohorts -> BGZF (1 thread)" \
    "bcftools merge -O z -o $OUT/b.vcf.gz ${COHORTS[*]}" \
    "$VCFR merge --threads 1 -O z -o $OUT/v.vcf.gz ${COHORTS[*]}"

# ------------------------------------------------ compression calibration --
# vcfr deflates with libdeflate, bcftools with zlib, and the two disagree about
# what a given level means: libdeflate's 6 is faster but weaker. Comparing both
# at "level 6" therefore compares different amounts of work. This sweep finds
# the vcfr level that reproduces bcftools' output size, which is the honest
# basis for the -O z rows above.
echo
echo "compression calibration (BGZF re-compress, $THREADS threads, 1 run):"
printf '%-28s %10s %14s\n' "encoder" "time" "output bytes"
printf '%-28s %10s %14s\n' "----------------------------" "----------" "--------------"
one() { local s e; s=$(date +%s.%N); eval "$1" >/dev/null 2>&1; e=$(date +%s.%N); printf '%.2f' "$(echo "$e - $s" | bc)"; }
d=$(one "bcftools view --threads $THREADS -O z -o $OUT/cal.gz $BIG")
printf '%-28s %9ss %14s\n' "bcftools -Oz (zlib 6)" "$d" "$(stat -c %s "$OUT/cal.gz")"
for L in 6 7 8; do
  d=$(one "$VCFR view --threads $THREADS -O z -l $L -o $OUT/cal.gz $BIG")
  printf '%-28s %9ss %14s\n' "vcfr -Oz -l $L" "$d" "$(stat -c %s "$OUT/cal.gz")"
done
