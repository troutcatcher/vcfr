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

# A full filesystem does not fail loudly here -- it silently inflates every
# timing through writeback throttling, and truncates outputs. Both tools get
# slower together, so the ratios still look plausible and nothing announces
# that the run is worthless. Refuse to start instead.
need=6
[[ -f $DATA/big.$SITES.$SAMPLES.vcf.gz ]] || need=$((need + 4))
for dir in "$DATA" "$OUT"; do
  avail=$(df -BG --output=avail "$dir" | tail -1 | tr -dc '0-9')
  if (( avail < need )); then
    echo "not enough free space on $(df --output=target "$dir" | tail -1): ${avail}G available, ~${need}G needed" >&2
    exit 1
  fi
done

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
#
# Reporting a speedup without resource usage is misleading: a tool can be
# faster because it is better, or just because it was handed more cores. So
# every run records wall time, CPU seconds and peak thread count.
#
# CPU comes from the shell's own `times` accounting, which is exact kernel
# rusage at reap time rather than a sampled estimate. It only works when
# `times` runs in THIS shell, so run_once reports through globals instead of
# stdout -- $(run_once ...) would put it in a subshell and read back zero.
# Timestamps use $EPOCHREALTIME so nothing forks between the two snapshots.

RUN_WALL=0 RUN_CPU=0 RUN_THREADS=1

cpu_snapshot() { # <file> -> child user+sys seconds
  awk 'NR==2{split($1,u,/[ms]/); split($2,s,/[ms]/);
             printf "%.3f", u[1]*60+u[2]+s[1]*60+s[2]}' "$1"
}

run_once() { # run_once <command string>; sets RUN_WALL RUN_CPU RUN_THREADS
  local start end pid poller rc
  # Flush dirty pages first. These commands write hundreds of megabytes, and
  # without this a row is timed while the previous row's writeback is still in
  # flight -- which showed up as the same bcftools invocation taking 22.8s in
  # one place and 28.2s in another.
  sync
  times > "$OUT/t0"
  start=$EPOCHREALTIME
  # `exec` so that $! is the benchmarked process itself. Plain `eval ... &`
  # backgrounds a wrapper subshell and $! names that instead, which makes the
  # thread count below read 1 for every tool.
  eval "exec $1" >/dev/null 2>"$OUT/cmd.err" &
  pid=$!
  # Peak thread count, sampled from /proc in a subshell so it cannot delay the
  # wall-clock measurement. Every read here is a bash builtin.
  (
    peak=1
    while [[ -d /proc/$pid/task ]]; do
      tasks=(/proc/"$pid"/task/*)
      # `if`, not `&&`: a false (( )) returns 1, which under `set -e` would
      # kill this subshell on the first sample that is not a new maximum.
      if (( ${#tasks[@]} > peak )); then peak=${#tasks[@]}; fi
      sleep 0.02
    done
    echo "$peak" > "$OUT/peak"
  ) &
  poller=$!
  rc=0
  wait "$pid" || rc=$?
  # Snapshot before reaping the poller, so only the benchmarked command's CPU
  # lands in the difference.
  times > "$OUT/t1"
  end=$EPOCHREALTIME
  wait "$poller" 2>/dev/null || true

  if [[ $rc -ne 0 ]]; then
    echo >&2
    echo "benchmark aborted: command failed (exit $rc)" >&2
    echo "  $1" >&2
    sed 's/^/  /' "$OUT/cmd.err" >&2
    exit 1
  fi
  RUN_WALL=$(echo "$end - $start" | bc)
  RUN_CPU=$(echo "$(cpu_snapshot "$OUT/t1") - $(cpu_snapshot "$OUT/t0")" | bc)
  RUN_THREADS=$(cat "$OUT/peak" 2>/dev/null || echo 1)
  rm -f "$OUT/peak"
}

# Best-of-REPS wall time, keeping the resource figures from the fastest run.
BEST_WALL=0 BEST_CPU=0 BEST_THREADS=1
best() { # best <command string>; sets BEST_*
  BEST_WALL=999999
  for _ in $(seq "$REPS"); do
    run_once "$1"
    if [[ $(echo "$RUN_WALL < $BEST_WALL" | bc) -eq 1 ]]; then
      BEST_WALL=$RUN_WALL BEST_CPU=$RUN_CPU BEST_THREADS=$RUN_THREADS
    fi
  done
}

hdr() {
  printf '%-40s %10s %10s %9s  %s\n' "$1" "$2" "$3" "$4" "$5"
}
hdr "operation" "bcftools" "vcfr" "speedup" "output size (bcftools / vcfr)"
hdr "----------------------------------------" "----------" "----------" "---------" "----------------------------"

# row <label> <bcftools cmd> <vcfr cmd>
#
# Prints wall time and speedup, then what each tool spent to get there: peak OS
# threads, CPU seconds and average cores in use. When both commands wrote
# $OUT/{b,v}.vcf.gz their sizes are shown too, since a BGZF speedup means
# nothing without knowing how hard each tool compressed.
row() {
  local sizes="" bt bc bp vt vc vp
  rm -f "$OUT/b.vcf.gz" "$OUT/v.vcf.gz"
  best "$2"; bt=$BEST_WALL bc=$BEST_CPU bp=$BEST_THREADS
  best "$3"; vt=$BEST_WALL vc=$BEST_CPU vp=$BEST_THREADS
  if [[ -s $OUT/b.vcf.gz && -s $OUT/v.vcf.gz ]]; then
    sizes="$(du -h "$OUT/b.vcf.gz" | cut -f1) / $(du -h "$OUT/v.vcf.gz" | cut -f1)"
  fi
  printf '%-40s %9.2fs %9.2fs %8sx  %s\n' \
    "$1" "$bt" "$vt" "$(echo "scale=2; $bt / $vt" | bc)" "$sizes"
  printf '%-40s   bcftools %2s thr %7ss cpu %5s cores  |  vcfr %2s thr %7ss cpu %5s cores\n' \
    "" "$bp" "$bc" "$(echo "scale=2; $bc / $bt" | bc)" \
       "$vp" "$vc" "$(echo "scale=2; $vc / $vt" | bc)"
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
echo "compression calibration (BGZF re-compress, $THREADS threads, best of $REPS):"
printf '%-28s %10s %14s\n' "encoder" "time" "output bytes"
printf '%-28s %10s %14s\n' "----------------------------" "----------" "--------------"
# Best-of-REPS like the rows above: a single run here lands right after the
# previous encoder flushed hundreds of megabytes, and measures writeback rather
# than compression.
one() { best "$1"; printf '%.2f' "$BEST_WALL"; }
d=$(one "bcftools view --threads $THREADS -O z -o $OUT/cal.gz $BIG")
printf '%-28s %9ss %14s\n' "bcftools -Oz (zlib 6)" "$d" "$(stat -c %s "$OUT/cal.gz")"
for L in 6 7 8; do
  d=$(one "$VCFR view --threads $THREADS -O z -l $L -o $OUT/cal.gz $BIG")
  printf '%-28s %9ss %14s\n' "vcfr -Oz -l $L" "$d" "$(stat -c %s "$OUT/cal.gz")"
done
