#!/usr/bin/env bash
# Profile-guided build of vcfr.
#
# PGO cannot be the default `cargo build` because it needs a training run:
# build instrumented, exercise the hot paths on representative data, then
# rebuild optimised against the collected profile. Measured on the benchmark
# machine it is worth ~5% on the level-0 encoder and smaller amounts on the
# parsing paths, at zero cost to output bytes.
#
#   scripts/build-pgo.sh            # leaves the binary at target/pgo/release/vcfr
set -euo pipefail
cd "$(dirname "$0")/.."

PROFDIR=$(mktemp -d)
WORK=$(mktemp -d)
trap 'rm -rf "$PROFDIR" "$WORK"' EXIT

LLVM_BIN=$(dirname "$(find "$(rustc --print sysroot)" -name llvm-profdata | head -1)")
[[ -x $LLVM_BIN/llvm-profdata ]] || {
  echo "llvm-profdata not found; run: rustup component add llvm-tools" >&2
  exit 1
}

echo "== instrumented build"
RUSTFLAGS="-Cprofile-generate=$PROFDIR" cargo build --release --target-dir target/pgo-gen
cargo build --release --example gen

echo "== training run"
GEN=target/release/examples/gen
V=target/pgo-gen/release/vcfr
$GEN --sites 20000 --samples 100 --contigs 2 > "$WORK/t.vcf"
$GEN --sites 20000 --samples 40 --contigs 2 --sample-prefix A --gt-seed 5 --drop-pct 15 > "$WORK/a.vcf"
$GEN --sites 20000 --samples 40 --contigs 2 --sample-prefix B --gt-seed 9 --drop-pct 15 > "$WORK/b.vcf"
$V view --threads 2 -O z -l 0 -o "$WORK/t0.gz" "$WORK/t.vcf"
$V view --threads 2 -O z -l 6 -o "$WORK/t6.gz" "$WORK/t.vcf"
$V view --threads 2 -O v -o /dev/null "$WORK/t0.gz"
$V view --threads 2 -s "$(head -1 <($V view -h "$WORK/t.vcf" | tail -1 | cut -f10-12 --output-delimiter=,))" \
      -O v -o /dev/null "$WORK/t.vcf" 2>/dev/null || true
$V merge --threads 2 -O v -o /dev/null "$WORK/a.vcf" "$WORK/b.vcf"

echo "== optimised rebuild"
"$LLVM_BIN/llvm-profdata" merge -o "$PROFDIR/merged.profdata" "$PROFDIR"/*.profraw
RUSTFLAGS="-Cprofile-use=$PROFDIR/merged.profdata" cargo build --release --target-dir target/pgo

echo "done: target/pgo/release/vcfr"
