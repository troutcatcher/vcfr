#!/usr/bin/env bash
# Differential test: every vcfr command is compared against the equivalent
# bcftools invocation and must match byte for byte (## header lines excluded,
# since bcftools stamps its own provenance into them).
set -uo pipefail

cd "$(dirname "$0")/.."
VCFR=${VCFR:-./target/release/vcfr}
GEN=${GEN:-./target/release/examples/gen}
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

command -v bcftools >/dev/null || { echo "bcftools not found; skipping"; exit 0; }
[[ -x $VCFR ]] || { echo "build first: cargo build --release --example gen && cargo build --release"; exit 1; }

pass=0; fail=0
# Strip provenance headers and sort the INFO keys: key order is not meaningful
# in VCF and bcftools reorders them as a side effect of updating AC/AN.
body() {
  grep -v '^##' "$1" | awk -F'\t' 'BEGIN{OFS="\t"}
    /^#/ {print; next}
    { n = split($8, a, ";"); for (i = 1; i < n; i++) for (j = 1; j <= n-i; j++)
        if (a[j] > a[j+1]) { t = a[j]; a[j] = a[j+1]; a[j+1] = t }
      s = a[1]; for (i = 2; i <= n; i++) s = s ";" a[i]; $8 = s; print }'
}

check() { # check <name> <vcfr-output-file> <bcftools-output-file>
  if diff -q <(body "$2") <(body "$3") >/dev/null; then
    printf 'ok    %-34s %s records\n' "$1" "$(grep -vc '^#' "$2")"; pass=$((pass+1)); return 0
  else
    printf 'FAIL  %s\n' "$1"; diff <(body "$2") <(body "$3") | head -6; fail=$((fail+1))
  fi
}

echo "== generating fixtures in $WORK"
$GEN --sites 4000 --samples 25 --contigs 3                                            > "$WORK/all.vcf"
$GEN --sites 4000 --samples 10 --contigs 3 --sample-prefix A --gt-seed 11 --drop-pct 20 --drop-seed 3 > "$WORK/a.vcf"
$GEN --sites 4000 --samples 13 --contigs 3 --sample-prefix B --gt-seed 22 --drop-pct 20 --drop-seed 9 > "$WORK/b.vcf"
$GEN --sites 4000 --samples  7 --contigs 3 --sample-prefix C --gt-seed 33 --drop-pct 35 --drop-seed 5 > "$WORK/c.vcf"
for f in all a b c; do
  $VCFR view -O z -o "$WORK/$f.vcf.gz" "$WORK/$f.vcf"
  bcftools index -f "$WORK/$f.vcf.gz"
done

# ---------------------------------------------------------------- view -----
# vcfr streams, so -r is compared against the streaming bcftools flag (-t).
run_view() { # run_view <name> <vcfr args...> -- <bcftools args...>
  local name="$1"; shift
  local va=() ba=() cur=va
  for x in "$@"; do
    if [[ $x == "--" ]]; then cur=ba; continue; fi
    if [[ $cur == va ]]; then va+=("$x"); else ba+=("$x"); fi
  done
  $VCFR view "${va[@]}" "$WORK/all.vcf.gz" > "$WORK/v.out" 2>"$WORK/v.err"
  bcftools view "${ba[@]}" "$WORK/all.vcf.gz" > "$WORK/b.out" 2>"$WORK/b.err"
  check "$name" "$WORK/v.out" "$WORK/b.out"
}

run_view "view passthrough"        --  
run_view "view -s subset"          -s S00003,S00010,S00021 -- -s S00003,S00010,S00021
run_view "view -s exclude"         -s ^S00000,S00024 -- -s ^S00000,S00024
run_view "view -s reordered"       -s S00024,S00000,S00012 -- -s S00024,S00000,S00012
run_view "view -s single"          -s S00007 -- -s S00007
run_view "view -s --no-update"     -I -s S00007 -- -I -s S00007
run_view "view -v snps"            -v snps -- -v snps
run_view "view -v indels"          -v indels -- -v indels
run_view "view -V snps"            -V snps -- -V snps
run_view "view -m3 multiallelic"   -m 3 -- -m 3
run_view "view -m2 -M2 biallelic"  -m 2 -M 2 -- -m 2 -M 2
run_view "view -f PASS"            -f PASS -- -f PASS
run_view "view -G"                 -G -- -G
run_view "view -H"                 -H -- -H
run_view "view --min-qual"         --min-qual 500 -- -i "QUAL>=500"
run_view "view -r single contig"   -r chr2 -- -t chr2
run_view "view -r interval"        -r chr1:40000-120000 -- -t chr1:40000-120000
run_view "view -r multi"           -r chr1:40000-60000,chr3:10000-90000 -- -t chr1:40000-60000,chr3:10000-90000
run_view "view -r + -s + -v"       -r chr2 -s S00001,S00002 -v snps -- -t chr2 -s S00001,S00002 -v snps

printf 'chr1\t40000\t60000\nchr3\t10000\t90000\n' > "$WORK/regions.txt"
$VCFR view -R "$WORK/regions.txt" "$WORK/all.vcf.gz" > "$WORK/v.out"
bcftools view -T "$WORK/regions.txt" "$WORK/all.vcf.gz" > "$WORK/b.out"
check "view -R regions file" "$WORK/v.out" "$WORK/b.out"

# reading uncompressed input and stdin
$VCFR view -v snps "$WORK/all.vcf" > "$WORK/v.out"
bcftools view -v snps "$WORK/all.vcf" > "$WORK/b.out"
check "view plain input" "$WORK/v.out" "$WORK/b.out"
cat "$WORK/all.vcf.gz" | $VCFR view -v snps - > "$WORK/v.out"
check "view from stdin" "$WORK/v.out" "$WORK/b.out"

# BGZF written by vcfr must be readable by htslib, and vice versa
$VCFR view -O z -l 9 -o "$WORK/rt.vcf.gz" "$WORK/all.vcf"
bcftools view "$WORK/rt.vcf.gz" > "$WORK/v.out"
bcftools view "$WORK/all.vcf" > "$WORK/b.out"
check "BGZF written by vcfr" "$WORK/v.out" "$WORK/b.out"
bcftools index -f "$WORK/rt.vcf.gz" && echo "ok    vcfr BGZF is indexable"
# level 0 is vcfr's own DEFLATE encoder: htslib must read it back exactly
$VCFR view -O z -l 0 -o "$WORK/rt0.vcf.gz" "$WORK/all.vcf"
bcftools view "$WORK/rt0.vcf.gz" > "$WORK/v.out"
check "BGZF written by vcfr -l 0" "$WORK/v.out" "$WORK/b.out"
bcftools index -f "$WORK/rt0.vcf.gz" && echo "ok    vcfr -l 0 BGZF is indexable"
# --codec rust routes levels 1-6 to the high-effort Rust encoder: htslib must
# read those streams back exactly too
for rl in 1 6; do
  $VCFR view -O z --codec rust -l $rl -o "$WORK/rtr.vcf.gz" "$WORK/all.vcf"
  bcftools view "$WORK/rtr.vcf.gz" > "$WORK/v.out"
  check "BGZF written by vcfr --codec rust -l $rl" "$WORK/v.out" "$WORK/b.out"
done
bcftools index -f "$WORK/rtr.vcf.gz" && echo "ok    vcfr --codec rust BGZF is indexable"

bgzip -c "$WORK/all.vcf" > "$WORK/hts.vcf.gz"
$VCFR view "$WORK/hts.vcf.gz" > "$WORK/v.out"
check "BGZF written by bgzip" "$WORK/v.out" "$WORK/b.out"
gzip -c "$WORK/all.vcf" > "$WORK/plain.vcf.gz"
$VCFR view "$WORK/plain.vcf.gz" > "$WORK/v.out"
check "plain gzip input" "$WORK/v.out" "$WORK/b.out"

# -------------------------------------------------------------- concat -----
for c in chr1 chr2 chr3; do
  $VCFR view -r $c -O z -o "$WORK/part.$c.vcf.gz" "$WORK/all.vcf.gz"
  bcftools index -f "$WORK/part.$c.vcf.gz"
done
PARTS="$WORK/part.chr1.vcf.gz $WORK/part.chr2.vcf.gz $WORK/part.chr3.vcf.gz"
$VCFR concat $PARTS > "$WORK/v.out"
bcftools concat $PARTS 2>/dev/null > "$WORK/b.out"
check "concat 3 parts" "$WORK/v.out" "$WORK/b.out"
diff -q <(body "$WORK/v.out") <(body "$WORK/all.vcf") >/dev/null \
  && { echo "ok    concat round-trips the original"; pass=$((pass+1)); } \
  || { echo "FAIL  concat round-trip"; fail=$((fail+1)); }

$VCFR concat --naive -O z -o "$WORK/naive.vcf.gz" $PARTS
bcftools view "$WORK/naive.vcf.gz" > "$WORK/v.out"
check "concat --naive" "$WORK/v.out" "$WORK/b.out"
bcftools index -f "$WORK/naive.vcf.gz" && echo "ok    --naive output is indexable"

printf '%s\n' $PARTS > "$WORK/list.txt"
$VCFR concat -F "$WORK/list.txt" > "$WORK/v.out"
check "concat -F file list" "$WORK/v.out" "$WORK/b.out"

# --------------------------------------------------------------- merge -----
# Sites whose ALT sets disagree, exercising allele renumbering.
cat > "$WORK/x1.vcf" <<'EOF'
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="pass">
##FILTER=<ID=LowQual,Description="lq">
##FILTER=<ID=q10,Description="q10">
##contig=<ID=chr1,length=1000>
##INFO=<ID=AC,Number=A,Type=Integer,Description="ac">
##INFO=<ID=AN,Number=1,Type=Integer,Description="an">
##INFO=<ID=DP,Number=1,Type=Integer,Description="dp">
##INFO=<ID=EFF,Number=A,Type=String,Description="per-alt">
##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="flag">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="ad">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	X1	X2
chr1	100	.	A	G	50	PASS	AC=1;AN=4;DP=10;EFF=mis	GT:AD	0/1:5,6	0/0:9,1
chr1	200	rsA	C	T,G	60	PASS	AC=1,1;AN=4;DP=20;EFF=a,b	GT:AD	1/2:1,2,3	0/0:8,1,1
chr1	300	.	G	GT	70	q10	AC=2;AN=4;DP=30;SOMATIC	GT:AD	1|1:0,7	0/0:6,0
chr1	500	.	A	AT	70	PASS	AC=2;AN=4;DP=30	GT:AD	1/1:0,7	0/0:6,0
EOF
cat > "$WORK/x2.vcf" <<'EOF'
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="pass">
##FILTER=<ID=LowQual,Description="lq">
##FILTER=<ID=q10,Description="q10">
##contig=<ID=chr1,length=1000>
##INFO=<ID=AC,Number=A,Type=Integer,Description="ac">
##INFO=<ID=AN,Number=1,Type=Integer,Description="an">
##INFO=<ID=DP,Number=1,Type=Integer,Description="dp">
##INFO=<ID=EFF,Number=A,Type=String,Description="per-alt">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="ad">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="dp">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Y1
chr1	100	rsB	A	T	80	PASS	AC=2;AN=2;DP=7;EFF=syn	GT:AD:DP	1/1:0,9:9
chr1	200	.	C	G	90	LowQual	AC=1;AN=2;DP=8;EFF=z	GT:AD:DP	0/1:4,4:8
chr1	400	.	T	A	10	PASS	AC=0;AN=2;DP=9;EFF=nc	GT:AD:DP	0/0:9,0:9
chr1	500	.	A	G	70	PASS	AC=1;AN=2;DP=5;EFF=q	GT:AD:DP	0/1:3,3:6
EOF
for f in x1 x2; do bgzip -f -c "$WORK/$f.vcf" > "$WORK/$f.vcf.gz"; bcftools index -f "$WORK/$f.vcf.gz"; done
$VCFR merge "$WORK/a.vcf.gz" "$WORK/b.vcf.gz" > "$WORK/v.out"
bcftools merge "$WORK/a.vcf.gz" "$WORK/b.vcf.gz" > "$WORK/b.out"
check "merge 2 cohorts" "$WORK/v.out" "$WORK/b.out"

$VCFR merge "$WORK/a.vcf.gz" "$WORK/b.vcf.gz" "$WORK/c.vcf.gz" > "$WORK/v.out"
bcftools merge "$WORK/a.vcf.gz" "$WORK/b.vcf.gz" "$WORK/c.vcf.gz" > "$WORK/b.out"
check "merge 3 cohorts" "$WORK/v.out" "$WORK/b.out"

$VCFR merge -0 "$WORK/a.vcf.gz" "$WORK/b.vcf.gz" > "$WORK/v.out"
bcftools merge -0 "$WORK/a.vcf.gz" "$WORK/b.vcf.gz" > "$WORK/b.out"
check "merge -0 missing-to-ref" "$WORK/v.out" "$WORK/b.out"

# bcftools merge always recalculates AC/AN, so -I is vcfr-only: check that it
# leaves the first contributing file's values untouched.
$VCFR merge -I "$WORK/x1.vcf.gz" "$WORK/x2.vcf.gz" 2>/dev/null | grep -v "^#" | head -1 | cut -f8 > "$WORK/v.out"
grep -q "AN=4" "$WORK/v.out" \
  && { echo "ok    merge --no-update keeps the input AN"; pass=$((pass+1)); } \
  || { echo "FAIL  merge --no-update (AN was recalculated): $(cat "$WORK/v.out")"; fail=$((fail+1)); }

for m in none snps indels both all; do
  $VCFR merge -m $m "$WORK/x1.vcf.gz" "$WORK/x2.vcf.gz" > "$WORK/v.out"
  bcftools merge -m $m "$WORK/x1.vcf.gz" "$WORK/x2.vcf.gz" > "$WORK/b.out"
  check "merge -m $m" "$WORK/v.out" "$WORK/b.out"
done

$VCFR merge --info-rules DP:max "$WORK/a.vcf.gz" "$WORK/b.vcf.gz" > "$WORK/v.out"
bcftools merge --info-rules DP:max "$WORK/a.vcf.gz" "$WORK/b.vcf.gz" > "$WORK/b.out"
check "merge --info-rules DP:max" "$WORK/v.out" "$WORK/b.out"

$VCFR merge -O z -o "$WORK/m.vcf.gz" "$WORK/a.vcf.gz" "$WORK/b.vcf.gz"
bcftools view "$WORK/m.vcf.gz" > "$WORK/v.out"
bcftools merge "$WORK/a.vcf.gz" "$WORK/b.vcf.gz" > "$WORK/b.out"
check "merge -O z" "$WORK/v.out" "$WORK/b.out"

$VCFR merge "$WORK/x1.vcf.gz" "$WORK/x2.vcf.gz" > "$WORK/v.out"
bcftools merge "$WORK/x1.vcf.gz" "$WORK/x2.vcf.gz" > "$WORK/b.out"
check "merge allele renumbering" "$WORK/v.out" "$WORK/b.out"

# bcftools adds INFO/AN and INFO/AC only when an input header declares them
# (Beagle output, for one, declares neither); vcfr must follow the same rule.
for variant in undeclared declared; do
  extra=""
  [[ $variant == declared ]] && extra='##INFO=<ID=AC,Number=A,Type=Integer,Description="a">
##INFO=<ID=AN,Number=1,Type=Integer,Description="n">'
  for sm in M N; do
    {
      echo '##fileformat=VCFv4.2'
      echo '##contig=<ID=chr1>'
      echo '##INFO=<ID=DP,Number=1,Type=Integer,Description="d">'
      [[ -n $extra ]] && echo "$extra"
      echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="g">'
      printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s1\n' "$sm"
      printf 'chr1\t100\t.\tA\tG\t.\tPASS\tDP=7\tGT\t0|1\n'
    } > "$WORK/$sm.vcf"
    bgzip -f "$WORK/$sm.vcf" && bcftools index -f "$WORK/$sm.vcf.gz"
  done
  $VCFR merge "$WORK/M.vcf.gz" "$WORK/N.vcf.gz" > "$WORK/v.out"
  bcftools merge "$WORK/M.vcf.gz" "$WORK/N.vcf.gz" > "$WORK/b.out"
  check "merge AC/AN rule ($variant)" "$WORK/v.out" "$WORK/b.out"
done

# A downstream reader that closes early (`| head`) must exit vcfr cleanly: no
# "Broken pipe" on stderr, no error exit code. merge and concat --naive each
# had a write site whose io::Error was stringified before ignore_broken_pipe
# could see its ErrorKind, so the pipe closing surfaced as a spurious failure.
pipe_is_silent() { # pipe_is_silent <label> <cmd...>
  local label="$1"; shift
  "$@" 2>"$WORK/pipe.err" | head -c1 >/dev/null
  local rc=$?
  if [[ -s "$WORK/pipe.err" ]]; then
    printf 'FAIL  %-34s stderr: %s\n' "$label" "$(cat "$WORK/pipe.err")"; fail=$((fail+1))
  elif [[ $rc -ne 0 && $rc -ne 141 ]]; then
    printf 'FAIL  %-34s pipeline exit %s\n' "$label" "$rc"; fail=$((fail+1))
  else
    printf 'ok    %s\n' "$label"; pass=$((pass+1))
  fi
}
pipe_is_silent "view | head is silent"   $VCFR view "$WORK/all.vcf.gz"
pipe_is_silent "concat | head is silent" $VCFR concat $PARTS
pipe_is_silent "concat --naive | head is silent" $VCFR concat --naive -O z $PARTS
pipe_is_silent "merge | head is silent"  $VCFR merge "$WORK/a.vcf.gz" "$WORK/b.vcf.gz"

echo
echo "passed: $pass   failed: $fail"
[[ $fail -eq 0 ]]
