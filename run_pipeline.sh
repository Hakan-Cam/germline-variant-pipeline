#!/usr/bin/env bash
# =============================================================================
# Germline small-variant pipeline: FASTQ -> BAM -> VCF -> annotated, tiered report
#
#   1. simulate   synthetic reference, gene model, reads and truth set
#   2. index      bwa / samtools reference indexes
#   3. align      bwa mem -> fixmate -> sort -> markdup -> index
#   4. qc         flagstat, contig coverage, per-base depth across coding exons
#   5. call       bcftools mpileup | call | norm | filter
#   6. annotate   consequence + HGVS + knowledge-base lookup + tiering (Python)
#   7. benchmark  precision / recall / genotype + annotation concordance vs truth
#   8. report     Markdown summary report
#
# Usage: ./run_pipeline.sh [-o OUTDIR] [-s SAMPLE] [-t THREADS] [-c COVERAGE] [-r SEED] [-m MIN_F1]
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTDIR="results"
SAMPLE="SAMPLE01"
THREADS=2
COVERAGE=40
SEED=42
MIN_F1=0.95
# Hard-filter thresholds (soft filters: failing records are kept but labelled)
MIN_QUAL=30
MIN_DP=10

usage() { sed -n '2,15p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }

while getopts ":o:s:t:c:r:m:h" opt; do
  case "$opt" in
    o) OUTDIR="$OPTARG" ;;
    s) SAMPLE="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    c) COVERAGE="$OPTARG" ;;
    r) SEED="$OPTARG" ;;
    m) MIN_F1="$OPTARG" ;;
    h) usage 0 ;;
    *) echo "Unknown option: -$OPTARG" >&2; usage 1 ;;
  esac
done

DATA="$OUTDIR/data"
REF="$DATA/reference/chrS.fa"
GTF="$DATA/reference/genes.gtf"
ALN="$OUTDIR/alignment"
QC="$OUTDIR/qc"
VAR="$OUTDIR/variants"
LOG="$OUTDIR/logs"
mkdir -p "$ALN" "$QC" "$VAR" "$LOG"

log()  { printf '[%s] %s\n' "$(date '+%H:%M:%S')" "$*" | tee -a "$LOG/pipeline.log" >&2; }
step() { log "== $* =="; }
trap 'log "ERROR: pipeline failed at line $LINENO (exit $?)"' ERR

# ---- 0. dependencies --------------------------------------------------------
step "0/8 checking dependencies"
for tool in bwa samtools bcftools bgzip tabix python3 awk; do
  command -v "$tool" >/dev/null 2>&1 || { log "missing dependency: $tool"; exit 127; }
done
{
  echo "bwa       $(bwa 2>&1 | awk '/Version/{print $2}')"
  echo "samtools  $(samtools --version | head -1 | cut -d' ' -f2)"
  echo "bcftools  $(bcftools --version | head -1 | cut -d' ' -f2)"
  echo "python    $(python3 --version | cut -d' ' -f2)"
} > "$LOG/software_versions.txt"
sed 's/^/    /' "$LOG/software_versions.txt" >&2

# ---- 1. simulate ------------------------------------------------------------
step "1/8 simulating data (seed=$SEED, ${COVERAGE}x)"
python3 "$SCRIPT_DIR/scripts/simulate_data.py" --outdir "$DATA" --sample "$SAMPLE" \
  --seed "$SEED" --coverage "$COVERAGE" 2>&1 | tee -a "$LOG/pipeline.log" >&2

# ---- 2. index ---------------------------------------------------------------
step "2/8 indexing reference"
bwa index "$REF" 2> "$LOG/bwa_index.log"
samtools faidx "$REF"

# ---- 3. align ---------------------------------------------------------------
step "3/8 aligning reads (bwa mem, $THREADS threads)"
BAM="$ALN/$SAMPLE.markdup.bam"
RG="@RG\tID:${SAMPLE}.L001\tSM:${SAMPLE}\tLB:${SAMPLE}.lib1\tPL:ILLUMINA"
bwa mem -t "$THREADS" -R "$RG" "$REF" \
    "$DATA/reads/${SAMPLE}_R1.fastq.gz" "$DATA/reads/${SAMPLE}_R2.fastq.gz" 2> "$LOG/bwa_mem.log" \
  | samtools fixmate -m -u - - \
  | samtools sort -u -@ "$THREADS" -T "$ALN/tmp.$SAMPLE" - \
  | samtools markdup -@ "$THREADS" -f "$QC/$SAMPLE.markdup_stats.txt" - "$BAM"
samtools index "$BAM"

# ---- 4. qc ------------------------------------------------------------------
step "4/8 alignment QC"
samtools flagstat -@ "$THREADS" "$BAM" > "$QC/$SAMPLE.flagstat.txt"
samtools coverage "$BAM" > "$QC/$SAMPLE.coverage.tsv"
# Coding exons as BED (GTF is 1-based closed, BED is 0-based half-open)
awk -F'\t' 'BEGIN{OFS="\t"} $3=="CDS"{
    match($9, /gene_name "[^"]+"/); g=substr($9, RSTART+11, RLENGTH-12);
    print $1, $4-1, $5, g }' "$GTF" | sort -k1,1 -k2,2n > "$QC/cds.bed"
samtools depth -a -b "$QC/cds.bed" "$BAM" > "$QC/$SAMPLE.cds_depth.tsv"
sed -n '1p;7p' "$QC/$SAMPLE.flagstat.txt" | sed 's/^/    /' >&2

# ---- 5. call ----------------------------------------------------------------
step "5/8 calling variants (bcftools)"
RAW="$VAR/$SAMPLE.raw.vcf.gz"
CALLS="$VAR/$SAMPLE.filtered.vcf.gz"
bcftools mpileup -Ou -f "$REF" -a FORMAT/AD,FORMAT/DP -q 20 -Q 20 --max-depth 1000 "$BAM" 2> "$LOG/mpileup.log" \
  | bcftools call -mv -Oz -o "$RAW"
bcftools norm -f "$REF" -m -any -Ou "$RAW" 2> "$LOG/norm.log" \
  | bcftools filter -s LowQual -e "QUAL<$MIN_QUAL || INFO/DP<$MIN_DP" -Oz -o "$CALLS"
tabix -f -p vcf "$CALLS"
bcftools stats "$CALLS" > "$QC/$SAMPLE.bcftools_stats.txt"
n_pass=$(bcftools view -H -f PASS "$CALLS" | wc -l)
n_all=$(bcftools view -H "$CALLS" | wc -l)
log "    $n_all records, $n_pass PASS"

# Normalise the truth set exactly like the calls so representations match
TRUTH="$VAR/$SAMPLE.truth.norm.vcf.gz"
bcftools norm -f "$REF" -m -any "$DATA/truth/$SAMPLE.truth.vcf" -Oz -o "$TRUTH" 2>> "$LOG/norm.log"

# ---- 6. annotate ------------------------------------------------------------
step "6/8 annotating and prioritising"
python3 "$SCRIPT_DIR/scripts/annotate_variants.py" --vcf "$CALLS" --fasta "$REF" --gtf "$GTF" \
  --kb "$DATA/kb/known_variants.tsv" --out "$VAR/$SAMPLE.annotated.tsv"

# ---- 7. benchmark -----------------------------------------------------------
step "7/8 benchmarking against truth (min F1 = $MIN_F1)"
python3 "$SCRIPT_DIR/scripts/benchmark.py" --truth-vcf "$TRUTH" \
  --truth-table "$DATA/truth/truth_variants.tsv" --calls-vcf "$CALLS" \
  --annotated "$VAR/$SAMPLE.annotated.tsv" \
  --out-json "$VAR/$SAMPLE.benchmark.json" --out-tsv "$VAR/$SAMPLE.benchmark.tsv" \
  --min-f1 "$MIN_F1"

# ---- 8. report --------------------------------------------------------------
step "8/8 writing report"
python3 "$SCRIPT_DIR/scripts/make_report.py" --sample "$SAMPLE" --outdir "$OUTDIR" \
  --out "$OUTDIR/$SAMPLE.report.md"

log "done -> $OUTDIR/$SAMPLE.report.md"
