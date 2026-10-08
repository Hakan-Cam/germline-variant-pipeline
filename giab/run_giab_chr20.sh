#!/usr/bin/env bash
# =============================================================================
# Benchmark the pipeline on real data: GIAB HG002 (NA24385), GRCh38 chr20,
# scored with hap.py against the NIST v4.2.1 small-variant benchmark.
#
#   1. fetch      chr20 reference, truth VCF + confident-region BED (cached)
#   2. reads      stream HG002 35x NovaSeq PCR-free reads for REGION -> paired FASTQ
#   3. align      same align_reads() as the synthetic pipeline (bwa mem ... markdup)
#   4. qc         flagstat + depth over REGION
#   5. call       same call_variants() (parallel scatter/gather bcftools)
#   6. hap.py     docker (jmcdani20/hap.py:v0.3.12, vcfeval engine) or a native install
#   7. report     Markdown + JSON summary of SNP / INDEL precision, recall, F1
#
# Usage: giab/run_giab_chr20.sh [-R REGION] [-t THREADS] [-o OUTDIR] [-c CHUNK_BP]
#                               [-H auto|docker|native|skip]
#                               [-b BAM] [-f REF_FASTA] [-v TRUTH_VCF] [-B TRUTH_BED]
#   -R  region to benchmark, default chr20 (all 64 Mb); e.g. chr20:10000000-15000000 for a quick run
#   -b/-f/-v/-B override the default URLs with local files (useful offline or for other samples)
# =============================================================================
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=../lib/steps.sh
source "$ROOT/lib/steps.sh"

# ---- defaults ---------------------------------------------------------------
GIAB_BASE="https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38"
TRUTH_VCF="$GIAB_BASE/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED="$GIAB_BASE/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
BAM="https://storage.googleapis.com/deepvariant/case-study-testdata/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam"
REF="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr20.fa.gz"
HAPPY_IMAGE="jmcdani20/hap.py:v0.3.12"
REGION="chr20"
OUTDIR="giab_results"
THREADS=$(( $(nproc) > 16 ? 16 : $(nproc) ))
CHUNK=5000000
HAPPY_MODE="auto"
SAMPLE="HG002"
MIN_QUAL=30
MIN_DP=10

usage() { sed -n '2,20p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }
while getopts ":R:t:o:c:H:b:f:v:B:h" opt; do
  case "$opt" in
    R) REGION="$OPTARG" ;;  t) THREADS="$OPTARG" ;;  o) OUTDIR="$OPTARG" ;;
    c) CHUNK="$OPTARG" ;;   H) HAPPY_MODE="$OPTARG" ;;
    b) BAM="$OPTARG" ;;     f) REF="$OPTARG" ;;  v) TRUTH_VCF="$OPTARG" ;;  B) TRUTH_BED="$OPTARG" ;;
    h) usage 0 ;;           *) echo "Unknown option: -$OPTARG" >&2; usage 1 ;;
  esac
done

CHROM="${REGION%%:*}"
is_url() { [[ "$1" =~ ^(https?|ftp|gs|s3):// ]]; }
is_url "$BAM" || BAM="$(cd "$(dirname "$BAM")" && pwd)/$(basename "$BAM")"   # stays valid after cd
mkdir -p "$OUTDIR"
OUTDIR="$(cd "$OUTDIR" && pwd)"
DL="$OUTDIR/downloads"; REFDIR="$OUTDIR/reference"; TRUTH="$OUTDIR/truth"
READS="$OUTDIR/reads"; ALN="$OUTDIR/alignment"; VAR="$OUTDIR/variants"
QC="$OUTDIR/qc"; LOG="$OUTDIR/logs"; HAPPY="$OUTDIR/happy"
mkdir -p "$DL" "$REFDIR" "$TRUTH" "$READS" "$ALN" "$VAR" "$QC" "$LOG" "$HAPPY"
TAG="$(echo "$REGION" | tr ':-' '__')"

log()  { printf '[%s] %s\n' "$(date '+%H:%M:%S')" "$*" | tee -a "$LOG/giab.log" >&2; }
step() { log "== $* =="; }
trap 'log "ERROR: failed at line $LINENO (exit $?)"' ERR

# fetch SRC DEST : download a URL once (resumable, cached) or link a local file
fetch() {
  local src=$1 dest=$2
  [[ -s "$dest" ]] && { log "    cached  $(basename "$dest")"; return; }
  if is_url "$src"; then
    log "    download $(basename "$src")"
    curl -fL --retry 3 --retry-delay 5 -C - -o "$dest.part" "$src"
    mv "$dest.part" "$dest"
  else
    [[ -s "$src" ]] || { log "missing input file: $src"; exit 1; }
    ln -sf "$(cd "$(dirname "$src")" && pwd)/$(basename "$src")" "$dest"
  fi
}

step "0/7 checking dependencies (region $REGION, $THREADS threads)"
for tool in bwa samtools bcftools tabix bgzip curl python3 awk; do
  command -v "$tool" >/dev/null 2>&1 || { log "missing dependency: $tool"; exit 127; }
done
START=$SECONDS

# ---- 1. reference, truth VCF, confident regions ----------------------------
step "1/7 fetching reference and GIAB v4.2.1 benchmark"
REF_FA="$REFDIR/$CHROM.fa"
if [[ ! -s "$REF_FA.bwt" ]]; then
  if is_url "$REF"; then
    [[ "$CHROM" == "chr20" ]] || { log "default reference is chr20 only; pass -f for $CHROM"; exit 1; }
    fetch "$REF" "$DL/$(basename "$REF")"
    gunzip -c "$DL/$(basename "$REF")" > "$REF_FA"
  else
    samtools faidx "$REF" "$CHROM" > "$REF_FA"
  fi
  samtools faidx "$REF_FA"
  log "    bwa index $CHROM"
  bwa index "$REF_FA" 2> "$LOG/bwa_index.log"
fi

fetch "$TRUTH_VCF" "$DL/truth.vcf.gz"
is_url "$TRUTH_VCF" && fetch "$TRUTH_VCF.tbi" "$DL/truth.vcf.gz.tbi"
[[ -s "$DL/truth.vcf.gz.tbi" ]] || tabix -f -p vcf "$DL/truth.vcf.gz"
fetch "$TRUTH_BED" "$DL/truth.bed"

TRUTH_REGION_VCF="$TRUTH/HG002.truth.$TAG.vcf.gz"
TRUTH_REGION_BED="$TRUTH/HG002.confident.$TAG.bed"
bcftools view -r "$REGION" -Oz -o "$TRUTH_REGION_VCF" "$DL/truth.vcf.gz"
tabix -f -p vcf "$TRUTH_REGION_VCF"
# clip the confident regions to REGION (BED is 0-based half-open)
awk -v region="$REGION" 'BEGIN{OFS="\t"; n=split(region,a,/[:-]/); c=a[1]; s=(n>=3?a[2]-1:0); e=(n>=3?a[3]:1e12)}
  $1==c && $3>s && $2<e { print $1, ($2<s?s:$2), ($3>e?e:$3) }' "$DL/truth.bed" > "$TRUTH_REGION_BED"
CONF_BP=$(awk '{s+=$3-$2} END{print s+0}' "$TRUTH_REGION_BED")
N_TRUTH=$(bcftools view -H "$TRUTH_REGION_VCF" | wc -l)
log "    truth: $N_TRUTH variants, $CONF_BP bp in confident regions"

# ---- 2. reads ---------------------------------------------------------------
step "2/7 extracting HG002 read pairs for $REGION"
R1="$READS/$SAMPLE.$TAG.R1.fastq.gz"; R2="$READS/$SAMPLE.$TAG.R2.fastq.gz"
if [[ ! -s "$R2" ]]; then
  # Streams only the requested region over HTTPS (htslib fetches the .bai);
  # drops secondary/supplementary records and pairs whose mate lies outside REGION.
  (cd "$READS" && samtools view -u -F 0x900 "$BAM" "$REGION") \
    | samtools collate -u -O -@ "$THREADS" - "$READS/collate.$TAG" \
    | samtools fastq -n -1 "$R1" -2 "$R2" -0 /dev/null -s /dev/null - 2> "$LOG/fastq.log"
fi
N_PAIRS=$(( $(gzip -cd "$R1" | wc -l) / 4 ))
log "    $N_PAIRS read pairs"

# ---- 3. align ---------------------------------------------------------------
step "3/7 aligning to $CHROM (bwa mem, $THREADS threads)"
OUT_BAM="$ALN/$SAMPLE.$TAG.markdup.bam"
align_reads "$REF_FA" "$R1" "$R2" "$SAMPLE" "$OUT_BAM" "$THREADS" "$QC/markdup_stats.txt" "$LOG/bwa_mem.log"

# ---- 4. qc ------------------------------------------------------------------
step "4/7 alignment QC"
samtools flagstat -@ "$THREADS" "$OUT_BAM" > "$QC/flagstat.txt"
samtools coverage -r "$REGION" "$OUT_BAM" > "$QC/coverage.tsv"
MEAN_DEPTH=$(awk 'NR==2{printf "%.1f", $7}' "$QC/coverage.tsv")
log "    mean depth over $REGION: ${MEAN_DEPTH}x"

# ---- 5. call ----------------------------------------------------------------
step "5/7 calling variants (bcftools, scatter over ${CHUNK} bp windows)"
RAW="$VAR/$SAMPLE.$TAG.raw.vcf.gz"; CALLS="$VAR/$SAMPLE.$TAG.filtered.vcf.gz"
call_variants "$REF_FA" "$OUT_BAM" "$RAW" "$CALLS" "$THREADS" "$MIN_QUAL" "$MIN_DP" "$LOG" "$CHUNK" "$REGION"
log "    $(bcftools view -H "$CALLS" | wc -l) records, $(bcftools view -H -f PASS "$CALLS" | wc -l) PASS"

# ---- 6. hap.py --------------------------------------------------------------
step "6/7 scoring with hap.py"
if [[ "$HAPPY_MODE" == "auto" ]]; then
  if command -v hap.py >/dev/null 2>&1; then HAPPY_MODE=native
  elif command -v docker >/dev/null 2>&1 && docker info >/dev/null 2>&1; then HAPPY_MODE=docker
  else HAPPY_MODE=skip; fi
fi
rel() { echo "${1#"$OUTDIR"/}"; }   # path relative to OUTDIR (mounted as /data in docker)
HAPPY_ARGS=(-f "$(rel "$TRUTH_REGION_BED")" -r "$(rel "$REF_FA")" -o "happy/$SAMPLE.$TAG"
            -l "$CHROM" --threads "$THREADS")
case "$HAPPY_MODE" in
  docker)
    log "    docker image $HAPPY_IMAGE (vcfeval engine)"
    docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v "$OUTDIR":/data -w /data "$HAPPY_IMAGE" \
      /opt/hap.py/bin/hap.py "$(rel "$TRUTH_REGION_VCF")" "$(rel "$CALLS")" "${HAPPY_ARGS[@]}" \
      --engine=vcfeval > "$LOG/happy.log" 2>&1 ;;
  native)
    engine=(); command -v rtg >/dev/null 2>&1 && engine=(--engine=vcfeval)
    log "    native hap.py (${engine[*]:-xcmp engine})"
    (cd "$OUTDIR" && hap.py "$(rel "$TRUTH_REGION_VCF")" "$(rel "$CALLS")" "${HAPPY_ARGS[@]}" \
      "${engine[@]}") > "$LOG/happy.log" 2>&1 ;;
  skip)
    log "    hap.py not available (install Docker or 'conda install -c bioconda hap.py'); stopping after calling"
    exit 0 ;;
  *) log "unknown -H mode: $HAPPY_MODE"; exit 1 ;;
esac

# ---- 7. report --------------------------------------------------------------
step "7/7 writing report"
python3 "$ROOT/scripts/summarize_happy.py" \
  --summary "$HAPPY/$SAMPLE.$TAG.summary.csv" \
  --out-md "$OUTDIR/$SAMPLE.$TAG.happy_report.md" --out-json "$OUTDIR/$SAMPLE.$TAG.happy_summary.json" \
  --meta "Sample=HG002 (NA24385), GIAB Ashkenazi son" \
  --meta "Region=$REGION" \
  --meta "Reads=NovaSeq PCR-free 2x150, 35x WGS (PrecisionFDA Truth v2); $N_PAIRS pairs extracted" \
  --meta "Mean depth after realignment=${MEAN_DEPTH}x" \
  --meta "Truth set=NIST v4.2.1 GRCh38 ($N_TRUTH variants in region, $CONF_BP bp confident)" \
  --meta "Caller=bcftools $(bcftools --version | head -1 | cut -d' ' -f2) mpileup/call -mv; soft filter QUAL<$MIN_QUAL or DP<$MIN_DP" \
  --meta "Scoring=hap.py ($HAPPY_MODE)" \
  --meta "Runtime=$(( (SECONDS - START) / 60 )) min on $THREADS threads"
log "done -> $OUTDIR/$SAMPLE.$TAG.happy_report.md"
