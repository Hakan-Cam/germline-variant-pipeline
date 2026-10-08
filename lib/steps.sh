#!/usr/bin/env bash
# Shared pipeline steps, sourced by run_pipeline.sh and giab/run_giab_chr20.sh.
# Keeping alignment and calling in one place guarantees the synthetic test and the
# GIAB benchmark exercise exactly the same commands and parameters.

# align_reads REF R1 R2 SAMPLE OUT_BAM THREADS MARKDUP_STATS LOGFILE
#   bwa mem (with read group) -> fixmate -m -> coordinate sort -> markdup -> index
align_reads() {
  local ref=$1 r1=$2 r2=$3 sample=$4 out_bam=$5 threads=$6 mdstats=$7 logfile=$8
  local rg="@RG\tID:${sample}.L001\tSM:${sample}\tLB:${sample}.lib1\tPL:ILLUMINA"
  bwa mem -t "$threads" -R "$rg" "$ref" "$r1" "$r2" 2> "$logfile" \
    | samtools fixmate -m -u - - \
    | samtools sort -u -@ "$threads" -m 1G -T "${out_bam%.bam}.sorttmp" - \
    | samtools markdup -@ "$threads" -f "$mdstats" - "$out_bam"
  samtools index -@ "$threads" "$out_bam"
}

# make_chunks FASTA_FAI CHUNK_BP [REGION]
#   Print 1-based "chrom:start-end" windows covering REGION ("chr", "chr:a-b"),
#   or every contig in the .fai when REGION is empty.
make_chunks() {
  local fai=$1 chunk=$2 region=${3:-}
  awk -v chunk="$chunk" -v region="$region" 'BEGIN{FS=OFS="\t"
      if (region != "") { n = split(region, a, /[:-]/); rc = a[1]; rs = (n >= 3 ? a[2] : 1); re = (n >= 3 ? a[3] : 0) } }
    {
      if (region != "" && $1 != rc) next
      s = (region != "" ? rs : 1); e = (region != "" && re > 0 ? re : $2)
      if (e > $2) e = $2
      for (p = s; p <= e; p += chunk) { q = p + chunk - 1; if (q > e) q = e; print $1 ":" p "-" q }
    }' "$fai"
}

# _call_chunk REF BAM TMPDIR INDEX REGION   (run in parallel by call_variants)
_call_chunk() {
  set -euo pipefail
  local ref=$1 bam=$2 tmp=$3 idx=$4 region=$5
  local out; out=$(printf '%s/chunk.%05d.vcf.gz' "$tmp" "$idx")
  bcftools mpileup -Ou -f "$ref" -r "$region" -a FORMAT/AD,FORMAT/DP -q 20 -Q 20 \
      --max-depth 1000 "$bam" 2>> "$tmp/mpileup.log" \
    | bcftools call -mv -Oz -o "$out" 2>> "$tmp/call.log"
}

# call_variants REF BAM RAW_VCF FILTERED_VCF THREADS MIN_QUAL MIN_DP LOGDIR [CHUNK_BP] [REGION]
#   Scatter mpileup|call over genomic windows in parallel, gather, then
#   normalise (left-align, split multi-allelics) and soft-filter.
call_variants() {
  local ref=$1 bam=$2 raw=$3 filtered=$4 threads=$5 min_qual=$6 min_dp=$7 logdir=$8
  local chunk=${9:-5000000} region=${10:-}
  local tmp; tmp=$(mktemp -d "${raw%.vcf.gz}.chunks.XXXXXX")
  make_chunks "$ref.fai" "$chunk" "$region" > "$tmp/regions.txt"
  [[ -s "$tmp/regions.txt" ]] || { echo "no regions to call (region='$region')" >&2; return 1; }

  export -f _call_chunk
  awk '{print NR, $1}' "$tmp/regions.txt" \
    | xargs -P "$threads" -L 1 bash -c '_call_chunk "$0" "$1" "$2" "$3" "$4"' "$ref" "$bam" "$tmp"

  bcftools concat -Oz -o "$raw" "$tmp"/chunk.*.vcf.gz 2> "$logdir/concat.log"
  bcftools norm -f "$ref" -m -any -Ou "$raw" 2> "$logdir/norm.log" \
    | bcftools filter -s LowQual -e "QUAL<$min_qual || INFO/DP<$min_dp" -Oz -o "$filtered"
  tabix -f -p vcf "$filtered"
  cat "$tmp/mpileup.log" > "$logdir/mpileup.log" 2>/dev/null || true
  rm -rf "$tmp"
}
