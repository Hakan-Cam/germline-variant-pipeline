# Real-data benchmark: GIAB HG002 chr20 + hap.py

The synthetic test proves the pipeline logic; this module measures how the **same alignment and calling steps** (`lib/steps.sh`) perform on a real human genome, using the field-standard approach: Genome in a Bottle truth set + Illumina's `hap.py`.

| Component | Source |
|---|---|
| Sample | HG002 / NA24385 (GIAB Ashkenazi son) |
| Reads | NovaSeq PCR-free 2x150, 35x WGS (PrecisionFDA Truth v2), chr20 BAM hosted by the DeepVariant team — streamed by region, converted back to FASTQ and **re-aligned from scratch** |
| Reference | GRCh38 chr20 (UCSC) |
| Truth | NIST v4.2.1 GRCh38 benchmark VCF + `noinconsistent` confident-region BED |
| Scoring | `hap.py` v0.3.12 with the `vcfeval` engine (Docker), or a native install |

## Run it

```bash
# quick test, 5 Mb (a few minutes)
giab/run_giab_chr20.sh -R chr20:10000000-15000000

# full chromosome 20 (64 Mb)
giab/run_giab_chr20.sh -R chr20 -t 16
```

Requirements: `bwa samtools bcftools tabix curl` plus **one** of

- Docker (on WSL2: Docker Desktop with WSL integration) — the script pulls `jmcdani20/hap.py:v0.3.12`, or
- a conda env with hap.py: `mamba create -n happy -c conda-forge -c bioconda hap.py rtg-tools` then `conda activate happy`.

Without either, use `-H skip` to stop after variant calling.

No local setup? Run it on GitHub: **Actions → giab-hg002-benchmark → Run workflow**. The report appears on the run's summary page and is uploaded as an artifact.

## Outputs (`giab_results/`)

```
HG002.<region>.happy_report.md     run details + SNP/INDEL precision, recall, F1 (ALL and PASS)
HG002.<region>.happy_summary.json
happy/                             raw hap.py output (summary.csv, extended.csv, annotated VCF)
qc/  logs/  variants/  alignment/
downloads/                         cached reference / truth files (reused on re-runs)
```

## Design notes

- **Re-alignment, not re-use:** the source BAM is converted to FASTQ so the benchmark covers `bwa mem → markdup → bcftools`, i.e. the whole pipeline.
- **Scatter/gather calling:** `call_variants` splits the region into windows (`-c`, default 5 Mb) and runs `mpileup | call` in parallel, then concatenates, normalises and filters.
- **ALL vs PASS rows:** hap.py is run without `--pass-only`, so the report shows what the `QUAL<30 or DP<10` soft filter gains in precision and costs in recall.
- **Known simplifications:** reads are aligned to chr20 alone (no decoys/other chromosomes, so some reads from paralogous regions can mis-map), mates outside the region are dropped, and no stratification by difficult regions is applied yet (`hap.py --stratification` with the GIAB v3 stratifications is the natural next step).
