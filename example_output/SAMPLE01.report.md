# Germline variant report — SAMPLE01

_Generated 2026-10-08 by `run_pipeline.sh`. Synthetic data for pipeline demonstration — not for clinical use._

## 1. Alignment QC

| Metric | Value |
|---|---|
| Total reads | 2,844 |
| Mapped (primary) | 100.00% |
| Properly paired | 100.00% |
| Duplicates flagged | 2 |
| Genome breadth (≥1x) | 99.58% |
| Mean depth | 39.9x |
| Mean MAPQ | 60 |

**Coding-exon coverage** (clinical reporting typically requires ≥20x across targets):

| Gene | CDS bp | Mean depth | Min depth | % bases ≥20x |
|---|---|---|---|---|
| GENE_A | 501 | 35.6 | 19 | 99.6% |
| GENE_B | 423 | 46.3 | 32 | 100.0% |

## 2. Variant calls

19 records after normalisation; **19 PASS** (14 SNVs, 5 indels; 16 het, 3 hom-alt). 0 failed filters (QUAL<30 or DP<10).

| Consequence | Count |
|---|---|
| intergenic_variant | 5 |
| missense_variant | 3 |
| stop_gained | 2 |
| frameshift_variant | 2 |
| synonymous_variant | 2 |
| start_lost | 1 |
| splice_donor_variant | 1 |
| splice_acceptor_variant | 1 |
| inframe_deletion | 1 |
| intron_variant | 1 |

## 3. Benchmark against truth set

| Class | TP | FP | FN | Precision | Recall | F1 | Genotype concordance | Annotation concordance |
|---|---|---|---|---|---|---|---|---|
| SNV | 14 | 0 | 0 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| INDEL | 5 | 0 | 0 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |

## 4. Prioritised variants

Rule-based triage (simplified, ACMG-inspired): **Tier 1** reported P/LP or rare predicted loss-of-function; **Tier 2** rare protein-altering (VUS-like); **Tier 3** benign/common/low impact.

### Tier 1 — review first (8)

| Gene | HGVS c. | HGVS p. | Consequence | Zygosity | VAF (DP) | KB | Reason |
|---|---|---|---|---|---|---|---|
| GENE_A | c.2T>C | p.Met1? | start_lost | het | 0.44 (41) | — | Predicted loss-of-function; rare/absent in population |
| GENE_A | c.118C>T | p.Arg40Ter | stop_gained | het | 0.56 (34) | Pathogenic | Reported Pathogenic in knowledge base |
| GENE_A | c.155+1G>A | p.? | splice_donor_variant | het | 0.55 (38) | — | Predicted loss-of-function; rare/absent in population |
| GENE_A | c.310del | p.Gly104Alafs | frameshift_variant | het | 0.44 (36) | — | Predicted loss-of-function; rare/absent in population |
| GENE_B | c.300_301insA | p.Leu101Thrfs | frameshift_variant | het | 0.45 (49) | — | Predicted loss-of-function; rare/absent in population |
| GENE_B | c.268C>T | p.Gln90Ter | stop_gained | hom_alt | 1.00 (50) | — | Predicted loss-of-function; rare/absent in population |
| GENE_B | c.191-2A>G | p.? | splice_acceptor_variant | het | 0.48 (46) | — | Predicted loss-of-function; rare/absent in population |
| GENE_B | c.164A>G | p.Lys55Arg | missense_variant | het | 0.57 (40) | Likely_pathogenic | Reported Likely_pathogenic in knowledge base |

### Tier 2 — variants of uncertain significance (2)

| Gene | HGVS c. | HGVS p. | Consequence | Zygosity | VAF (DP) | KB | Reason |
|---|---|---|---|---|---|---|---|
| GENE_A | c.178C>T | p.Arg60Trp | missense_variant | hom_alt | 1.00 (42) | — | Rare protein-altering variant; needs review (VUS-like) |
| GENE_A | c.399_401del | p.Ser134del | inframe_deletion | het | 0.49 (37) | — | Rare protein-altering variant; needs review (VUS-like) |

### Tier 3 — not prioritised (9)

- Low predicted impact: 7
- Reported Benign in knowledge base: 2

Full annotation table: `variants/SAMPLE01.annotated.tsv`

## 5. Software

```
bwa       0.7.17-r1188
samtools  1.19.2
bcftools  1.19
python    3.13.16
```
