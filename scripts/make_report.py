#!/usr/bin/env python3
"""Assemble a Markdown run report from pipeline outputs (QC, calls, benchmark, tiers)."""
import argparse
import csv
import json
import os
import re
from collections import Counter, defaultdict
from datetime import date


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def md_table(headers, rows):
    out = ["| " + " | ".join(headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


def flagstat_metrics(path):
    m = {}
    with open(path) as fh:
        for line in fh:
            n = int(line.split()[0])
            if line.split()[3:5] == ["in", "total"]:
                m["total"] = n
            elif " primary mapped (" in line:
                m["mapped_pct"] = re.search(r"\(([\d.]+)%", line).group(1)
            elif " properly paired (" in line:
                m["paired_pct"] = re.search(r"\(([\d.]+)%", line).group(1)
            elif " primary duplicates" in line:
                m["dups"] = n
    return m


def cds_depth(path, bed):
    gene_of = []
    with open(bed) as fh:
        for line in fh:
            c, s, e, g = line.split()[:4]
            gene_of.append((c, int(s) + 1, int(e), g))
    depths = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            c, p, d = line.split()
            p = int(p)
            for gc, s, e, g in gene_of:
                if gc == c and s <= p <= e:
                    depths[g].append(int(d))
                    break
    rows = []
    for g, ds in sorted(depths.items()):
        rows.append([g, len(ds), f"{sum(ds) / len(ds):.1f}", min(ds),
                     f"{100 * sum(d >= 20 for d in ds) / len(ds):.1f}%"])
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    d, s = a.outdir, a.sample

    fs = flagstat_metrics(f"{d}/qc/{s}.flagstat.txt")
    cov = read_tsv(f"{d}/qc/{s}.coverage.tsv")
    bench = json.load(open(f"{d}/variants/{s}.benchmark.json"))
    ann = read_tsv(f"{d}/variants/{s}.annotated.tsv")
    versions = open(f"{d}/logs/software_versions.txt").read().strip()

    L = [f"# Germline variant report — {s}", "",
         f"_Generated {date.today().isoformat()} by `run_pipeline.sh`. "
         "Synthetic data for pipeline demonstration — not for clinical use._", ""]

    # --- QC
    L += ["## 1. Alignment QC", ""]
    c = cov[0]
    L.append(md_table(["Metric", "Value"], [
        ["Total reads", f"{fs['total']:,}"],
        ["Mapped (primary)", f"{fs['mapped_pct']}%"],
        ["Properly paired", f"{fs['paired_pct']}%"],
        ["Duplicates flagged", f"{fs.get('dups', 0):,}"],
        ["Genome breadth (≥1x)", f"{float(c['coverage']):.2f}%"],
        ["Mean depth", f"{float(c['meandepth']):.1f}x"],
        ["Mean MAPQ", c["meanmapq"]],
    ]))
    L += ["", "**Coding-exon coverage** (clinical reporting typically requires ≥20x across targets):", ""]
    L.append(md_table(["Gene", "CDS bp", "Mean depth", "Min depth", "% bases ≥20x"],
                      cds_depth(f"{d}/qc/{s}.cds_depth.tsv", f"{d}/qc/cds.bed")))

    # --- calls
    passed = [r for r in ann if r["filter"] == "PASS"]
    types = Counter("SNV" if len(r["ref"]) == len(r["alt"]) == 1 else "INDEL" for r in passed)
    zyg = Counter(r["zygosity"] for r in passed)
    L += ["", "## 2. Variant calls", "",
          f"{len(ann)} records after normalisation; **{len(passed)} PASS** "
          f"({types['SNV']} SNVs, {types['INDEL']} indels; {zyg['het']} het, {zyg['hom_alt']} hom-alt). "
          f"{len(ann) - len(passed)} failed filters (QUAL<30 or DP<10).", ""]
    csq = Counter(r["consequence"] for r in passed)
    L.append(md_table(["Consequence", "Count"], sorted(csq.items(), key=lambda x: -x[1])))

    # --- benchmark
    L += ["", "## 3. Benchmark against truth set", ""]
    L.append(md_table(
        ["Class", "TP", "FP", "FN", "Precision", "Recall", "F1", "Genotype concordance", "Annotation concordance"],
        [[t, b["TP"], b["FP"], b["FN"], f"{b['precision']:.3f}", f"{b['recall']:.3f}", f"{b['f1']:.3f}",
          f"{b['genotype_concordance']:.3f}", f"{b['annotation_concordance']:.3f}"] for t, b in bench.items()]))

    # --- prioritised
    cols = ["Gene", "HGVS c.", "HGVS p.", "Consequence", "Zygosity", "VAF (DP)", "KB", "Reason"]

    def rows_for(tier):
        return [[r["gene"], r["hgvs_c"], r["hgvs_p"], r["consequence"], r["zygosity"],
                 f"{r['vaf']} ({r['dp']})",
                 r["kb_clinsig"] if r["kb_clinsig"] != "." else "—", r["tier_reason"]]
                for r in ann if r["tier"] == tier]

    L += ["", "## 4. Prioritised variants", "",
          "Rule-based triage (simplified, ACMG-inspired): **Tier 1** reported P/LP or rare predicted "
          "loss-of-function; **Tier 2** rare protein-altering (VUS-like); **Tier 3** benign/common/low impact.", ""]
    for tier, title in (("Tier1", "Tier 1 — review first"), ("Tier2", "Tier 2 — variants of uncertain significance")):
        rows = rows_for(tier)
        L += [f"### {title} ({len(rows)})", ""]
        L += [md_table(cols, rows) if rows else "_None_", ""]
    t3 = [r for r in ann if r["tier"] == "Tier3"]
    reasons = Counter(re.sub(r"\(.*\)", "", r["tier_reason"]).strip() for r in t3)
    L += [f"### Tier 3 — not prioritised ({len(t3)})", ""]
    L += [f"- {k}: {v}" for k, v in reasons.most_common()]
    L += ["", f"Full annotation table: `variants/{s}.annotated.tsv`", ""]

    L += ["## 5. Software", "", "```", versions, "```", ""]
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w") as fh:
        fh.write("\n".join(L))


if __name__ == "__main__":
    main()
