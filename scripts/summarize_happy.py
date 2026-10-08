#!/usr/bin/env python3
"""Turn a hap.py <prefix>.summary.csv into a Markdown report and a compact JSON summary."""
import argparse
import csv
import json

FIELDS = [("TRUTH.TOTAL", "Truth"), ("TRUTH.TP", "TP"), ("TRUTH.FN", "FN"), ("QUERY.FP", "FP"),
          ("FP.gt", "FP (genotype)"), ("METRIC.Recall", "Recall"), ("METRIC.Precision", "Precision"),
          ("METRIC.F1_Score", "F1")]
RATES = {"METRIC.Recall", "METRIC.Precision", "METRIC.F1_Score"}


def read_summary(path):
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    out = {}
    for r in rows:
        out.setdefault(r["Type"], {})[r["Filter"]] = r
    return out


def fmt(key, val):
    if val in (None, "", "nan", "NaN"):
        return "—"
    if key in RATES:
        return f"{float(val):.4f}"
    return f"{int(float(val)):,}"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--summary", required=True, help="hap.py *.summary.csv")
    ap.add_argument("--out-md", required=True)
    ap.add_argument("--out-json", required=True)
    ap.add_argument("--meta", action="append", default=[], help="KEY=VALUE line for the run-details table")
    a = ap.parse_args()

    data = read_summary(a.summary)
    md = ["# GIAB HG002 benchmark (hap.py)", "", "## Run details", "", "| | |", "|---|---|"]
    md += [f"| {k} | {v} |" for k, v in (m.split("=", 1) for m in a.meta)]
    md += ["", "## Accuracy against NIST v4.2.1 confident regions", "",
           "| Type | Filter | " + " | ".join(lbl for _, lbl in FIELDS) + " |",
           "|---|---|" + "---|" * len(FIELDS)]
    summary = {}
    for vtype in ("SNP", "INDEL"):
        for filt in ("PASS", "ALL"):
            r = data.get(vtype, {}).get(filt)
            if not r:
                continue
            md.append(f"| {vtype} | {filt} | " + " | ".join(fmt(k, r.get(k)) for k, _ in FIELDS) + " |")
            summary.setdefault(vtype, {})[filt] = {
                k.split(".")[-1].lower(): (float(r[k]) if k in RATES else int(float(r[k])))
                for k, _ in FIELDS if r.get(k) not in (None, "", "nan", "NaN")}
    md += ["", "**ALL** = every call; **PASS** = after the `LowQual` soft filter. "
           "FP (genotype) counts calls at a true site with the wrong zygosity.", ""]
    with open(a.out_md, "w") as fh:
        fh.write("\n".join(md))
    with open(a.out_json, "w") as fh:
        json.dump(summary, fh, indent=2)
    for vtype, d in summary.items():
        p = d.get("PASS", d.get("ALL"))
        print(f"[hap.py] {vtype:5s} recall={p['recall']:.4f} precision={p['precision']:.4f} F1={p['f1_score']:.4f}")


if __name__ == "__main__":
    main()
