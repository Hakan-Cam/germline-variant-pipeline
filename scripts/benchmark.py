#!/usr/bin/env python3
"""
Benchmark variant calls and annotations against the simulated truth set.

Reports, separately for SNVs and indels:
  * TP / FP / FN, precision, recall and F1 (PASS calls vs truth, exact allele match
    after both VCFs are left-aligned and normalised with bcftools norm),
  * genotype concordance for true positives (het vs hom-alt),
  * annotation concordance: predicted consequence vs the designed consequence.

Writes a JSON summary and a per-variant TSV. With --min-f1 the script exits
non-zero if any variant class falls below the threshold, so it can gate CI.
"""
import argparse
import csv
import gzip
import json
import re
import sys


def zyg(gt):
    a = re.split(r"[/|]", gt)
    if len(set(a)) == 1 and a[0] not in ("0", "."):
        return "hom_alt"
    return "het" if "0" in a and len(set(a)) > 1 else "other"


def vtype(ref, alt):
    return "SNV" if len(ref) == len(alt) == 1 else "INDEL"


def read_vcf(path, pass_only=False):
    out = {}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if pass_only and f[6] not in ("PASS", "."):
                continue
            gt = dict(zip(f[8].split(":"), f[9].split(":"))).get("GT", "./.")
            for alt in f[4].split(","):
                out[(f[0], int(f[1]), f[3], alt)] = dict(id=f[2], gt=gt)
    return out


def safe_div(a, b):
    return round(a / b, 4) if b else 0.0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth-vcf", required=True, help="normalised truth VCF")
    ap.add_argument("--truth-table", required=True, help="truth_variants.tsv (expected consequences)")
    ap.add_argument("--calls-vcf", required=True, help="normalised, filtered calls")
    ap.add_argument("--annotated", required=True, help="annotate_variants.py output")
    ap.add_argument("--out-json", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--min-f1", type=float, default=None)
    args = ap.parse_args()

    truth = read_vcf(args.truth_vcf)
    calls = read_vcf(args.calls_vcf, pass_only=True)
    with open(args.truth_table) as fh:
        expected = {r["id"]: r["expected_consequence"] for r in csv.DictReader(fh, delimiter="\t")}
    with open(args.annotated) as fh:
        observed = {(r["chrom"], int(r["pos"]), r["ref"], r["alt"]): r["consequence"]
                    for r in csv.DictReader(fh, delimiter="\t")}

    stats = {t: dict(TP=0, FP=0, FN=0, gt_match=0, csq_match=0) for t in ("SNV", "INDEL")}
    rows = []
    for key in sorted(set(truth) | set(calls), key=lambda k: (k[0], k[1])):
        t = vtype(key[2], key[3])
        tv, cv = truth.get(key), calls.get(key)
        status = "TP" if tv and cv else ("FN" if tv else "FP")
        stats[t][status] += 1
        gt_ok = csq_ok = ""
        exp_csq = expected.get(tv["id"], ".") if tv else "."
        obs_csq = observed.get(key, ".")
        if status == "TP":
            gt_ok = zyg(tv["gt"]) == zyg(cv["gt"])
            csq_ok = exp_csq == obs_csq
            stats[t]["gt_match"] += gt_ok
            stats[t]["csq_match"] += csq_ok
        rows.append(dict(id=tv["id"] if tv else ".", chrom=key[0], pos=key[1], ref=key[2], alt=key[3],
                         type=t, status=status, truth_gt=tv["gt"] if tv else ".",
                         call_gt=cv["gt"] if cv else ".", genotype_match=gt_ok,
                         expected_consequence=exp_csq, observed_consequence=obs_csq,
                         consequence_match=csq_ok))

    summary, failed = {}, False
    for t, s in stats.items():
        p, r = safe_div(s["TP"], s["TP"] + s["FP"]), safe_div(s["TP"], s["TP"] + s["FN"])
        f1 = safe_div(2 * p * r, p + r)
        summary[t] = dict(TP=s["TP"], FP=s["FP"], FN=s["FN"], precision=p, recall=r, f1=f1,
                          genotype_concordance=safe_div(s["gt_match"], s["TP"]),
                          annotation_concordance=safe_div(s["csq_match"], s["TP"]))
        if args.min_f1 is not None and f1 < args.min_f1:
            failed = True

    with open(args.out_json, "w") as fh:
        json.dump(summary, fh, indent=2)
    with open(args.out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    for t, s in summary.items():
        print(f"[benchmark] {t:5s} TP={s['TP']:<3} FP={s['FP']:<3} FN={s['FN']:<3} "
              f"precision={s['precision']:.3f} recall={s['recall']:.3f} F1={s['f1']:.3f} "
              f"GT-concordance={s['genotype_concordance']:.3f} "
              f"annotation-concordance={s['annotation_concordance']:.3f}", file=sys.stderr)
    if failed:
        print(f"[benchmark] FAIL: F1 below --min-f1 {args.min_f1}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
