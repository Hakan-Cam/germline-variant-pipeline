#!/usr/bin/env python3
"""
Annotate and prioritise variants from a normalised VCF.

For every variant this script:
  1. finds overlapping transcripts in a GTF (CDS features),
  2. predicts the consequence from sequence (VEP/Sequence Ontology terms):
     synonymous, missense, stop_gained, stop_lost, start_lost, frameshift,
     in-frame insertion/deletion, splice donor/acceptor/region, intronic, intergenic,
  3. writes simplified HGVS c./p. notation,
  4. joins a knowledge base (clinical significance + population allele frequency),
  5. assigns a review tier with a short, human-readable reason.

Tiering is a simplified, ACMG-inspired triage for demonstration only:
  Tier 1  reported P/LP, or rare predicted loss-of-function (HIGH impact)
  Tier 2  rare protein-altering variant (MODERATE impact) - VUS-like, needs review
  Tier 3  reported B/LB, population AF >= 1%, or low/modifier impact
  Filtered  failed caller QC filters

Assumptions (documented limitations): annotated transcripts are coding-only
(exons == CDS, no UTRs); HGVS is simplified (no 3' shifting, insertions are not
converted to dup). Standard library only.
"""
import argparse
import csv
import gzip
import re
import sys
from collections import defaultdict

COMP = str.maketrans("ACGTN", "TGCAN")
BASES = "TCAG"
AA1 = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {a + b + c: AA1[16 * i + 4 * j + k]
               for i, a in enumerate(BASES) for j, b in enumerate(BASES) for k, c in enumerate(BASES)}
AA3 = dict(A="Ala", R="Arg", N="Asn", D="Asp", C="Cys", Q="Gln", E="Glu", G="Gly", H="His",
           I="Ile", L="Leu", K="Lys", M="Met", F="Phe", P="Pro", S="Ser", T="Thr", W="Trp",
           Y="Tyr", V="Val", X="Xaa")
AA3["*"] = "Ter"

IMPACT = {
    "stop_gained": "HIGH", "frameshift_variant": "HIGH", "splice_donor_variant": "HIGH",
    "splice_acceptor_variant": "HIGH", "start_lost": "HIGH", "stop_lost": "HIGH",
    "missense_variant": "MODERATE", "inframe_deletion": "MODERATE",
    "inframe_insertion": "MODERATE", "protein_altering_variant": "MODERATE",
    "coding_sequence_variant": "MODERATE",
    "synonymous_variant": "LOW", "stop_retained_variant": "LOW", "splice_region_variant": "LOW",
    "intron_variant": "MODIFIER", "intergenic_variant": "MODIFIER",
}
SEVERITY = ["HIGH", "MODERATE", "LOW", "MODIFIER"]
PLP = {"Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"}
BLB = {"Benign", "Likely_benign", "Benign/Likely_benign"}
COMMON_AF = 0.01


def revcomp(s):
    return s.translate(COMP)[::-1]


def translate(cds, to_stop=False):
    prot = []
    for i in range(0, len(cds) - 2, 3):
        aa = CODON_TABLE.get(cds[i:i + 3], "X")
        prot.append(aa)
        if to_stop and aa == "*":
            break
    return "".join(prot)


def aa3(s):
    return "".join(AA3[a] for a in s)


# --------------------------------------------------------------------------
# Inputs
# --------------------------------------------------------------------------
def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            elif line:
                buf.append(line)
    if name:
        seqs[name] = "".join(buf).upper()
    return seqs


class Transcript:
    def __init__(self, tid, gene, chrom, strand):
        self.id, self.gene, self.chrom, self.strand = tid, gene, chrom, strand
        self.exons = []  # CDS segments, genomic (start, end), 1-based inclusive

    def finalise(self, genome):
        self.exons.sort()
        self.start, self.end = self.exons[0][0], self.exons[-1][1]
        gpos = [p for s, e in self.exons for p in range(s, e + 1)]
        if self.strand == "-":
            gpos = gpos[::-1]
        self.cds_gpos = gpos                                   # transcript order
        self.cds_index = {g: i + 1 for i, g in enumerate(gpos)}  # genomic -> c. position
        self.genomic_cds = "".join(genome[s - 1:e] for s, e in self.exons)
        self.cds = self.genomic_cds if self.strand == "+" else revcomp(self.genomic_cds)
        self.protein = translate(self.cds)
        self.introns = [(a[1] + 1, b[0] - 1) for a, b in zip(self.exons, self.exons[1:])]


def read_gtf(path, genome):
    txs = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] != "CDS":
                continue
            attrs = dict(re.findall(r'(\S+) "([^"]*)"', f[8]))
            tid = attrs["transcript_id"]
            if tid not in txs:
                txs[tid] = Transcript(tid, attrs.get("gene_name", attrs["gene_id"]), f[0], f[6])
            txs[tid].exons.append((int(f[3]), int(f[4])))
    for t in txs.values():
        t.finalise(genome[t.chrom])
    return list(txs.values())


def read_kb(path):
    kb = {}
    if not path:
        return kb
    with open(path) as fh:
        rows = (l for l in fh if not l.startswith("#"))
        for r in csv.DictReader(rows, delimiter="\t"):
            kb[(r["chrom"], int(r["pos"]), r["ref"], r["alt"])] = r
    return kb


def read_vcf(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            fmt = dict(zip(f[8].split(":"), f[9].split(":"))) if len(f) > 9 else {}
            for alt in f[4].split(","):
                yield dict(chrom=f[0], pos=int(f[1]), ref=f[3].upper(), alt=alt.upper(),
                           qual=f[5], filter=f[6], fmt=fmt)


# --------------------------------------------------------------------------
# Consequence prediction
# --------------------------------------------------------------------------
def variant_shape(pos, ref, alt):
    """Return (kind, affected genomic positions). Insertions affect the gap pos|pos+1."""
    if len(ref) == len(alt) == 1:
        return "snv", [pos]
    if len(ref) > len(alt) and len(alt) == 1 and ref[0] == alt[0]:
        return "del", list(range(pos + 1, pos + len(ref)))
    if len(alt) > len(ref) and len(ref) == 1 and ref[0] == alt[0]:
        return "ins", [pos, pos + 1]
    return "complex", list(range(pos, pos + len(ref)))


def splice_effect(tx, kind, pos, affected):
    """Check canonical 2-bp splice sites and the 3-8 bp intronic splice region."""
    hit_region = False
    for s, e in tx.introns:
        left_site, right_site = {s, s + 1}, {e - 1, e}
        if kind == "ins":   # insertion disrupts a site only if it falls *between* its two bases
            hits_left, hits_right = pos == s, pos == e - 1
        else:
            hits_left, hits_right = bool(left_site & set(affected)), bool(right_site & set(affected))
        if hits_left or hits_right:
            donor = (hits_left and tx.strand == "+") or (hits_right and tx.strand == "-")
            return "splice_donor_variant" if donor else "splice_acceptor_variant"
        region = set(range(s + 2, s + 8)) | set(range(e - 7, e - 1))
        if region & set(affected):
            hit_region = True
    return "splice_region_variant" if hit_region else None


def intronic_hgvs(tx, pos, ref, alt):
    """c.N+d / c.N-d notation for an intronic SNV."""
    for s, e in tx.introns:
        if s <= pos <= e:
            if tx.strand == "+":
                up, d_up = tx.cds_index[s - 1], pos - (s - 1)
                down, d_down = tx.cds_index[e + 1], (e + 1) - pos
                r, a = ref, alt
            else:
                up, d_up = tx.cds_index[e + 1], (e + 1) - pos
                down, d_down = tx.cds_index[s - 1], pos - (s - 1)
                r, a = ref.translate(COMP), alt.translate(COMP)
            loc = f"{up}+{d_up}" if d_up <= d_down else f"{down}-{d_down}"
            return f"c.{loc}{r}>{a}"
    return "."


def coding_effect(tx, kind, pos, ref, alt):
    """Apply the variant to the CDS and compare translations."""
    gp = tx.cds_index
    # Build alt CDS in genomic orientation (index into tx.genomic_cds)
    g_sorted = tx.cds_gpos if tx.strand == "+" else tx.cds_gpos[::-1]
    gidx = {g: i for i, g in enumerate(g_sorted)}
    if kind == "ins":
        i = gidx[pos]
        alt_g = tx.genomic_cds[:i + 1] + alt[1:] + tx.genomic_cds[i + 1:]
    else:
        if kind == "del":       # VCF anchor base is not deleted; it may lie in the intron
            first, last, repl = pos + 1, pos + len(ref) - 1, ""
        else:
            first, last, repl = pos, pos + len(ref) - 1, alt
        if not (first in gidx and last in gidx and gidx[last] - gidx[first] == last - first):
            return "coding_sequence_variant", ".", "."
        i, j = gidx[first], gidx[last]
        alt_g = tx.genomic_cds[:i] + repl + tx.genomic_cds[j + 1:]
    alt_cds = alt_g if tx.strand == "+" else revcomp(alt_g)

    # ---- HGVS c.
    if kind == "snv":
        c = gp[pos]
        r, a = (ref, alt) if tx.strand == "+" else (ref.translate(COMP), alt.translate(COMP))
        hgvs_c = f"c.{c}{r}>{a}"
    elif kind == "del":
        cs = sorted(gp[p] for p in range(pos + 1, pos + len(ref)))
        hgvs_c = f"c.{cs[0]}del" if len(cs) == 1 else f"c.{cs[0]}_{cs[-1]}del"
    elif kind == "ins":
        ins = alt[1:] if tx.strand == "+" else revcomp(alt[1:])
        a, b = sorted((gp[pos], gp[pos + 1]))
        hgvs_c = f"c.{a}_{b}ins{ins}"
    else:
        cs = sorted(gp[p] for p in range(pos, pos + len(ref)))
        seq = alt if tx.strand == "+" else revcomp(alt)
        hgvs_c = f"c.{cs[0]}_{cs[-1]}delins{seq}"

    ref_p = tx.protein
    delta = len(alt_cds) - len(tx.cds)

    # ---- frameshift
    if delta % 3 != 0:
        alt_p = translate(alt_cds, to_stop=True)
        i = next((k for k in range(min(len(ref_p), len(alt_p))) if ref_p[k] != alt_p[k]),
                 min(len(ref_p), len(alt_p)))
        if i < len(alt_p) and alt_p[i] == "*":
            return "frameshift_variant", hgvs_c, f"p.{AA3[ref_p[i]]}{i + 1}Ter"
        new = AA3[alt_p[i]] if i < len(alt_p) else "Xaa"
        return "frameshift_variant", hgvs_c, f"p.{AA3[ref_p[i]]}{i + 1}{new}fs"

    alt_p = translate(alt_cds)

    # ---- substitutions (same length)
    if delta == 0:
        diffs = [k for k in range(len(ref_p)) if ref_p[k] != alt_p[k]]
        if not diffs:
            k = (gp[pos] - 1) // 3
            csq = "stop_retained_variant" if ref_p[k] == "*" else "synonymous_variant"
            return csq, hgvs_c, f"p.{AA3[ref_p[k]]}{k + 1}="
        k = diffs[0]
        r, a = ref_p[k], alt_p[k]
        if k == 0 and r == "M":
            return "start_lost", hgvs_c, "p.Met1?"
        if r == "*":
            return "stop_lost", hgvs_c, f"p.Ter{k + 1}{AA3[a]}ext*?"
        if a == "*":
            return "stop_gained", hgvs_c, f"p.{AA3[r]}{k + 1}Ter"
        if len(diffs) == 1:
            return "missense_variant", hgvs_c, f"p.{AA3[r]}{k + 1}{AA3[a]}"
        return "protein_altering_variant", hgvs_c, f"p.{AA3[r]}{k + 1}delins{aa3(alt_p[k:diffs[-1] + 1])}"

    # ---- in-frame indels: trim common prefix/suffix of the proteins
    i = 0
    while i < min(len(ref_p), len(alt_p)) and ref_p[i] == alt_p[i]:
        i += 1
    s = 0
    while (s < min(len(ref_p), len(alt_p)) - i and ref_p[len(ref_p) - 1 - s] == alt_p[len(alt_p) - 1 - s]):
        s += 1
    deleted, inserted = ref_p[i:len(ref_p) - s], alt_p[i:len(alt_p) - s]
    if i == 0 and deleted.startswith("M"):
        return "start_lost", hgvs_c, "p.Met1?"
    if "*" in inserted:
        return "stop_gained", hgvs_c, f"p.{AA3[ref_p[i]]}{i + 1}Ter"
    if "*" in deleted:
        return "stop_lost", hgvs_c, "p.?"
    if deleted and not inserted:
        span = f"{AA3[deleted[0]]}{i + 1}" + (f"_{AA3[deleted[-1]]}{i + len(deleted)}" if len(deleted) > 1 else "")
        return "inframe_deletion", hgvs_c, f"p.{span}del"
    if inserted and not deleted:
        return "inframe_insertion", hgvs_c, f"p.{AA3[ref_p[i - 1]]}{i}_{AA3[ref_p[i]]}{i + 1}ins{aa3(inserted)}"
    csq = "inframe_deletion" if delta < 0 else "inframe_insertion"
    return csq, hgvs_c, f"p.{AA3[ref_p[i]]}{i + 1}delins{aa3(inserted)}"


def annotate_against(tx, var):
    pos, ref, alt = var["pos"], var["ref"], var["alt"]
    kind, affected = variant_shape(pos, ref, alt)
    in_cds = [p in tx.cds_index for p in affected]

    splice = splice_effect(tx, kind, pos, affected)
    if splice in ("splice_donor_variant", "splice_acceptor_variant"):
        hgvs_c = intronic_hgvs(tx, pos, ref, alt) if kind == "snv" else "."
        return splice, hgvs_c, "p.?"
    if all(in_cds):
        return coding_effect(tx, kind, pos, ref, alt)
    if any(in_cds):
        return "coding_sequence_variant", ".", "p.?"   # indel spanning an exon/intron boundary
    hgvs_c = intronic_hgvs(tx, pos, ref, alt) if kind == "snv" else "."
    return (splice or "intron_variant"), hgvs_c, "."


def annotate(var, transcripts):
    best = None
    for tx in transcripts:
        if tx.chrom != var["chrom"] or var["pos"] + len(var["ref"]) < tx.start or var["pos"] > tx.end:
            continue
        csq, c, p = annotate_against(tx, var)
        rank = SEVERITY.index(IMPACT[csq])
        if best is None or rank < best[0]:
            best = (rank, tx, csq, c, p)
    if best is None:
        return dict(gene=".", transcript=".", strand=".", consequence="intergenic_variant",
                    impact="MODIFIER", hgvs_c=".", hgvs_p=".")
    _, tx, csq, c, p = best
    return dict(gene=tx.gene, transcript=tx.id, strand=tx.strand, consequence=csq,
                impact=IMPACT[csq], hgvs_c=c, hgvs_p=p)


# --------------------------------------------------------------------------
# Prioritisation
# --------------------------------------------------------------------------
def assign_tier(row):
    if row["filter"] not in ("PASS", "."):
        return "Filtered", f"Failed caller filter ({row['filter']})"
    sig, af = row["kb_clinsig"], row["kb_pop_af"]
    af = float(af) if af not in (".", "") else None
    if sig in PLP:
        return "Tier1", f"Reported {sig} in knowledge base"
    if sig in BLB:
        return "Tier3", f"Reported {sig} in knowledge base"
    if af is not None and af >= COMMON_AF:
        return "Tier3", f"Common in population (AF={af:g})"
    if row["impact"] == "HIGH":
        return "Tier1", "Predicted loss-of-function; rare/absent in population"
    if row["impact"] == "MODERATE":
        return "Tier2", "Rare protein-altering variant; needs review (VUS-like)"
    return "Tier3", f"Low predicted impact ({row['consequence']})"


def genotype_fields(fmt):
    gt = fmt.get("GT", "./.")
    alleles = re.split(r"[/|]", gt)
    zyg = ("hom_alt" if len(set(alleles)) == 1 and alleles[0] not in ("0", ".") else
           "het" if "0" in alleles and len(set(alleles)) > 1 else
           "hom_ref" if set(alleles) == {"0"} else "other")
    ad = fmt.get("AD", ".").split(",")
    dp = fmt.get("DP", ".")
    if len(ad) >= 2 and all(x.isdigit() for x in ad[:2]):
        r, a = int(ad[0]), int(ad[1])
        vaf = f"{a / (r + a):.2f}" if r + a else "."
        return gt, zyg, dp, ad[0], ad[1], vaf
    return gt, zyg, dp, ".", ".", "."


COLUMNS = ["tier", "tier_reason", "chrom", "pos", "ref", "alt", "gene", "transcript", "strand",
           "consequence", "impact", "hgvs_c", "hgvs_p", "gt", "zygosity", "dp", "ad_ref", "ad_alt",
           "vaf", "qual", "filter", "kb_clinsig", "kb_pop_af"]
TIER_ORDER = {"Tier1": 0, "Tier2": 1, "Tier3": 2, "Filtered": 3}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", required=True, help="normalised, biallelic VCF (.vcf or .vcf.gz)")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--kb", help="known-variant TSV (chrom pos ref alt clinical_significance population_af)")
    ap.add_argument("--out", required=True, help="output TSV")
    args = ap.parse_args()

    genome = read_fasta(args.fasta)
    transcripts = read_gtf(args.gtf, genome)
    kb = read_kb(args.kb)

    rows = []
    for var in read_vcf(args.vcf):
        if var["alt"] in ("*", "."):
            continue
        ann = annotate(var, transcripts)
        hit = kb.get((var["chrom"], var["pos"], var["ref"], var["alt"]), {})
        gt, zyg, dp, adr, ada, vaf = genotype_fields(var["fmt"])
        row = dict(chrom=var["chrom"], pos=var["pos"], ref=var["ref"], alt=var["alt"],
                   qual=var["qual"], filter=var["filter"], gt=gt, zygosity=zyg, dp=dp,
                   ad_ref=adr, ad_alt=ada, vaf=vaf,
                   kb_clinsig=hit.get("clinical_significance", "."),
                   kb_pop_af=hit.get("population_af", "."), **ann)
        row["tier"], row["tier_reason"] = assign_tier(row)
        rows.append(row)

    rows.sort(key=lambda r: (TIER_ORDER[r["tier"]], r["chrom"], r["pos"]))
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    counts = defaultdict(int)
    for r in rows:
        counts[r["tier"]] += 1
    summary = ", ".join(f"{t}={counts[t]}" for t in TIER_ORDER if counts[t])
    print(f"[annotate] {len(rows)} variants annotated ({summary}) -> {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
