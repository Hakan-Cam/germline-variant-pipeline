#!/usr/bin/env python3
"""
Simulate a small, fully reproducible germline test case with a known truth set.

Outputs (under --outdir):
  reference/chrS.fa              synthetic ~10.7 kb chromosome
  reference/genes.gtf            gene model: GENE_A (+ strand, 3 exons), GENE_B (- strand, 2 exons)
  reads/<sample>_R1.fastq.gz     paired-end 2x150 Illumina-like reads from a diploid genome
  reads/<sample>_R2.fastq.gz
  truth/<sample>.truth.vcf       19 planted variants (SNVs + indels, het + hom)
  truth/truth_variants.tsv       each variant's designed (expected) consequence
  kb/known_variants.tsv          toy "ClinVar/gnomAD-like" knowledge base for prioritisation

Variants are designed in transcript (c.) space so that the expected consequence
(synonymous, missense, stop_gained, start_lost, frameshift, in-frame deletion,
splice donor/acceptor, intronic, intergenic) is known in advance. The annotation
step must rediscover these consequences independently from VCF + GTF + FASTA.

All data are synthetic. Only the Python standard library is used.
"""
import argparse
import gzip
import os
import random

COMP = str.maketrans("ACGTN", "TGCAN")
STOPS = {"TAA", "TAG", "TGA"}
SENSE_CODONS = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"
                if a + b + c not in STOPS]
TRANSITION = {"A": "G", "G": "A", "C": "T", "T": "C"}


def revcomp(s):
    return s.translate(COMP)[::-1]


def rand_seq(rng, n, gc=0.45):
    w = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]  # A C G T
    return "".join(rng.choices("ACGT", weights=w, k=n))


# --------------------------------------------------------------------------
# Design: genes, planted codons and truth variants
# --------------------------------------------------------------------------
# planted: 1-based codon number -> codon, so designed variants have known effects
GENES = [
    dict(gene="GENE_A", strand="+", exon_lens=[155, 205, 141], intron_lens=[420, 380],
         stop="TAA", planted={25: "GCC", 40: "CGA", 60: "CGG", 80: "GAT"}),
    dict(gene="GENE_B", strand="-", exon_lens=[190, 233], intron_lens=[450],
         stop="TGA", planted={30: "TCA", 55: "AAA", 90: "CAG"}),
]

# Transcript-space SNVs: (id, gene, location, ref, alt, GT, expected consequence)
# location: ("c", pos) coding position, or ("intron", k, offset) with offset > 0
# counted from the intron's 5' end, or offset < 0 counted back from its 3' end.
TX_SNVS = [
    ("V01", "GENE_A", ("c", 75), "C", "T", "0/1", "synonymous_variant"),     # Ala25Ala
    ("V02", "GENE_A", ("c", 118), "C", "T", "0/1", "stop_gained"),           # Arg40Ter
    ("V03", "GENE_A", ("c", 178), "C", "T", "1/1", "missense_variant"),      # Arg60Trp
    ("V04", "GENE_A", ("c", 238), "G", "A", "0/1", "missense_variant"),      # Asp80Asn
    ("V05", "GENE_A", ("c", 2), "T", "C", "0/1", "start_lost"),              # Met1?
    ("V06", "GENE_A", ("intron", 1, 1), "G", "A", "0/1", "splice_donor_variant"),
    ("V07", "GENE_A", ("intron", 1, 200), None, None, "0/1", "intron_variant"),
    ("V10", "GENE_B", ("c", 90), "A", "G", "0/1", "synonymous_variant"),     # Ser30Ser
    ("V11", "GENE_B", ("c", 164), "A", "G", "0/1", "missense_variant"),      # Lys55Arg
    ("V12", "GENE_B", ("c", 268), "C", "T", "1/1", "stop_gained"),           # Gln90Ter
    ("V14", "GENE_B", ("intron", 1, -2), "A", "G", "0/1", "splice_acceptor_variant"),
]
# Transcript-space indels: (id, gene, kind, cds_start, cds_end_or_insseq, GT, expected)
TX_INDELS = [
    ("V08", "GENE_A", "del", 310, 310, "0/1", "frameshift_variant"),
    ("V09", "GENE_A", "del", 400, 402, "0/1", "inframe_deletion"),
    ("V13", "GENE_B", "ins", 300, "A", "0/1", "frameshift_variant"),
]
# Intergenic variants: (id, anchor, offset, kind, payload, GT)
# anchor: "chr" (absolute), "afterA" (bp after GENE_A end), "afterB" (bp after GENE_B end)
INTERGENIC = [
    ("V15", "chr", 1000, "snv", None, "1/1"),
    ("V16", "chr", 2000, "snv", None, "0/1"),
    ("V17", "afterA", 1000, "snv", None, "0/1"),
    ("V19", "afterA", 1800, "ins", "GC", "0/1"),
    ("V18", "afterB", 1500, "del", 1, "0/1"),
]
# Toy knowledge base entries: id -> (clinical significance, population AF)
KB = {
    "V02": ("Pathogenic", 0.0),
    "V04": ("Benign", 0.12),
    "V11": ("Likely_pathogenic", 0.00004),
    "V16": ("Benign", 0.35),
}


# --------------------------------------------------------------------------
# Genome construction
# --------------------------------------------------------------------------
def build_cds(rng, spec):
    n_codons = sum(spec["exon_lens"]) // 3
    assert sum(spec["exon_lens"]) % 3 == 0, "CDS length must be a multiple of 3"
    codons = ["ATG"] + rng.choices(SENSE_CODONS, k=n_codons - 2) + [spec["stop"]]
    for num, codon in spec["planted"].items():
        codons[num - 1] = codon
    return "".join(codons)


def build_gene(rng, spec):
    """Return the unspliced transcript (transcript orientation) and exon/intron offsets."""
    cds = build_cds(rng, spec)
    pre, exons, introns = [], [], []   # exons: (pre_start, cds_start); introns: (pre_start, length)
    cds_pos, pre_pos = 1, 0
    for i, elen in enumerate(spec["exon_lens"]):
        exons.append((pre_pos, cds_pos, elen))
        pre.append(cds[cds_pos - 1: cds_pos - 1 + elen])
        cds_pos += elen
        pre_pos += elen
        if i < len(spec["intron_lens"]):
            ilen = spec["intron_lens"][i]
            introns.append((pre_pos, ilen))
            pre.append("GT" + rand_seq(rng, ilen - 4) + "AG")
            pre_pos += ilen
    return dict(spec, cds=cds, pre="".join(pre), exons=exons, introns=introns)


def pre_index_for_cds(gene, c):
    for pre_start, cds_start, elen in gene["exons"]:
        if cds_start <= c < cds_start + elen:
            return pre_start + (c - cds_start)
    raise ValueError(f"c.{c} outside CDS of {gene['gene']}")


def genomic_pos(gene, pre_idx):
    """1-based genomic coordinate for a 0-based index into the unspliced transcript."""
    if gene["strand"] == "+":
        return gene["gstart"] + pre_idx
    return gene["gend"] - pre_idx


def to_genome_base(gene, base):
    return base if gene["strand"] == "+" else base.translate(COMP)


def build_genome(rng):
    genes = [build_gene(rng, s) for s in GENES]
    parts, pos = [], 1
    layout = [("flank", 3000), ("gene", 0), ("flank", 2500), ("gene", 1), ("flank", 3000)]
    for kind, val in layout:
        if kind == "flank":
            seq = rand_seq(rng, val)
        else:
            g = genes[val]
            seq = g["pre"] if g["strand"] == "+" else revcomp(g["pre"])
            g["gstart"], g["gend"] = pos, pos + len(seq) - 1
        parts.append(seq)
        pos += len(seq)
    return "".join(parts), {g["gene"]: g for g in genes}


# --------------------------------------------------------------------------
# Variants
# --------------------------------------------------------------------------
def make_variants(genome, genes, rng):
    variants = []  # dicts: id chrom pos ref alt gt expected gene

    def add(vid, pos, ref, alt, gt, expected, gene):
        assert genome[pos - 1: pos - 1 + len(ref)] == ref, f"{vid}: REF mismatch at {pos}"
        variants.append(dict(id=vid, pos=pos, ref=ref, alt=alt, gt=gt,
                             expected=expected, gene=gene))

    for vid, gname, loc, ref_t, alt_t, gt, expected in TX_SNVS:
        g = genes[gname]
        if loc[0] == "c":
            idx = pre_index_for_cds(g, loc[1])
        else:
            pre_start, ilen = g["introns"][loc[1] - 1]
            off = loc[2]
            idx = pre_start + off - 1 if off > 0 else pre_start + ilen + off
        tx_base = g["pre"][idx]
        if ref_t is None:
            ref_t, alt_t = tx_base, TRANSITION[tx_base]
        assert tx_base == ref_t, f"{vid}: designed ref {ref_t} != sequence {tx_base}"
        add(vid, genomic_pos(g, idx), to_genome_base(g, ref_t), to_genome_base(g, alt_t),
            gt, expected, gname)

    for vid, gname, kind, c_start, payload, gt, expected in TX_INDELS:
        g = genes[gname]
        if kind == "del":
            gpos = sorted(genomic_pos(g, pre_index_for_cds(g, c)) for c in range(c_start, payload + 1))
            anchor = gpos[0] - 1
            ref = genome[anchor - 1: gpos[-1]]
            add(vid, anchor, ref, ref[0], gt, expected, gname)
        else:  # insertion after c_start (transcript orientation)
            left_c = c_start if g["strand"] == "+" else c_start + 1
            anchor = genomic_pos(g, pre_index_for_cds(g, left_c))
            ins = payload if g["strand"] == "+" else revcomp(payload)
            ref = genome[anchor - 1]
            add(vid, anchor, ref, ref + ins, gt, expected, gname)

    ends = {"chr": 0, "afterA": genes["GENE_A"]["gend"], "afterB": genes["GENE_B"]["gend"]}
    for vid, anchor_name, off, kind, payload, gt in INTERGENIC:
        pos = ends[anchor_name] + off
        ref = genome[pos - 1]
        if kind == "snv":
            add(vid, pos, ref, TRANSITION[ref], gt, "intergenic_variant", ".")
        elif kind == "ins":
            add(vid, pos, ref, ref + payload, gt, "intergenic_variant", ".")
        else:
            ref = genome[pos - 1: pos + payload]
            add(vid, pos, ref, ref[0], gt, "intergenic_variant", ".")

    variants.sort(key=lambda v: v["pos"])
    for a, b in zip(variants, variants[1:]):
        assert b["pos"] - a["pos"] > 30, f"{a['id']} and {b['id']} are too close"
    return variants


def make_haplotypes(genome, variants, rng):
    haps = [list(genome), list(genome)]
    for v in sorted(variants, key=lambda v: -v["pos"]):  # right-to-left keeps coordinates valid
        targets = [0, 1] if v["gt"] == "1/1" else [rng.randint(0, 1)]
        for h in targets:
            s = v["pos"] - 1
            haps[h][s: s + len(v["ref"])] = list(v["alt"])
    return ["".join(h) for h in haps]


# --------------------------------------------------------------------------
# Reads
# --------------------------------------------------------------------------
def add_errors(read, rate, rng):
    if rate <= 0:
        return read
    out = list(read)
    for i in range(len(out)):
        if rng.random() < rate:
            out[i] = rng.choice([b for b in "ACGT" if b != out[i]])
    return "".join(out)


def simulate_reads(haps, args, rng, r1_path, r2_path):
    genome_len = len(haps[0])
    n_pairs = int(args.coverage * genome_len / (2 * args.read_len))
    qual = "F" * args.read_len
    with gzip.open(r1_path, "wt") as r1, gzip.open(r2_path, "wt") as r2:
        for i in range(n_pairs):
            h = i % 2
            seq = haps[h]
            flen = max(args.read_len + 20, int(rng.gauss(args.frag_mean, args.frag_sd)))
            start = rng.randint(0, len(seq) - flen)
            frag = seq[start: start + flen]
            if rng.random() < 0.5:
                frag = revcomp(frag)
            a = add_errors(frag[: args.read_len], args.error_rate, rng)
            b = add_errors(revcomp(frag)[: args.read_len], args.error_rate, rng)
            name = f"{args.sample}:{i}:hap{h + 1}"
            r1.write(f"@{name}/1\n{a}\n+\n{qual}\n")
            r2.write(f"@{name}/2\n{b}\n+\n{qual}\n")
    return n_pairs


# --------------------------------------------------------------------------
# Writers
# --------------------------------------------------------------------------
def write_fasta(path, name, seq, width=60):
    with open(path, "w") as fh:
        fh.write(f">{name}\n")
        for i in range(0, len(seq), width):
            fh.write(seq[i: i + width] + "\n")


def write_gtf(path, chrom, genes):
    def attrs(g, extra=""):
        return f'gene_id "{g["gene"]}"; transcript_id "{g["gene"]}-201"; gene_name "{g["gene"]}";{extra}'

    with open(path, "w") as fh:
        for g in genes.values():
            s, e, st = g["gstart"], g["gend"], g["strand"]
            fh.write(f'{chrom}\tsimulated\tgene\t{s}\t{e}\t.\t{st}\t.\tgene_id "{g["gene"]}"; gene_name "{g["gene"]}";\n')
            fh.write(f"{chrom}\tsimulated\ttranscript\t{s}\t{e}\t.\t{st}\t.\t{attrs(g)}\n")
            for n, (pre_start, cds_start, elen) in enumerate(g["exons"], 1):
                p1 = genomic_pos(g, pre_start)
                p2 = genomic_pos(g, pre_start + elen - 1)
                lo, hi = min(p1, p2), max(p1, p2)
                frame = (3 - (cds_start - 1) % 3) % 3
                ex = f' exon_number "{n}";'
                fh.write(f"{chrom}\tsimulated\texon\t{lo}\t{hi}\t.\t{st}\t.\t{attrs(g, ex)}\n")
                fh.write(f"{chrom}\tsimulated\tCDS\t{lo}\t{hi}\t.\t{st}\t{frame}\t{attrs(g, ex)}\n")


def write_truth(outdir, chrom, genome_len, variants, sample):
    with open(os.path.join(outdir, f"{sample}.truth.vcf"), "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write(f"##contig=<ID={chrom},length={genome_len}>\n")
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fh.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n")
        for v in variants:
            fh.write(f"{chrom}\t{v['pos']}\t{v['id']}\t{v['ref']}\t{v['alt']}\t.\tPASS\t.\tGT\t{v['gt']}\n")
    with open(os.path.join(outdir, "truth_variants.tsv"), "w") as fh:
        fh.write("id\tchrom\tpos\tref\talt\tgt\tgene\texpected_consequence\n")
        for v in variants:
            fh.write(f"{v['id']}\t{chrom}\t{v['pos']}\t{v['ref']}\t{v['alt']}\t{v['gt']}\t{v['gene']}\t{v['expected']}\n")


def write_kb(path, chrom, variants):
    with open(path, "w") as fh:
        fh.write("# Toy knowledge base for a synthetic genome - NOT real clinical data\n")
        fh.write("chrom\tpos\tref\talt\tgene\tclinical_significance\tpopulation_af\n")
        for v in variants:
            if v["id"] in KB:
                sig, af = KB[v["id"]]
                fh.write(f"{chrom}\t{v['pos']}\t{v['ref']}\t{v['alt']}\t{v['gene']}\t{sig}\t{af}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--sample", default="SAMPLE01")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--coverage", type=float, default=40.0)
    ap.add_argument("--read-len", type=int, default=150)
    ap.add_argument("--frag-mean", type=int, default=350)
    ap.add_argument("--frag-sd", type=int, default=35)
    ap.add_argument("--error-rate", type=float, default=0.002)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    chrom = "chrS"
    genome, genes = build_genome(rng)
    variants = make_variants(genome, genes, rng)
    haps = make_haplotypes(genome, variants, rng)

    dirs = {d: os.path.join(args.outdir, d) for d in ("reference", "reads", "truth", "kb")}
    for d in dirs.values():
        os.makedirs(d, exist_ok=True)

    write_fasta(os.path.join(dirs["reference"], f"{chrom}.fa"), chrom, genome)
    write_gtf(os.path.join(dirs["reference"], "genes.gtf"), chrom, genes)
    write_truth(dirs["truth"], chrom, len(genome), variants, args.sample)
    write_kb(os.path.join(dirs["kb"], "known_variants.tsv"), chrom, variants)
    n = simulate_reads(haps, args, rng,
                       os.path.join(dirs["reads"], f"{args.sample}_R1.fastq.gz"),
                       os.path.join(dirs["reads"], f"{args.sample}_R2.fastq.gz"))
    print(f"[simulate] genome {len(genome):,} bp | {len(genes)} genes | "
          f"{len(variants)} truth variants | {n:,} read pairs (~{args.coverage:g}x)")


if __name__ == "__main__":
    main()
