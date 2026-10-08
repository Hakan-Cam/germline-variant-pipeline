"""Unit tests for consequence prediction on hand-built + and - strand transcripts."""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import annotate_variants as av  # noqa: E402

# Toy contig (1-based coordinates):
#   1-10   flank
#   11-19  exon 1  ATG GCC CG   (CDS c.1-9 for the + strand gene)
#   20-29  intron  GT....AG
#   30-38  exon 2  A TGG CAG TAA ... (c.10-18)
EXON1 = "ATGGCCCGA"          # Met Ala Arg
INTRON = "GTAAAAAAAG"        # 10 bp
EXON2 = "TGGCAGTAA"          # Trp Gln Ter
FLANK = "CCCCCCCCCC"
GENOME = FLANK + EXON1 + INTRON + EXON2 + FLANK
E1, E2 = (11, 19), (30, 38)


def make_tx(strand, genome):
    tx = av.Transcript("TX1", "GENE", "chr1", strand)
    tx.exons = [E1, E2]
    tx.finalise(genome)
    return tx


@pytest.fixture
def plus():
    return make_tx("+", GENOME)


@pytest.fixture
def minus():
    # Same exons, but the gene is read on the reverse strand of a reverse-complemented contig.
    g = av.revcomp(GENOME)
    n = len(g)
    tx = av.Transcript("TX2", "GENEM", "chr1", "-")
    tx.exons = [(n - E2[1] + 1, n - E2[0] + 1), (n - E1[1] + 1, n - E1[0] + 1)]
    tx.finalise(g)
    return tx, g


def var(pos, ref, alt):
    return dict(chrom="chr1", pos=pos, ref=ref, alt=alt)


def test_protein(plus):
    assert plus.cds == EXON1 + EXON2
    assert plus.protein == "MARWQ*"


@pytest.mark.parametrize("pos,ref,alt,csq,hgvs_c,hgvs_p", [
    (16, "C", "T", "synonymous_variant", "c.6C>T", "p.Ala2="),
    (17, "C", "T", "stop_gained", "c.7C>T", "p.Arg3Ter"),      # CGA -> TGA
    (18, "G", "C", "missense_variant", "c.8G>C", "p.Arg3Pro"),
    (12, "T", "C", "start_lost", "c.2T>C", "p.Met1?"),
    (36, "T", "C", "stop_lost", "c.16T>C", "p.Ter6Glnext*?"),
    (20, "G", "A", "splice_donor_variant", "c.9+1G>A", "p.?"),
    (29, "G", "C", "splice_acceptor_variant", "c.10-1G>C", "p.?"),
    (23, "A", "G", "splice_region_variant", "c.9+4A>G", "."),
    (5, "C", "T", "intergenic_variant", ".", "."),
])
def test_plus_strand_snvs(plus, pos, ref, alt, csq, hgvs_c, hgvs_p):
    ann = av.annotate(var(pos, ref, alt), [plus])
    assert (ann["consequence"], ann["hgvs_c"], ann["hgvs_p"]) == (csq, hgvs_c, hgvs_p)


def test_plus_frameshift_and_inframe(plus):
    fs = av.annotate(var(13, "GG", "G"), [plus])          # delete c.4 (G)
    assert fs["consequence"] == "frameshift_variant" and fs["hgvs_c"] == "c.4del"
    inframe = av.annotate(var(13, "GGCC", "G"), [plus])   # delete c.4_6 (GCC = Ala2)
    assert inframe["consequence"] == "inframe_deletion"
    assert inframe["hgvs_p"] == "p.Ala2del"
    ins = av.annotate(var(14, "G", "GAAA"), [plus])       # in-frame Lys insertion
    assert ins["consequence"] == "inframe_insertion"


def test_minus_strand_matches_plus(minus):
    tx, g = minus
    assert tx.cds == EXON1 + EXON2
    n = len(g)
    # c.7C>T on the + gene sits at genomic 17; on the mirrored contig it is n-17+1 with G>A
    ann = av.annotate(dict(chrom="chr1", pos=n - 17 + 1, ref="G", alt="A"), [tx])
    assert (ann["consequence"], ann["hgvs_c"], ann["hgvs_p"]) == ("stop_gained", "c.7C>T", "p.Arg3Ter")
    donor = av.annotate(dict(chrom="chr1", pos=n - 20 + 1, ref="C", alt="T"), [tx])
    assert (donor["consequence"], donor["hgvs_c"]) == ("splice_donor_variant", "c.9+1G>A")


def test_tiering():
    base = dict(filter="PASS", kb_clinsig=".", kb_pop_af=".", impact="HIGH", consequence="stop_gained")
    assert av.assign_tier(base)[0] == "Tier1"
    assert av.assign_tier({**base, "kb_pop_af": "0.2"})[0] == "Tier3"
    assert av.assign_tier({**base, "impact": "MODERATE"})[0] == "Tier2"
    assert av.assign_tier({**base, "impact": "LOW", "kb_clinsig": "Pathogenic"})[0] == "Tier1"
    assert av.assign_tier({**base, "filter": "LowQual"})[0] == "Filtered"
