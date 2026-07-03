"""The canonical synthetic cohort: PXX (XX), PXY (XY, WGD+diploid), PXXY (XXY).

Minimal constructed truth — just enough to fire every relevant Python-script code
path and edge case. See python/TESTING.md.
"""

from __future__ import annotations

from .model import (
    Assay, Clone, Cohort, GermlineSite, Patient, Sample, Segment, SegmentCN,
    SequencingRun, SomaticVariant, neutral_split, wgd_neutral_split,
)


def _fill_neutral(patient: Patient) -> None:
    """Fill in a neutral SegmentCN for every (segment, tumor sample) left
    unspecified, using the contig's germline ploidy and the sample's WGD count.
    Keeps the cohort definitions terse."""
    wgd_by_sample = {s.name: s.wgd for s in patient.tumor_samples}
    for seg in patient.segments:
        g1, g2 = neutral_split(patient.chromosomal_ploidy(seg.contig))
        for name, w in wgd_by_sample.items():
            if name in seg.per_sample:
                continue
            a1, a2 = wgd_neutral_split(g1, g2, w)
            seg.per_sample[name] = SegmentCN(cn_a1=a1, cn_a2=a2)


def _het(contig, position, seg_contig=None, haplotype="A", pop_af=0.5):
    return GermlineSite(contig=contig, position=position, ref="A", alt="G",
                        pop_af=pop_af, genotype="0/1", haplotype=haplotype,
                        segment_contig=seg_contig or contig)


def _hom(contig, position, genotype="1/1", pop_af=0.9):
    return GermlineSite(contig=contig, position=position, ref="C", alt="T",
                        pop_af=pop_af, genotype=genotype)


def _junk(contig, position, ref="AC", alt="G"):
    # multi-base ref/alt => not a biallelic SNV => must be filtered by genotype.VCF
    return GermlineSite(contig=contig, position=position, ref=ref, alt=alt,
                        pop_af=0.3, genotype="0/1", is_snv=False)


# --------------------------------------------------------------------------- #
# PXX — female (XX), no WGD
# --------------------------------------------------------------------------- #
def _build_pxx() -> Patient:
    N = Sample("PXX-N", is_normal=True, timepoint=0, runs=[SequencingRun(Assay.WES, "WES_A")])
    T1 = Sample("PXX-T1", False, 1, [SequencingRun(Assay.WES, "WES_A")], purity=0.8, ploidy=2.0)
    T2 = Sample("PXX-T2", False, 2,
                [SequencingRun(Assay.ULP_WGS), SequencingRun(Assay.DEEP_PANEL, "PANEL")],
                purity=0.6, ploidy=2.0)
    T3 = Sample("PXX-T3", False, 3,
                [SequencingRun(Assay.ULP_WGS), SequencingRun(Assay.WES, "WES_B")],
                purity=0.7, ploidy=2.0)
    # out-of-patient contamination: T2 a touch contaminated, others clean
    T2.contamination_true = 0.03

    clones = [
        Clone("C1", None, {"PXX-T1": 1.0, "PXX-T2": 1.0, "PXX-T3": 1.0}),   # trunk
        Clone("C2", "C1", {"PXX-T1": 0.6, "PXX-T2": 0.0, "PXX-T3": 0.5}),   # branch
        Clone("C3", "C1", {"PXX-T1": 0.0, "PXX-T2": 0.5, "PXX-T3": 0.0}),   # private to T2
    ]

    segments = [
        Segment("1", 100_000, 400_000, "NEUTRAL", {}),
        Segment("1", 500_000, 900_000, "LOH",
                {"PXX-T1": SegmentCN(2, 0), "PXX-T3": SegmentCN(2, 0)}, is_cnloh=True, n_hets=4),
        Segment("1", 1_000_000, 1_200_000, "FOCAL_BALANCED_HIGHCN",
                {"PXX-T1": SegmentCN(5, 5)}, n_hets=2),  # f≈0.5 high CN, few hets -> resplit
        Segment("2", 100_000, 600_000, "HZ_DEL", {"PXX-T1": SegmentCN(0, 0)}, n_hets=0),
        Segment("2", 700_000, 1_300_000, "IMBALANCED_GAIN",
                {"PXX-T1": SegmentCN(3, 1), "PXX-T3": SegmentCN(3, 1)}, n_hets=5),
        Segment("3", 100_000, 500_000, "HIGH_AMP", {"PXX-T1": SegmentCN(8, 1)}, n_hets=3),
        # Subclonal gain: event (3,1) over baseline (1,1) at ccf 0.5. NOTE: a subclonal
        # CN event is generally NOT identifiable from copy ratio alone — its mixed
        # observed allelic CN (here ~2,1) aliases to a *clonal* lower-integer state
        # (2,1). map_to_absolute's per-segment CCF therefore reports such a segment as
        # clonal; subclonal *variant* CCFs are instead recovered by calculate_ccf from
        # VAF (which is identifiable). The observed copy ratio IS the ccf-mix, though.
        Segment("3", 600_000, 900_000, "SUBCLONAL_GAIN",
                {"PXX-T2": SegmentCN(3, 1, ccf=0.5, baseline_a1=1, baseline_a2=1)}, n_hets=4),
        Segment("X", 100_000, 800_000, "LOH",
                {"PXX-T1": SegmentCN(2, 0)}, is_cnloh=True, n_hets=3),  # XX: ploidy 2
    ]

    germline = [
        # pop_af deliberately spans common<->rare so the genotyping pop_af prior is
        # exercised (a common SNP should yield a higher variant genotype likelihood
        # than a rare one with the same read evidence). Asserted directly in
        # test_genotype.py::test_pop_af_prior_*.
        _het("1", 600_000, "1", pop_af=0.45), _het("1", 650_000, "1", pop_af=0.08),
        _het("1", 700_000, "1", haplotype="B", pop_af=0.30),
        _het("2", 800_000, "2", pop_af=0.45), _het("2", 900_000, "2", haplotype="B", pop_af=0.20),
        _het("2", 1_000_000, "2", pop_af=0.10),
        _hom("1", 300_000), _hom("3", 200_000, genotype="0/0"),
        _het("X", 400_000, "X", pop_af=0.40), _het("X", 500_000, "X", pop_af=0.15),
        _junk("1", 350_000), _junk("2", 650_000, ref="A", alt="."),
    ]

    somatic = [
        SomaticVariant("1", 250_000, "C", "T", "SNV", "C1", multiplicity=1, in_abs_maf=True),
        SomaticVariant("3", 300_000, "G", "A", "SNV", "C1", multiplicity=2, in_abs_maf=True),  # on HIGH_AMP
        SomaticVariant("2", 800_000, "A", "G", "SNV", "C2", multiplicity=1, in_abs_maf=False),  # rescue, subclonal
        SomaticVariant("3", 700_000, "T", "C", "SNV", "C3", multiplicity=1, in_abs_maf=False),  # subclonal T2
        SomaticVariant("1", 200_000, "AT", "A", "INDEL", "C1", multiplicity=1, in_abs_maf=False),
        SomaticVariant("MT", 8_000, "G", "A", "SNV", "C1", multiplicity=1, note="MT"),         # dropped
        SomaticVariant("Y", 100_000, "C", "G", "SNV", "C1", multiplicity=1, note="unmapped"),   # Y absent from female segtab
        SomaticVariant("2", 1_100_000, "A", "T", "SNV", "C1", multiplicity=1, note="zero_cov"),
    ]

    p = Patient("PXX", "XX", [N, T1, T2, T3], clones, segments, germline, somatic)
    _fill_neutral(p)
    return p


# --------------------------------------------------------------------------- #
# PXY — male (XY); WGD in T1, diploid in T3
# --------------------------------------------------------------------------- #
def _build_pxy() -> Patient:
    N = Sample("PXY-N", True, 0, [SequencingRun(Assay.WES, "WES_A")])
    # two WGD events -> quadruploid (ploidy 4)
    T1 = Sample("PXY-T1", False, 1,
                [SequencingRun(Assay.WES, "WES_A"), SequencingRun(Assay.WES, "WES_A")],
                purity=0.85, ploidy=4.0, wgd=2)
    T2 = Sample("PXY-T2", False, 2, [SequencingRun(Assay.DEEP_PANEL, "PANEL")],
                purity=0.7, ploidy=2.0, wgd=0)
    # diploid male tumor: exercises haploid-chrX / chrY events in a non-WGD background
    T3 = Sample("PXY-T3", False, 3, [SequencingRun(Assay.WES, "WES_A")],
                purity=0.75, ploidy=2.0, wgd=0)
    # single WGD event -> triploid (ploidy 3); WGD=1 vs WGD=2 coverage
    T4 = Sample("PXY-T4", False, 4, [SequencingRun(Assay.WES, "WES_A")],
                purity=0.8, ploidy=3.0, wgd=1)

    clones = [
        Clone("C1", None, {"PXY-T1": 1.0, "PXY-T2": 1.0, "PXY-T3": 1.0, "PXY-T4": 1.0}),
        Clone("C2", "C1", {"PXY-T1": 0.7, "PXY-T2": 0.4, "PXY-T3": 0.0, "PXY-T4": 0.6}),
    ]

    segments = [
        Segment("1", 100_000, 600_000, "NEUTRAL", {}),
        Segment("1", 700_000, 1_200_000, "IMBALANCED_GAIN",
                {"PXY-T1": SegmentCN(4, 2), "PXY-T3": SegmentCN(3, 1)}, n_hets=5),
        Segment("2", 100_000, 800_000, "LOH",
                {"PXY-T1": SegmentCN(4, 0)}, is_cnloh=True, n_hets=4),
        # diploid-male haploid-chrX events (would be masked under WGD everywhere):
        Segment("X", 100_000, 700_000, "X_LOH", {"PXY-T3": SegmentCN(0, 0)}, n_hets=0),   # 1 -> 0
        Segment("X", 800_000, 1_400_000, "X_GAIN", {"PXY-T3": SegmentCN(0, 2)}, n_hets=2),  # 1 -> 2
        Segment("Y", 50_000, 400_000, "Y_LOSS", {"PXY-T3": SegmentCN(0, 0)}, n_hets=0),     # 1 -> 0
    ]

    germline = [
        _het("1", 800_000, "1"), _het("1", 900_000, "1", haplotype="B"), _het("1", 1_000_000, "1"),
        _het("2", 300_000, "2"), _het("2", 400_000, "2"),
        _hom("1", 200_000), _het("X", 900_000, "X"),  # X het only meaningful where ploidy>=2
    ]

    somatic = [
        SomaticVariant("1", 300_000, "C", "T", "SNV", "C1", multiplicity=1, in_abs_maf=True),
        SomaticVariant("1", 800_000, "G", "A", "SNV", "C2", multiplicity=1, in_abs_maf=False),
        SomaticVariant("X", 200_000, "A", "G", "SNV", "C1", multiplicity=1, in_abs_maf=False),  # haploid X
    ]

    p = Patient("PXY", "XY", [N, T1, T2, T3, T4], clones, segments, germline, somatic)
    _fill_neutral(p)
    return p


# --------------------------------------------------------------------------- #
# PXXY — XXY (non-standard karyotype)
# --------------------------------------------------------------------------- #
def _build_pxxy() -> Patient:
    N = Sample("PXXY-N", True, 0, [SequencingRun(Assay.WES, "WES_A")])
    T1 = Sample("PXXY-T1", False, 1, [SequencingRun(Assay.WES, "WES_A")],
                purity=0.8, ploidy=2.0)

    clones = [Clone("C1", None, {"PXXY-T1": 1.0})]

    segments = [
        Segment("1", 100_000, 800_000, "NEUTRAL", {}),
        Segment("X", 100_000, 700_000, "LOH", {"PXXY-T1": SegmentCN(2, 0)}, is_cnloh=True, n_hets=3),  # X ploidy 2
        Segment("Y", 50_000, 400_000, "NEUTRAL", {}),  # Y ploidy 1
    ]

    germline = [
        _het("1", 300_000, "1"), _het("1", 400_000, "1", haplotype="B"),
        _het("X", 300_000, "X"), _het("X", 400_000, "X", haplotype="B"),
    ]

    somatic = [SomaticVariant("1", 200_000, "C", "T", "SNV", "C1", multiplicity=1, in_abs_maf=True)]

    p = Patient("PXXY", "XXY", [N, T1], clones, segments, germline, somatic)
    _fill_neutral(p)
    return p


def build_cohort() -> Cohort:
    return Cohort([_build_pxx(), _build_pxy(), _build_pxxy()])
