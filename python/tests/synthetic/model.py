"""Ground-truth model for the synthetic golden cohort.

These dataclasses describe an *invented* truth (clones, copy-number events,
germline sites, somatic variants, purity/ploidy/sex). The simulation layer
(simulate.py / emit.py) forward-simulates the upstream tool-output files from
this truth; the truth tables themselves are the answer key against which the
Python scripts are validated.

The cohort is deliberately minimal and NOT biologically realistic — it is the
smallest construction that fires every relevant code path and edge case.
See python/TESTING.md.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum


# --------------------------------------------------------------------------- #
# Reference
# --------------------------------------------------------------------------- #
# A tiny synthetic reference: a few autosomes, both allosomes, and MT. Lengths
# are small but large enough to place a handful of segments and loci per contig.
CONTIGS: list[tuple[str, int]] = [
    ("1", 2_000_000),
    ("2", 2_000_000),
    ("3", 1_000_000),
    ("X", 1_500_000),
    ("Y", 500_000),
    ("MT", 16_000),
]
CONTIG_NAMES = [c for c, _ in CONTIGS]
CONTIG_LENGTHS = dict(CONTIGS)


def chr_name(contig: str, chr_prefixed: bool) -> str:
    """Map a canonical contig name to the chr-prefixed convention if requested."""
    if not chr_prefixed:
        return contig
    return "chrM" if contig == "MT" else f"chr{contig}"


# --------------------------------------------------------------------------- #
# Assays
# --------------------------------------------------------------------------- #
class Assay(str, Enum):
    WES = "WES"            # decent genome coverage (CR) + decent depth (allelic)
    ULP_WGS = "ULP_WGS"    # genome-wide CR, but ~no allelic depth
    DEEP_PANEL = "DEEP_PANEL"  # deep allelic at targets, no genome segmentation


# default (use_for_tCR, use_for_aCR) per assay
ASSAY_USE = {
    Assay.WES: (True, True),
    Assay.ULP_WGS: (True, False),
    Assay.DEEP_PANEL: (False, True),
}


@dataclass(frozen=True)
class SequencingRun:
    assay: Assay
    # which target set the run's data lives on. WES panels can differ between runs
    # to drive harmonize intersection cases ("WES_A"/"WES_B"/"WES_C"); ULP_WGS is
    # genome-wide (None); DEEP_PANEL uses "PANEL".
    target_set: str | None = None
    use_for_tcr: bool | None = None
    use_for_acr: bool | None = None
    is_paired_end: bool = True

    def resolved_use(self) -> tuple[bool, bool]:
        t, a = ASSAY_USE[self.assay]
        return (t if self.use_for_tcr is None else self.use_for_tcr,
                a if self.use_for_acr is None else self.use_for_acr)


@dataclass
class Sample:
    name: str
    is_normal: bool
    timepoint: int
    runs: list[SequencingRun]
    purity: float = 0.0          # tumor purity; 0 for normals (models normal admixture)
    ploidy: float = 2.0          # genome ploidy (= normal_ploidy + wgd for a diploid genome)
    wgd: int = 0                 # number of whole-genome-doubling events (0/1/2); see wgd_neutral_split
    contamination_true: float = 0.0  # out-of-patient contamination, constant per sample


def neutral_split(ploidy: int) -> tuple[int, int]:
    """Balanced germline allelic split for a contig of the given ploidy."""
    return ploidy // 2, ploidy - ploidy // 2


def wgd_neutral_split(g1: int, g2: int, w: int) -> tuple[int, int]:
    """Neutral allelic copies after ``w`` whole-genome-doubling events.

    A WGD here is the *gain of one copy* (genome-wide, all contigs), so one WGD
    takes a diploid 1+1 to 3 copies (1+2), and two WGD events reach quadruploid
    2+2 = 4 — not a single ×2 doubling. Extra copies are added major-allele-first,
    then minor, keeping the split as balanced as possible. Applies to allosomes
    too (e.g. male X 0+1 -> 0+2 after one WGD)."""
    a1, a2 = sorted((g1, g2))
    for k in range(w):
        if k % 2 == 0:
            a2 += 1
        else:
            a1 += 1
    return a1, a2


# --------------------------------------------------------------------------- #
# Ground-truth biology
# --------------------------------------------------------------------------- #
@dataclass
class Clone:
    clone_id: str
    parent_id: str | None
    ccf: dict[str, float]        # sample_name -> cancer cell fraction in [0, 1]


@dataclass
class SegmentCN:
    """Allelic copy-number state of a segment in one sample.

    For a subclonal event, ``ccf`` < 1 and the remaining (1 - ccf) fraction of
    tumour cells carries the ``baseline`` state (defaults to balanced diploid for
    the contig's ploidy). total_cn is derived.
    """
    cn_a1: int
    cn_a2: int
    ccf: float = 1.0
    baseline_a1: int | None = None
    baseline_a2: int | None = None

    @property
    def total_cn(self) -> int:
        return self.cn_a1 + self.cn_a2


@dataclass
class Segment:
    contig: str
    start: int
    end: int
    event_type: str              # NEUTRAL / LOH / HZ_DEL / BALANCED_GAIN / ...
    per_sample: dict[str, SegmentCN]
    is_cnloh: bool = False
    n_hets: int = 5              # informative het count on the segment (few is fine)


@dataclass
class GermlineSite:
    contig: str
    position: int
    ref: str
    alt: str
    pop_af: float
    genotype: str                # "0/0" | "0/1" | "1/1" (germline, shared across samples)
    haplotype: str | None = None  # "A" | "B" for hets: which parental haplotype carries alt
    is_snv: bool = True          # False => should be filtered by the VCF SNV mask
    segment_contig: str | None = None  # which CR segment this het informs (for phasing)


@dataclass
class SomaticVariant:
    contig: str
    position: int
    ref: str
    alt: str
    vtype: str                   # "SNV" | "INDEL"
    clone_id: str                # determines per-sample CCF via the clone tree
    multiplicity: int = 1        # copies carrying the alt allele
    in_abs_maf: bool = False     # already present in ABSOLUTE ABS_MAF (vs rescue-only)
    note: str = ""               # e.g. "MT", "unmapped", "zero_cov"


@dataclass
class Patient:
    patient_id: str
    sex: str                     # "XX" | "XY" | "XXY"
    samples: list[Sample]
    clones: list[Clone] = field(default_factory=list)
    segments: list[Segment] = field(default_factory=list)
    germline_sites: list[GermlineSite] = field(default_factory=list)
    somatic_variants: list[SomaticVariant] = field(default_factory=list)
    normal_ploidy: int = 2

    @property
    def nX(self) -> int:
        return self.sex.count("X")

    @property
    def nY(self) -> int:
        return self.sex.count("Y")

    def chromosomal_ploidy(self, contig: str) -> int:
        c = contig[3:] if contig.lower().startswith("chr") else contig
        c = "MT" if c == "M" else c
        if c == "X":
            return self.nX
        if c == "Y":
            return self.nY
        if c == "MT":
            return 1
        return self.normal_ploidy

    @property
    def tumor_samples(self) -> list[Sample]:
        return [s for s in self.samples if not s.is_normal]

    @property
    def normal_samples(self) -> list[Sample]:
        return [s for s in self.samples if s.is_normal]


@dataclass
class Cohort:
    patients: list[Patient]

    def patient(self, patient_id: str) -> Patient:
        return next(p for p in self.patients if p.patient_id == patient_id)
