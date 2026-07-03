"""Forward simulation: ground truth + RNG -> per-locus counts and copy ratios.

The copy-ratio formula is the exact inverse of acs_conversion + map_to_cn, so the
downstream pipeline recovers the truth CN (modulo clustering). Allele fractions
follow the standard allelic-CN -> BAF mixture with tumour purity and a constant
out-of-patient contamination term.
"""

from __future__ import annotations

import numpy as np

from .model import Patient, Sample, SegmentCN, neutral_split, wgd_neutral_split

# Map a germline het's parental haplotype label to an allele index (a1/a2).
HAP_TO_ALLELE = {"A": 1, "B": 2}


def _sample(patient: Patient, name: str) -> Sample:
    return next(s for s in patient.samples if s.name == name)


def cn_state_at(patient: Patient, sample_name: str, contig: str, position: int) -> SegmentCN:
    """Allelic CN state at a locus in a sample (segment overlap, else neutral)."""
    for seg in patient.segments:
        if seg.contig == contig and seg.start <= position <= seg.end and sample_name in seg.per_sample:
            return seg.per_sample[sample_name]
    g1, g2 = neutral_split(patient.chromosomal_ploidy(contig))
    a1, a2 = wgd_neutral_split(g1, g2, _sample(patient, sample_name).wgd)
    return SegmentCN(a1, a2)


def _effective_copies(cn: SegmentCN, allele: int) -> tuple[float, float]:
    """(alt-allele copies, total copies) in tumour cells, mixing the (1-ccf)
    baseline population for subclonal segments."""
    event_alt = cn.cn_a1 if allele == 1 else cn.cn_a2
    event_tot = cn.total_cn
    if cn.ccf >= 1.0 or (cn.baseline_a1 is None and cn.baseline_a2 is None):
        return event_alt, event_tot
    b1 = 0 if cn.baseline_a1 is None else cn.baseline_a1
    b2 = 0 if cn.baseline_a2 is None else cn.baseline_a2
    base_alt = b1 if allele == 1 else b2
    base_tot = b1 + b2
    alt = cn.ccf * event_alt + (1 - cn.ccf) * base_alt
    tot = cn.ccf * event_tot + (1 - cn.ccf) * base_tot
    return alt, tot


def expected_vaf(alt_copies: float, total_copies: float, purity: float,
                 germline_alt: int = 1, germline_total: int = 2) -> float:
    """Allelic-CN -> observed alt fraction with normal-cell admixture (purity)."""
    num = purity * alt_copies + (1 - purity) * germline_alt
    den = purity * total_copies + (1 - purity) * germline_total
    return float(num / den) if den > 0 else 0.0


def het_alt_fraction(patient: Patient, sample: Sample, contig: str, position: int,
                     haplotype: str) -> float:
    cn = cn_state_at(patient, sample.name, contig, position)
    allele = HAP_TO_ALLELE.get(haplotype, 1)
    alt, tot = _effective_copies(cn, allele)
    return expected_vaf(alt, tot, sample.purity)


def germline_site_vaf(patient: Patient, sample: Sample, contig: str, position: int,
                      genotype: str, haplotype: str | None) -> float:
    """Expected alt fraction at a germline SNP site for GetPileupSummaries."""
    cn = cn_state_at(patient, sample.name, contig, position)
    if genotype == "0/0":
        return expected_vaf(0, cn.total_cn, sample.purity, germline_alt=0)
    if genotype == "1/1":
        return expected_vaf(cn.total_cn, cn.total_cn, sample.purity, germline_alt=2)
    return het_alt_fraction(patient, sample, contig, position, haplotype or "A")


def somatic_alt_fraction(patient: Patient, sample: Sample, contig: str, position: int,
                         ccf: float, multiplicity: int) -> float:
    """Somatic variant VAF: ccf*purity fraction of cells carry `multiplicity`
    alt copies out of the local total copies."""
    cn = cn_state_at(patient, sample.name, contig, position)
    total = cn.total_cn
    num = sample.purity * ccf * multiplicity
    den = sample.purity * total + (1 - sample.purity) * 2
    return float(num / den) if den > 0 else 0.0


def sample_segmentation(patient: Patient, sample_name: str):
    """Genome-wide segmentation for a sample: the explicit event segments plus
    neutral filler covering the gaps, so the length-weighted mean CN ≈ ploidy
    (what real ModelSegments emits, and what the de-novo map_to_cn normalization
    assumes). Yields dicts with contig/start/end/cn/n_hets/is_cnloh/event_type."""
    from .model import CONTIGS

    out = []
    for contig, length in CONTIGS:
        if contig == "MT":
            continue
        events = sorted((s for s in patient.segments
                         if s.contig == contig and sample_name in s.per_sample),
                        key=lambda s: s.start)
        pos = 1
        for ev in events:
            if ev.start > pos:
                mid = (pos + ev.start - 1) // 2
                out.append(dict(contig=contig, start=pos, end=ev.start - 1,
                                cn=cn_state_at(patient, sample_name, contig, mid),
                                n_hets=5, is_cnloh=False, event_type="NEUTRAL"))
            out.append(dict(contig=contig, start=ev.start, end=ev.end,
                            cn=ev.per_sample[sample_name], n_hets=ev.n_hets,
                            is_cnloh=ev.is_cnloh, event_type=ev.event_type))
            pos = ev.end + 1
        if pos <= length:
            mid = (pos + length) // 2
            out.append(dict(contig=contig, start=pos, end=length,
                            cn=cn_state_at(patient, sample_name, contig, mid),
                            n_hets=5, is_cnloh=False, event_type="NEUTRAL"))
    return out


def expected_log2cr(cn_total: float, chr_ploidy: int, purity: float, ploidy: float,
                    normal_ploidy: int = 2) -> float:
    """Inverse of acs_conversion(tau)->map_to_cn(CN): generating CR with this makes
    the pipeline recover ``cn_total``."""
    if chr_ploidy <= 0:
        return float("nan")
    D = chr_ploidy * ((1 - purity) + purity * ploidy / normal_ploidy)
    val = (purity * cn_total + (1 - purity)) / D
    return float(np.log2(max(val, 1e-6)))


def delta_b(chr_ploidy: int, purity: float, ploidy: float, normal_ploidy: int = 2):
    """The (delta, b) of map_to_cn's affine CN<->copy-ratio map, so allelic copy
    numbers can be turned into haplotype-specific copy ratios (hscr)."""
    D = (1 - purity) * chr_ploidy + purity * ploidy * chr_ploidy / normal_ploidy
    if D <= 0:
        return float("nan"), float("nan")
    return purity / D, (1 - purity) * chr_ploidy / D


def observed_allelic_cn(cn: SegmentCN) -> tuple[float, float]:
    """The CCF-mixed *observed* allelic copy number (what the copy-ratio signal
    reflects): ccf*event + (1-ccf)*baseline. For clonal segments (ccf>=1) this is
    just the event state. Subclonal segments therefore carry a fractional observed
    copy number, which is what lets a CCF estimator recover ccf<1 from the data."""
    if cn.ccf >= 1.0 or (cn.baseline_a1 is None and cn.baseline_a2 is None):
        return float(cn.cn_a1), float(cn.cn_a2)
    b1 = 0.0 if cn.baseline_a1 is None else float(cn.baseline_a1)
    b2 = 0.0 if cn.baseline_a2 is None else float(cn.baseline_a2)
    return (cn.ccf * cn.cn_a1 + (1 - cn.ccf) * b1,
            cn.ccf * cn.cn_a2 + (1 - cn.ccf) * b2)


def observed_total_cn(cn: SegmentCN) -> float:
    a1, a2 = observed_allelic_cn(cn)
    return a1 + a2


def expected_maf(cn: SegmentCN, purity: float) -> float:
    """Minor allele fraction observed for a segment (for ModelSegments MAF),
    using the CCF-mixed observed allelic copy number."""
    a1, a2 = observed_allelic_cn(cn)
    minor, total = min(a1, a2), a1 + a2
    f = expected_vaf(minor, total, purity, germline_alt=1, germline_total=2)
    return min(f, 1 - f)


def draw_counts(rng, exp_vaf: float, depth: int, error: float,
                contamination: float, contam_vaf: float) -> tuple[int, int, int] | None:
    """Draw (ref, alt, other_alt) at a locus, or None if no coverage (ULP)."""
    if depth <= 0:
        return None
    n = int(rng.poisson(depth))
    if n == 0:
        return None
    eff = (1 - contamination) * exp_vaf + contamination * contam_vaf
    eff = min(max(eff, 0.0), 1.0)
    alt = int(rng.binomial(n, eff))
    other = int(rng.binomial(n, error))
    ref = max(n - alt - other, 0)
    return ref, alt, other


def draw_log2cr(rng, exp_log2cr: float, cr_depth: int) -> float:
    """Add depth-scaled Gaussian noise to the expected log2 copy ratio."""
    if np.isnan(exp_log2cr):
        return exp_log2cr
    sigma = 0.05 + 2.0 / max(cr_depth, 1)
    return float(exp_log2cr + rng.normal(0.0, sigma))
