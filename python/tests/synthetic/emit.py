"""Emit upstream tool-output files from the simulated cohort.

G1 scope: reference .dict, GetPileupSummaries pileups, DenoiseReadCounts copy
ratios, CalculateContamination contamination + .segments, ModelSegments
modelFinal.seg + af.param, the germline reference VCF, and a genotyped gvcf.
All gzip-friendly TSVs/text; deterministic given the RNG seed.
"""

from __future__ import annotations

import gzip
import pathlib

import numpy as np
import pandas as pd

from .assays import PROFILES, target_bins
from .model import CONTIGS, Assay, Patient, Sample, chr_name
from .simulate import (
    cn_state_at, delta_b, draw_counts, draw_log2cr, expected_log2cr, expected_maf,
    germline_site_vaf, observed_allelic_cn, observed_total_cn, sample_segmentation,
    somatic_alt_fraction,
)

ABS_COLUMNS = [
    "sample", "Chromosome", "Start.bp", "End.bp", "n_probes", "length", "seg_sigma",
    "W", "total_copy_ratio", "modal_total_cn", "expected_total_cn", "total_HZ",
    "total_amp", "corrected_total_cn", "rescaled_total_cn", "bi.allelic", "copy.ratio",
    "hscr.a1", "hscr.a2", "modal.a1", "modal.a2", "expected.a1", "expected.a2",
    "subclonal.a1", "subclonal.a2", "cancer.cell.frac.a1", "ccf.ci95.low.a1",
    "ccf.ci95.high.a1", "cancer.cell.frac.a2", "ccf.ci95.low.a2", "ccf.ci95.high.a2",
    "LOH", "HZ", "SC_HZ", "amp.a1", "amp.a2", "rescaled.cn.a1", "rescaled.cn.a2",
]

LN2 = float(np.log(2))

ERROR_RATE = 0.005
PILEUP_HEADER = "contig\tposition\tref_count\talt_count\tother_alt_count\tallele_frequency\n"


# --------------------------------------------------------------------------- #
# reference dict
# --------------------------------------------------------------------------- #
def write_dict(path, chr_prefixed: bool):
    lines = ["@HD\tVN:1.6"]
    for contig, length in CONTIGS:
        lines.append(f"@SQ\tSN:{chr_name(contig, chr_prefixed)}\tLN:{length}")
    pathlib.Path(path).write_text("\n".join(lines) + "\n")


def _open(path):
    path = str(path)
    return gzip.open(path, "wt") if path.endswith(".gz") else open(path, "w")


# --------------------------------------------------------------------------- #
# pileups (GetPileupSummaries) — one per (sample, run) with allelic data
# --------------------------------------------------------------------------- #
def write_pileup(path, patient: Patient, sample: Sample, assay: Assay, rng,
                 chr_prefixed: bool, panel_only: bool):
    prof = PROFILES[assay]
    rows = []
    panel = {(c, s, e) for c, s, e, _ in target_bins("PANEL")} if panel_only else None
    for site in patient.germline_sites:
        if not site.is_snv:
            continue
        if panel_only and not _in_bins(site.contig, site.position, panel):
            continue
        vaf = germline_site_vaf(patient, sample, site.contig, site.position,
                                site.genotype, site.haplotype)
        counts = draw_counts(rng, vaf, prof.het_depth, ERROR_RATE,
                             sample.contamination_true, contam_vaf=0.5)
        if counts is None:
            continue
        ref, alt, other = counts
        rows.append((chr_name(site.contig, chr_prefixed), site.position, ref, alt, other,
                     round(site.pop_af, 4)))
    with _open(path) as fh:
        fh.write(f"#<METADATA>SAMPLE={sample.name}\n")
        fh.write(PILEUP_HEADER)
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    return len(rows)


def _in_bins(contig, position, bins):
    return any(c == contig and s <= position <= e for c, s, e in bins)


# --------------------------------------------------------------------------- #
# denoised copy ratios (DenoiseReadCounts) — one per (sample, run) with CR
# --------------------------------------------------------------------------- #
def write_denoised_cr(path, patient: Patient, sample: Sample, target_set, rng,
                      chr_prefixed: bool, cr_depth: int):
    rows = []
    for contig, start, end, _idx in target_bins(target_set):
        mid = (start + end) // 2
        cn = cn_state_at(patient, sample.name, contig, mid)
        cp = patient.chromosomal_ploidy(contig)
        log2cr = expected_log2cr(cn.total_cn, cp, sample.purity, sample.ploidy,
                                 patient.normal_ploidy)
        log2cr = draw_log2cr(rng, log2cr, cr_depth)
        if np.isnan(log2cr):
            continue
        rows.append((chr_name(contig, chr_prefixed), start, end, round(log2cr, 6)))
    with _open(path) as fh:
        for contig, length in CONTIGS:
            fh.write(f"@SQ\tSN:{chr_name(contig, chr_prefixed)}\tLN:{length}\n")
        fh.write("CONTIG\tSTART\tEND\tLOG2_COPY_RATIO\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    return len(rows)


# --------------------------------------------------------------------------- #
# contamination + .segments (CalculateContamination)
# --------------------------------------------------------------------------- #
def write_contamination(path, sample: Sample):
    with open(path, "w") as fh:
        fh.write("sample\tcontamination\tcontamination_error\n")
        err = max(sample.contamination_true * 0.1, 1e-3)
        fh.write(f"{sample.name}\t{sample.contamination_true}\t{err}\n")


def write_segments(path, patient: Patient, sample: Sample, chr_prefixed: bool):
    rows = []
    for s in sample_segmentation(patient, sample.name):
        maf = expected_maf(s["cn"], sample.purity)
        rows.append((chr_name(s["contig"], chr_prefixed), s["start"], s["end"], round(maf, 5)))
    with open(path, "w") as fh:
        fh.write(f"#<METADATA>SAMPLE={sample.name}\n")
        fh.write("contig\tstart\tend\tminor_allele_fraction\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")


# --------------------------------------------------------------------------- #
# ModelSegments modelFinal.seg + af.param
# --------------------------------------------------------------------------- #
MS_COLS = ["CONTIG", "START", "END", "NUM_POINTS_COPY_RATIO", "NUM_POINTS_ALLELE_FRACTION",
           "LOG2_COPY_RATIO_POSTERIOR_10", "LOG2_COPY_RATIO_POSTERIOR_50", "LOG2_COPY_RATIO_POSTERIOR_90",
           "MINOR_ALLELE_FRACTION_POSTERIOR_10", "MINOR_ALLELE_FRACTION_POSTERIOR_50",
           "MINOR_ALLELE_FRACTION_POSTERIOR_90", "CALL"]


def write_model_seg(path, patient: Patient, sample: Sample, chr_prefixed: bool,
                    include_call: bool = True):
    rows = []
    for s in sample_segmentation(patient, sample.name):
        cn = s["cn"]
        cp = patient.chromosomal_ploidy(s["contig"])
        obs_total = observed_total_cn(cn)  # CCF-mixed observed copy ratio
        log2_50 = expected_log2cr(obs_total, cp, sample.purity, sample.ploidy, patient.normal_ploidy)
        if np.isnan(log2_50):
            continue
        maf_50 = expected_maf(cn, sample.purity)
        call = "0"
        if include_call:
            call = "+" if obs_total > cp else "-" if obs_total < cp else "0"
        rows.append([chr_name(s["contig"], chr_prefixed), s["start"], s["end"],
                     max(s["n_hets"], 1), s["n_hets"],
                     round(log2_50 - 0.05, 6), round(log2_50, 6), round(log2_50 + 0.05, 6),
                     round(max(maf_50 - 0.02, 0.0), 5), round(maf_50, 5),
                     round(min(maf_50 + 0.02, 0.5), 5), call])
    cols = MS_COLS if include_call else MS_COLS[:-1]
    rows = [r[:len(cols)] for r in rows]
    with open(path, "w") as fh:
        fh.write("@HD\tVN:1.6\n")
        fh.write(f"@RG\tID:GATKCopyNumber\tSM:{sample.name}\n")
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")


def write_af_param(path, mean_bias: float = 1.05):
    with open(path, "w") as fh:
        fh.write("@param\n")
        fh.write("PARAMETER_NAME\tPOSTERIOR_10\tPOSTERIOR_50\tPOSTERIOR_90\n")
        fh.write(f"MEAN_BIAS\t{mean_bias - 0.05}\t{mean_bias}\t{mean_bias + 0.05}\n")
        fh.write("BIAS_VARIANCE\t0.005\t0.01\t0.02\n")
        fh.write("OUTLIER_PROBABILITY\t0.01\t0.02\t0.03\n")


# --------------------------------------------------------------------------- #
# germline reference VCF + genotyped gvcf
# --------------------------------------------------------------------------- #
def write_reference_vcf(path, patient: Patient, chr_prefixed: bool):
    lines = ["##fileformat=VCFv4.2"]
    for contig, length in CONTIGS:
        lines.append(f"##contig=<ID={chr_name(contig, chr_prefixed)},length={length}>")
    lines.append("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]))
    for g in patient.germline_sites:
        info = f"AF={g.pop_af}"
        lines.append("\t".join([chr_name(g.contig, chr_prefixed), str(g.position), ".",
                                g.ref, g.alt, ".", "PASS", info, "GT"]))
    pathlib.Path(path).write_text("\n".join(lines) + "\n")


# --------------------------------------------------------------------------- #
# caller MAFs (Mutect/annotation) — somatic SNV / INDEL
# --------------------------------------------------------------------------- #
MAF_COLS = ["Hugo_Symbol", "Chromosome", "Start_position", "End_position",
            "Reference_Allele", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode",
            "t_alt_count", "t_ref_count", "n_alt_count", "n_ref_count", "tumor_f"]


def write_caller_maf(path, patient: Patient, sample: Sample, rng, vtype: str,
                     chr_prefixed: bool, depth: int = 120):
    clone_ccf = {c.clone_id: c.ccf for c in patient.clones}
    rows = []
    for v in patient.somatic_variants:
        if v.vtype != vtype:
            continue
        ccf = clone_ccf.get(v.clone_id, {}).get(sample.name, 0.0)
        vaf = somatic_alt_fraction(patient, sample, v.contig, v.position, ccf, v.multiplicity)
        d = 0 if v.note == "zero_cov" else depth
        if d <= 0:
            t_alt, t_ref = 0, 0
        else:
            n = int(rng.poisson(d))
            t_alt = int(rng.binomial(n, min(max(vaf, 0.0), 1.0)))
            t_ref = n - t_alt
        n_ref = int(rng.poisson(depth)) if depth > 0 else 0
        tumor_f = (t_alt / (t_alt + t_ref)) if (t_alt + t_ref) > 0 else ""
        end = v.position + max(len(v.ref) - 1, 0)
        rows.append(["GENE", chr_name(v.contig, chr_prefixed), v.position, end,
                     v.ref, v.alt, sample.name, t_alt, t_ref, 0, n_ref,
                     round(tumor_f, 5) if tumor_f != "" else ""])
    with open(path, "w") as fh:
        fh.write("\t".join(MAF_COLS) + "\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")


# --------------------------------------------------------------------------- #
# ABSOLUTE segtab (the has_absolute path of map_to_absolute)
# --------------------------------------------------------------------------- #
def write_absolute_segtab(path, patient: Patient, sample: Sample, chr_prefixed: bool):
    """Emit an ABSOLUTE-style segtab from truth, on the same genome-wide
    segmentation as the acs/somix seg so the two concat cleanly. Provides the
    absolute CN calls so map_to_absolute maps them through (no de-novo rescue)."""
    segs = sample_segmentation(patient, sample.name)
    total_len = sum(s["end"] - s["start"] for s in segs) or 1
    rows = []
    for s in segs:
        cp = patient.chromosomal_ploidy(s["contig"])
        if cp <= 0:
            continue  # match acs seg, which drops ploidy-0 contigs
        cn = s["cn"]
        delta, b = delta_b(cp, sample.purity, sample.ploidy, patient.normal_ploidy)
        a1, a2, tot = cn.cn_a1, cn.cn_a2, cn.total_cn   # integer (clonal) call = event state
        obs1, obs2 = observed_allelic_cn(cn)            # CCF-mixed observed copy ratio
        obs_tot = obs1 + obs2
        length = s["end"] - s["start"]
        row = {
            "sample": sample.name, "Chromosome": chr_name(s["contig"], chr_prefixed),
            "Start.bp": s["start"], "End.bp": s["end"], "n_probes": 10, "length": length,
            "seg_sigma": 0.05, "W": round(length / total_len, 6),
            "total_copy_ratio": round(obs_tot, 6), "modal_total_cn": tot, "expected_total_cn": float(tot),
            "total_HZ": int(tot == 0), "total_amp": int(tot >= 7),
            "corrected_total_cn": float(tot), "rescaled_total_cn": float(tot),
            "bi.allelic": int(a1 > 0 and a2 > 0), "copy.ratio": round(obs_tot, 6),
            "hscr.a1": round(obs1 * delta * cp + b, 6), "hscr.a2": round(obs2 * delta * cp + b, 6),
            "modal.a1": a1, "modal.a2": a2, "expected.a1": float(a1), "expected.a2": float(a2),
            "subclonal.a1": int(cn.ccf < 1), "subclonal.a2": int(cn.ccf < 1),
            "cancer.cell.frac.a1": round(cn.ccf, 4), "ccf.ci95.low.a1": round(max(cn.ccf - 0.1, 0), 4),
            "ccf.ci95.high.a1": round(min(cn.ccf + 0.1, 1), 4),
            "cancer.cell.frac.a2": round(cn.ccf, 4), "ccf.ci95.low.a2": round(max(cn.ccf - 0.1, 0), 4),
            "ccf.ci95.high.a2": round(min(cn.ccf + 0.1, 1), 4),
            "LOH": int((a1 == 0 or a2 == 0) and tot > 0), "HZ": int(tot == 0), "SC_HZ": 0,
            "amp.a1": int(a1 >= 4), "amp.a2": int(a2 >= 4),
            "rescaled.cn.a1": a1, "rescaled.cn.a2": a2,
        }
        rows.append([row[c] for c in ABS_COLUMNS])
    with open(path, "w") as fh:
        fh.write("# ABSOLUTE segtab (synthetic)\n")
        fh.write("\t".join(ABS_COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")


# --------------------------------------------------------------------------- #
# SOMIX segmentation (alternative CR-seg input to map_to_absolute)
# --------------------------------------------------------------------------- #
SOMIX_COLS = ["contig", "Start.bp", "End.bp", "n_markers", "n_snps", "f_MAP",
              "log_tCR", "sem_log_tCR", "f_hessian", "sample_id"]


def write_somix_seg(path, patient: Patient, sample: Sample, chr_prefixed: bool):
    rows = []
    for s in sample_segmentation(patient, sample.name):
        cn = s["cn"]
        cp = patient.chromosomal_ploidy(s["contig"])
        log2cr = expected_log2cr(observed_total_cn(cn), cp, sample.purity, sample.ploidy, patient.normal_ploidy)
        if np.isnan(log2cr):
            continue
        # map_to_cn somix path: tau = exp(log_tCR) * chr_ploidy = chr_ploidy * 2^log2cr
        log_tcr = log2cr * LN2
        f = expected_maf(cn, sample.purity)
        rows.append([chr_name(s["contig"], chr_prefixed), s["start"], s["end"],
                     max(s["n_hets"], 1), s["n_hets"], round(f, 5),
                     round(log_tcr, 6), 0.05, -1000.0, sample.name])
    with open(path, "w") as fh:
        fh.write("\t".join(SOMIX_COLS) + "\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")


def write_genotyped_gvcf(path, patient: Patient, chr_prefixed: bool):
    """Genotyped germline gvcf with phased GT (0|1 / 1|0) at hets — for
    pileup_to_allelic_counts. Phase from the truth haplotype."""
    lines = ["##fileformat=VCFv4.2"]
    for contig, length in CONTIGS:
        lines.append(f"##contig=<ID={chr_name(contig, chr_prefixed)},length={length}>")
    lines.append("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                            "INFO", "FORMAT", patient.patient_id]))
    for g in patient.germline_sites:
        if not g.is_snv:
            continue
        if g.genotype == "0/1":
            gt = "0|1" if g.haplotype == "B" else "1|0"
        else:
            gt = g.genotype
        lines.append("\t".join([chr_name(g.contig, chr_prefixed), str(g.position), ".",
                                g.ref, g.alt, ".", "PASS", ".", "GT", gt]))
    with _open(path) as fh:
        fh.write("\n".join(lines) + "\n")
