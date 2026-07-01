"""Drive emission of the full synthetic data tree from the cohort.

Deterministic given ``seed``: a single RNG is threaded through emission in a
fixed (patient, sample, run) order. Returns a files index (also written to
files.json) that tests use to locate inputs.
"""

from __future__ import annotations

import json
import pathlib

import numpy as np

from . import emit
from .assays import PROFILES
from .model import Assay, Cohort


def generate_cohort(cohort: Cohort, outdir, seed: int = 0,
                    chr_prefixed: bool = False) -> dict:
    outdir = pathlib.Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)

    emit.write_dict(outdir / "reference.dict", chr_prefixed=False)
    emit.write_dict(outdir / "reference.chr.dict", chr_prefixed=True)

    index: dict = {"chr_prefixed": chr_prefixed, "patients": {}}

    for patient in cohort.patients:
        pdir = outdir / patient.patient_id
        pdir.mkdir(parents=True, exist_ok=True)
        pentry: dict = {"reference_vcf": None, "gvcf": None, "samples": {}}

        ref_vcf = pdir / f"{patient.patient_id}.reference.vcf"
        emit.write_reference_vcf(ref_vcf, patient, chr_prefixed)
        pentry["reference_vcf"] = str(ref_vcf)

        # plain text (not .gz): deterministic bytes for the drift check + inspectable.
        # pileup_to_allelic_counts.read_gvcf handles plain and gzipped alike.
        gvcf = pdir / f"{patient.patient_id}.genotyped.g.vcf"
        emit.write_genotyped_gvcf(gvcf, patient, chr_prefixed)
        pentry["gvcf"] = str(gvcf)

        for sample in patient.samples:
            sentry: dict = {"runs": [], "contamination": None, "segments": None,
                            "model_seg": None, "af_param": None}

            contam = pdir / f"{sample.name}.contamination"
            emit.write_contamination(contam, sample)
            sentry["contamination"] = str(contam)

            seg = pdir / f"{sample.name}.segments"
            emit.write_segments(seg, patient, sample, chr_prefixed)
            sentry["segments"] = str(seg)

            mseg = pdir / f"{sample.name}.modelFinal.seg"
            emit.write_model_seg(mseg, patient, sample, chr_prefixed)
            sentry["model_seg"] = str(mseg)

            afp = pdir / f"{sample.name}.modelFinal.af.param"
            emit.write_af_param(afp)
            sentry["af_param"] = str(afp)

            if not sample.is_normal:
                abs_segtab = pdir / f"{sample.name}.absolute.segtab.tsv"
                emit.write_absolute_segtab(abs_segtab, patient, sample, chr_prefixed)
                sentry["absolute_segtab"] = str(abs_segtab)

                somix = pdir / f"{sample.name}.somix.seg"
                emit.write_somix_seg(somix, patient, sample, chr_prefixed)
                sentry["somix_seg"] = str(somix)

                snv = pdir / f"{sample.name}.snv.maf"
                emit.write_caller_maf(snv, patient, sample, rng, "SNV", chr_prefixed)
                sentry["snv_maf"] = str(snv)

                indel = pdir / f"{sample.name}.indel.maf"
                emit.write_caller_maf(indel, patient, sample, rng, "INDEL", chr_prefixed)
                sentry["indel_maf"] = str(indel)

            for ridx, run in enumerate(sample.runs):
                prof = PROFILES[run.assay]
                rentry = {"assay": run.assay.value, "target_set": run.target_set,
                          "pileup": None, "denoised_cr": None}
                if prof.has_allelic:
                    pu = pdir / f"{sample.name}.{run.assay.value}.{ridx}.pileup"
                    emit.write_pileup(pu, patient, sample, run.assay, rng,
                                      chr_prefixed, panel_only=(run.assay == Assay.DEEP_PANEL))
                    rentry["pileup"] = str(pu)
                if prof.has_cr:
                    cr = pdir / f"{sample.name}.{run.assay.value}.{ridx}.denoised_CR.tsv"
                    target = run.target_set if run.assay == Assay.WES else None
                    emit.write_denoised_cr(cr, patient, sample, target, rng,
                                           chr_prefixed, prof.cr_depth)
                    rentry["denoised_cr"] = str(cr)
                sentry["runs"].append(rentry)

            pentry["samples"][sample.name] = sentry
        index["patients"][patient.patient_id] = pentry

    (outdir / "files.json").write_text(json.dumps(index, indent=2) + "\n")
    return index
