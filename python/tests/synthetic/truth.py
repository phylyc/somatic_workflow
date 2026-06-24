"""Emit the ground-truth answer key (and cohort manifest) for the synthetic cohort.

Truth tables are written as gzip-compressed TSVs; the manifest as JSON. These are
the answer key the Python scripts are validated against (semantic golden) and the
input description that drives validate_inputs scenarios.
"""

from __future__ import annotations

import json
import pathlib

import pandas as pd

from .model import CONTIGS, Cohort, Patient


def _patient_rows(p: Patient) -> list[dict]:
    rows = []
    for s in p.samples:
        rows.append(dict(
            patient_id=p.patient_id, sample_id=s.name, is_normal=s.is_normal,
            timepoint=s.timepoint, sex=p.sex, purity=s.purity, ploidy=s.ploidy,
            wgd=s.wgd, contamination_true=s.contamination_true,
        ))
    return rows


def _clone_rows(p: Patient) -> list[dict]:
    rows = []
    for c in p.clones:
        for sample_id, ccf in c.ccf.items():
            rows.append(dict(patient_id=p.patient_id, clone_id=c.clone_id,
                             parent_clone_id=c.parent_id, sample_id=sample_id, ccf=ccf))
    return rows


def _segment_rows(p: Patient) -> list[dict]:
    rows = []
    for seg in p.segments:
        for sample_id, cn in seg.per_sample.items():
            rows.append(dict(
                patient_id=p.patient_id, contig=seg.contig, start=seg.start, end=seg.end,
                event_type=seg.event_type, is_cnloh=seg.is_cnloh, n_hets=seg.n_hets,
                sample_id=sample_id, cn_a1=cn.cn_a1, cn_a2=cn.cn_a2, total_cn=cn.total_cn,
                ccf=cn.ccf, baseline_a1=cn.baseline_a1, baseline_a2=cn.baseline_a2,
            ))
    return rows


def _germline_rows(p: Patient) -> list[dict]:
    return [dict(patient_id=p.patient_id, contig=g.contig, position=g.position,
                 ref=g.ref, alt=g.alt, pop_af=g.pop_af, genotype=g.genotype,
                 haplotype=g.haplotype, is_snv=g.is_snv, segment_contig=g.segment_contig)
            for g in p.germline_sites]


def _somatic_rows(p: Patient) -> list[dict]:
    clone_ccf = {c.clone_id: c.ccf for c in p.clones}
    rows = []
    for v in p.somatic_variants:
        for s in p.tumor_samples:
            rows.append(dict(
                patient_id=p.patient_id, contig=v.contig, position=v.position,
                ref=v.ref, alt=v.alt, vtype=v.vtype, clone_id=v.clone_id,
                multiplicity=v.multiplicity, in_abs_maf=v.in_abs_maf, note=v.note,
                sample_id=s.name, ccf=clone_ccf.get(v.clone_id, {}).get(s.name, 0.0),
            ))
    return rows


def _manifest(cohort: Cohort) -> dict:
    patients = []
    for p in cohort.patients:
        samples = []
        for s in p.samples:
            runs = []
            for r in s.runs:
                tcr, acr = r.resolved_use()
                runs.append(dict(assay=r.assay.value, target_set=r.target_set,
                                 use_for_tCR=tcr, use_for_aCR=acr,
                                 is_paired_end=r.is_paired_end))
            samples.append(dict(sample_id=s.name, is_normal=s.is_normal,
                                timepoint=s.timepoint, purity=s.purity,
                                ploidy=s.ploidy, wgd=s.wgd, runs=runs))
        patients.append(dict(
            patient_id=p.patient_id, sex=p.sex, normal_ploidy=p.normal_ploidy,
            normal_sample_names=[s.name for s in p.normal_samples],
            samples=samples,
        ))
    return dict(reference={"contigs": [{"name": c, "length": ln} for c, ln in CONTIGS]},
                patients=patients)


def write_truth(cohort: Cohort, outdir: str | pathlib.Path) -> dict[str, pathlib.Path]:
    outdir = pathlib.Path(outdir)
    truth_dir = outdir / "truth"
    truth_dir.mkdir(parents=True, exist_ok=True)

    tables = {
        "patient": [r for p in cohort.patients for r in _patient_rows(p)],
        "clones": [r for p in cohort.patients for r in _clone_rows(p)],
        "segments": [r for p in cohort.patients for r in _segment_rows(p)],
        "germline_sites": [r for p in cohort.patients for r in _germline_rows(p)],
        "somatic_variants": [r for p in cohort.patients for r in _somatic_rows(p)],
    }
    written = {}
    for name, rows in tables.items():
        # plain TSV (not gzip) so the committed ground truth is easy to inspect
        path = truth_dir / f"{name}.truth.tsv"
        pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
        written[name] = path

    manifest_path = outdir / "manifest.json"
    manifest_path.write_text(json.dumps(_manifest(cohort), indent=2) + "\n")
    written["manifest"] = manifest_path
    return written
