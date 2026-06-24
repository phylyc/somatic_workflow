"""Meta-tests for the synthetic cohort generator (G0).

Validate that the invented ground truth is internally consistent and that truth
emission is deterministic. These guard the generator itself before any script is
validated against it.
"""

import json

import pandas as pd
import pytest

from tests.synthetic import build_cohort, write_truth
from tests.synthetic.model import CONTIG_LENGTHS


@pytest.fixture(scope="module")
def cohort():
    return build_cohort()


@pytest.mark.unit
def test_three_patients_with_expected_sexes(cohort):
    sexes = {p.patient_id: p.sex for p in cohort.patients}
    assert sexes == {"PXX": "XX", "PXY": "XY", "PXXY": "XXY"}


@pytest.mark.unit
def test_each_patient_has_normal_and_two_tumors(cohort):
    for p in cohort.patients:
        assert len(p.normal_samples) >= 1, p.patient_id
        assert len(p.tumor_samples) >= 1, p.patient_id
    # PXX/PXY engineered with >= 2 tumour timepoints
    assert len(cohort.patient("PXX").tumor_samples) >= 2
    assert len(cohort.patient("PXY").tumor_samples) >= 2


@pytest.mark.unit
def test_allelic_cn_sums_to_total(cohort):
    for p in cohort.patients:
        for seg in p.segments:
            for sample_id, cn in seg.per_sample.items():
                assert cn.cn_a1 + cn.cn_a2 == cn.total_cn
                assert cn.cn_a1 >= 0 and cn.cn_a2 >= 0


@pytest.mark.unit
def test_every_tumor_sample_covered_on_every_segment(cohort):
    # _fill_neutral must leave no (segment, tumor sample) unspecified
    for p in cohort.patients:
        tumor = {s.name for s in p.tumor_samples}
        for seg in p.segments:
            assert tumor.issubset(set(seg.per_sample)), (p.patient_id, seg.contig, seg.start)


@pytest.mark.unit
def test_ccf_in_unit_interval(cohort):
    for p in cohort.patients:
        for c in p.clones:
            for ccf in c.ccf.values():
                assert 0.0 <= ccf <= 1.0
        for seg in p.segments:
            for cn in seg.per_sample.values():
                assert 0.0 <= cn.ccf <= 1.0


@pytest.mark.unit
def test_segments_within_contig_bounds_and_nonoverlapping(cohort):
    for p in cohort.patients:
        by_contig = {}
        for seg in p.segments:
            assert seg.start < seg.end
            assert seg.end <= CONTIG_LENGTHS[seg.contig], (p.patient_id, seg.contig)
            by_contig.setdefault(seg.contig, []).append((seg.start, seg.end))
        for contig, ivals in by_contig.items():
            ivals.sort()
            for (s1, e1), (s2, e2) in zip(ivals, ivals[1:]):
                assert e1 <= s2, f"{p.patient_id} {contig} overlap {e1}>{s2}"


@pytest.mark.unit
def test_allosome_ploidy_matches_karyotype(cohort):
    pxx, pxy, pxxy = (cohort.patient(x) for x in ("PXX", "PXY", "PXXY"))
    assert pxx.chromosomal_ploidy("X") == 2 and pxx.chromosomal_ploidy("Y") == 0
    assert pxy.chromosomal_ploidy("X") == 1 and pxy.chromosomal_ploidy("Y") == 1
    assert pxxy.chromosomal_ploidy("X") == 2 and pxxy.chromosomal_ploidy("Y") == 1


@pytest.mark.unit
def test_pxy_has_wgd_and_diploid_tumors(cohort):
    pxy = cohort.patient("PXY")
    wgd = {s.name: s.wgd for s in pxy.tumor_samples}
    assert any(wgd.values()) and not all(wgd.values()), "need both WGD and diploid tumors"


@pytest.mark.unit
def test_pxy_covers_wgd_counts_0_1_2(cohort):
    levels = {s.wgd for s in cohort.patient("PXY").tumor_samples}
    assert {0, 1, 2}.issubset(levels)


@pytest.mark.unit
def test_wgd_neutral_cn_scales_with_event_count(cohort):
    pxy = cohort.patient("PXY")
    neutral = next(s for s in pxy.segments
                   if s.contig == "1" and s.event_type == "NEUTRAL")
    # diploid 2 -> +1 per WGD event, genome-wide
    assert neutral.per_sample["PXY-T3"].total_cn == 2   # wgd=0
    assert neutral.per_sample["PXY-T4"].total_cn == 3   # wgd=1
    assert neutral.per_sample["PXY-T1"].total_cn == 4   # wgd=2
    # WGD reaches the allosomes too: male X 0+1 -> 0+2 after one WGD
    xseg = next((s for s in pxy.segments if s.contig == "X"
                 and "PXY-T4" in s.per_sample), None)
    if xseg is not None:
        assert xseg.per_sample["PXY-T4"].total_cn == 2


@pytest.mark.unit
def test_multi_run_and_assay_coverage(cohort):
    assays = {r.assay.value for p in cohort.patients for s in p.samples for r in s.runs}
    assert assays == {"WES", "ULP_WGS", "DEEP_PANEL"}
    # at least one sample carries multiple sequencing runs
    multi = [s.name for p in cohort.patients for s in p.samples if len(s.runs) > 1]
    assert multi


@pytest.mark.unit
def test_event_types_and_edge_variants_present(cohort):
    events = {seg.event_type for p in cohort.patients for seg in p.segments}
    for required in ["NEUTRAL", "LOH", "HZ_DEL", "IMBALANCED_GAIN", "HIGH_AMP",
                     "SUBCLONAL_GAIN", "FOCAL_BALANCED_HIGHCN"]:
        assert required in events, required
    notes = {v.note for p in cohort.patients for v in p.somatic_variants}
    assert {"MT", "unmapped", "zero_cov"}.issubset(notes)
    # both rescue-only and already-in-ABS_MAF somatic variants exist
    flags = {v.in_abs_maf for p in cohort.patients for v in p.somatic_variants}
    assert flags == {True, False}
    # a non-SNV germline row exists to exercise the VCF filter
    assert any(not g.is_snv for p in cohort.patients for g in p.germline_sites)


@pytest.mark.io
def test_write_truth_emits_tables_and_manifest(tmp_path):
    written = write_truth(build_cohort(), tmp_path)
    for name in ["patient", "clones", "segments", "germline_sites",
                 "somatic_variants", "manifest"]:
        assert written[name].exists()
    seg = pd.read_csv(written["segments"], sep="\t")
    assert (seg["cn_a1"] + seg["cn_a2"] == seg["total_cn"]).all()
    manifest = json.loads(written["manifest"].read_text())
    assert {p["patient_id"] for p in manifest["patients"]} == {"PXX", "PXY", "PXXY"}


@pytest.mark.regression
def test_committed_data_matches_regeneration(synth):
    """The inspectable committed tree (python/tests/synthetic/data/) must equal a
    fresh regeneration byte-for-byte, so it can't silently go stale and the user
    inspects exactly what the suite exercises. Regenerate with:
        python -m tests.synthetic.make_synthetic_testing_data --outdir tests/synthetic/data
    """
    import pathlib

    committed = pathlib.Path(__file__).resolve().parent / "synthetic" / "data"
    if not committed.exists():
        pytest.skip("committed data tree absent; regenerate into tests/synthetic/data")

    def tree(base):
        return {p.relative_to(base) for p in base.rglob("*")
                if p.is_file() and p.name != "files.json"}

    ct, rt = tree(committed), tree(synth.dir)
    assert ct == rt, f"file set differs: only-committed={ct - rt}, only-regen={rt - ct}"
    for rel in sorted(ct):
        assert (committed / rel).read_bytes() == (synth.dir / rel).read_bytes(), \
            f"committed data is stale at {rel}; regenerate tests/synthetic/data"


@pytest.mark.io
def test_truth_emission_deterministic(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    write_truth(build_cohort(), a)
    write_truth(build_cohort(), b)
    # manifest is byte-identical; truth tables equal as frames
    assert (a / "manifest.json").read_text() == (b / "manifest.json").read_text()
    for name in ["patient", "clones", "segments", "germline_sites", "somatic_variants"]:
        fa = pd.read_csv(a / "truth" / f"{name}.truth.tsv", sep="\t")
        fb = pd.read_csv(b / "truth" / f"{name}.truth.tsv", sep="\t")
        pd.testing.assert_frame_equal(fa, fb)
