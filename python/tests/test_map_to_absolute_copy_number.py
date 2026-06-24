"""Tests for map_to_absolute_copy_number.py.

Two layers:
  * CHARACTERIZATION (regression-style): run map_to_cn end-to-end on a committed
    synthetic ACS seg and assert the three output files match committed snapshots.
    These snapshots encode *current behavior* (generated from the code at the time
    the suite was written) — they are the safety net that proves the nested-helper
    extraction is behavior-preserving. They are NOT a semantically-validated golden
    baseline — semantic recovery is checked instead by the black-box section below.
  * UNIT: the pure helpers extracted to module level.
"""

import pathlib

import numpy as np
import pandas as pd
import pytest

import map_to_absolute_copy_number as m2a
from helpers import run_script, assert_frame_close, run_map_to_absolute, synth_truth_segs

FIX = pathlib.Path(__file__).resolve().parent / "fixtures" / "map_to_absolute"
SNAP = pathlib.Path(__file__).resolve().parent / "snapshots" / "map_to_absolute"
ACS_SEG = FIX / "synthetic.acs.seg"

OUTPUT_FILES = {
    "allelic": [
        "TestS.segtab.allelic.completed.txt",
        "TestS.IGV.seg.allelic.completed.txt",
        "TestS.rescued_intervals.allelic.txt",
    ],
    "total": [
        "TestS.segtab.total.completed.txt",
        "TestS.IGV.seg.total.completed.txt",
        "TestS.rescued_intervals.total.txt",
    ],
}


@pytest.mark.regression
@pytest.mark.parametrize("mode", ["allelic", "total"])
def test_map_to_cn_characterization(tmp_path, mode):
    out = tmp_path / mode
    run_script("map_to_absolute_copy_number.py", [
        "--outdir", str(out), "--purity", "0.8", "--ploidy", "2.0",
        "--sample", "TestS", "--sex", "XXY",
        "--acs_cr_seg", str(ACS_SEG), "--copy_num_type", mode,
    ])
    for fname in OUTPUT_FILES[mode]:
        got = pd.read_csv(out / fname, sep="\t", low_memory=False)
        exp = pd.read_csv(SNAP / mode / fname, sep="\t", low_memory=False)
        assert_frame_close(got, exp, atol=1e-6, rtol=1e-6)


# --------------------------------------------------------------------------- #
# extracted pure helpers
# --------------------------------------------------------------------------- #
@pytest.mark.unit
@pytest.mark.parametrize("chrom,nX,nY,normal,expected", [
    ("X", 2, 1, 2, 2), ("chrX", 1, 1, 2, 1), ("Y", 2, 1, 2, 1),
    ("chrY", 2, 0, 2, 0), ("1", 2, 1, 2, 2), ("7", 1, 1, 3, 3),
])
def test_chromosomal_ploidy(chrom, nX, nY, normal, expected):
    assert m2a.chromosomal_ploidy(chrom, nX, nY, normal) == expected


@pytest.mark.unit
@pytest.mark.parametrize("values,expected", [
    ([2, 2, 4, 4], 2),         # tie -> nearest to 2
    ([3, 3, 3, 1], 3),         # clear mode
    ([0, 0, 5], 0),
])
def test_wmode_tiebreak(values, expected):
    w = np.ones(len(values))
    assert m2a.wmode(np.array(values, dtype=float), w) == expected


@pytest.mark.unit
def test_wmode_weights_matter():
    # value 4 has more weight than value 2 -> mode is 4
    assert m2a.wmode(np.array([2.0, 4.0]), np.array([0.1, 10.0])) == 4


@pytest.mark.stat
def test_map_cn_to_cluster_snaps_to_nearest():
    cv = np.array([1.0, 2.0, 2.5, 3.0])
    is_int = np.modf(cv)[0] == 0
    assert m2a.map_cn_to_cluster(2.02, 0.1, cv, is_int) == pytest.approx(2.0)


@pytest.mark.stat
def test_map_cn_to_cluster_no_valid_cluster_returns_input():
    cv = np.array([1.0, 2.0, 3.0])
    is_int = np.modf(cv)[0] == 0
    # cn far from every cluster with tiny sigma -> nothing within tail threshold
    assert m2a.map_cn_to_cluster(100.0, 1e-6, cv, is_int) == pytest.approx(100.0)


@pytest.mark.stat
def test_map_cn_to_cluster_integer_preference():
    # integer (2.0) and fractional (2.5) both plausible; integer preferred
    cv = np.array([2.0, 2.5])
    is_int = np.modf(cv)[0] == 0
    assert m2a.map_cn_to_cluster(2.25, 0.5, cv, is_int) == pytest.approx(2.0)


def _seg_row(**kw):
    base = dict(rescaled_total_cn=2.0, Chromosome="1", is_parental_haploid=False,
                SegLabelCNLOH=1, **{"mu.minor.abs": np.nan})
    base.update(kw)
    return pd.Series(base)


@pytest.mark.unit
def test_split_alleles_loh_when_cn_below_one():
    assert m2a.split_alleles_for_segment(_seg_row(rescaled_total_cn=0.4), {}) == (0, 0.4)


@pytest.mark.unit
def test_split_alleles_parental_haploid():
    row = _seg_row(rescaled_total_cn=2.0, is_parental_haploid=True)
    assert m2a.split_alleles_for_segment(row, {}) == (0, 2.0)


@pytest.mark.unit
def test_split_alleles_no_cnloh_label():
    row = _seg_row(rescaled_total_cn=3.0, SegLabelCNLOH=0)
    assert m2a.split_alleles_for_segment(row, {}) == (0, 3.0)


@pytest.mark.unit
def test_split_alleles_balanced_plateau():
    row = _seg_row(rescaled_total_cn=4.0, Chromosome="1")
    assert m2a.split_alleles_for_segment(row, {"1": 4}) == (2.0, 2.0)


@pytest.mark.unit
def test_split_alleles_uses_mu_minor():
    row = _seg_row(rescaled_total_cn=3.0, Chromosome="1", **{"mu.minor.abs": 1.0})
    # c0=3 (odd) -> not balanced; split by mu.minor.abs clipped to [0, CN/2]
    assert m2a.split_alleles_for_segment(row, {"1": 3}) == (1.0, 2.0)


@pytest.mark.unit
def test_split_alleles_fallback():
    row = _seg_row(rescaled_total_cn=3.0, Chromosome="1")  # mu.minor.abs NaN
    assert m2a.split_alleles_for_segment(row, {"1": 3}) == (1, 2)


# --------------------------------------------------------------------------- #
# extracted CCF functions: allele_ccf / segment_ccf (M5b)
# --------------------------------------------------------------------------- #
def _ccf_seg(purity=0.8, ploidy=2.0, cp=2, normal_ploidy=2):
    """A minimal valid seg + (delta, chr_ploidy, b) for the CCF functions.
    HSCR is placed at integer allelic-comb states (clonal); the last segment is
    high-CN to exercise the high-copy regime."""
    D = (1 - purity) * cp + purity * ploidy * cp / normal_ploidy
    delta, b = purity / D, (1 - purity) * cp / D
    rows = [(1, 1), (0, 2), (1, 3), (7, 8)]  # neutral, LOH, gain, high-CN (total 15)
    a1 = np.array([r[0] for r in rows], dtype=float)
    a2 = np.array([r[1] for r in rows], dtype=float)
    n = len(rows)
    seg = pd.DataFrame({
        "Chromosome": ["1"] * n,
        "hscr.a1": a1 * delta * cp + b, "hscr.a2": a2 * delta * cp + b,
        "sigma.a1": [0.05] * n, "sigma.a2": [0.05] * n,
        "modal.a1": a1.astype(int), "modal.a2": a2.astype(int),
        "rescaled.cn.a1": a1, "rescaled.cn.a2": a2,
        "CN": a1 + a2, "CN.sigma": [0.1] * n,
        "is_parental_haploid": [False] * n, "W": [1.0 / n] * n,
    })
    return seg, pd.Series([delta] * n), pd.Series([float(cp)] * n), pd.Series([b] * n)


@pytest.mark.stat
def test_allele_ccf_shapes_and_bounds():
    seg, delta, cp, b = _ccf_seg()
    hat, low, high = m2a.allele_ccf(seg, delta, cp, b, "a1")
    assert hat.shape == low.shape == high.shape == (len(seg),)
    for arr in (hat, low, high):
        finite = arr[np.isfinite(arr)]
        assert ((finite >= 0) & (finite <= 1)).all()
    ok = np.isfinite(hat) & np.isfinite(low) & np.isfinite(high)
    assert (low[ok] <= hat[ok] + 1e-9).all() and (hat[ok] <= high[ok] + 1e-9).all()


@pytest.mark.stat
def test_allele_ccf_deterministic():
    seg, delta, cp, b = _ccf_seg()
    a = m2a.allele_ccf(seg, delta, cp, b, "a2")
    c = m2a.allele_ccf(seg, delta, cp, b, "a2")
    for x, y in zip(a, c):
        np.testing.assert_array_equal(x, y)


@pytest.mark.stat
def test_segment_ccf_columns_and_bounds():
    seg, delta, cp, b = _ccf_seg()
    out = m2a.segment_ccf(seg, delta, cp, b, {"1": 2}, nX=2, nY=0, normal_ploidy=2)
    expected_cols = {
        "ccf_mean", "ccf_median", "ccf_mode", "ccf_ci95_low", "ccf_ci95_high",
        "ccf_ci95_width", "ccf_entropy", "ccf_max_p", "ccf_support_50",
        "ccf_support_90", "current_pair_entropy", "current_pair_effn",
        "tail_mass_high_cn", "high_copy_flag", "ambiguous_ccf_flag",
        "baseline_q1", "baseline_q2",
    }
    assert set(out.columns) == expected_cols
    assert list(out.index) == list(seg.index)
    med = out["ccf_median"].to_numpy()
    finite = med[np.isfinite(med)]
    assert ((finite >= 0) & (finite <= 1)).all()


@pytest.mark.stat
def test_segment_ccf_flags_high_copy_segment():
    seg, delta, cp, b = _ccf_seg()
    out = m2a.segment_ccf(seg, delta, cp, b, {"1": 2}, nX=2, nY=0, normal_ploidy=2)
    # the last segment (total CN 15) is in the high-copy regime
    assert out["high_copy_flag"].iloc[-1] == 1
    assert out["high_copy_flag"].iloc[0] == 0  # neutral segment is not


# --------------------------------------------------------------------------- #
# Black-box expectations (synthetic cohort). The has-absolute path is used
# (ABSOLUTE calls authoritative); CN/allelic/flags are checked against truth.
# --------------------------------------------------------------------------- #
def _seg_at(seg, contig, start):
    r = seg[(seg.Chromosome.astype(str) == contig) & (seg["Start.bp"] == start)]
    return r.iloc[0] if not r.empty else None


@pytest.mark.regression
def test_map_segment_ccf_bounds_and_aliasing(synth, tmp_path):
    # ccf_median is COMPUTED by segment_ccf (not the ABSOLUTE passthrough). A
    # near-integer subclonal CN aliases to clonal (documented identifiability limit);
    # we assert the valid invariant (bounded) rather than ill-posed recovery.
    seg = run_map_to_absolute(synth, tmp_path / "m", "PXX", "PXX-T2", use_absolute=True)
    med = seg["ccf_median"].to_numpy(dtype=float)
    fin = med[np.isfinite(med)]
    assert ((fin >= 0) & (fin <= 1)).all()
    assert "ccf_median" in seg.columns and "cancer.cell.frac.a1" in seg.columns


@pytest.mark.regression
def test_map_amp_flags(synth, tmp_path):
    seg = run_map_to_absolute(synth, tmp_path / "m", "PXX", "PXX-T1", use_absolute=True)
    amp = _seg_at(seg, "3", 100000)        # HIGH_AMP {8,1}: amp on the 8-copy allele only
    assert int(amp["amp.a1"]) + int(amp["amp.a2"]) == 1
    focal = _seg_at(seg, "1", 1000000)     # FOCAL {5,5}: both alleles amplified
    assert int(focal["amp.a1"]) == 1 and int(focal["amp.a2"]) == 1


@pytest.mark.regression
def test_map_male_sex_chromosome_recovery(synth, tmp_path):
    seg = run_map_to_absolute(synth, tmp_path / "m", "PXY", "PXY-T3", use_absolute=True)
    x_loh = _seg_at(seg, "X", 100000)      # X_LOH, truth total 0 -> HZ
    x_gain = _seg_at(seg, "X", 800000)     # X_GAIN {0,2} -> total 2, LOH
    y_loss = _seg_at(seg, "Y", 50000)      # Y_LOSS, truth total 0 -> HZ
    assert round(x_loh["rescaled_total_cn"]) == 0 and x_loh["HZ"] == 1
    assert round(x_gain["rescaled_total_cn"]) == 2 and x_gain["LOH"] == 1
    assert round(y_loss["rescaled_total_cn"]) == 0 and y_loss["HZ"] == 1


@pytest.mark.regression
def test_map_xxy_chrx_diploid_chry_haploid(synth, tmp_path):
    seg = run_map_to_absolute(synth, tmp_path / "m", "PXXY", "PXXY-T1", use_absolute=True)
    x_loh = _seg_at(seg, "X", 100000)      # LOH {2,0} -> total 2
    assert round(x_loh["rescaled_total_cn"]) == 2 and x_loh["LOH"] == 1
    y_neut = _seg_at(seg, "Y", 50000)      # NEUTRAL haploid {0,1} -> total 1
    assert round(y_neut["rescaled_total_cn"]) == 1


@pytest.mark.regression
def test_map_female_chrx_diploid_no_chry_and_flanks_neutral(synth, tmp_path):
    seg = run_map_to_absolute(synth, tmp_path / "m", "PXX", "PXX-T1", use_absolute=True)
    assert not (seg["Chromosome"].astype(str) == "Y").any()    # female: no chrY
    x_loh = _seg_at(seg, "X", 100000)                          # chrX LOH {2,0} -> total 2
    assert round(x_loh["rescaled_total_cn"]) == 2 and x_loh["LOH"] == 1
    flank = _seg_at(seg, "1", 1)                               # chr1:1-99999 filler
    if flank is not None:
        assert round(flank["rescaled_total_cn"]) == 2
        assert flank["LOH"] == 0 and flank["HZ"] == 0


@pytest.mark.regression
def test_map_allelic_sum_identity_and_nonnegativity(synth, tmp_path):
    seg = run_map_to_absolute(synth, tmp_path / "m", "PXX", "PXX-T1", use_absolute=True)
    assert (seg["rescaled.cn.a1"] >= 0).all() and (seg["rescaled.cn.a2"] >= 0).all()
    assert (seg["modal.a1"] + seg["modal.a2"] == seg["modal_total_cn"]).all()


@pytest.mark.regression
def test_map_every_absolute_segment_present(synth, tmp_path):
    seg = run_map_to_absolute(synth, tmp_path / "m", "PXX", "PXX-T1", use_absolute=True)
    abs_in = pd.read_csv(synth.files["patients"]["PXX"]["samples"]["PXX-T1"]["absolute_segtab"],
                         sep="\t", comment="#")
    out_iv = set(zip(seg["Chromosome"].astype(str), seg["Start.bp"], seg["End.bp"]))
    in_iv = set(zip(abs_in["Chromosome"].astype(str), abs_in["Start.bp"], abs_in["End.bp"]))
    assert in_iv <= out_iv


# --------------------------------------------------------------------------- #
# Truth-based CN recovery over the full acs -> map chain (ABSOLUTE authoritative,
# the WGD/somix/de-novo variants), migrated from the synthetic-pipeline suite.
# --------------------------------------------------------------------------- #
@pytest.mark.regression
def test_map_to_absolute_has_absolute_recovers_cn_and_flags(synth, tmp_path):
    # primary path: map_to_absolute post-processes ABSOLUTE -> exact CN/allelic + flags
    seg = run_map_to_absolute(synth, tmp_path / "m", "PXX", "PXX-T1", use_absolute=True)
    truth = synth_truth_segs(synth, "PXX", "PXX-T1")
    seen = 0
    for _, t in truth.iterrows():
        row = seg[(seg.Chromosome.astype(str) == str(t.contig)) & (seg["Start.bp"] == t.start)]
        if row.empty:
            continue
        seen += 1
        r = row.iloc[0]
        assert r["rescaled_total_cn"] == pytest.approx(t.total_cn, abs=1e-6)
        assert round(r["modal_total_cn"]) == t.total_cn
        assert {round(r["rescaled.cn.a1"]), round(r["rescaled.cn.a2"])} == {t.cn_a1, t.cn_a2}
        if t.total_cn == 0:
            assert r["HZ"] == 1
        if (t.cn_a1 == 0) ^ (t.cn_a2 == 0):  # one allele lost, total > 0
            assert r["LOH"] == 1
    assert seen >= 6


@pytest.mark.regression
@pytest.mark.parametrize("sample", ["PXY-T4", "PXY-T1"])  # WGD=1 (ploidy3), WGD=2 (ploidy4)
def test_map_to_absolute_wgd_recovery(synth, tmp_path, sample):
    seg = run_map_to_absolute(synth, tmp_path / sample, "PXY", sample, use_absolute=True)
    truth = synth_truth_segs(synth, "PXY", sample)
    for _, t in truth.iterrows():
        row = seg[(seg.Chromosome.astype(str) == str(t.contig)) & (seg["Start.bp"] == t.start)]
        if row.empty:
            continue
        assert round(seg.loc[row.index[0], "modal_total_cn"]) == t.total_cn


@pytest.mark.regression
def test_map_to_absolute_somix_path(synth, tmp_path):
    # somix CR-seg parsing branch, with ABSOLUTE calls for tight recovery
    seg = run_map_to_absolute(synth, tmp_path / "sx", "PXX", "PXX-T1",
                              somix=True, use_absolute=True)
    truth = synth_truth_segs(synth, "PXX", "PXX-T1")
    for _, t in truth.iterrows():
        row = seg[(seg.Chromosome.astype(str) == str(t.contig)) & (seg["Start.bp"] == t.start)]
        if row.empty:
            continue
        assert round(row.iloc[0]["modal_total_cn"]) == t.total_cn


@pytest.mark.regression
def test_map_to_absolute_denovo_runs_and_flags_homozygous_del(synth, tmp_path):
    # de-novo (rescue) path has no ABSOLUTE: CN is normalized to sum(W*CN)=ploidy,
    # so absolute CN recovery is approximate by design. Check it runs and that the
    # homozygous-deletion segment survives as the lowest-CN segment.
    seg = run_map_to_absolute(synth, tmp_path / "dn", "PXX", "PXX-T1")
    assert len(seg) > 0
    assert seg["rescaled_total_cn"].min() < 0.5
