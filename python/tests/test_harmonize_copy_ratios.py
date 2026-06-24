"""Unit / stat tests for harmonize_copy_ratios.py."""

import numpy as np
import pandas as pd
import pytest

import harmonize_copy_ratios as hcr
from helpers import run_script, synth_runs
from tests.synthetic.assays import target_bins
from tests.synthetic.scenarios import write_cr_for_target_set, write_empty_cr


def make_h(**kw):
    defaults = dict(
        copy_ratio_files=[],
        column_names=["CONTIG", "START", "END", "LOG2_COPY_RATIO"],
        column_types=[str, int, int, float],
        agg_col=["LOG2_COPY_RATIO"], agg_func=["mean"],
        min_target_length=2, verbose=False,
    )
    defaults.update(kw)
    h = hcr.Harmonizer(**defaults)
    h.contigs = ["1", "2", "X"]
    return h


def cr_rows(sample, rows):
    """rows: (contig, start, end, log2cr[, call])."""
    ncol = len(rows[0])
    cols = ["CONTIG", "START", "END", "LOG2_COPY_RATIO", "CALL"][:ncol]
    df = pd.DataFrame(rows, columns=cols)
    df["SAMPLE"] = sample
    return df


# --------------------------------------------------------------------------- #
# interval algebra
# --------------------------------------------------------------------------- #
@pytest.mark.unit
def test_non_overlapping_breakpoints():
    h = make_h()
    df = pd.DataFrame({"CONTIG": ["1", "1"], "START": [1, 5], "END": [10, 15],
                       "LOG2_COPY_RATIO": [0.1, 0.2]})
    out = h.non_overlapping(df)
    breaks = sorted(set(out["START"]).union(set(out["END"])))
    assert breaks == [1, 5, 10, 15]


@pytest.mark.unit
def test_merge_abutting_same_value_merges():
    h = make_h()
    idx = pd.MultiIndex.from_tuples([("1", 1, 10), ("1", 10, 20)],
                                    names=["CONTIG", "START", "END"])
    df = pd.DataFrame({"LOG2_COPY_RATIO": [0.30, 0.30]}, index=idx)
    out = h.merge_abutting_intervals(df)
    assert len(out) == 1
    assert out.index[0] == ("1", 1, 20)


@pytest.mark.unit
def test_merge_abutting_different_value_kept():
    h = make_h()
    idx = pd.MultiIndex.from_tuples([("1", 1, 10), ("1", 10, 20)],
                                    names=["CONTIG", "START", "END"])
    df = pd.DataFrame({"LOG2_COPY_RATIO": [0.30, 0.90]}, index=idx)
    out = h.merge_abutting_intervals(df)
    assert len(out) == 2


@pytest.mark.unit
def test_select_intervals_min_length():
    h = make_h()
    idx = pd.MultiIndex.from_tuples([("1", 1, 2), ("1", 10, 100)],
                                    names=["CONTIG", "START", "END"])
    df = pd.DataFrame({"x": [1, 2]}, index=idx)
    out = h.select_intervals(df, min_length=10)
    assert len(out) == 1 and out.index[0] == ("1", 10, 100)


@pytest.mark.unit
def test_sort_genomic_positions():
    h = make_h()
    idx = pd.MultiIndex.from_tuples(
        [("2", 1, 5), ("1", 10, 20), ("1", 1, 9)],
        names=["CONTIG", "START", "END"])
    out = h.sort_genomic_positions(idx)
    assert list(out) == [("1", 1, 9), ("1", 10, 20), ("2", 1, 5)]


@pytest.mark.unit
@pytest.mark.parametrize("func,values,expected", [
    ("first", [3, 1, 2], 3),
    ("last", [3, 1, 2], 2),
    ("strongest_signal", [0.1, -0.9, 0.3], -0.9),
])
def test_agg_func_map(func, values, expected):
    f = hcr.Harmonizer.agg_func_map(func)
    assert f(pd.Series(values)) == expected


@pytest.mark.unit
def test_agg_func_map_error_mean():
    f = hcr.Harmonizer.agg_func_map("error_mean")
    assert f(pd.Series([3.0, 4.0])) == pytest.approx(np.sqrt((9 + 16) / 2))


# --------------------------------------------------------------------------- #
# _harmonize integration
# --------------------------------------------------------------------------- #
@pytest.mark.stat
def test_harmonize_intersection_common_grid():
    h = make_h()
    h.copy_ratios = [
        cr_rows("A", [("1", 1, 10, 0.1)]),
        cr_rows("B", [("1", 1, 10, 0.5)]),
    ]
    out = h._harmonize()
    assert out.shape[0] == 1
    assert out[("A", "LOG2_COPY_RATIO")].iloc[0] == pytest.approx(0.1)
    assert out[("B", "LOG2_COPY_RATIO")].iloc[0] == pytest.approx(0.5)


@pytest.mark.stat
def test_harmonize_empty_intersection():
    # samples on disjoint contigs -> no shared interval under INTERSECTION
    h = make_h()
    h.copy_ratios = [
        cr_rows("A", [("1", 1, 100, 0.1)]),
        cr_rows("B", [("2", 1, 100, 0.5)]),
    ]
    out = h._harmonize()
    assert out.empty


@pytest.mark.stat
def test_harmonize_intersection_drops_contig_in_subset():
    # a contig present in only some samples is dropped under INTERSECTION
    h = make_h()
    h.copy_ratios = [
        cr_rows("A", [("1", 1, 100, 0.1), ("2", 1, 100, 0.2)]),
        cr_rows("B", [("1", 1, 100, 0.5)]),   # no contig 2
    ]
    out = h._harmonize()
    assert set(out.index.get_level_values("CONTIG")) == {"1"}


@pytest.mark.stat
def test_harmonize_union_keeps_all():
    h = make_h(interval_set_rule="UNION")
    h.copy_ratios = [
        cr_rows("A", [("1", 1, 100, 0.1)]),
        cr_rows("B", [("2", 1, 100, 0.5)]),
    ]
    out = h._harmonize()
    assert out.shape[0] == 2


@pytest.mark.stat
def test_harmonize_filters_germline_cnv_calls():
    h = make_h(
        column_names=["CONTIG", "START", "END", "LOG2_COPY_RATIO", "CALL"],
        column_types=[str, int, int, float, str],
        agg_col=["LOG2_COPY_RATIO", "CALL"], agg_func=["mean", "first"],
        normal_sample="N", filter_germline_calls=True,
    )
    h.copy_ratios = [
        cr_rows("N", [("1", 1, 100, 0.0, "+"), ("1", 200, 300, 0.0, "0")]),
        cr_rows("T", [("1", 1, 100, 0.9, "0"), ("1", 200, 300, 0.1, "0")]),
    ]
    out = h._harmonize()
    # the [1,100] interval is "+" in the normal -> dropped; [200,300] survives
    assert out.shape[0] == 1
    starts = out.index.get_level_values("START")
    assert list(starts) == [200]


# --------------------------------------------------------------------------- #
# Black-box expectations (synthetic cohort, contract-derived).
# --------------------------------------------------------------------------- #
_HARM = ["--column_names", "CONTIG", "START", "END", "LOG2_COPY_RATIO",
         "--column_types", "str", "int", "int", "float",
         "--agg_col", "LOG2_COPY_RATIO", "--agg_func", "mean", "--min_target_length", "1"]


def _pxx_t1_cr(synth):
    return next(r["denoised_cr"] for r in synth_runs(synth, "PXX", "PXX-T1") if r["denoised_cr"])


@pytest.mark.regression
def test_harmonize_preserves_value_on_matching_grid(synth, tmp_path):
    cr = _pxx_t1_cr(synth)
    out = tmp_path / "o"
    out.mkdir()
    run_script("harmonize_copy_ratios.py",
               ["-I", cr, "-S", "PXX-T1", "-O", str(out), "--suffix", ".h.tsv",
                "--interval_set_rule", "INTERSECTION", *_HARM])
    src = pd.read_csv(cr, sep="\t", comment="@")
    harm = pd.read_csv(out / "PXX-T1.h.tsv", sep="\t", comment="@")
    src_val = {(str(c), int(s)): v for c, s, v in
               zip(src.CONTIG, src.START, src.LOG2_COPY_RATIO)}
    assert len(harm) == len(src)
    for r in harm.itertuples():
        assert r.LOG2_COPY_RATIO == pytest.approx(src_val[(str(r.CONTIG), int(r.START))], abs=1e-6)


@pytest.mark.regression
def test_harmonize_min_target_length_boundary(synth, tmp_path):
    # every synthetic CR interval is exactly 50000 bp -> min_target_length 50001 drops all
    cr = _pxx_t1_cr(synth)
    out = tmp_path / "o"
    out.mkdir()
    base = ["-I", cr, "-S", "PXX-T1", "-O", str(out), "--suffix", ".h.tsv",
            "--interval_set_rule", "INTERSECTION", "--column_names", "CONTIG", "START",
            "END", "LOG2_COPY_RATIO", "--column_types", "str", "int", "int", "float",
            "--agg_col", "LOG2_COPY_RATIO", "--agg_func", "mean"]
    run_script("harmonize_copy_ratios.py", base + ["--min_target_length", "50001"])
    assert pd.read_csv(out / "PXX-T1.h.tsv", sep="\t", comment="@").empty


@pytest.mark.regression
def test_harmonize_output_sorted_disjoint_unique(synth, tmp_path):
    # output intervals are sorted & non-overlapping per contig, with no duplicates
    cr = _pxx_t1_cr(synth)
    ulp = next(r["denoised_cr"] for r in synth_runs(synth, "PXX", "PXX-T2")
               if r["assay"] == "ULP_WGS" and r["denoised_cr"])
    out = tmp_path / "o"
    out.mkdir()
    run_script("harmonize_copy_ratios.py",
               ["-I", cr, "-I", ulp, "-S", "PXX-T1", "-S", "PXX-T2", "-O", str(out),
                "--suffix", ".h.tsv", "--interval_set_rule", "UNION", *_HARM])
    harm = pd.read_csv(out / "PXX-T1.h.tsv", sep="\t", comment="@")
    assert len(harm) == len(harm.drop_duplicates(["CONTIG", "START", "END"]))
    for _c, g in harm.groupby("CONTIG"):
        starts = g["START"].to_numpy()
        ends = g["END"].to_numpy()
        assert (starts[1:] >= ends[:-1]).all()   # sorted & non-overlapping


# --------------------------------------------------------------------------- #
# Interval intersection / union on the real synthetic target grids, plus the
# degenerate empty/single-sample cases (migrated from the scenario suite).
# --------------------------------------------------------------------------- #
def _pxx_t1(synth):
    p = synth.cohort.patient("PXX")
    return p, next(s for s in p.samples if s.name == "PXX-T1")


def _bin_starts(target_set, patient):
    return {(c, s) for c, s, _e, _i in target_bins(target_set)
            if patient.chromosomal_ploidy(c) > 0}


@pytest.mark.regression
def test_harmonize_intersection_grid(synth, tmp_path):
    # PXX-T1 = WES_A (even bins); PXX-T2 ULP run = genome-wide bins.
    t1_cr = next(r["denoised_cr"] for r in synth_runs(synth, "PXX", "PXX-T1") if r["denoised_cr"])
    t2_cr = next(r["denoised_cr"] for r in synth_runs(synth, "PXX", "PXX-T2")
                 if r["assay"] == "ULP_WGS" and r["denoised_cr"])
    out = tmp_path / "harm"
    out.mkdir()
    run_script("harmonize_copy_ratios.py", [
        "-I", t1_cr, "-I", t2_cr, "-S", "PXX-T1", "-S", "PXX-T2",
        "-O", str(out), "--suffix", ".harmonized.CR.tsv",
        "--interval_set_rule", "INTERSECTION", *_HARM])
    harm = pd.read_csv(out / "PXX-T1.harmonized.CR.tsv", sep="\t", comment="@")
    # intersection of WES_A (even bins) and genome-wide == WES_A bins, minus
    # contigs with ploidy 0 (female chrY has no copy ratio -> not emitted).
    patient = synth.cohort.patient("PXX")
    exp_starts = _bin_starts("WES_A", patient)
    got_starts = set(zip(harm["CONTIG"].astype(str), harm["START"].astype(int)))
    assert got_starts == exp_starts
    assert len(harm) == len(exp_starts)


@pytest.mark.regression
def test_empty_cr_harmonizes_without_crash(synth, tmp_path):
    e1 = write_empty_cr(tmp_path / "a.cr")
    e2 = write_empty_cr(tmp_path / "b.cr")
    out = tmp_path / "o"
    out.mkdir()
    run_script("harmonize_copy_ratios.py",
               ["-I", e1, "-I", e2, "-S", "A", "-S", "B", "-O", str(out),
                "--suffix", ".h.tsv", "--interval_set_rule", "INTERSECTION", *_HARM])
    assert (out / "A.h.tsv").exists()


@pytest.mark.regression
def test_single_sample_harmonize_keeps_own_intervals(synth, tmp_path):
    patient, sample = _pxx_t1(synth)
    rng = np.random.default_rng(1)
    cr = write_cr_for_target_set(tmp_path / "wesA.cr", patient, sample, "WES_A", rng)
    out = tmp_path / "o"
    out.mkdir()
    run_script("harmonize_copy_ratios.py",
               ["-I", cr, "-S", "S1", "-O", str(out), "--suffix", ".h.tsv",
                "--interval_set_rule", "INTERSECTION", *_HARM])
    harm = pd.read_csv(out / "S1.h.tsv", sep="\t", comment="@")
    assert len(harm) == len(_bin_starts("WES_A", patient))


@pytest.mark.regression
def test_empty_intersection_disjoint_wes(synth, tmp_path):
    patient, sample = _pxx_t1(synth)
    rng = np.random.default_rng(2)
    a = write_cr_for_target_set(tmp_path / "wesA.cr", patient, sample, "WES_A", rng)
    c = write_cr_for_target_set(tmp_path / "wesC.cr", patient, sample, "WES_C", rng)
    out = tmp_path / "o"
    out.mkdir()
    run_script("harmonize_copy_ratios.py",
               ["-I", a, "-I", c, "-S", "A", "-S", "C", "-O", str(out),
                "--suffix", ".h.tsv", "--interval_set_rule", "INTERSECTION", *_HARM])
    harm = pd.read_csv(out / "A.h.tsv", sep="\t", comment="@")
    assert harm.empty  # WES_A (even bins) and WES_C (odd bins) are disjoint


@pytest.mark.regression
def test_partial_intersection_wes(synth, tmp_path):
    patient, sample = _pxx_t1(synth)
    rng = np.random.default_rng(3)
    a = write_cr_for_target_set(tmp_path / "wesA.cr", patient, sample, "WES_A", rng)
    b = write_cr_for_target_set(tmp_path / "wesB.cr", patient, sample, "WES_B", rng)
    out = tmp_path / "o"
    out.mkdir()
    run_script("harmonize_copy_ratios.py",
               ["-I", a, "-I", b, "-S", "A", "-S", "B", "-O", str(out),
                "--suffix", ".h.tsv", "--interval_set_rule", "INTERSECTION", *_HARM])
    harm = pd.read_csv(out / "A.h.tsv", sep="\t", comment="@")
    expected = _bin_starts("WES_A", patient) & _bin_starts("WES_B", patient)
    got = set(zip(harm["CONTIG"].astype(str), harm["START"].astype(int)))
    assert got == expected
    assert 0 < len(expected) < len(_bin_starts("WES_A", patient))
