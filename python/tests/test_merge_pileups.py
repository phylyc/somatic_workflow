"""Unit / IO / stat tests for merge_pileups.py."""

import gzip

import pandas as pd
import pytest

import merge_pileups as mp
from helpers import args_ns, pileup_df, write_pileup, PILEUP_COLUMNS, run_script, synth_runs
from tests.synthetic.scenarios import write_empty_pileup


# --------------------------------------------------------------------------- #
# get_contigs
# --------------------------------------------------------------------------- #
@pytest.mark.unit
def test_get_contigs_default_when_none():
    contigs = mp.get_contigs(None)
    assert contigs[:3] == ["1", "2", "3"]
    assert "X" in contigs and "Y" in contigs and "MT" in contigs
    assert "chr1" in contigs and "chrM" in contigs


@pytest.mark.unit
def test_get_contigs_from_dict(ref_dict):
    assert mp.get_contigs(ref_dict) == ["1", "2", "X"]


# --------------------------------------------------------------------------- #
# sort_genomic_positions
# --------------------------------------------------------------------------- #
@pytest.mark.unit
def test_sort_genomic_positions_orders_by_contig_then_position():
    idx = pd.MultiIndex.from_tuples(
        [("2", 5), ("1", 10), ("1", 2), ("2", 1)], names=["contig", "position"]
    )
    out = mp.sort_genomic_positions(idx, contig_order=["1", "2"])
    assert list(out) == [("1", 2), ("1", 10), ("2", 1), ("2", 5)]


@pytest.mark.unit
def test_sort_genomic_positions_appends_unknown_contigs():
    idx = pd.MultiIndex.from_tuples(
        [("99", 1), ("1", 1)], names=["contig", "position"]
    )
    out = mp.sort_genomic_positions(idx, contig_order=["1"])
    # known contig first, unknown appended after
    assert list(out) == [("1", 1), ("99", 1)]


# --------------------------------------------------------------------------- #
# get_header_and_df / write_header_and_df
# --------------------------------------------------------------------------- #
@pytest.mark.io
def test_get_header_and_df_reads_header_and_table(tmp_path):
    df = pileup_df([("1", 100, 10, 5, 0, 0.33)])
    path = write_pileup(tmp_path / "s.pileup", df, sample="Test-X")
    header, got = mp.get_header_and_df(path, columns=PILEUP_COLUMNS)
    assert "SAMPLE=Test-X" in header
    assert list(got.columns) == PILEUP_COLUMNS
    assert got.shape == (1, 6)


@pytest.mark.io
def test_get_header_and_df_missing_file_returns_empty(tmp_path):
    with pytest.warns(UserWarning):
        header, df = mp.get_header_and_df(str(tmp_path / "nope.pileup"),
                                          columns=PILEUP_COLUMNS)
    assert header is None
    assert df.empty
    assert list(df.columns) == PILEUP_COLUMNS


@pytest.mark.io
def test_write_header_and_df_roundtrip_gzip(tmp_path):
    df = pileup_df([("1", 1, 3, 4, 0, 0.5)])
    out = tmp_path / "s.pileup.gz"
    mp.write_header_and_df("#h\n", df, str(out))
    with gzip.open(out, "rt") as fh:
        text = fh.read()
    assert text.startswith("#h\n")
    assert "ref_count" in text


# --------------------------------------------------------------------------- #
# merge_pileups core
# --------------------------------------------------------------------------- #
def _read_out(out_dir, sample):
    return pd.read_csv(out_dir / f"{sample}.pileup", sep="\t", comment="#")


@pytest.mark.stat
def test_merge_sums_counts_and_takes_max_af(tmp_path, tmp_out):
    # Same sample name across two runs -> summed counts, max AF.
    p1 = write_pileup(tmp_path / "r1.pileup",
                      pileup_df([("1", 100, 10, 5, 1, 0.2)]), sample="S")
    p2 = write_pileup(tmp_path / "r2.pileup",
                      pileup_df([("1", 100, 4, 6, 2, 0.7)]), sample="S")
    args = args_ns(sample=["S", "S"], pileup=[p1, p2], ref_dict=None,
                   output_dir=str(tmp_out), min_read_depth=0,
                   compress_output=False, verbose=False)
    mp.merge_pileups(args)
    out = _read_out(tmp_out, "S")
    row = out.iloc[0]
    assert row["ref_count"] == 14 and row["alt_count"] == 11
    assert row["other_alt_count"] == 3
    assert row["allele_frequency"] == pytest.approx(0.7)


@pytest.mark.stat
def test_merge_dedupes_within_run(tmp_path, tmp_out):
    # duplicate (contig, position) within one run is dropped before merging
    p = write_pileup(tmp_path / "r.pileup",
                     pileup_df([("1", 100, 10, 5, 0, 0.2),
                                ("1", 100, 99, 99, 9, 0.9)]), sample="S")
    args = args_ns(sample=["S"], pileup=[p], ref_dict=None,
                   output_dir=str(tmp_out), min_read_depth=0,
                   compress_output=False, verbose=False)
    mp.merge_pileups(args)
    out = _read_out(tmp_out, "S")
    assert len(out) == 1
    assert out.iloc[0]["ref_count"] == 10  # first kept


@pytest.mark.stat
def test_merge_min_read_depth_filter(tmp_path, tmp_out):
    p = write_pileup(tmp_path / "r.pileup",
                     pileup_df([("1", 1, 1, 1, 0, 0.5),     # depth 2
                                ("1", 2, 20, 20, 0, 0.5)]),  # depth 40
                     sample="S")
    args = args_ns(sample=["S"], pileup=[p], ref_dict=None,
                   output_dir=str(tmp_out), min_read_depth=10,
                   compress_output=False, verbose=False)
    mp.merge_pileups(args)
    out = _read_out(tmp_out, "S")
    assert list(out["position"]) == [2]


@pytest.mark.stat
def test_merge_empty_input_no_crash(tmp_path, tmp_out):
    p = write_pileup(tmp_path / "empty.pileup",
                     pileup_df([]), sample="S")
    args = args_ns(sample=["S"], pileup=[p], ref_dict=None,
                   output_dir=str(tmp_out), min_read_depth=0,
                   compress_output=False, verbose=False)
    mp.merge_pileups(args)
    out = _read_out(tmp_out, "S")
    assert out.empty


@pytest.mark.stat
def test_merge_two_distinct_samples_written_separately(tmp_path, tmp_out):
    pn = write_pileup(tmp_path / "n.pileup",
                      pileup_df([("1", 1, 5, 5, 0, 0.5)]), sample="N")
    pt = write_pileup(tmp_path / "t.pileup",
                      pileup_df([("1", 1, 8, 2, 0, 0.2)]), sample="T")
    args = args_ns(sample=["N", "T"], pileup=[pn, pt], ref_dict=None,
                   output_dir=str(tmp_out), min_read_depth=0,
                   compress_output=False, verbose=False)
    mp.merge_pileups(args)
    assert (tmp_out / "N.pileup").exists()
    assert (tmp_out / "T.pileup").exists()
    assert _read_out(tmp_out, "T").iloc[0]["ref_count"] == 8


@pytest.mark.stat
def test_merge_output_sorted_by_contig_order(tmp_path, tmp_out, ref_dict):
    p = write_pileup(tmp_path / "r.pileup",
                     pileup_df([("X", 5, 1, 1, 0, 0.5),
                                ("1", 9, 1, 1, 0, 0.5),
                                ("2", 3, 1, 1, 0, 0.5)]), sample="S")
    args = args_ns(sample=["S"], pileup=[p], ref_dict=ref_dict,
                   output_dir=str(tmp_out), min_read_depth=0,
                   compress_output=False, verbose=False)
    mp.merge_pileups(args)
    out = _read_out(tmp_out, "S")
    assert list(out["contig"].astype(str)) == ["1", "2", "X"]


# --------------------------------------------------------------------------- #
# Black-box expectations (synthetic cohort): drive the real script on the
# generated cohort and check the output against the known ground truth.
# --------------------------------------------------------------------------- #
@pytest.mark.regression
def test_merge_carries_population_af_not_recomputed(synth, tmp_path):
    # allele_frequency is the population AF carried as max across runs, NOT the
    # observed VAF recomputed from merged counts.
    runs = [r["pileup"] for r in synth_runs(synth, "PXY", "PXY-T1") if r["pileup"]]
    out = tmp_path / "o"
    out.mkdir()
    argv = []
    for p in runs:
        argv += ["--sample", "PXY-T1", "--pileup", p]
    argv += ["--output_dir", str(out)]
    run_script("merge_pileups.py", argv)
    merged = pd.read_csv(out / "PXY-T1.pileup", sep="\t", comment="#")
    ins = [pd.read_csv(p, sep="\t", comment="#") for p in runs]
    exp_af = pd.concat(ins).groupby(["contig", "position"])["allele_frequency"].max()
    for r in merged.itertuples():
        assert r.allele_frequency == pytest.approx(exp_af[(str(r.contig), int(r.position))])
    # 1:200000 is hom-alt: observed VAF ~1.0 but population AF is < 1
    row = merged[(merged.contig.astype(str) == "1") & (merged.position == 200000)].iloc[0]
    observed = row.alt_count / (row.ref_count + row.alt_count)
    assert abs(row.allele_frequency - observed) > 0.05


@pytest.mark.regression
def test_merge_pileups_sums_runs(synth, tmp_path):
    # PXY-T1 has two WES runs -> two pileups; merged counts == exact sum.
    runs = synth_runs(synth, "PXY", "PXY-T1")
    pileups = [r["pileup"] for r in runs if r["pileup"]]
    assert len(pileups) == 2

    out = tmp_path / "merge"
    out.mkdir()
    # merge_pileups zips --sample with --pileup, so repeat the sample name per run
    argv = []
    for p in pileups:
        argv += ["--sample", "PXY-T1", "--pileup", p]
    argv += ["--output_dir", str(out)]
    run_script("merge_pileups.py", argv)

    got = pd.read_csv(out / "PXY-T1.pileup", sep="\t", comment="#")
    ins = [pd.read_csv(p, sep="\t", comment="#") for p in pileups]
    exp = pd.concat(ins).groupby(["contig", "position"]).agg(
        ref_count=("ref_count", "sum"), alt_count=("alt_count", "sum"),
        other_alt_count=("other_alt_count", "sum"),
        allele_frequency=("allele_frequency", "max")).reset_index()
    m = got.merge(exp, on=["contig", "position"], suffixes=("", "_exp"))
    assert len(m) == len(exp)
    for c in ["ref_count", "alt_count", "other_alt_count"]:
        assert (m[c] == m[f"{c}_exp"]).all(), c


@pytest.mark.regression
def test_empty_pileup_merges_to_empty(synth, tmp_path):
    # header-only pileup -> empty merged output, exit 0 (legacy empty-input case)
    p = write_empty_pileup(tmp_path / "e.pileup", "EMPTY")
    out = tmp_path / "o"
    out.mkdir()
    run_script("merge_pileups.py", ["--sample", "EMPTY", "--pileup", p,
                                    "--output_dir", str(out)])
    df = pd.read_csv(out / "EMPTY.pileup", sep="\t", comment="#")
    assert df.empty


@pytest.mark.regression
def test_chr_prefixed_merge_roundtrips_contigs(synth_chr, tmp_path):
    # merge two WES runs of PXY-T1 (chr-named) and confirm chr contigs round-trip
    pileups = [r["pileup"] for r in synth_runs(synth_chr, "PXY", "PXY-T1") if r["pileup"]]
    out = tmp_path / "o"
    out.mkdir()
    argv = []
    for p in pileups:
        argv += ["--sample", "PXY-T1", "--pileup", p]
    argv += ["--output_dir", str(out)]
    run_script("merge_pileups.py", argv)
    merged = pd.read_csv(out / "PXY-T1.pileup", sep="\t", comment="#")
    assert merged["contig"].astype(str).str.startswith("chr").all()


@pytest.mark.regression
def test_merge_min_read_depth_filters_on_merged_depth(synth, tmp_path):
    # --min_read_depth is applied to the MERGED (summed-across-runs) depth, so a
    # locus that is shallow in each run individually can still pass once combined.
    runs = [r["pileup"] for r in synth_runs(synth, "PXY", "PXY-T1") if r["pileup"]]
    ins = [pd.read_csv(p, sep="\t", comment="#") for p in runs]
    merged = pd.concat(ins).groupby(["contig", "position"]).agg(
        ref_count=("ref_count", "sum"), alt_count=("alt_count", "sum"),
        other_alt_count=("other_alt_count", "sum")).reset_index()
    # merge_pileups thresholds on the total merged depth (ref + alt + other_alt)
    depth = {(str(r.contig), int(r.position)): r.ref_count + r.alt_count + r.other_alt_count
             for r in merged.itertuples()}
    lo, hi = min(depth.values()), max(depth.values())
    assert lo < hi, "need a spread of merged depths to exercise the boundary"
    thr = lo + 1  # drops at least the shallowest locus, keeps the deepest

    out = tmp_path / "o"
    out.mkdir()
    argv = []
    for p in runs:
        argv += ["--sample", "PXY-T1", "--pileup", p]
    argv += ["--output_dir", str(out), "--min_read_depth", str(thr)]
    run_script("merge_pileups.py", argv)
    got = pd.read_csv(out / "PXY-T1.pileup", sep="\t", comment="#")
    kept = set(zip(got["contig"].astype(str), got["position"].astype(int)))
    assert kept == {k for k, d in depth.items() if d >= thr}
    assert len(kept) < len(depth)  # filtering actually removed something
