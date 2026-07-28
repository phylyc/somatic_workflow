"""Unit / stat tests for calculate_cancer_cell_fraction.py."""

import numpy as np
import pandas as pd
import pytest
import scipy.stats as st

import calculate_cancer_cell_fraction as ccf
from helpers import args_ns, run_ccf


# --------------------------------------------------------------------------- #
# pure utilities
# --------------------------------------------------------------------------- #
@pytest.mark.unit
@pytest.mark.parametrize("sex,expected", [
    ("FEMALE", "XX"), ("female", "XX"), ("Male", "XY"), ("MALE", "XY"),
    ("unknown", "XXY"), ("nan", "XXY"), ("none", "XXY"),
    ("XX", "XX"), ("XY", "XY"), ("XXY", "XXY"), ("xyy", "XYY"),
])
def test_standardize_sex(sex, expected):
    assert ccf.standardize_sex(sex) == expected


@pytest.mark.unit
@pytest.mark.parametrize("chrom,expected", [
    ("chr1", "1"), ("1", "1"), ("CHR2", "2"), ("chrX", "X"), ("X", "X"),
    ("chrM", "MT"), ("M", "MT"), ("MT", "MT"), (5, "5"),
])
def test_normalize_chromosome(chrom, expected):
    assert ccf.normalize_chromosome(chrom) == expected


@pytest.mark.unit
@pytest.mark.parametrize("chrom,sex,ploidy", [
    ("X", "XX", 2), ("X", "XY", 1), ("Y", "XY", 1), ("Y", "XX", 0),
    ("1", "XX", 2), ("chr1", "XY", 2), ("X", "XXY", 2), ("Y", "XXY", 1),
])
def test_get_chromosomal_ploidy(chrom, sex, ploidy):
    assert ccf.get_chromosomal_ploidy(chrom, sex=sex, normal_ploidy=2) == ploidy


@pytest.mark.unit
def test_sort_genomic_df_canonical_order_stable():
    df = pd.DataFrame({
        "Chromosome": ["X", "2", "1", "1", "MT"],
        "Start_position": [1, 1, 50, 10, 1],
        "End_position": [2, 2, 60, 20, 2],
    })
    out = ccf.sort_genomic_df(df)
    assert list(out["Chromosome"]) == ["1", "1", "2", "X", "MT"]
    # within chr 1, sorted by start (stable)
    assert list(out.loc[out["Chromosome"] == "1", "Start_position"]) == [10, 50]


@pytest.mark.unit
def test_sort_genomic_df_empty_missing_columns_noop():
    out = ccf.sort_genomic_df(pd.DataFrame())
    assert out.empty
    assert list(out.columns) == []


@pytest.mark.unit
def test_normalize_maf_columns():
    df = pd.DataFrame(columns=["a(b)", "c", "d(e)f"])
    out = ccf.normalize_maf_columns(df)
    assert list(out.columns) == ["a.b.", "c", "d.e.f"]


@pytest.mark.unit
def test_to_numeric_coerces_and_skips_missing():
    df = pd.DataFrame({"a": ["1", "2", "x"], "b": ["0.5", "NA", "0.1"]})
    out = ccf.to_numeric(df, ["a", "b", "missing"])
    assert out["a"].tolist()[:2] == [1.0, 2.0]
    assert np.isnan(out["a"].iloc[2])
    assert np.isnan(out["b"].iloc[1])


@pytest.mark.unit
def test_ensure_columns_adds_missing():
    df = pd.DataFrame({"a": [1, 2]})
    out = ccf.ensure_columns(df, ["a", "b", "c"], fill_value=0)
    assert (out["b"] == 0).all() and (out["c"] == 0).all()
    assert out["a"].tolist() == [1, 2]


@pytest.mark.unit
def test_choose_output_columns_drops_internal_and_appends_extras():
    # note: pandas .empty is True when there are no rows, so include a row.
    caller = pd.DataFrame({"Chromosome": ["1"], "Start_position": [1],
                           "_chrom": ["1"], "_variant_key": ["1:1:1:A:T"]})
    out = ccf.choose_output_columns(caller, pd.DataFrame())
    assert "_chrom" not in out and "_variant_key" not in out
    assert "Chromosome" in out
    assert "ccf_hat" in out  # an ABSOLUTE extra column appended


@pytest.mark.unit
def test_make_variant_key():
    df = pd.DataFrame({
        "_chrom": ["1", "2"], "Start_position": [100, 200],
        "End_position": [100, 200], "Reference_Allele": ["A", "G"],
        "Tumor_Seq_Allele2": ["T", "A"],
    })
    assert ccf.make_variant_key(df).tolist() == ["1:100:100:A:T", "2:200:200:G:A"]


# --------------------------------------------------------------------------- #
# IO / preparation
# --------------------------------------------------------------------------- #
@pytest.mark.io
def test_read_table_missing_returns_empty(tmp_path):
    out = ccf.read_table(str(tmp_path / "nope.tsv"))
    assert out.empty


@pytest.mark.io
def test_read_optional_maf_none():
    assert ccf.read_optional_maf(None).empty


# --------------------------------------------------------------------------- #
# quote / comment handling
#
# MAF is a plain tab-delimited format: `"` is a literal data character and `#`
# only marks a leading header block. Reading with pandas' defaults treats `"` as
# a quote character, so a line carrying an odd number of quotes opens a quoted
# region that swallows every following tab AND newline until the next quote.
# That merges whole records into one field: the surviving row loses its field
# boundaries (values shift left, tail padded with NaN) and the absorbed rows
# disappear. Observed in production as ABS_MAF rows whose Mutect2 QC block was
# silently all-NaN and whose CCF grid landed 72 columns early.
# --------------------------------------------------------------------------- #
_QUOTED_MAF = (
    "Hugo_Symbol\tDrugBank\tref_context\tHGNC_Alias_names\tt_alt_count\tMPOS\n"
    "AAA\tInsulin(DB00071)\tACGTACGT\tplain alias\t11\t20\n"
    'FEV\tInsulin(DB00071)\tGGGGCCCC\tFEV (fifth Ewing variant)", "FEV, ETS"\t22\t21\n'
    "BBB\tGentamicin(DB00798)\tTTTTAAAA\tanother alias\t33\t22\n"
    "CCC\tUrokinase(DB00013)\tCCCCGGGG\tthird alias\t44\t23\n"
)


@pytest.mark.io
def test_read_table_odd_quotes_do_not_merge_records(tmp_path):
    path = tmp_path / "quoted.maf"
    path.write_text(_QUOTED_MAF)

    out = ccf.read_table(str(path))

    assert len(out) == 4, "no record may be absorbed by an unbalanced quote"
    assert out["t_alt_count"].tolist() == [11, 22, 33, 44]
    assert out["MPOS"].tolist() == [20, 21, 22, 23], "QC block must survive intact"
    for col in out.columns:
        if out[col].dtype == object:
            assert not out[col].astype(str).str.contains("\t", regex=False).any(), (
                f"column {col} contains a merged field"
            )


@pytest.mark.io
def test_read_table_does_not_truncate_at_inline_hash(tmp_path):
    path = tmp_path / "hash.maf"
    path.write_text(
        "#version 2.4\n"
        "Hugo_Symbol\tDescription\tt_alt_count\n"
        "AAA\tno hash here\t11\n"
        "BBB\tprotein #2 subunit\t22\n"
    )

    out = ccf.read_table(str(path))

    assert len(out) == 2, "leading comment block must be skipped, data kept"
    assert out["Description"].tolist() == ["no hash here", "protein #2 subunit"]
    assert out["t_alt_count"].tolist() == [11, 22], "inline '#' must not amputate the row"


@pytest.mark.unit
@pytest.mark.parametrize("raw,expected", [
    # complete QUOTE_MINIMAL encodings -- these get decoded
    ('"quoted"', "quoted"),
    ('"""X"', '"X'),
    ('""', ""),
    ('"a ""b"" c"', 'a "b" c'),
    # not quoted fields -- MAF treats " as literal data, so these pass through
    ("unquoted", "unquoted"),
    ('has ""inner"" quotes', 'has ""inner"" quotes'),
    ('"', '"'),
    ("", ""),
    # Funcotator HGNC alias/previous-name lists: begin and end with a quote but are
    # NOT a quoted field. Stripping the outer pair would leave 'plexin 2", "plexin-A2'.
    ('"plexin 2", "plexin-A2"', '"plexin 2", "plexin-A2"'),
    ('"MDS1/EVI1-like", "transcription factor MEL1"', '"MDS1/EVI1-like", "transcription factor MEL1"'),
    # truncated alias (the &beta; entity cut at the ';') -- unbalanced, left alone
    ('"transforming growth factor-&', '"transforming growth factor-&'),
])
def test_unquote_fields(raw, expected):
    out = ccf.unquote_fields(pd.DataFrame({"a": [raw]}))
    assert out["a"].iloc[0] == expected


@pytest.mark.unit
def test_unquote_fields_leaves_non_object_columns():
    df = pd.DataFrame({"n": [1, 2], "f": [1.5, 2.5], "s": ['"x"', "y"]})
    out = ccf.unquote_fields(df)
    assert out["n"].tolist() == [1, 2]
    assert out["f"].tolist() == [1.5, 2.5]
    assert out["s"].tolist() == ["x", "y"]


@pytest.mark.io
def test_write_table_is_rectangular_and_unquoted(tmp_path):
    path = tmp_path / "out.tsv"
    ccf.write_table(
        pd.DataFrame({"a": ["x", "y"], "b": ["has\ttab", "has\nnewline"], "c": [1, 2]}),
        str(path),
    )

    text = path.read_text()
    assert '"' not in text, "output must never be quoted"
    widths = {len(line.split("\t")) for line in text.rstrip("\n").split("\n")}
    assert widths == {3}, f"every line must carry the same field count, got {widths}"


@pytest.mark.io
def test_read_table_write_table_round_trip_preserves_quoted_values(tmp_path):
    src = tmp_path / "src.maf"
    src.write_text(_QUOTED_MAF)
    first = ccf.read_table(str(src))

    dst = tmp_path / "dst.maf"
    ccf.write_table(first, str(dst))
    second = ccf.read_table(str(dst))

    pd.testing.assert_frame_equal(first, second)


@pytest.mark.io
def test_calculate_ccf_all_empty_mafs_writes_header(tmp_path):
    segtab = tmp_path / "S.segtab.allelic.completed.txt"
    pd.DataFrame({
        "sample": ["S"],
        "Chromosome": ["1"],
        "Start.bp": [1],
        "End.bp": [100],
        "modal_total_cn": [2],
        "rescaled_total_cn": [2.0],
        "rescaled.cn.a1": [1.0],
        "rescaled.cn.a2": [1.0],
        "modal.a1": [1],
        "modal.a2": [1],
    }).to_csv(segtab, sep="\t", index=False)
    empty_abs = tmp_path / "empty_abs.maf"
    empty_snv = tmp_path / "empty_snv.maf"
    empty_indel = tmp_path / "empty_indel.maf"
    for path in [empty_abs, empty_snv, empty_indel]:
        path.write_text("")

    ccf.calculate_ccf(args_ns(
        outdir=str(tmp_path), purity=1.0, ploidy=2.0, ssnv_skew=1.0,
        sample="S", sex="XY", absolute_maf=str(empty_abs), absolute_segtab=str(segtab),
        snv_maf=str(empty_snv), indel_maf=str(empty_indel), normal_ploidy=2,
        copy_num_type="allelic",
    ))

    out = pd.read_csv(tmp_path / "S.ABS_MAF.allelic.completed.txt", sep="\t", low_memory=False)
    assert out.empty
    assert "Chromosome" in out.columns
    assert "Start_position" in out.columns
    assert "ccf_hat" in out.columns


@pytest.mark.unit
def test_infer_sample_name_priority():
    seg = pd.DataFrame({"sample": ["SEG"]})
    caller = pd.DataFrame({"Tumor_Sample_Barcode": ["CALLER"]})
    abs_maf = pd.DataFrame({"Tumor_Sample_Barcode": ["ABS"]})
    assert ccf.infer_sample_name(args_ns(sample="CLI"), seg, caller, abs_maf) == "CLI"
    assert ccf.infer_sample_name(args_ns(sample=None), seg, caller, abs_maf) == "CALLER"
    assert ccf.infer_sample_name(args_ns(sample=None), seg, pd.DataFrame(), abs_maf) == "ABS"
    assert ccf.infer_sample_name(args_ns(sample=None), seg, pd.DataFrame(), pd.DataFrame()) == "SEG"
    assert ccf.infer_sample_name(args_ns(sample=None), pd.DataFrame({"sample": []}),
                                 pd.DataFrame(), pd.DataFrame()) == "sample"


@pytest.mark.unit
def test_concat_callers_dedupes_on_variant_key():
    snv = pd.DataFrame({"_variant_key": ["1:1:1:A:T", "1:2:2:C:G"], "x": [1, 2]})
    indel = pd.DataFrame({"_variant_key": ["1:1:1:A:T"], "x": [99]})
    out = ccf.concat_callers(snv, indel)
    assert len(out) == 2
    assert out.loc[out["_variant_key"] == "1:1:1:A:T", "x"].iloc[0] == 1  # keep first


@pytest.mark.unit
def test_prepare_segtab_indexes_and_sorts():
    seg = pd.DataFrame({
        "Chromosome": ["2", "1"], "Start.bp": ["1", "1"], "End.bp": ["100", "100"],
        "modal_total_cn": ["3", "2"], "sample": ["S", "S"],
    })
    out = ccf.prepare_segtab(seg, sex="XXY", normal_ploidy=2)
    assert list(out["_chrom"]) == ["1", "2"]            # sorted
    assert list(out["_seg_idx1"]) == [1, 2]             # 1-indexed after sort
    assert out["_chr_ploidy"].tolist() == [2, 2]


@pytest.mark.unit
def test_prepare_segtab_renames_tcr_columns():
    seg = pd.DataFrame({"Chromosome": ["1"], "Start.bp": [1], "End.bp": [100],
                        "modal_cn": [3], "sample": ["S"]})
    out = ccf.prepare_segtab(seg, sex="XXY", normal_ploidy=2)
    assert "modal_total_cn" in out.columns
    assert out["modal_total_cn"].iloc[0] == 3


@pytest.mark.io
def test_map_variants_to_segments_overlap_and_ploidy():
    variants = pd.DataFrame({
        "_chrom": ["1", "1", "2"],
        "Start_position": [50, 5000, 10],
        "End_position": [50, 5000, 10],
    })
    seg = ccf.prepare_segtab(pd.DataFrame({
        "Chromosome": ["1", "2"], "Start.bp": [1, 1], "End.bp": [100, 100],
        "modal_total_cn": [3, 4], "rescaled.cn.a1": [1, 2],
        "rescaled.cn.a2": [2, 2], "sample": ["S", "S"],
    }), sex="XY", normal_ploidy=2)
    out = ccf.map_variants_to_segments(variants, seg)
    # variant at chr1:50 maps into chr1 segment; chr1:5000 is past the segment end
    assert out["q_hat"].iloc[0] == 3
    assert np.isnan(out["q_hat"].iloc[1])
    assert out["local_cn_a1"].iloc[0] == 1 and out["local_cn_a2"].iloc[0] == 2
    # chr2 ploidy under XY is normal (autosome) = 2
    assert out["normal_allele_count"].iloc[2] == 2


@pytest.mark.io
def test_map_variants_to_segments_total_mode_no_allelic_columns():
    variants = pd.DataFrame({"_chrom": ["1"], "Start_position": [50],
                             "End_position": [50]})
    seg = ccf.prepare_segtab(pd.DataFrame({
        "Chromosome": ["1"], "Start.bp": [1], "End.bp": [100],
        "modal_total_cn": [2], "sample": ["S"]}), sex="XXY", normal_ploidy=2)
    out = ccf.map_variants_to_segments(variants, seg)
    assert out["q_hat"].iloc[0] == 2
    assert np.isnan(out["local_cn_a1"].iloc[0])  # gracefully NaN


# --------------------------------------------------------------------------- #
# statistical cores
# --------------------------------------------------------------------------- #
@pytest.mark.stat
@pytest.mark.parametrize("k,a,b,n", [(0, 2, 3, 10), (5, 1, 1, 10), (10, 4, 2, 10)])
def test_d_beta_binom_matches_scipy(k, a, b, n):
    got = ccf.d_beta_binom(k, a, b, n)
    exp = st.betabinom.pmf(k, n, a, b)
    assert got == pytest.approx(exp, rel=1e-9)
    # log path
    assert ccf.d_beta_binom(k, a, b, n, log=True) == pytest.approx(np.log(exp), rel=1e-9)


@pytest.mark.stat
def test_d_beta_binom_uniform_when_a_b_one():
    n = 7
    vals = np.array([ccf.d_beta_binom(k, 1, 1, n) for k in range(n + 1)])
    np.testing.assert_allclose(vals, np.full(n + 1, 1 / (n + 1)), rtol=1e-9)
    assert vals.sum() == pytest.approx(1.0)


@pytest.mark.stat
def test_log_density_f_zero_and_one_are_binomial():
    eps = 1e-3
    d = ccf.log_density_alt_cond_coverage(alt=0, cov=10, f=[0.0, 0.5, 1.0],
                                          ssnv_skew=1.0, rho=100, epsilon=eps)
    assert d[0] == pytest.approx(st.binom.logpmf(0, 10, eps / 3))
    assert d[2] == pytest.approx(st.binom.logpmf(0, 10, 1 - eps / 3))
    assert np.isfinite(d[1])


@pytest.mark.stat
def test_log_density_nan_to_neg_inf():
    d = ccf.log_density_alt_cond_coverage(alt=5, cov=np.nan, f=[0.5],
                                          ssnv_skew=1.0, rho=100, epsilon=1e-3)
    assert d[0] == -np.inf


@pytest.mark.stat
def test_calc_ccf_posterior_ll_grid_length():
    out = ccf.calc_ccf_posterior_ll_grid(alt=5, ref=5, alpha=0.5, q=2,
                                         normal_allele_count=2, ssnv_skew=1.0,
                                         rho=100, epsilon=1e-3)
    assert out.shape == (101,)


def _mut(alt, ref, q=2.0, nac=2.0, total_q=2.0):
    return pd.DataFrame({"t_alt_count": [alt], "t_ref_count": [ref],
                         "q_hat": [q], "normal_allele_count": [nac],
                         "total_q_hat": [total_q]})


@pytest.mark.stat
def test_calc_ccf_dens_rows_normalize():
    dens = ccf.calc_ccf_dens(_mut(20, 0), alpha=1.0, ssnv_skew=1.0, rho=100,
                             epsilon=1e-3)
    assert dens.shape == (1, 101)
    assert np.nansum(dens[0]) == pytest.approx(1.0, rel=1e-6)


@pytest.mark.stat
def test_calc_ccf_dens_qhat_zero_floors_to_one():
    # An observed variant on a segment whose modal total CN rounds to 0 must still get a
    # CCF posterior: the copy number is floored at 1 rather than dropped with an empty grid.
    kw = dict(alpha=1.0, ssnv_skew=1.0, rho=100, epsilon=1e-3)
    dens0 = ccf.calc_ccf_dens(_mut(10, 10, q=0.0), **kw)
    dens1 = ccf.calc_ccf_dens(_mut(10, 10, q=1.0), **kw)
    assert not np.all(np.isnan(dens0[0]))
    assert np.nansum(dens0[0]) == pytest.approx(1.0, rel=1e-6)
    np.testing.assert_allclose(dens0[0], dens1[0], equal_nan=True)


@pytest.mark.stat
def test_calc_ccf_dens_qhat_nan_is_nan():
    # No copy number at all (e.g. total-CR segtab lacking total CN) leaves an empty grid.
    dens = ccf.calc_ccf_dens(_mut(10, 10, q=np.nan), alpha=1.0, ssnv_skew=1.0,
                             rho=100, epsilon=1e-3)
    assert np.all(np.isnan(dens[0]))


@pytest.mark.stat
def test_calc_ccf_95ci_all_nan_row_stays_nan():
    # A row with no posterior must report ccf_hat/low/high all as NaN, not [0, 1].
    dens = np.full((1, ccf.CCF_GRID.size), np.nan)
    ci = ccf.calc_ccf_95ci(dens)
    assert np.isnan(ci["ccf_hat"].iloc[0])
    assert np.isnan(ci["ccf_CI95_low"].iloc[0])
    assert np.isnan(ci["ccf_CI95_high"].iloc[0])


@pytest.mark.stat
def test_calc_ccf_dens_clonal_peaks_high():
    # high VAF, full purity, diploid -> CCF posterior peaks near 1
    dens = ccf.calc_ccf_dens(_mut(40, 0), alpha=1.0, ssnv_skew=1.0, rho=100,
                             epsilon=1e-3)
    peak = ccf.CCF_GRID[np.nanargmax(dens[0])]
    assert peak >= 0.9


@pytest.mark.stat
def test_calc_ccf_95ci_orders_low_hat_high():
    dens = ccf.calc_ccf_dens(_mut(10, 10), alpha=1.0, ssnv_skew=1.0, rho=100,
                             epsilon=1e-3)
    ci = ccf.calc_ccf_95ci(dens)
    lo, hat, hi = ci["ccf_CI95_low"].iloc[0], ci["ccf_hat"].iloc[0], ci["ccf_CI95_high"].iloc[0]
    assert lo <= hat <= hi


@pytest.mark.stat
def test_calc_ccf_ci_width_shrinks_with_coverage():
    # same VAF (50%), deeper coverage -> tighter 95% CI on the CCF posterior
    kw = dict(alpha=1.0, ssnv_skew=1.0, rho=100, epsilon=1e-3)
    shallow = ccf.calc_ccf_95ci(ccf.calc_ccf_dens(_mut(10, 10), **kw))
    deep = ccf.calc_ccf_95ci(ccf.calc_ccf_dens(_mut(150, 150), **kw))
    w_shallow = shallow["ccf_CI95_high"].iloc[0] - shallow["ccf_CI95_low"].iloc[0]
    w_deep = deep["ccf_CI95_high"].iloc[0] - deep["ccf_CI95_low"].iloc[0]
    assert w_deep < w_shallow


@pytest.mark.stat
def test_power_calc_monotone_in_coverage():
    kw = dict(eps=1e-3 / 3, fdr=0.5e-7, delta=0.1, ssnv_skew=1.0, rho=100)
    low = ccf.power_calc(n=10, **kw)
    high = ccf.power_calc(n=500, **kw)
    assert 0.0 <= low <= high <= 1.0
    assert high > low


@pytest.mark.stat
def test_power_calc_monotone_in_delta():
    kw = dict(n=100, eps=1e-3 / 3, fdr=0.5e-7, ssnv_skew=1.0, rho=100)
    assert ccf.power_calc(delta=0.01, **kw) <= ccf.power_calc(delta=0.5, **kw)


@pytest.mark.stat
def test_power_calc_nan_guard():
    assert np.isnan(ccf.power_calc(n=np.nan, eps=1e-3, fdr=0.5e-7, delta=0.1,
                                   ssnv_skew=1.0, rho=100))


@pytest.mark.stat
def test_power_calc_single_read_monotone():
    a = ccf.power_calc_for_single_read(n=5, delta=0.1, ssnv_skew=1.0, rho=100)
    b = ccf.power_calc_for_single_read(n=200, delta=0.1, ssnv_skew=1.0, rho=100)
    assert 0.0 <= a <= b <= 1.0


@pytest.mark.stat
def test_mode_ssnv_power_calc_falls_back_to_total_q():
    mut = _mut(10, 10, q=np.nan, total_q=2.0)
    power = ccf.mode_ssnv_power_calc(mut, alpha=0.5, ssnv_skew=1.0, rho=100)
    assert power.shape == (1,)
    assert 0.0 <= power[0] <= 1.0


# --------------------------------------------------------------------------- #
# Black-box expectations (synthetic cohort, contract-derived; acs -> map -> ccf).
# --------------------------------------------------------------------------- #
def _ccf_at(df, contig, pos):
    r = df[(df.Chromosome.astype(str) == contig) & (df.Start_position == pos)]
    return float(r["ccf_hat"].iloc[0]) if not r.empty else None


@pytest.mark.regression
def test_ccf_per_sample_isolation(synth, tmp_path):
    # 2:800000 is subclonal (ccf 0.6) in PXX-T1 but absent (ccf 0) in PXX-T2
    t1 = run_ccf(synth, tmp_path / "t1", "PXX", "PXX-T1")
    t2 = run_ccf(synth, tmp_path / "t2", "PXX", "PXX-T2")
    c1 = _ccf_at(t1, "2", 800000)
    c2 = _ccf_at(t2, "2", 800000)
    assert c1 is not None and c1 >= 0.4
    assert c2 is None or c2 <= 0.2


@pytest.mark.regression
def test_ccf_bounds_multiplicity_and_posterior_sum(synth, tmp_path):
    out = run_ccf(synth, tmp_path / "t1", "PXX", "PXX-T1")
    hat = out["ccf_hat"].to_numpy(dtype=float)
    lo = out["ccf_CI95_low"].to_numpy(dtype=float)
    hi = out["ccf_CI95_high"].to_numpy(dtype=float)
    ok = np.isfinite(hat)
    assert ((hat[ok] >= 0) & (hat[ok] <= 1)).all()
    b = np.isfinite(lo) & np.isfinite(hi) & ok
    assert (lo[b] <= hat[b] + 1e-9).all() and (hat[b] <= hi[b] + 1e-9).all()
    # multiplicity-2 clonal variant -> ccf ~1, never > 1
    r = out[(out.Chromosome.astype(str) == "3") & (out.Start_position == 300000)]
    if not r.empty:
        assert float(r["ccf_hat"].iloc[0]) >= 0.8
    # per-variant CCF posterior over the grid sums to ~1 (finite rows)
    grid = [c for c in ccf.CCF_COLUMNS if c in out.columns]
    if grid:
        sums = out[grid].sum(axis=1).to_numpy(dtype=float)
        fin = sums[np.isfinite(sums) & (sums > 0)]
        np.testing.assert_allclose(fin, 1.0, atol=1e-6)


@pytest.mark.regression
def test_ccf_detection_power_and_calibration(synth, tmp_path):
    out = run_ccf(synth, tmp_path / "t1", "PXX", "PXX-T1")
    dp = out["detection_power"].to_numpy(dtype=float)
    fin = dp[np.isfinite(dp)]
    assert ((fin >= 0) & (fin <= 1)).all()
    r = out[(out.Chromosome.astype(str) == "1") & (out.Start_position == 250000)]  # clonal
    if not r.empty:
        assert float(r["detection_power"].iloc[0]) >= 0.9
    # 95% CI brackets truth ccf for the large majority of retained variants
    truth = synth.truth["somatic_variants"].query(
        "patient_id=='PXX' and sample_id=='PXX-T1'").copy()
    truth["note"] = truth["note"].fillna("")
    tmap = {(str(t.contig), int(t.position)): float(t.ccf)
            for _, t in truth.iterrows() if t.note == ""}
    covered = total = 0
    for r in out.itertuples():
        key = (str(r.Chromosome), int(r.Start_position))
        if key in tmap and np.isfinite(r.ccf_CI95_low) and np.isfinite(r.ccf_CI95_high):
            total += 1
            covered += int(r.ccf_CI95_low - 1e-9 <= tmap[key] <= r.ccf_CI95_high + 1e-9)
    assert total >= 3 and covered / total >= 0.8


@pytest.mark.regression
def test_ccf_cross_sample_clonal_consistency(synth, tmp_path):
    # a clonal variant (PXY 1:300000) recovers ccf ~1 across samples of differing
    # purity/ploidy/WGD (T1 = WGD2/ploidy4, T4 = WGD1/ploidy3)
    for sample in ("PXY-T1", "PXY-T4"):
        out = run_ccf(synth, tmp_path / sample, "PXY", sample)
        c = _ccf_at(out, "1", 300000)
        assert c is not None and c >= 0.8, (sample, c)


@pytest.mark.regression
def test_calculate_ccf_recovers_ccf(synth, tmp_path):
    # ccf_hat recovered vs truth ccf over the full acs -> map -> ccf chain;
    # MT / unmapped / zero-coverage variants are dropped.
    out = run_ccf(synth, tmp_path, "PXX", "PXX-T1")
    out_pos = {(str(c), int(p)) for c, p in zip(out["Chromosome"], out["Start_position"])}

    truth = synth.truth["somatic_variants"].query(
        "patient_id=='PXX' and sample_id=='PXX-T1'").copy()
    truth["note"] = truth["note"].fillna("")  # empty note reads back as NaN
    by_pos = {(str(t.contig), int(t.position)): t for _, t in truth.iterrows()}
    ccf_hat = {(str(c), int(p)): h for c, p, h in
               zip(out["Chromosome"], out["Start_position"], out["ccf_hat"])}

    for note in ["MT", "unmapped", "zero_cov"]:
        for (c, p), t in by_pos.items():
            if t.note == note:
                assert (c, p) not in out_pos, (note, c, p)

    checked = 0
    for (c, p), t in by_pos.items():
        if t.note or (c, p) not in ccf_hat:
            continue
        checked += 1
        if t.ccf >= 0.99:
            assert ccf_hat[(c, p)] >= 0.8, (c, p, ccf_hat[(c, p)])
        elif t.ccf == 0.0:
            assert ccf_hat[(c, p)] <= 0.3, (c, p, ccf_hat[(c, p)])
        else:  # subclonal (e.g. 0.6)
            assert abs(ccf_hat[(c, p)] - t.ccf) <= 0.25, (c, p, ccf_hat[(c, p)], t.ccf)
    assert checked >= 3
