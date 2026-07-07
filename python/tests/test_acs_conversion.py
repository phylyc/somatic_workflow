"""Unit / stat tests for acs_conversion.py (ModelSegments -> AllelicCapSeg).

The conversion is a single monolithic function reading/writing files, so we
drive it via ``convert_model_segments_to_alleliccapseg(args)`` with synthetic
ModelSegments .seg inputs and assert the computed AllelicCapSeg columns.
"""

import numpy as np
import pandas as pd
import pytest

import acs_conversion as acs
from helpers import args_ns, write_tsv, run_script
from tests.synthetic.simulate import expected_maf
from tests.synthetic.model import SegmentCN

MS_COLS = [
    "CONTIG", "START", "END",
    "NUM_POINTS_COPY_RATIO", "NUM_POINTS_ALLELE_FRACTION",
    "LOG2_COPY_RATIO_POSTERIOR_10", "LOG2_COPY_RATIO_POSTERIOR_50",
    "LOG2_COPY_RATIO_POSTERIOR_90",
    "MINOR_ALLELE_FRACTION_POSTERIOR_10", "MINOR_ALLELE_FRACTION_POSTERIOR_50",
    "MINOR_ALLELE_FRACTION_POSTERIOR_90",
]


def seg_row(contig="1", start=1, end=1001, n_cr=10, n_af=10,
            cr10=-0.05, cr50=0.0, cr90=0.05,
            maf10=0.4, maf50=0.45, maf90=0.48):
    return dict(zip(MS_COLS, [contig, start, end, n_cr, n_af,
                              cr10, cr50, cr90, maf10, maf50, maf90]))


def write_modelseg(path, rows):
    df = pd.DataFrame(rows, columns=MS_COLS)
    return write_tsv(path, df, header_lines=["@HD\tVN:1.6", "@SQ\tSN:1\tLN:100000"])


def run_acs(tmp_path, rows, **kw):
    seg = write_modelseg(tmp_path / "test.seg", rows)
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    args = args_ns(
        seg=seg, sample_id=kw.pop("sample_id", None),
        af_parameters=kw.pop("af_parameters", None),
        output_dir=str(out_dir), min_hets=kw.pop("min_hets", 0),
        min_probes=kw.pop("min_probes", 0),
        maf90_threshold=kw.pop("maf90_threshold", 0.49),
        sex=kw.pop("sex", "XXY"), normal_ploidy=kw.pop("normal_ploidy", 2),
        verbose=False,
    )
    acs.convert_model_segments_to_alleliccapseg(args)
    out = pd.read_csv(out_dir / "test.acs.seg", sep="\t")
    return out, out_dir


# --------------------------------------------------------------------------- #
# sex -> ploidy
# --------------------------------------------------------------------------- #
@pytest.mark.unit
@pytest.mark.parametrize("sex,contig,ploidy", [
    ("XX", "1", 2), ("XY", "1", 2),
    ("XX", "X", 2), ("XY", "X", 1),
    ("XXY", "X", 2), ("XXY", "Y", 1),
    ("Female", "X", 2), ("Male", "X", 1),
])
def test_tau_uses_sex_ploidy(tmp_path, sex, contig, ploidy):
    # tau = chr_ploidy * 2**LOG2_CR_50; with CR50=0 -> tau == chr_ploidy
    out, _ = run_acs(tmp_path, [seg_row(contig=contig, cr50=0.0, maf90=0.40)],
                     sex=sex)
    if ploidy == 0:
        assert np.isnan(out["tau"].iloc[0])
    else:
        assert out["tau"].iloc[0] == pytest.approx(ploidy)


@pytest.mark.stat
def test_tau_scales_with_log2_cr(tmp_path):
    out, _ = run_acs(tmp_path, [seg_row(cr50=1.0, maf90=0.40)], sex="XXY")
    # autosome ploidy 2, 2**1 = 2 -> tau = 4
    assert out["tau"].iloc[0] == pytest.approx(4.0)


# --------------------------------------------------------------------------- #
# f handling
# --------------------------------------------------------------------------- #
@pytest.mark.stat
def test_f_copied_from_maf50(tmp_path):
    out, _ = run_acs(tmp_path, [seg_row(maf50=0.31, maf90=0.40)])
    assert out["f"].iloc[0] == pytest.approx(0.31)


@pytest.mark.stat
def test_f_set_to_half_above_maf90_threshold(tmp_path):
    out, _ = run_acs(tmp_path, [seg_row(maf50=0.31, maf90=0.495)],
                     maf90_threshold=0.49)
    assert out["f"].iloc[0] == pytest.approx(0.5)


@pytest.mark.stat
def test_f_nan_on_x_for_male(tmp_path):
    out, _ = run_acs(tmp_path, [seg_row(contig="X", maf50=0.4, maf90=0.4)],
                     sex="XY")
    assert np.isnan(out["f"].iloc[0])


@pytest.mark.stat
def test_f_always_nan_on_y(tmp_path):
    out, _ = run_acs(tmp_path, [seg_row(contig="Y", maf50=0.4, maf90=0.4)],
                     sex="XY")
    assert np.isnan(out["f"].iloc[0])


@pytest.mark.stat
def test_ploidy_zero_masks_tau(tmp_path):
    # Y in an XX sample: chr_ploidy 0 -> tau & sigma.tau NaN (dropped downstream)
    out, _ = run_acs(tmp_path, [seg_row(contig="Y", maf90=0.4)], sex="XX")
    assert np.isnan(out["tau"].iloc[0])
    assert np.isnan(out["sigma.tau"].iloc[0])


# --------------------------------------------------------------------------- #
# error propagation
# --------------------------------------------------------------------------- #
@pytest.mark.stat
def test_sigma_tau_error_propagation(tmp_path):
    cr10, cr50, cr90 = -0.1, 0.0, 0.1
    out, _ = run_acs(tmp_path, [seg_row(cr10=cr10, cr50=cr50, cr90=cr90,
                                        maf90=0.40)], sex="XXY")
    ploidy = 2
    sigma_log2_tau = (cr90 - cr10) / 2.563
    var_log_tau = np.log(2) ** 2 * sigma_log2_tau ** 2
    expected = ploidy * np.exp(np.log(2) * cr50 + var_log_tau / 2) * \
        np.sqrt(np.exp(var_log_tau) - 1)
    assert out["sigma.tau"].iloc[0] == pytest.approx(expected, rel=1e-9)


@pytest.mark.stat
def test_mu_minor_major_and_sigmas(tmp_path):
    out, _ = run_acs(tmp_path, [seg_row(cr50=0.0, maf50=0.3, maf90=0.40,
                                        maf10=0.25)], sex="XXY")
    row = out.iloc[0]
    tau, f = row["tau"], row["f"]
    assert row["mu.minor"] == pytest.approx(f * tau)
    assert row["mu.major"] == pytest.approx((1 - f) * tau)
    sigma_f = (0.40 - 0.25) / 2.563
    sigma_tau = row["sigma.tau"]
    assert row["sigma.minor"] == pytest.approx(
        np.sqrt(tau ** 2 * sigma_f ** 2 + f ** 2 * sigma_tau ** 2), rel=1e-9)
    assert row["sigma.major"] == pytest.approx(
        np.sqrt(tau ** 2 * sigma_f ** 2 + (1 - f) ** 2 * sigma_tau ** 2),
        rel=1e-9)


# --------------------------------------------------------------------------- #
# CNLOH labelling (label_cnloh, exercised through output)
# --------------------------------------------------------------------------- #
@pytest.mark.stat
def test_cnloh_default_label_two(tmp_path):
    out, _ = run_acs(tmp_path, [seg_row(contig="1", maf50=0.45, maf90=0.40)])
    assert out["SegLabelCNLOH"].iloc[0] == 2


@pytest.mark.stat
def test_cnloh_label_one_for_adjacent_pair(tmp_path):
    # identical CR -> tau equal (pass_t); f differs by >0.01 (not pass_f).
    # smaller-f segment gets label 1.
    rows = [
        seg_row(contig="1", start=1, end=1001, maf50=0.45, maf90=0.46),
        seg_row(contig="1", start=1001, end=2001, maf50=0.30, maf90=0.40),
    ]
    out, _ = run_acs(tmp_path, rows)
    assert list(out["SegLabelCNLOH"]) == [2, 1]


@pytest.mark.stat
def test_cnloh_label_zero_for_triplet(tmp_path):
    rows = [
        seg_row(contig="1", start=1, end=1001, maf50=0.40, maf90=0.45),
        seg_row(contig="1", start=1001, end=2001, maf50=0.20, maf90=0.30),
        seg_row(contig="1", start=2001, end=3001, maf50=0.45, maf90=0.47),
    ]
    out, _ = run_acs(tmp_path, rows)
    assert list(out["SegLabelCNLOH"]) == [2, 0, 2]


@pytest.mark.stat
def test_cnloh_no_labeling_across_chromosomes(tmp_path):
    rows = [
        seg_row(contig="1", maf50=0.45, maf90=0.46),
        seg_row(contig="2", maf50=0.30, maf90=0.40),
    ]
    out, _ = run_acs(tmp_path, rows)
    assert list(out["SegLabelCNLOH"]) == [2, 2]


# --------------------------------------------------------------------------- #
# skew output
# --------------------------------------------------------------------------- #
@pytest.mark.io
def test_skew_from_mean_bias(tmp_path):
    af = pd.DataFrame({"PARAMETER_NAME": ["MEAN_BIAS"], "POSTERIOR_50": [1.0]})
    af_path = write_tsv(tmp_path / "af.param", af, header_lines=["@comment"])
    out, out_dir = run_acs(tmp_path, [seg_row(maf90=0.40)], af_parameters=af_path)
    skew = np.loadtxt(out_dir / "test.acs.seg.skew")
    assert float(skew) == pytest.approx(2.0 / (1.0 + 1.0))  # == 1.0


@pytest.mark.io
def test_skew_sentinel_without_af_parameters(tmp_path):
    out, out_dir = run_acs(tmp_path, [seg_row(maf90=0.40)], af_parameters=None)
    skew = np.loadtxt(out_dir / "test.acs.seg.skew")
    assert float(skew) == pytest.approx(-1.0)


# --------------------------------------------------------------------------- #
# filtering
# --------------------------------------------------------------------------- #
@pytest.mark.stat
def test_min_hets_and_min_probes_drop_rows(tmp_path):
    rows = [
        seg_row(contig="1", start=1, end=1001, n_cr=100, n_af=100, maf90=0.40),
        seg_row(contig="1", start=1001, end=2001, n_cr=1, n_af=1, maf90=0.40),
    ]
    out, _ = run_acs(tmp_path, rows, min_hets=5, min_probes=5)
    assert len(out) == 1
    assert out["Start.bp"].iloc[0] == 1


# --------------------------------------------------------------------------- #
# somix segmentation input
# --------------------------------------------------------------------------- #
SOMIX_COLS = ["contig", "Start.bp", "End.bp", "n_markers", "n_snps", "f_MAP",
              "log_tCR", "sem_log_tCR", "f_hessian", "sample_id"]


def somix_row(contig="1", start=1, end=1001, n_markers=10, n_snps=10,
              f_map=0.42, log_tcr=0.1, sem_log_tcr=0.05, f_hessian=-1000.0,
              sample_id="S1"):
    return dict(zip(SOMIX_COLS, [contig, start, end, n_markers, n_snps, f_map,
                                 log_tcr, sem_log_tcr, f_hessian, sample_id]))


def run_acs_somix(tmp_path, rows, **kw):
    seg = write_tsv(tmp_path / "test.somix.seg", pd.DataFrame(rows, columns=SOMIX_COLS))
    out_dir = tmp_path / "out"
    out_dir.mkdir(exist_ok=True)
    sample_id = kw.pop("sample_id", None)
    args = args_ns(
        seg=seg, sample_id=sample_id,
        af_parameters=kw.pop("af_parameters", None),
        output_dir=str(out_dir), min_hets=kw.pop("min_hets", 0),
        min_probes=kw.pop("min_probes", 0),
        maf90_threshold=kw.pop("maf90_threshold", 0.49),
        sex=kw.pop("sex", "XXY"), normal_ploidy=kw.pop("normal_ploidy", 2),
        verbose=False,
    )
    acs.convert_model_segments_to_alleliccapseg(args)
    # --sample_id, when given, names the output after the sample; otherwise the input basename
    prefix = sample_id if sample_id is not None else "test.somix"
    out = pd.read_csv(out_dir / f"{prefix}.acs.seg", sep="\t")
    return out, out_dir


@pytest.mark.stat
def test_somix_tau_and_sigma_tau(tmp_path):
    # tau = chr_ploidy * exp(log_tCR); log-normal error propagation for sigma.tau
    log_tcr, sem = 0.2, 0.05
    out, _ = run_acs_somix(tmp_path, [somix_row(log_tcr=log_tcr, sem_log_tcr=sem)])
    ploidy = 2
    assert out["tau"].iloc[0] == pytest.approx(ploidy * np.exp(log_tcr))
    expected_sigma = ploidy * np.exp(log_tcr + sem ** 2 / 2) * np.sqrt(np.exp(sem ** 2) - 1)
    assert out["sigma.tau"].iloc[0] == pytest.approx(expected_sigma, rel=1e-9)


@pytest.mark.stat
def test_somix_f_passthrough_when_clearly_imbalanced(tmp_path):
    # tight posterior (large |f_hessian|) well below 0.5 -> f kept, not snapped
    out, _ = run_acs_somix(tmp_path, [somix_row(f_map=0.30, f_hessian=-1e6)],
                           maf90_threshold=0.49)
    assert out["f"].iloc[0] == pytest.approx(0.30)


@pytest.mark.stat
def test_somix_f_snapped_to_half_when_ci90_reaches_half(tmp_path):
    # near-balanced MAP whose CI90 upper bound reaches 0.5 -> snapped to 0.5
    out, _ = run_acs_somix(tmp_path, [somix_row(f_map=0.47, f_hessian=-1000.0)],
                           maf90_threshold=0.49)
    assert out["f"].iloc[0] == pytest.approx(0.5)


@pytest.mark.stat
def test_somix_f_snap_gated_by_variance(tmp_path):
    # same f_MAP: a wide posterior (small |f_hessian|) snaps to 0.5, a tight one stays
    # imbalanced -> the CI90 gate, not f_MAP alone, decides
    f_map = 0.45
    wide, _ = run_acs_somix(tmp_path, [somix_row(f_map=f_map, f_hessian=-50.0)],
                            maf90_threshold=0.49)
    tight, _ = run_acs_somix(tmp_path, [somix_row(f_map=f_map, f_hessian=-1e6)],
                             maf90_threshold=0.49)
    assert wide["f"].iloc[0] == pytest.approx(0.5)
    assert tight["f"].iloc[0] == pytest.approx(f_map)


@pytest.mark.stat
def test_somix_snap_uses_90th_percentile_not_two_sigma(tmp_path):
    # f_MAP + z90*sd < threshold < f_MAP + 2*sd: the CI90 (z90 = 1.2816) gate must NOT
    # snap, confirming a 90th-percentile bound rather than a 2-sigma bound.
    f_map, sd = 0.44, 0.03
    assert f_map + 1.2816 * sd < 0.49 < f_map + 2 * sd   # the test's own precondition
    out, _ = run_acs_somix(tmp_path, [somix_row(f_map=f_map, f_hessian=-1.0 / sd ** 2)],
                           maf90_threshold=0.49)
    assert out["f"].iloc[0] == pytest.approx(f_map)


@pytest.mark.stat
def test_somix_mu_and_var_f_from_hessian(tmp_path):
    # var_f = -1 / f_hessian feeds the minor/major sigmas
    f_hessian = -400.0
    out, _ = run_acs_somix(tmp_path, [somix_row(f_map=0.3, f_hessian=f_hessian)])
    row = out.iloc[0]
    tau, f, var_f = row["tau"], row["f"], -1.0 / f_hessian
    assert row["mu.minor"] == pytest.approx(f * tau)
    assert row["mu.major"] == pytest.approx((1 - f) * tau)
    assert row["sigma.minor"] == pytest.approx(
        np.sqrt(tau ** 2 * var_f + f ** 2 * row["sigma.tau"] ** 2), rel=1e-9)
    assert row["sigma.major"] == pytest.approx(
        np.sqrt(tau ** 2 * var_f + (1 - f) ** 2 * row["sigma.tau"] ** 2), rel=1e-9)


@pytest.mark.stat
def test_somix_sex_handling(tmp_path):
    # male: X/Y scaled by haploid ploidy 1, f NaN on both
    out, _ = run_acs_somix(tmp_path, [somix_row(contig="X", log_tcr=0.1),
                                      somix_row(contig="Y", log_tcr=0.1)], sex="XY")
    by_chr = out.set_index(out["Chromosome"].astype(str))
    assert by_chr.loc["X", "tau"] == pytest.approx(np.exp(0.1))   # ploidy 1
    assert np.isnan(by_chr.loc["X", "f"])
    assert np.isnan(by_chr.loc["Y", "f"])


@pytest.mark.stat
def test_somix_subset_by_sample_id(tmp_path):
    rows = [somix_row(contig="1", start=1, end=1001, sample_id="S1"),
            somix_row(contig="1", start=1, end=1001, sample_id="S2"),
            somix_row(contig="2", start=1, end=1001, sample_id="S2")]
    out, _ = run_acs_somix(tmp_path, rows, sample_id="S2")
    assert len(out) == 2
    assert set(out["Chromosome"].astype(str)) == {"1", "2"}


@pytest.mark.io
def test_somix_skew_from_lambda_map(tmp_path):
    # MEAN_BIAS is lambda_MAP for the row matching --sample_id
    af = pd.DataFrame({"sample_id": ["S1", "S2"], "s_MAP": [70.0, 65.0],
                       "lambda_MAP": [1.5, 1.1]})
    af_path = write_tsv(tmp_path / "af.param", af)
    out, out_dir = run_acs_somix(tmp_path, [somix_row(sample_id="S1")],
                                 sample_id="S1", af_parameters=af_path)
    skew = float(np.loadtxt(out_dir / "S1.acs.seg.skew"))   # --sample_id names the output
    assert skew == pytest.approx(2.0 / (1.0 + 1.5))   # lambda_MAP of S1


@pytest.mark.io
def test_somix_skew_requires_sample_id_when_multi_row(tmp_path):
    af = pd.DataFrame({"sample_id": ["S1", "S2"], "s_MAP": [70.0, 65.0],
                       "lambda_MAP": [1.5, 1.1]})
    af_path = write_tsv(tmp_path / "af.param", af)
    with pytest.raises(Exception):
        run_acs_somix(tmp_path, [somix_row(sample_id="S1")], af_parameters=af_path)


@pytest.mark.unit
def test_acs_missing_seg_errors(tmp_path):
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    base = dict(af_parameters=None, sample_id=None, output_dir=str(out_dir),
                min_hets=0, min_probes=0, maf90_threshold=0.49, sex="XXY",
                normal_ploidy=2, verbose=False)
    with pytest.raises(Exception):
        acs.convert_model_segments_to_alleliccapseg(args_ns(seg=None, **base))


@pytest.mark.unit
def test_acs_unknown_seg_format_errors(tmp_path):
    # a segmentation with neither ModelSegments nor somix marker columns is rejected
    seg = write_tsv(tmp_path / "test.seg", pd.DataFrame({"contig": ["1"], "value": [0.1]}))
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    args = args_ns(seg=seg, sample_id=None, af_parameters=None, output_dir=str(out_dir),
                   min_hets=0, min_probes=0, maf90_threshold=0.49, sex="XXY",
                   normal_ploidy=2, verbose=False)
    with pytest.raises(Exception):
        acs.convert_model_segments_to_alleliccapseg(args)


# --------------------------------------------------------------------------- #
# Black-box expectations (synthetic cohort, contract-derived).
# --------------------------------------------------------------------------- #
ACS_COLUMNS = ["Chromosome", "Start.bp", "End.bp", "n_probes", "length", "n_hets",
               "f", "tau", "sigma.tau", "mu.minor", "sigma.minor", "mu.major",
               "sigma.major", "SegLabelCNLOH"]


def _run_acs_synth(synth, out, patient, sample, sex, **extra):
    P = synth.files["patients"][patient]["samples"][sample]
    out.mkdir(parents=True, exist_ok=True)
    argv = ["--output_dir", str(out), "--seg", P["model_seg"],
            "--af_parameters", P["af_param"], "--sex", sex]
    for k, v in extra.items():
        argv += [f"--{k}", str(v)]
    run_script("acs_conversion.py", argv)
    return pd.read_csv(out / f"{sample}.modelFinal.acs.seg", sep="\t")


@pytest.mark.regression
def test_acs_output_column_schema(synth, tmp_path):
    acs = _run_acs_synth(synth, tmp_path / "s", "PXX", "PXX-T1", "XX")
    assert list(acs.columns) == ACS_COLUMNS


@pytest.mark.regression
def test_acs_male_x_tau_uses_haploid_ploidy(synth, tmp_path):
    P = synth.files["patients"]["PXY"]["samples"]["PXY-T1"]
    acs = _run_acs_synth(synth, tmp_path / "s", "PXY", "PXY-T1", "XY")
    seg = pd.read_csv(P["model_seg"], sep="\t", comment="@")
    log2 = dict(zip(zip(seg.CONTIG.astype(str), seg.START),
                    seg.LOG2_COPY_RATIO_POSTERIOR_50))
    for r in acs.itertuples():
        cr = 2.0 ** log2[(str(r.Chromosome), r._2)]   # _2 == Start.bp
        ploidy = 1 if str(r.Chromosome) in ("X", "Y") else 2
        if np.isnan(r.tau):
            continue
        assert r.tau == pytest.approx(ploidy * cr, rel=1e-4), (r.Chromosome, r._2)


@pytest.mark.regression
def test_acs_sex_chromosome_presence(synth, tmp_path):
    fem = _run_acs_synth(synth, tmp_path / "fem", "PXX", "PXX-T1", "XX")
    assert not (fem["Chromosome"].astype(str) == "Y").any()   # female: chrY ploidy 0
    xxy = _run_acs_synth(synth, tmp_path / "xxy", "PXXY", "PXXY-T1", "XXY")
    assert (xxy["Chromosome"].astype(str) == "Y").any()       # XXY: chrY ploidy 1


@pytest.mark.regression
def test_acs_sigma_na_iff_f_na(synth, tmp_path):
    # contract (per maintainer): where f is NA, sigma.minor and sigma.major are NA;
    # where f is present they are positive. PXY (male) has NA f on X/Y.
    acs = _run_acs_synth(synth, tmp_path / "s", "PXY", "PXY-T1", "XY")
    for f, smin, smaj, chrom in zip(acs["f"], acs["sigma.minor"], acs["sigma.major"],
                                    acs["Chromosome"]):
        if pd.isna(f):
            assert pd.isna(smin) and pd.isna(smaj), chrom
        else:
            assert smin > 0 and smaj > 0, chrom


@pytest.mark.regression
def test_acs_maf90_threshold_snaps_f_to_half(synth, tmp_path):
    seg = pd.read_csv(synth.files["patients"]["PXX"]["samples"]["PXX-T1"]["model_seg"],
                      sep="\t", comment="@")
    row = seg[(seg.CONTIG.astype(str) == "2") & (seg.START == 700000)].iloc[0]
    m90 = float(row.MINOR_ALLELE_FRACTION_POSTERIOR_90)
    assert 0.0 < m90 < 0.5
    below = _run_acs_synth(synth, tmp_path / "lo", "PXX", "PXX-T1", "XX", maf90_threshold=m90 - 0.05)
    above = _run_acs_synth(synth, tmp_path / "hi", "PXX", "PXX-T1", "XX", maf90_threshold=m90 + 0.05)

    def f_at(df):
        return float(df[(df.Chromosome.astype(str) == "2") & (df["Start.bp"] == 700000)]["f"].iloc[0])

    assert f_at(below) == pytest.approx(0.5)   # MAF_90 > threshold -> snapped to 0.5
    assert f_at(above) < 0.5                    # MAF_90 < threshold -> kept imbalanced


@pytest.mark.regression
def test_acs_conversion_recovers_f(synth, tmp_path):
    # minor allele fraction f recovered from truth allelic CN (autosomes)
    acs_df = _run_acs_synth(synth, tmp_path / "s", "PXX", "PXX-T1", "XX",
                            maf90_threshold=0.49)
    seg_truth = synth.truth["segments"]
    seg_truth = seg_truth[(seg_truth.patient_id == "PXX") & (seg_truth.sample_id == "PXX-T1")]
    purity = float(synth.truth["patient"].query(
        "patient_id=='PXX' and sample_id=='PXX-T1'")["purity"].iloc[0])
    for _, row in acs_df.iterrows():
        match = seg_truth[(seg_truth.contig.astype(str) == str(row["Chromosome"])) &
                          (seg_truth.start == row["Start.bp"])]
        if match.empty:
            continue
        cn = SegmentCN(int(match.cn_a1.iloc[0]), int(match.cn_a2.iloc[0]))
        exp_f = expected_maf(cn, purity)
        # autosomes only (X handling differs by sex); skip rows snapped to 0.5
        if str(row["Chromosome"]) in {"X", "Y"} or pd.isna(row["f"]) or row["f"] == 0.5:
            continue
        assert abs(row["f"] - exp_f) < 0.05, (row["Chromosome"], row["Start.bp"], row["f"], exp_f)


@pytest.mark.regression
def test_acs_passthrough_probe_het_and_length(synth, tmp_path):
    # n_probes / n_hets carry the ModelSegments point counts unchanged; length is
    # the segment span. These pass through the conversion untouched.
    acs_df = _run_acs_synth(synth, tmp_path / "s", "PXX", "PXX-T1", "XX")
    seg = pd.read_csv(synth.files["patients"]["PXX"]["samples"]["PXX-T1"]["model_seg"],
                      sep="\t", comment="@")
    probes = dict(zip(zip(seg.CONTIG.astype(str), seg.START), seg.NUM_POINTS_COPY_RATIO))
    hets = dict(zip(zip(seg.CONTIG.astype(str), seg.START), seg.NUM_POINTS_ALLELE_FRACTION))
    checked = 0
    for r in acs_df.itertuples():
        key = (str(r.Chromosome), r._2)   # _2 == Start.bp
        assert int(r.n_probes) == int(probes[key]), key
        assert int(r.n_hets) == int(hets[key]), key
        assert r.length == pytest.approx(r._3 - r._2), key   # _3 == End.bp
        checked += 1
    assert checked >= 6


@pytest.mark.regression
def test_acs_conversion_skew(synth, tmp_path):
    # skew = 2 / (1 + MEAN_BIAS); the cohort's af_parameters carry MEAN_BIAS = 1.05
    out = tmp_path / "s"
    _run_acs_synth(synth, out, "PXX", "PXX-T1", "XX")
    skew = float(np.loadtxt(out / "PXX-T1.modelFinal.acs.seg.skew"))
    assert skew == pytest.approx(2.0 / (1.0 + 1.05), rel=1e-6)
