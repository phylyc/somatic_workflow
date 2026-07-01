"""Unit / IO / stat tests for genotype.py."""

import types

import numpy as np
import pandas as pd
import pytest

import genotype as gt
from helpers import (args_ns, pileup_df, write_pileup, run_genotype, run_script,
                     gt_class, synth_wes_pileup, synth_runs)
from tests.synthetic.scenarios import write_empty_pileup


# --------------------------------------------------------------------------- #
# helpers local to this module
# --------------------------------------------------------------------------- #
def write_vcf9(path, records, contigs=("1", "2", "X")):
    """Write a 9-column (no-sample) VCF that genotype.VCF can parse."""
    lines = ["##fileformat=VCFv4.2"]
    for c in contigs:
        lines.append(f"##contig=<ID={c},length=1000>")
    lines.append("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                             "FILTER", "INFO", "FORMAT"]))
    for chrom, pos, ref, alt in records:
        lines.append("\t".join([str(chrom), str(pos), ".", ref, alt, ".",
                                 "PASS", ".", "GT"]))
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def named(bam_sample_name, file_path="f"):
    return types.SimpleNamespace(bam_sample_name=bam_sample_name,
                                 file_path=file_path)


# --------------------------------------------------------------------------- #
# Genotyper pure functions
# --------------------------------------------------------------------------- #
@pytest.mark.unit
@pytest.mark.parametrize("sex,x,y", [
    ("XX", 2, 0), ("XY", 1, 1), ("XXY", 2, 1), ("Female", 2, 0),
    ("Male", 1, 1), ("XYY", 1, 1), ("XXX", 2, 0),
])
def test_get_contig_ploidies(sex, x, y):
    p = gt.Genotyper.get_contig_ploidies(sex)
    assert p["1"] == 2 and p["MT"] == 1
    assert p["X"] == x and p["Y"] == y
    assert p["chrX"] == x  # chr-prefixed mirror


@pytest.mark.unit
def test_get_error_prob_formula_and_clip():
    g = gt.Genotyper(min_error_rate=1e-4, max_error_rate=1e-1)
    # 3/2 * other/total = 1.5 * 3 / 30 = 0.15 -> clipped to max 0.1
    p = pileup_df([("1", 1, 17, 10, 3, 0.5)])
    assert g.get_error_prob(p) == pytest.approx(0.1)
    # tiny error -> clipped to min
    p2 = pileup_df([("1", 1, 1000, 0, 0, 0.5)])
    assert g.get_error_prob(p2) == pytest.approx(1e-4)


@pytest.mark.unit
def test_select_confident_calls_threshold():
    g = gt.Genotyper(min_genotype_likelihood=0.95)
    df = pd.DataFrame({
        "0/0": [0.99, 0.5], "0/1": [0.0, 0.3], "1/1": [0.0, 0.1],
        "./.": [0.01, 0.1],
    })
    out = g.select_confident_calls(df, normalized=True)
    assert len(out) == 1 and out.index[0] == 0


@pytest.mark.unit
def test_cross_check_sample_name_agreement():
    assert gt.cross_check_sample_name(named("S"), named("S"), named(None)) == "S"


@pytest.mark.unit
def test_cross_check_sample_name_disagreement_warns():
    with pytest.warns(UserWarning):
        name = gt.cross_check_sample_name(named("A"), named("B"))
    assert name in ("A", "B")


@pytest.mark.unit
def test_cross_check_sample_name_all_none():
    assert gt.cross_check_sample_name(named(None), named(None)) is None


@pytest.mark.unit
def test_join_pileups_sums_and_max_af():
    p1, p2 = gt.Pileup(), gt.Pileup()
    p1.df = pileup_df([("1", 1, 10, 5, 1, 0.2)])
    p2.df = pileup_df([("1", 1, 4, 6, 2, 0.7)])
    p1.bam_sample_name = p2.bam_sample_name = "S"
    j = gt.join_pileups([p1, p2])
    row = j.df.iloc[0]
    assert row["ref_count"] == 14 and row["alt_count"] == 11
    assert row["other_alt_count"] == 3 and row["allele_frequency"] == 0.7


@pytest.mark.unit
def test_join_contaminations_inverse_variance_weighted():
    c1, c2 = gt.Contamination(), gt.Contamination()
    c1.value, c1.error, c1.bam_sample_name = 0.10, 0.01, "S"
    c2.value, c2.error, c2.bam_sample_name = 0.20, 0.02, "S"
    j = gt.join_contaminations([c1, c2])
    # weights 1/err^2 -> heavily favors the lower-error estimate (0.10)
    assert 0.10 <= j.value <= 0.12


# --------------------------------------------------------------------------- #
# interval algebra
# --------------------------------------------------------------------------- #
@pytest.mark.unit
def test_non_overlapping_decomposes_two_intervals():
    df = pd.DataFrame({"start": [1, 5], "end": [10, 15],
                       "minor_allele_fraction": [0.4, 0.3]})
    out = gt.non_overlapping(df)
    breaks = sorted(set(out["start"]).union(set(out["end"])))
    assert breaks == [1, 5, 10, 15]


@pytest.mark.unit
def test_merge_abutting_intervals_idempotent():
    idx = pd.MultiIndex.from_tuples(
        [("1", 1, 10), ("1", 10, 20)], names=["contig", "start", "end"])
    df = pd.DataFrame({"minor_allele_fraction": [0.3, 0.3]}, index=idx)
    once = gt.merge_abutting_intervals(df)
    twice = gt.merge_abutting_intervals(once)
    assert len(once) == len(twice)


# --------------------------------------------------------------------------- #
# Phaser static helpers
# --------------------------------------------------------------------------- #
@pytest.mark.unit
def test_exact_sign_max_gray_finds_opposite_signs():
    Vc = np.array([[1.0, 1.0, 1.0], [-1.0, -1.0, -1.0]])
    x = gt.Phaser.exact_sign_max_gray(Vc)
    assert x[0] != x[1]  # anti-correlated rows get opposite signs


@pytest.mark.unit
def test_best_improvement_monotone_objective():
    Vc = np.array([[1.0, 1.0, 1.0], [-1.0, -1.0, -1.0]])
    x = gt.Phaser.best_improvement(Vc, np.array([1, 1]))
    assert x[0] != x[1]


@pytest.mark.unit
def test_spectral_init_deterministic():
    Vc = np.array([[1.0, 0.2, -0.3], [0.9, 0.1, -0.4], [-1.0, 0.0, 0.5]])
    a = gt.Phaser.spectral_init(Vc)
    b = gt.Phaser.spectral_init(Vc)
    np.testing.assert_array_equal(a, b)


@pytest.mark.unit
def test_avg_pairwise_corr_identical_rows():
    rng = np.random.default_rng(0)
    base = rng.standard_normal(30)
    M = np.vstack([base, base])
    mask = np.ones(30, dtype=bool)
    assert gt.Phaser.avg_pairwise_corr(M, mask) == pytest.approx(1.0, abs=1e-6)


@pytest.mark.unit
def test_connected_components_shared_vs_disjoint():
    phaser = gt.Phaser()
    shared = phaser.connected_components_from_tracks(
        [np.array([1, 0, 0]), np.array([1, 0, 0])])
    assert len(shared) == 1 and sorted(shared[0]) == [0, 1]
    disjoint = phaser.connected_components_from_tracks(
        [np.array([1, 0, 0]), np.array([0, 0, 1])])
    assert len(disjoint) == 2


# --------------------------------------------------------------------------- #
# IO classes
# --------------------------------------------------------------------------- #
@pytest.mark.io
def test_pileup_parses_header_and_filters(tmp_path):
    df = pileup_df([("1", 1, 1, 0, 0, 0.1),     # depth 1
                    ("1", 2, 10, 10, 0, 0.5)])   # depth 20
    path = write_pileup(tmp_path / "p.pileup", df, sample="Test-Y")
    pu = gt.Pileup(path, min_read_depth=5)
    assert pu.bam_sample_name == "Test-Y"
    assert list(pu.df["position"]) == [2]


@pytest.mark.io
def test_contamination_defaults_when_no_file():
    c = gt.Contamination(None)
    assert c.value == 0
    assert c.error == gt.Contamination._min_error


@pytest.mark.io
def test_contamination_parses(tmp_path):
    p = tmp_path / "c.tsv"
    p.write_text("sample\tcontamination\tcontamination_error\nS\t0.05\t0.01\n")
    c = gt.Contamination(str(p))
    assert c.value == pytest.approx(0.05)
    assert c.bam_sample_name == "S"


@pytest.mark.io
def test_segments_dot_segments_format(tmp_path):
    p = tmp_path / "s.segments"
    p.write_text("contig\tstart\tend\tminor_allele_fraction\n1\t1\t100\t0.4\n")
    s = gt.Segments(str(p))
    assert list(s.df.columns) == ["contig", "start", "end", "minor_allele_fraction"]
    assert s.df["minor_allele_fraction"].iloc[0] == pytest.approx(0.4)


@pytest.mark.io
def test_afparameters_defaults_and_parse(tmp_path):
    d = gt.AFParameters(None, default_ref_bias=1.3)
    assert d.ref_bias == 1.3
    p = tmp_path / "af.param"
    p.write_text("PARAMETER_NAME\tPOSTERIOR_50\nMEAN_BIAS\t1.1\n"
                 "BIAS_VARIANCE\t0.01\nOUTLIER_PROBABILITY\t0.02\n")
    a = gt.AFParameters(str(p))
    assert a.ref_bias == pytest.approx(1.1)
    assert a.outlier_probability == pytest.approx(0.02)


@pytest.mark.io
def test_vcf_filters_to_biallelic_snvs_and_dedupes(tmp_path):
    vcf = write_vcf9(tmp_path / "v.vcf", [
        ("1", 100, "A", "G"),    # SNV keep
        ("1", 200, "AC", "G"),   # indel drop
        ("1", 300, "A", "ATG"),  # indel drop
        ("2", 50, "C", "T"),     # SNV keep
        ("2", 60, "C", "."),     # spanning/empty drop
        ("1", 100, "A", "T"),    # duplicate pos -> both dropped (keep=False)
    ])
    v = gt.VCF(vcf)
    positions = set(v.df.index.get_level_values("position"))
    assert positions == {50}  # 100 dropped as dup, indels/'.' dropped
    assert v.contigs[:3] == ["1", "2", "X"]


@pytest.mark.io
def test_join_vcfs_single_and_empty(tmp_path):
    assert gt.join_vcfs([]) is None
    v = gt.VCF(write_vcf9(tmp_path / "v.vcf", [("1", 100, "A", "G")]))
    assert gt.join_vcfs([v]) is v


@pytest.mark.io
def test_join_vcfs_dedupes_shared_locus(tmp_path):
    # two inputs share ("1", 100); each also has a unique locus
    v1 = gt.VCF(write_vcf9(tmp_path / "v1.vcf",
                           [("1", 100, "A", "G"), ("2", 50, "C", "T")]))
    v2 = gt.VCF(write_vcf9(tmp_path / "v2.vcf",
                           [("1", 100, "A", "T"), ("2", 60, "C", "T")]))
    joint = gt.join_vcfs([v1, v2])
    # the shared locus collapses to a single row
    assert joint.df.loc[[("1", 100)]].shape[0] == 1
    # unique loci from both inputs are preserved (3 distinct loci total)
    assert joint.df.index.is_unique
    assert set(joint.df.index) == {("1", 100), ("2", 50), ("2", 60)}


# --------------------------------------------------------------------------- #
# statistical cores
# --------------------------------------------------------------------------- #
@pytest.mark.stat
@pytest.mark.parametrize("ref,alt,expected_gt", [
    (50, 0, "0/0"), (0, 50, "1/1"), (25, 25, "0/1"),
])
def test_genotype_likelihoods_argmax(ref, alt, expected_gt):
    g = gt.Genotyper(sex="XXY")
    pu = pileup_df([("1", 100, ref, alt, 0, 0.2)])
    out = g.calculate_genotype_likelihoods_per_segment(pu, minor_af=0.47)
    row = out.iloc[0]
    best = row[["0/0", "0/1", "1/1"]].astype(float).idxmax()
    assert best == expected_gt


@pytest.mark.stat
def test_het_likelihood_is_sum_of_phased():
    g = gt.Genotyper(sex="XXY")
    pu = pileup_df([("1", 100, 25, 25, 0, 0.2)])
    out = g.calculate_genotype_likelihoods_per_segment(pu, minor_af=0.47)
    row = out.iloc[0]
    assert row["0/1"] == pytest.approx(row["0|1"] + row["1|0"], rel=1e-9)


@pytest.mark.stat
def test_ploidy_zero_only_outlier():
    # Y in an XX sample has ploidy 0 -> 0/0, 0/1, 1/1 all ~0, outlier dominates
    g = gt.Genotyper(sex="XX")
    pu = pileup_df([("Y", 100, 25, 25, 0, 0.2)])
    out = g.calculate_genotype_likelihoods_per_segment(pu, minor_af=0.47)
    row = out.iloc[0]
    assert row["./."] > row[["0/0", "0/1", "1/1"]].astype(float).max()


@pytest.mark.stat
def test_pop_af_prior_raises_het_likelihood():
    # Same het read evidence; a common SNP (high popAF) must get a higher variant
    # (het) genotype likelihood than a rare one — the population-AF prior.
    g = gt.Genotyper(sex="XXY")
    common = g.calculate_genotype_likelihoods_per_segment(
        pileup_df([("1", 100, 25, 25, 0, 0.45)]), minor_af=0.47).iloc[0]
    rare = g.calculate_genotype_likelihoods_per_segment(
        pileup_df([("1", 100, 25, 25, 0, 0.05)]), minor_af=0.47).iloc[0]
    assert common["0/1"] > rare["0/1"]


@pytest.mark.stat
def test_pop_af_prior_raises_hom_alt_likelihood():
    g = gt.Genotyper(sex="XXY")
    common = g.calculate_genotype_likelihoods_per_segment(
        pileup_df([("1", 100, 0, 50, 0, 0.45)])).iloc[0]
    rare = g.calculate_genotype_likelihoods_per_segment(
        pileup_df([("1", 100, 0, 50, 0, 0.05)])).iloc[0]
    assert common["1/1"] > rare["1/1"]


@pytest.mark.stat
def test_joint_genotype_likelihood_normalizes_and_calls():
    g = gt.Genotyper(sex="XXY", min_genotype_likelihood=0.0)
    pu1 = pileup_df([("1", 100, 0, 50, 0, 0.2), ("1", 200, 25, 25, 0, 0.2)])
    pu2 = pileup_df([("1", 100, 0, 40, 0, 0.2), ("1", 200, 20, 20, 0, 0.2)])
    pl1 = gt.PileupLikelihood(g.calculate_genotype_likelihoods(pu1), "T", "T")
    pl2 = gt.PileupLikelihood(g.calculate_genotype_likelihoods(pu2), "N", "N")
    joint = g.get_joint_genotype_likelihood([pl1, pl2])
    # the joint table carries only the unphased genotypes
    probs = joint[["0/0", "0/1", "1/1", "./."]].sum(axis=1)
    np.testing.assert_allclose(probs.to_numpy(), 1.0, atol=1e-6)
    assert set(joint["GT"]).issubset({"0/0", "0/1", "1/1", "./."})
    # locus 100 (all alt) should be called 1/1
    assert joint.loc[("1", 100), "GT"] == "1/1"


@pytest.mark.stat
def test_joint_normal_weight_upweights_normal():
    g = gt.Genotyper(sex="XXY", min_genotype_likelihood=0.0)
    # tumor leans 1/1 (all alt), normal leans 0/1 (balanced); heavy normal weight
    pu_t = pileup_df([("1", 100, 0, 50, 0, 0.2)])
    pu_n = pileup_df([("1", 100, 25, 25, 0, 0.2)])
    plt_ = gt.PileupLikelihood(g.calculate_genotype_likelihoods(pu_t), "T", "T")
    pln = gt.PileupLikelihood(g.calculate_genotype_likelihoods(pu_n), "N", "N")
    joint = g.get_joint_genotype_likelihood([plt_, pln], normal_samples=["N"],
                                            normal_to_tumor_weight=100.0)
    assert joint.loc[("1", 100), "GT"] == "0/1"


# --------------------------------------------------------------------------- #
# Black-box expectations (synthetic cohort, contract-derived).
# --------------------------------------------------------------------------- #
@pytest.mark.regression
def test_genotype_invariant_to_sample_subset(synth, tmp_path):
    # germline genotype is a patient-level property: same call at every shared site
    full = run_genotype(synth, tmp_path / "full", "PXX", ["PXX-N", "PXX-T1", "PXX-T3"])
    subset = run_genotype(synth, tmp_path / "sub", "PXX", ["PXX-N", "PXX-T1"])
    fmap = {(str(r.CHROM), int(r.POS)): gt_class(r.PXX) for r in full.itertuples()}
    smap = {(str(r.CHROM), int(r.POS)): gt_class(r.PXX) for r in subset.itertuples()}
    shared = set(fmap) & set(smap)
    assert len(shared) >= 8
    for k in shared:
        assert fmap[k] == smap[k], (k, fmap[k], smap[k])


@pytest.mark.regression
def test_genotype_invariant_to_normal_flag(synth, tmp_path):
    # germline genotype class is intrinsic to the individual: passing --normal_sample
    # must not change the call set / classes for clean data
    with_n = run_genotype(synth, tmp_path / "n", "PXX", ["PXX-N", "PXX-T1", "PXX-T3"],
                          normal="PXX-N")
    without = run_genotype(synth, tmp_path / "x", "PXX", ["PXX-N", "PXX-T1", "PXX-T3"])
    a = {(str(r.CHROM), int(r.POS)): gt_class(r.PXX) for r in with_n.itertuples()}
    b = {(str(r.CHROM), int(r.POS)): gt_class(r.PXX) for r in without.itertuples()}
    assert set(a) == set(b)
    for k in a:
        assert a[k] == b[k], (k, a[k], b[k])


@pytest.mark.regression
def test_genotype_phasing_matches_truth_haplotype(synth, tmp_path):
    # within each CR segment, the het phase partition matches truth haplotype up to a flip
    vcf = run_genotype(synth, tmp_path / "g", "PXX", ["PXX-N", "PXX-T1", "PXX-T3"])
    gt = {(str(r.CHROM), int(r.POS)): str(r.PXX) for r in vcf.itertuples()}
    truth = synth.truth["germline_sites"].query("patient_id=='PXX'")
    truth = truth[(truth.genotype == "0/1") & truth.is_snv]
    checked = 0
    for _seg, grp in truth.groupby("segment_contig"):
        pairs = [(gt.get((str(g.contig), int(g.position))) == "1|0", g.haplotype == "A")
                 for _, g in grp.iterrows()
                 if gt.get((str(g.contig), int(g.position))) in ("0|1", "1|0")]
        if len(pairs) < 2:
            continue
        checked += 1
        agree = sum(o == h for o, h in pairs)
        assert agree == len(pairs) or agree == 0, (_seg, pairs)  # consistent up to flip
    assert checked >= 1


@pytest.mark.regression
def test_genotype_count_tables_correlation_and_provenance(synth, tmp_path):
    out = tmp_path / "g"
    vcf = run_genotype(synth, out, "PXX", ["PXX-N", "PXX-T1", "PXX-T3"])

    # count tables are row-aligned with the VCF (only genotyped SNV sites) and
    # reproduce the input pileup counts for PXX-T1
    ref = pd.read_csv(out / "PXX.germline.ref_count.tsv", sep="\t")
    alt = pd.read_csv(out / "PXX.germline.alt_count.tsv", sep="\t")
    assert len(ref) == len(vcf) == len(alt)
    pu = pd.read_csv(synth_wes_pileup(synth, "PXX", "PXX-T1"), sep="\t", comment="#")
    pmap = {(str(r.contig), int(r.position)): (r.ref_count, r.alt_count) for r in pu.itertuples()}
    checked = 0
    for i, v in enumerate(vcf.itertuples()):
        key = (str(v.CHROM), int(v.POS))
        if key in pmap:
            assert int(ref["PXX-T1"].iloc[i]) == pmap[key][0], key
            assert int(alt["PXX-T1"].iloc[i]) == pmap[key][1], key
            checked += 1
    assert checked >= 8

    # same-patient samples are highly correlated; diagonal ~1
    corr = pd.read_csv(out / "PXX.sample_correlation.tsv", sep="\t", index_col=0)
    for s in corr.index:
        assert corr.loc[s, s] == pytest.approx(1.0, abs=1e-6)
    off = [corr.loc[a, b] for a in corr.index for b in corr.columns if a != b]
    assert min(off) >= 0.9

    # called sites ⊆ candidate VCF; no Y / MT invented
    ref_vcf = pd.read_csv(synth.files["patients"]["PXX"]["reference_vcf"], sep="\t",
                          comment="#", header=None,
                          names=["CHROM", "POS", "ID", "REF", "ALT", "Q", "F", "I", "FMT"])
    candidates = set(zip(ref_vcf.CHROM.astype(str), ref_vcf.POS.astype(int)))
    called = set(zip(vcf.CHROM.astype(str), vcf.POS.astype(int)))
    assert called <= candidates
    assert not any(c in ("Y", "MT") for c, _ in called)


@pytest.mark.regression
def test_genotype_recovers_germline_calls(synth, tmp_path):
    # germline genotype class recovered at every SNV; non-SNV truth sites filtered
    P = synth.files["patients"]["PXX"]["samples"]
    out = tmp_path / "gt"
    out.mkdir()
    argv = ["--output_dir", str(out), "--patient", "PXX",
            "--variant", synth.files["patients"]["PXX"]["reference_vcf"]]
    for s in ["PXX-N", "PXX-T1", "PXX-T3"]:
        argv += ["--sample", s, "--pileup", synth_wes_pileup(synth, "PXX", s),
                 "--segments", P[s]["segments"], "--contamination", P[s]["contamination"]]
    argv += ["--normal_sample", "PXX-N", "--min_read_depth", "10",
             "--min_genotype_likelihood", "0.95", "--threads", "1"]
    run_script("genotype.py", argv)

    vcf = pd.read_csv(out / "PXX.germline.vcf", sep="\t", comment="#", header=None,
                      names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                             "INFO", "FORMAT", "PXX"])
    called = {(str(r.CHROM), int(r.POS)): gt_class(r.PXX) for r in vcf.itertuples()}

    truth = synth.truth["germline_sites"]
    truth = truth[truth.patient_id == "PXX"]
    want = {"0/0": "hom_ref", "1/1": "hom_alt", "0/1": "het"}
    for _, g in truth.iterrows():
        key = (str(g.contig), int(g.position))
        if not g.is_snv:
            assert key not in called, f"non-SNV {key} should be filtered"
            continue
        assert key in called, f"missing germline site {key}"
        assert called[key] == want[g.genotype], (key, called[key], g.genotype)


@pytest.mark.regression
def test_genotype_duplicate_sample_names(synth, tmp_path):
    # PXY-T1 has two WES runs; passing the same assigned name twice exercises
    # GenotypeData.deduplicate_samples (the legacy duplicate-sample case).
    pileups = [r["pileup"] for r in synth_runs(synth, "PXY", "PXY-T1") if r["pileup"]]
    assert len(pileups) == 2
    P = synth.files["patients"]["PXY"]["samples"]["PXY-T1"]
    out = tmp_path / "o"
    out.mkdir()
    argv = ["--output_dir", str(out), "--patient", "PXY",
            "--variant", synth.files["patients"]["PXY"]["reference_vcf"]]
    for p in pileups:
        argv += ["--sample", "PXY-T1", "--pileup", p,
                 "--segments", P["segments"], "--contamination", P["contamination"]]
    argv += ["--threads", "1"]
    run_script("genotype.py", argv)
    assert (out / "PXY.germline.vcf").exists()


@pytest.mark.regression
def test_genotype_without_normal(synth, tmp_path):
    # tumors only, no --normal_sample (legacy no-normal case)
    vcf = run_genotype(synth, tmp_path / "o", "PXX", ["PXX-T1", "PXX-T3"])
    assert len(vcf) > 0


@pytest.mark.regression
def test_genotype_empty_pileup_sample_is_degenerate_not_fatal(synth, tmp_path):
    # a sample whose pileup is empty must not break the patient genotype: the run
    # succeeds, the real samples still get called, and the empty sample contributes
    # no counts (all-zero/NA column) with a degenerate (NaN) self-correlation.
    P = synth.files["patients"]["PXX"]["samples"]
    empty = write_empty_pileup(tmp_path / "empty.pileup", "PXX-EMPTY")
    out = tmp_path / "g"
    out.mkdir()
    argv = ["--output_dir", str(out), "--patient", "PXX",
            "--variant", synth.files["patients"]["PXX"]["reference_vcf"]]
    for s in ["PXX-N", "PXX-T1"]:
        argv += ["--sample", s, "--pileup", synth_wes_pileup(synth, "PXX", s),
                 "--segments", P[s]["segments"], "--contamination", P[s]["contamination"]]
    # empty-pileup sample reuses T1's segments/contamination
    argv += ["--sample", "PXX-EMPTY", "--pileup", empty,
             "--segments", P["PXX-T1"]["segments"],
             "--contamination", P["PXX-T1"]["contamination"], "--threads", "1"]
    run_script("genotype.py", argv)   # exit 0 expected (not fatal)

    vcf = pd.read_csv(out / "PXX.germline.vcf", sep="\t", comment="#", header=None,
                      names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                             "INFO", "FORMAT", "PXX"])
    assert len(vcf) > 0   # the real samples still drive a non-empty call set

    alt = pd.read_csv(out / "PXX.germline.alt_count.tsv", sep="\t")
    # the empty sample's count column carries no signal
    col = pd.to_numeric(alt["PXX-EMPTY"], errors="coerce").fillna(0)
    assert (col == 0).all()
    # real samples retain a self-correlation of 1; the empty sample is degenerate
    corr = pd.read_csv(out / "PXX.sample_correlation.tsv", sep="\t", index_col=0)
    for s in ["PXX-N", "PXX-T1"]:
        assert corr.loc[s, s] == pytest.approx(1.0, abs=1e-6)
    assert pd.isna(corr.loc["PXX-EMPTY", "PXX-T1"]) or corr.loc["PXX-EMPTY", "PXX-T1"] == 0


@pytest.mark.regression
def test_genotype_calls_only_biallelic_snvs_aligned_to_count_tables(synth, tmp_path):
    # the genotyped call set spans only biallelic SNV sites, and the per-sample
    # count tables are row-aligned to exactly that set (one row per called SNV)
    out = tmp_path / "g"
    vcf = run_genotype(synth, out, "PXX", ["PXX-N", "PXX-T1", "PXX-T3"])
    assert (vcf.REF.str.len() == 1).all()
    assert (vcf.ALT.str.len() == 1).all()
    for name in ("ref_count", "alt_count"):
        tbl = pd.read_csv(out / f"PXX.germline.{name}.tsv", sep="\t")
        assert len(tbl) == len(vcf)
