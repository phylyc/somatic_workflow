"""Unit tests for validate_inputs.py: pure helpers + the rule engine in main()."""

import pytest

import validate_inputs as vi
from helpers import run_script


# --------------------------------------------------------------------------- #
# pure helpers
# --------------------------------------------------------------------------- #
@pytest.mark.unit
@pytest.mark.parametrize("s,expected", [
    ("true", True), ("True", True), ("TRUE", True),
    ("false", False), ("", False), ("1", False), ("yes", False), (None, False),
])
def test_parse_bool(s, expected):
    assert vi.parse_bool(s) is expected


@pytest.mark.unit
def test_sq_contigs_extracts_sn_in_order():
    lines = [
        "@HD\tVN:1.6\n",
        "@SQ\tSN:1\tLN:100\n",
        "@SQ\tSN:2\tLN:100\n",
        "@RG\tID:x\n",
        "@SQ\tSN:X\tLN:100\n",
    ]
    assert vi.sq_contigs(lines) == ["1", "2", "X"]


@pytest.mark.unit
@pytest.mark.parametrize("sub,full,expected", [
    (["1", "2"], ["1", "2", "3"], True),
    (["1", "3"], ["1", "2", "3"], True),
    (["3", "1"], ["1", "2", "3"], False),   # wrong order
    ([], ["1", "2"], True),                  # empty subseq
    (["4"], ["1", "2", "3"], False),         # missing element
])
def test_is_ordered_subsequence(sub, full, expected):
    assert vi.is_ordered_subsequence(sub, full) is expected


# --------------------------------------------------------------------------- #
# rule engine via main()
# --------------------------------------------------------------------------- #
def run_main(monkeypatch, capsys, argv):
    monkeypatch.setattr("sys.argv", ["validate_inputs.py", *argv])
    code = 0
    try:
        vi.main()
    except SystemExit as e:
        code = e.code or 0
    out = capsys.readouterr().out
    return code, out


def base_argv(n=1):
    """A minimal valid invocation: n bams, n bais, n sample names, all runs off."""
    argv = ["--n_bams", str(n), "--n_bais", str(n)]
    for i in range(n):
        argv += ["--sample_names", f"S{i}"]
    return argv


@pytest.mark.unit
def test_valid_minimal_passes(monkeypatch, capsys):
    code, out = run_main(monkeypatch, capsys, base_argv(1))
    assert code == 0
    assert "Input validation passed" in out


@pytest.mark.unit
def test_no_bams_errors(monkeypatch, capsys):
    code, out = run_main(monkeypatch, capsys,
                         ["--n_bams", "0", "--n_bais", "0"])
    assert code == 1
    assert "no bams provided" in out


@pytest.mark.unit
def test_bai_count_mismatch_errors(monkeypatch, capsys):
    code, out = run_main(monkeypatch, capsys,
                         ["--n_bams", "2", "--n_bais", "1",
                          "--sample_names", "A", "--sample_names", "B"])
    assert code == 1
    assert "bais has 1 entries but there are 2 bams" in out


@pytest.mark.unit
def test_per_run_array_length_mismatch_errors(monkeypatch, capsys):
    argv = base_argv(2)  # 2 bams, 2 sample names
    argv += ["--is_paired_end", "true"]  # only 1 entry for 2 bams
    code, out = run_main(monkeypatch, capsys, argv)
    assert code == 1
    assert "is_paired_end has 1 entries but there are 2 bams" in out


@pytest.mark.unit
def test_normal_sample_names_length_exempt(monkeypatch, capsys):
    # normal_sample_names is a subset, not 1:1 with bams -> no length error
    argv = base_argv(2)
    argv += ["--normal_sample_names", "S0"]
    code, out = run_main(monkeypatch, capsys, argv)
    assert code == 0


@pytest.mark.unit
def test_unknown_normal_sample_errors(monkeypatch, capsys):
    argv = base_argv(2)
    argv += ["--normal_sample_names", "TYPO"]
    code, out = run_main(monkeypatch, capsys, argv)
    assert code == 1
    assert "normal_sample_names not found among sample_names" in out


@pytest.mark.unit
def test_tcr_requires_pon_or_annotated_intervals(monkeypatch, capsys):
    argv = base_argv(1)
    argv += ["--run_collect_total_read_counts", "true"]
    code, out = run_main(monkeypatch, capsys, argv)
    assert code == 1
    assert "total copy ratio" in out and "DenoiseReadCounts" in out


@pytest.mark.unit
def test_tcr_satisfied_by_annotated_intervals(monkeypatch, capsys):
    argv = base_argv(1)
    argv += ["--run_collect_total_read_counts", "true",
             "--annotated_target_intervals", "annot.tsv"]
    code, out = run_main(monkeypatch, capsys, argv)
    assert code == 0


@pytest.mark.unit
def test_tcr_satisfied_by_nonempty_pon(monkeypatch, capsys, tmp_path):
    pon = tmp_path / "pon.hdf5"
    pon.write_text("not empty")
    argv = base_argv(1)
    argv += ["--run_collect_total_read_counts", "true",
             "--cnv_panel_of_normals", str(pon)]
    code, out = run_main(monkeypatch, capsys, argv)
    assert code == 0


@pytest.mark.unit
def test_empty_pon_does_not_satisfy_tcr(monkeypatch, capsys, tmp_path):
    pon = tmp_path / "pon.hdf5"
    pon.write_text("")  # 0 bytes
    argv = base_argv(1)
    argv += ["--run_collect_total_read_counts", "true",
             "--cnv_panel_of_normals", str(pon)]
    code, out = run_main(monkeypatch, capsys, argv)
    assert code == 1
    assert "0-byte" in out or "neither" in out


@pytest.mark.unit
def test_allelic_counts_requires_germline_alleles(monkeypatch, capsys):
    argv = base_argv(1)
    argv += ["--run_collect_allelic_read_counts", "true"]
    code, out = run_main(monkeypatch, capsys, argv)
    assert code == 1
    assert "common_germline_alleles" in out


@pytest.mark.unit
def test_allelic_counts_requires_index(monkeypatch, capsys):
    argv = base_argv(1)
    argv += ["--run_collect_allelic_read_counts", "true",
             "--has_common_germline_alleles", "true"]
    code, out = run_main(monkeypatch, capsys, argv)
    assert code == 1
    assert "without its index" in out


@pytest.mark.unit
def test_realignment_requires_image(monkeypatch, capsys):
    argv = base_argv(1)
    argv += ["--run_realignment_filter", "true"]
    code, out = run_main(monkeypatch, capsys, argv)
    assert code == 1
    assert "realignment" in out.lower()


@pytest.mark.unit
def test_clonal_decomp_without_model_segments_warns(monkeypatch, capsys):
    argv = base_argv(1)
    argv += ["--run_clonal_decomposition", "true"]
    code, out = run_main(monkeypatch, capsys, argv)
    assert code == 0  # warning, not error
    assert "WARNING" in out and "run_clonal_decomposition" in out


@pytest.mark.unit
def test_missing_sex_warns_not_errors(monkeypatch, capsys):
    code, out = run_main(monkeypatch, capsys, base_argv(1))
    assert code == 0
    assert "sex is not set" in out


@pytest.mark.unit
def test_run_used_for_neither_tcr_nor_acr_warns(monkeypatch, capsys):
    argv = base_argv(1)
    argv += ["--use_for_tCR", "false", "--use_for_aCR", "false"]
    code, out = run_main(monkeypatch, capsys, argv)
    assert code == 0
    assert "neither tCR nor aCR" in out


@pytest.mark.unit
def test_multiple_errors_reported_together(monkeypatch, capsys):
    # both bai mismatch and unknown normal -> count >= 2, single failure summary
    argv = ["--n_bams", "2", "--n_bais", "1",
            "--sample_names", "A", "--sample_names", "B",
            "--normal_sample_names", "TYPO"]
    code, out = run_main(monkeypatch, capsys, argv)
    assert code == 1
    assert "error(s)" in out
    assert out.count("ERROR:") >= 2


# --------------------------------------------------------------------------- #
# Synthetic-cohort manifest (end-to-end via subprocess): the generated manifest
# yields a well-formed config, and a corrupted bai count is rejected.
# --------------------------------------------------------------------------- #
def _manifest_argv(synth, patient_id):
    pat = next(p for p in synth.manifest["patients"] if p["patient_id"] == patient_id)
    sample_names, tcr, acr = [], [], []
    for s in pat["samples"]:
        for run in s["runs"]:
            sample_names.append(s["sample_id"])
            tcr.append(str(run["use_for_tCR"]).lower())
            acr.append(str(run["use_for_aCR"]).lower())
    n = len(sample_names)
    argv = ["--n_bams", str(n), "--n_bais", str(n)]
    for s in sample_names:
        argv += ["--sample_names", s]
    for nm in pat["normal_sample_names"]:
        argv += ["--normal_sample_names", nm]
    for t in tcr:
        argv += ["--use_for_tCR", t]
    for a in acr:
        argv += ["--use_for_aCR", a]
    return argv, n


@pytest.mark.regression
def test_validate_inputs_manifest_passes(synth):
    argv, _ = _manifest_argv(synth, "PXX")
    run_script("validate_inputs.py", argv, expect_returncode=0)


@pytest.mark.regression
def test_validate_inputs_manifest_detects_bai_mismatch(synth):
    argv, n = _manifest_argv(synth, "PXX")
    i = argv.index("--n_bais")
    argv[i + 1] = str(n + 1)
    run_script("validate_inputs.py", argv, expect_returncode=1)
