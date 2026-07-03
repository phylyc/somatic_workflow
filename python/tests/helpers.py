"""Test helpers: synthetic-fixture builders, frame comparison, CLI runners.

Kept dependency-light (numpy/pandas only) and deterministic.
"""

import argparse
import pathlib
import subprocess
import sys

import numpy as np
import pandas as pd

PYDIR = pathlib.Path(__file__).resolve().parent.parent

PILEUP_COLUMNS = ["contig", "position", "ref_count", "alt_count",
                  "other_alt_count", "allele_frequency"]


# --------------------------------------------------------------------------- #
# argparse.Namespace builder
# --------------------------------------------------------------------------- #
def args_ns(**kwargs):
    """Build an argparse.Namespace for calling a script's ``core(args)`` fn."""
    return argparse.Namespace(**kwargs)


# --------------------------------------------------------------------------- #
# DataFrame builders
# --------------------------------------------------------------------------- #
def pileup_df(rows):
    """rows: iterable of (contig, position, ref, alt, other_alt, af) tuples
    or dicts with the PILEUP_COLUMNS keys."""
    if rows and isinstance(rows[0], dict):
        df = pd.DataFrame(rows, columns=PILEUP_COLUMNS)
    else:
        df = pd.DataFrame(rows, columns=PILEUP_COLUMNS)
    return df.astype({
        "contig": str, "position": int, "ref_count": int, "alt_count": int,
        "other_alt_count": int, "allele_frequency": float,
    })


def write_pileup(path, df, sample="Test-S"):
    """Write a GATK GetPileupSummaries-style pileup with a metadata header."""
    path = pathlib.Path(path)
    with open(path, "w") as fh:
        fh.write(f"#<METADATA>SAMPLE={sample}\n")
        df.to_csv(fh, sep="\t", index=False)
    return str(path)


def write_tsv(path, df, header_lines=None):
    """Write a TSV optionally prefixed with raw header lines (e.g. @-comments)."""
    path = pathlib.Path(path)
    with open(path, "w") as fh:
        for line in header_lines or []:
            fh.write(line if line.endswith("\n") else line + "\n")
        df.to_csv(fh, sep="\t", index=False)
    return str(path)


def vcf_file(path, records, contigs=("1", "2", "X")):
    """Write a minimal VCFv4.2.

    records: iterable of (chrom, pos, ref, alt) or
             (chrom, pos, ref, alt, fmt, sample) tuples.
    """
    path = pathlib.Path(path)
    lines = ["##fileformat=VCFv4.2"]
    for c in contigs:
        lines.append(f"##contig=<ID={c},length=1000>")
    has_gt = any(len(r) > 4 for r in records)
    header = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
    if has_gt:
        header += ["FORMAT", "SAMPLE"]
    lines.append("\t".join(header))
    for r in records:
        chrom, pos, ref, alt = r[0], r[1], r[2], r[3]
        row = [str(chrom), str(pos), ".", ref, alt, ".", ".", "."]
        if has_gt:
            fmt = r[4] if len(r) > 4 else "GT"
            sample = r[5] if len(r) > 5 else "0/1"
            row += [fmt, sample]
        lines.append("\t".join(row))
    path.write_text("\n".join(lines) + "\n")
    return str(path)


# --------------------------------------------------------------------------- #
# Frame comparison
# --------------------------------------------------------------------------- #
def assert_frame_close(got, exp, float_cols=None, atol=1e-3, rtol=1e-5,
                       check_row_order=True):
    """Compare two DataFrames: exact on non-float columns, tolerant on floats.

    NaNs compare equal. Column set must match; row order optional.
    """
    got = got.reset_index(drop=True).copy()
    exp = exp.reset_index(drop=True).copy()
    assert set(got.columns) == set(exp.columns), (
        f"column mismatch: got {sorted(got.columns)} vs exp {sorted(exp.columns)}"
    )
    got = got[list(exp.columns)]
    if not check_row_order:
        got = got.sort_values(list(exp.columns)).reset_index(drop=True)
        exp = exp.sort_values(list(exp.columns)).reset_index(drop=True)
    assert len(got) == len(exp), f"row count {len(got)} != {len(exp)}"

    float_cols = set(float_cols or [])
    for col in exp.columns:
        g, e = got[col], exp[col]
        if col in float_cols or (e.dtype.kind == "f" and col not in float_cols):
            np.testing.assert_allclose(
                g.to_numpy(dtype=float), e.to_numpy(dtype=float),
                atol=atol, rtol=rtol, equal_nan=True,
                err_msg=f"float column {col!r} differs",
            )
        else:
            mismatch = ~((g == e) | (g.isna() & e.isna()))
            assert not mismatch.any(), (
                f"column {col!r} differs at rows "
                f"{list(np.flatnonzero(mismatch.to_numpy()))}"
            )


# --------------------------------------------------------------------------- #
# Golden-file harness (regression tier)
# --------------------------------------------------------------------------- #
def golden_check(update, golden_path, df, float_cols=None, **close_kw):
    """Compare ``df`` against a golden TSV, or (re)write it under --update-golden.

    Returns True if a comparison was performed, False if the baseline was just
    written. Raises pytest.skip via the caller's guard when missing and not
    updating (handled by ``require_golden``).
    """
    golden_path = pathlib.Path(golden_path)
    if update:
        golden_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(golden_path, sep="\t", index=False)
        return False
    exp = pd.read_csv(golden_path, sep="\t", low_memory=False)
    assert_frame_close(df, exp, float_cols=float_cols, **close_kw)
    return True


# --------------------------------------------------------------------------- #
# CLI runners
# --------------------------------------------------------------------------- #
def run_script(script, argv, expect_returncode=0):
    """Run ``python <python/script>`` as a subprocess (true E2E / __main__ wiring)."""
    cmd = [sys.executable, str(PYDIR / script), *map(str, argv)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if expect_returncode is not None:
        assert proc.returncode == expect_returncode, (
            f"{script} exited {proc.returncode}\n"
            f"STDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
        )
    return proc


# --------------------------------------------------------------------------- #
# Synthetic-cohort access + script-chain runners (shared by the per-script tests)
# --------------------------------------------------------------------------- #
def synth_runs(synth, patient, sample):
    return synth.files["patients"][patient]["samples"][sample]["runs"]


def synth_wes_pileup(synth, patient, sample):
    return next(r["pileup"] for r in synth_runs(synth, patient, sample)
               if r["pileup"] and r["assay"] == "WES")


def synth_scalar(synth, patient, sample, col):
    return synth.truth["patient"].query(
        "patient_id==@patient and sample_id==@sample")[col].iloc[0]


def synth_truth_segs(synth, patient, sample):
    return synth.truth["segments"].query("patient_id==@patient and sample_id==@sample")


def gt_class(gt):
    """Collapse a VCF GT string to {hom_ref, het, hom_alt}."""
    gt = str(gt).split(":", 1)[0]
    a = set(str(gt).replace("|", "/").split("/"))
    return "hom_ref" if a == {"0"} else "hom_alt" if a == {"1"} else "het"


def vcf_sample_field(row, sample_col, field="GT"):
    """Extract one FORMAT field from a VCF sample column in a pandas row."""
    fmt = str(row["FORMAT"]).split(":")
    values = str(row[sample_col]).split(":")
    try:
        idx = fmt.index(field)
    except ValueError:
        idx = 0
    return values[idx] if idx < len(values) else values[0]


def run_genotype(synth, out, patient, sample_names, normal=None):
    """Run genotype.py over the given samples; return the germline VCF as a DataFrame."""
    out = pathlib.Path(out)
    P = synth.files["patients"][patient]["samples"]
    argv = ["--output_dir", str(out), "--patient", patient,
            "--variant", synth.files["patients"][patient]["reference_vcf"]]
    for s in sample_names:
        argv += ["--sample", s, "--pileup", synth_wes_pileup(synth, patient, s),
                 "--segments", P[s]["segments"], "--contamination", P[s]["contamination"]]
    if normal:
        argv += ["--normal_sample", normal]
    argv += ["--threads", "1"]
    run_script("genotype.py", argv)
    vcf = pd.read_csv(out / f"{patient}.germline.vcf", sep="\t", comment="#", header=None,
                      names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                             "INFO", "FORMAT", patient])
    if not vcf.empty:
        vcf[patient] = vcf.apply(lambda row: vcf_sample_field(row, patient, "GT"), axis=1)
    return vcf


def run_map_to_absolute(synth, out, patient, sample, mode="allelic", somix=False,
                        use_absolute=False):
    """Run acs_conversion -> map_to_absolute; return the completed segtab DataFrame."""
    out = pathlib.Path(out)
    P = synth.files["patients"][patient]["samples"][sample]
    purity = float(synth_scalar(synth, patient, sample, "purity"))
    ploidy = float(synth_scalar(synth, patient, sample, "ploidy"))
    sex = str(synth_scalar(synth, patient, sample, "sex"))
    out.mkdir(parents=True, exist_ok=True)
    if somix:
        cr_arg = ["--somix_cr_seg", P["somix_seg"]]
    else:
        run_script("acs_conversion.py", ["--output_dir", str(out), "--seg",
                   P["model_seg"], "--af_parameters", P["af_param"], "--sex", sex])
        cr_arg = ["--acs_cr_seg", str(out / f"{sample}.modelFinal.acs.seg")]
    if use_absolute:
        cr_arg += ["--absolute_segtab", P["absolute_segtab"]]
    run_script("map_to_absolute_copy_number.py", [
        "--outdir", str(out), "--purity", str(purity), "--ploidy", str(ploidy),
        "--sample", sample, "--sex", sex, *cr_arg, "--copy_num_type", mode])
    return pd.read_csv(out / f"{sample}.segtab.{mode}.completed.txt", sep="\t")


def run_ccf(synth, out, patient, sample):
    """Run the acs -> map_to_absolute -> calculate_ccf chain; return the ABS_MAF."""
    out = pathlib.Path(out)
    run_map_to_absolute(synth, out, patient, sample, use_absolute=True)
    segtab = out / f"{sample}.segtab.allelic.completed.txt"
    P = synth.files["patients"][patient]["samples"][sample]
    purity = float(synth_scalar(synth, patient, sample, "purity"))
    ploidy = float(synth_scalar(synth, patient, sample, "ploidy"))
    cdir = out / "ccf"
    cdir.mkdir(parents=True, exist_ok=True)
    run_script("calculate_cancer_cell_fraction.py", [
        "--outdir", str(cdir), "--purity", str(purity), "--ploidy", str(ploidy),
        "--sex", str(synth_scalar(synth, patient, sample, "sex")),
        "--absolute_segtab", str(segtab), "--snv_maf", P["snv_maf"],
        "--indel_maf", P["indel_maf"]])
    return pd.read_csv(cdir / f"{sample}.ABS_MAF.allelic.completed.txt", sep="\t")
