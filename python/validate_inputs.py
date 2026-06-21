"""Pre-flight validation of MultiSampleSomaticWorkflow inputs.

Run as a gate task (ParseInput) before any expensive work. It collects ALL
violations and reports them together, then exits non-zero if any are errors, so
a user fixes one submission instead of discovering problems one failed task at a
time. WDL `version development` has no assert/fail builtin, which is why this
lives in a task.

Three classes of check:
  * shape / concordance   - every provided per-run array has length == #bams;
                            normal sample names are a subset of all sample names.
  * logical requirements  - implicit downstream needs driven by the run_* flags,
                            e.g. total-copy-ratio denoising needs either a (non
                            0-byte) CNV panel of normals or annotated intervals.
  * deep bam integrity    - opt-in (localizes bams): samtools quickcheck, plus
                            each bam's @SQ contig order must match the reference
                            .dict (a mismatch reliably trips up GATK downstream).

Everything that is merely sub-optimal (missing-but-optional resources) is a
warning, not an error.
"""

import argparse
import os
import subprocess
import sys


def read_lines(path):
    if not path:
        return None  # distinguishes "not provided" from "provided empty"
    with open(path) as f:
        return [line.rstrip("\n") for line in f if line.rstrip("\n") != ""]


def parse_bool(s):
    return str(s).lower() == "true"


def sq_contigs(header_lines):
    """Ordered SN values from @SQ lines of a SAM/dict header."""
    contigs = []
    for line in header_lines:
        if not line.startswith("@SQ"):
            continue
        for field in line.split("\t"):
            if field.startswith("SN:"):
                contigs.append(field[3:])
    return contigs


def is_ordered_subsequence(sub, full):
    """True if `sub` appears within `full` in the same relative order."""
    it = iter(full)
    return all(contig in it for contig in sub)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--n_bams", type=int, required=True)
    p.add_argument("--n_bais", type=int, required=True)

    # per-run arrays (a file of one value per line; absent flag == not provided)
    for name in ("sample_names", "normal_sample_names", "is_paired_end",
                 "use_for_tCR", "use_for_aCR", "timepoints",
                 "target_intervals", "annotated_target_intervals",
                 "cnv_panel_of_normals"):
        p.add_argument("--" + name + "_file")

    # resource presence (booleans; we never localize the resources just to test
    # for existence)
    for name in ("common_germline_alleles", "common_germline_alleles_idx",
                 "realignment_image", "germline_resource", "snv_panel_of_normals",
                 "germline_resource_v4_1", "snv_panel_of_normals_v4_1",
                 "funcotator_sources", "sex"):
        p.add_argument("--has_" + name, default="false")

    # run_* flags that imply downstream requirements
    for name in ("run_collect_total_read_counts", "run_model_segments",
                 "run_collect_allelic_read_counts", "run_contamination_model",
                 "run_realignment_filter", "run_variant_calling_mutect1",
                 "run_variant_calling", "run_variant_annotation",
                 "run_clonal_decomposition"):
        p.add_argument("--" + name, default="false")

    # deep tier
    p.add_argument("--deep", action="store_true")
    p.add_argument("--bams_file")
    p.add_argument("--ref_dict")

    args = p.parse_args()

    errors = []
    warnings = []
    n = args.n_bams

    has = {k[4:]: parse_bool(v) for k, v in vars(args).items() if k.startswith("has_")}
    run = {k: parse_bool(v) for k, v in vars(args).items() if k.startswith("run_")}

    # ---- shape / concordance ------------------------------------------------
    if n == 0:
        errors.append("no bams provided; at least one is required.")
    if args.n_bais != n:
        errors.append(
            "bais has {} entries but there are {} bams.".format(args.n_bais, n))

    arrays = {name: read_lines(getattr(args, name + "_file"))
              for name in ("sample_names", "normal_sample_names", "is_paired_end",
                           "use_for_tCR", "use_for_aCR", "timepoints",
                           "target_intervals", "annotated_target_intervals",
                           "cnv_panel_of_normals")}

    for name, vals in arrays.items():
        if name == "normal_sample_names":
            continue  # length is independent of #bams (subset of samples)
        if vals is not None and len(vals) != n:
            errors.append(
                "{} has {} entries but there are {} bams; every provided per-run "
                "array must align 1:1 with the bams.".format(name, len(vals), n))

    sample_names = arrays["sample_names"]
    normals = arrays["normal_sample_names"]
    if sample_names is not None and normals is not None:
        unknown = sorted(set(normals) - set(sample_names))
        if unknown:
            errors.append(
                "normal_sample_names not found among sample_names: {}. A typo here "
                "silently demotes a normal to a tumor.".format(", ".join(unknown)))

    # ---- logical requirements (driven by run_* flags) -----------------------
    # Which runs feed total copy ratio: explicit use_for_tCR, else all of them.
    tcr_flags = arrays["use_for_tCR"]
    annotated = arrays["annotated_target_intervals"]
    cnv_pon = arrays["cnv_panel_of_normals"]
    if run["run_collect_total_read_counts"] or run["run_model_segments"]:
        for i in range(n):
            wants_tcr = (tcr_flags is None) or (i < len(tcr_flags) and parse_bool(tcr_flags[i]))
            if not wants_tcr:
                continue
            pon_ok = cnv_pon is not None and i < len(cnv_pon) and os.path.getsize(cnv_pon[i]) > 0
            annot_ok = annotated is not None and i < len(annotated) and annotated[i] != ""
            if not (pon_ok or annot_ok):
                errors.append(
                    "run {} is used for total copy ratio but has neither a (non "
                    "0-byte) CNV panel of normals nor annotated target intervals; "
                    "DenoiseReadCounts needs one of them.".format(i))

    if run["run_collect_allelic_read_counts"] or run["run_contamination_model"]:
        if not has["common_germline_alleles"]:
            errors.append(
                "allelic counts / contamination are enabled but no "
                "common_germline_alleles resource was provided.")
        elif not has["common_germline_alleles_idx"]:
            errors.append("common_germline_alleles is provided without its index (.idx/.tbi).")

    if run["run_realignment_filter"] and not has["realignment_image"]:
        errors.append(
            "run_realignment_filter is on but no realignment_bwa_mem_index_image "
            "was provided.")

    if run["run_clonal_decomposition"] and not run["run_model_segments"]:
        warnings.append(
            "run_clonal_decomposition is on while run_model_segments is off; this "
            "only works if the samples already carry cached copy-ratio "
            "segmentation and af model parameters (input_patient).")

    # ---- soft warnings (sub-optimal, not fatal) -----------------------------
    if run["run_variant_calling_mutect1"] and not (
            has["germline_resource_v4_1"] and has["snv_panel_of_normals_v4_1"]):
        warnings.append(
            "Mutect1 runs without a v4.1 germline resource and/or SNV panel of "
            "normals; calling proceeds and Mutect2 cleans up, but supplying them "
            "improves Mutect1.")
    if run["run_variant_calling"] and not (
            has["germline_resource"] and has["snv_panel_of_normals"]):
        warnings.append(
            "Mutect2 runs without a germline resource and/or SNV panel of normals; "
            "this is possible but discouraged.")
    if run["run_variant_annotation"] and not has["funcotator_sources"]:
        warnings.append(
            "no Funcotator data sources provided; they will be downloaded from the "
            "GATK bundle at runtime, which is slow.")
    if not has["sex"]:
        warnings.append("sex is not set; allosome copy-number handling falls back to a default.")

    use_aCR = arrays["use_for_aCR"]
    for i in range(n):
        t = (tcr_flags is None) or (i < len(tcr_flags) and parse_bool(tcr_flags[i]))
        a = (use_aCR is None) or (i < len(use_aCR) and parse_bool(use_aCR[i]))
        if not t and not a:
            warnings.append("run {} is marked for neither tCR nor aCR; it will not "
                            "contribute to copy ratios.".format(i))

    # ---- deep bam integrity (opt-in) ----------------------------------------
    if args.deep:
        bams = read_lines(args.bams_file) or []
        for bam in bams:
            rc = subprocess.run(["samtools", "quickcheck", bam]).returncode
            if rc != 0:
                errors.append("samtools quickcheck failed for {} (truncated or "
                              "corrupt bam).".format(bam))
        if args.ref_dict and bams:
            with open(args.ref_dict) as f:
                ref_contigs = sq_contigs(f.readlines())
            ref_set = set(ref_contigs)
            for bam in bams:
                hdr = subprocess.run(["samtools", "view", "-H", bam],
                                     capture_output=True, text=True)
                if hdr.returncode != 0:
                    errors.append("could not read header of {}.".format(bam))
                    continue
                bam_contigs = sq_contigs(hdr.stdout.splitlines())
                missing = [c for c in bam_contigs if c not in ref_set]
                if missing:
                    errors.append("{} has contigs absent from the reference .dict: "
                                  "{}.".format(bam, ", ".join(missing[:5])))
                elif not is_ordered_subsequence(bam_contigs, ref_contigs):
                    errors.append(
                        "{} contig order does not match the reference .dict; this "
                        "trips up GATK. Set run_reorder_bam_contigs = true.".format(bam))

    # ---- report -------------------------------------------------------------
    for w in warnings:
        print("WARNING: " + w)
    for e in errors:
        print("ERROR: " + e)

    if errors:
        print("\nInput validation failed with {} error(s).".format(len(errors)))
        sys.exit(1)
    print("\nInput validation passed ({} warning(s)).".format(len(warnings)))


if __name__ == "__main__":
    main()
