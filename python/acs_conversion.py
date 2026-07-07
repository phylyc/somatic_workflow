import argparse
import numpy as np
import os
import pandas as pd
import time


def message(*args, **kwargs) -> None:
    print(f"{time.strftime('%H:%M:%S')} ", *args, **kwargs)
    return None


def parse_args():
    parser = argparse.ArgumentParser(
        prog="ModelSegmentsToAllelicCapSegConversion",
        description="""
            Convert a GATK ModelSegments segmentation (modelFinal.seg) or a somix
            segmentation into AllelicCapSeg (ACS) format for ABSOLUTE. Both are supplied
            via --seg and told apart automatically by their columns. Per segment it emits
            total copy ratio scaled by chromosomal ploidy
            (ModelSegments: tau = chr_ploidy * 2^LOG2_COPY_RATIO_POSTERIOR_50; somix:
            tau = chr_ploidy * exp(log_tCR)), the minor allele fraction f in [0, 0.5]
            (ModelSegments: MINOR_ALLELE_FRACTION_POSTERIOR_50, snapped to 0.5 when
            MINOR_ALLELE_FRACTION_POSTERIOR_90 > --maf90_threshold; somix: f_MAP), the
            allelic copy ratios mu.minor = f*tau and mu.major = (1-f)*tau with propagated
            sigmas, and a CNLOH label. Sex is handled per karyotype: a contig with
            chromosomal ploidy 0 (e.g. chrY in a female) carries no copy-ratio signal and
            its tau/f are NaN; haploid male X/Y are scaled by ploidy 1. Segments below
            --min_hets / --min_probes are dropped. A multi-sample somix file (and its
            --af_parameters) is subset to --sample_id.

            A companion '<prefix>.acs.seg.skew' file is written for ABSOLUTE: skew =
            2 / (1 + MEAN_BIAS) from the allele-fraction parameters when --af_parameters is
            given, otherwise the sentinel -1.0 (meaning "no skew / disabled" in
            total-copy-ratio mode). MEAN_BIAS is the ModelSegments MEAN_BIAS POSTERIOR_50,
            or, for somix af parameters (columns sample_id, s_MAP, lambda_MAP), the
            reference-bias multiplier lambda_MAP of the row matching --sample_id.
        """,
        epilog="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.usage = "acs_conversion.py --output_dir <output_dir> --seg <seg> [--af_parameters <af_parameters>] [--sample_id <sample_id>] [--min_hets <min_hets>] [--min_probes <min_probes>] [--maf90_threshold <maf90_threshold>] [--sex <sex>] [--verbose]"
    parser.add_argument("--output_dir",     type=str,   required=True,  help="Path to the output directory.")
    parser.add_argument("--sample_id",      type=str,   help="Sample identifier; selects the matching row from a multi-sample somix segmentation and --af_parameters file.")
    parser.add_argument("--seg",            type=str,   required=True,  help="Path to a segmentation file: either a GATK ModelSegments modelFinal.seg or a somix segmentation (auto-detected from its columns).")
    parser.add_argument("--af_parameters",  type=str,                   help="Path to the GATK ModelSegments modelFinal.af.param output file.")
    parser.add_argument("--min_hets",       type=int,   default=0,      help="Minimum number of heterozygous sites for AllelicCapSeg to call a segment.")
    parser.add_argument("--min_probes",     type=int,   default=0,      help="Minimum number of target intervals for AllelicCapSeg to call a segment.")
    parser.add_argument("--maf90_threshold",type=float, default=0.49,   help="Threshold of 90%% quantile for setting minor allele fraction to 0.5.")
    parser.add_argument("--sex",            type=str,   default="XXY",  help="Genotype sex of the patient for ploidy priors on X and Y chromosomes: {Female, Male, female, male, XX, XY, XXY, XYY, XXX, etc.}")
    parser.add_argument("--normal_ploidy",  type=int,   required=False, default=2, help="Normal/germline ploidy of that organism.")
    parser.add_argument("--verbose",        default=False,  action="store_true", help="Print information to stdout during execution.")
    return parser.parse_args()


def main():
    args = parse_args()
    print_args(args)
    convert_model_segments_to_alleliccapseg(args=args)


def print_args(args):
    if args.verbose:
        message("Calling ModelSegmentsToAllelicCapSegConversion")
        print("Arguments:")
        for key, value in vars(args).items():
            print(f"  {key}: {value}")
        print()


def convert_model_segments_to_alleliccapseg(args):
    """ Adapted from https://portal.firecloud.org/?return=terra#methods/lichtens/Model_Segments_PostProcessing/3
    """
    if args.sex in ["Female", "female"]:
        args.sex = "XX"
    if args.sex in ["Male", "male"]:
        args.sex = "XY"

    nX = args.sex.count("X")
    nY = args.sex.count("Y")

    # --seg is overloaded: it accepts either a GATK ModelSegments segmentation or a somix
    # segmentation, told apart by their columns below (like --af_parameters).
    if args.seg is None:
        raise Exception("--seg must be provided.")
    input_seg = args.seg

    # get the output file names
    prefix = os.path.basename(input_seg).removesuffix(".gz").removesuffix(".seg").removesuffix(".tsv") if args.sample_id is None else args.sample_id
    output_filename = os.path.join(args.output_dir, f"{prefix}.acs.seg")
    output_skew_filename = output_filename + ".skew"

    # define AllelicCapSeg columns
    alleliccapseg_seg_columns = [
        "Chromosome",
        "Start.bp",
        "End.bp",
        "n_probes",
        "length",
        "n_hets",
        "f",
        "tau",
        "sigma.tau",
        "mu.minor",
        "sigma.minor",
        "mu.major",
        "sigma.major",
        "SegLabelCNLOH"
    ]

    # Main function to assign CNLOH (copy-neutral loss of heterozygosity) labels
    def label_cnloh(data):
        """Annotates CNLOH labels on data based on conditions defined in the LabelSegsCNLOH logic from AllelicCapSeg."""
        n_seg = len(data)
        labels = np.full(n_seg, 2)  # Start with label 2 (no CNLOH) for all segments

        # Define helper functions for conditions
        def pass_t(tau, sigma, ix1, ix2, tau_thresh=0.08):
            """Checks if the difference in tau between two segments is below the threshold."""
            return ((tau[ix1] - tau[ix2]) ** 2 / (sigma[ix1] ** 2 + sigma[ix2] ** 2)) < tau_thresh

        def pass_f(f, ix1, ix2, f_thresh=0.01):
            """Checks if the difference in allele frequency (f) between two segments is below the threshold."""
            return abs(f[ix1] - f[ix2]) < f_thresh

        ## Unused check since we don't have access to alt or ref counts or outlier probabilities anymore.
        # def pass_pair(alt, ref, outlier_prob, ix):
        #     """Checks for heterozygosity balance in the segment."""
        #     cov = alt + ref
        #     het_out = outlier_prob[ix]
        #     out_ix = het_out > 0.5  # Consider SNPs with outlier probability > 0.5
        #     return (np.any((alt[~out_ix] / cov[~out_ix]) > 0.5) and
        #             np.any((alt[~out_ix] / cov[~out_ix]) < 0.5))

        # Pass 1: Assign label 1 based on CNLOH conditions between adjacent segments
        for i in range(n_seg - 1):
            if data["Chromosome"][i] != data["Chromosome"][i + 1]:
                continue
            if pd.isna(data["f"][i]) or pd.isna(data["f"][i + 1]):
                continue
            if pass_t(data["tau"], data["sigma.tau"], i, i + 1) and not pass_f(data["f"], i, i + 1):
                if data["f"][i] > data["f"][i + 1]:
                    labels[i + 1] = 1
                else:
                    labels[i] = 1

        # Pass 2: Assign label 0 based on stricter CNLOH conditions for triplets
        for j in range(n_seg - 2):
            if (data["Chromosome"][j] != data["Chromosome"][j + 1] or
                data["Chromosome"][j + 1] != data["Chromosome"][j + 2]
            ):
                continue
            if pd.isna(data["f"][j]) or pd.isna(data["f"][j + 1]) or pd.isna(data["f"][j + 2]):
                continue
            if (pass_t(data["tau"], data["sigma.tau"], j, j + 2) and
                pass_t(data["tau"], data["sigma.tau"], j, j + 1) and
                pass_t(data["tau"], data["sigma.tau"], j + 1, j + 2) and
                not pass_f(data["f"], j, j + 1, 0.1) and
                not pass_f(data["f"], j + 1, j + 2, 0.1) and
                # pass_pair(data["alt"], data["ref"], data["outlier_prob"], j + 1) and
                data["f"][j + 1] < data["f"][j] and
                data["f"][j + 1] < data["f"][j + 2]
            ):
                labels[j + 1] = 0

        return labels

    alleliccapseg_seg_pd = pd.DataFrame(columns=alleliccapseg_seg_columns)

    def chromosomal_ploidy(chromosome: pd.Series) -> pd.Series:
        return chromosome.apply(lambda chr: nX if chr in ["X", "chrX"] else nY if chr in ["Y", "chrY"] else args.normal_ploidy)

    # Detect the --seg format from its columns: ModelSegments carries
    # LOG2_COPY_RATIO_POSTERIOR_50, somix carries log_tCR.
    seg_pd = pd.read_csv(input_seg, sep="\t", comment="@", na_values="NA", low_memory=False)
    is_somix = "log_tCR" in seg_pd.columns
    if not is_somix and "LOG2_COPY_RATIO_POSTERIOR_50" not in seg_pd.columns:
        raise Exception("Unknown --seg input file format; expected GATK ModelSegments or somix columns.")

    if is_somix:
        # ------------------------------------------------------------------- #
        # somix segmentation input
        # ------------------------------------------------------------------- #
        somix_seg_pd = seg_pd
        if args.sample_id is not None and "sample_id" in somix_seg_pd.columns:
            somix_seg_pd = somix_seg_pd.loc[somix_seg_pd["sample_id"] == args.sample_id]
        somix_seg_pd = somix_seg_pd.reset_index(drop=True)

        alleliccapseg_seg_pd["Chromosome"] = somix_seg_pd["contig"]
        alleliccapseg_seg_pd["Start.bp"] = somix_seg_pd["Start.bp"]
        alleliccapseg_seg_pd["End.bp"] = somix_seg_pd["End.bp"]
        alleliccapseg_seg_pd["n_probes"] = somix_seg_pd["n_markers"]
        alleliccapseg_seg_pd["length"] = alleliccapseg_seg_pd["End.bp"] - alleliccapseg_seg_pd["Start.bp"]
        alleliccapseg_seg_pd["n_hets"] = somix_seg_pd["n_snps"]

        chr_ploidy = chromosomal_ploidy(alleliccapseg_seg_pd["Chromosome"])

        # somix reports the natural-log total copy ratio (log_tCR) and its standard error
        # (sem_log_tCR). This mirrors the somix parsing in map_to_absolute_copy_number.py
        # so the two paths agree: tau = exp(log_tCR) * chr_ploidy, with log-normal error
        # propagation for sigma.tau.
        alleliccapseg_seg_pd["tau"] = np.exp(somix_seg_pd["log_tCR"]) * chr_ploidy
        alleliccapseg_seg_pd["sigma.tau"] = chr_ploidy * np.exp(somix_seg_pd["log_tCR"] + somix_seg_pd["sem_log_tCR"] ** 2 / 2) * np.sqrt(np.exp(somix_seg_pd["sem_log_tCR"] ** 2) - 1)

        # f_MAP is already the minor allele fraction in [0, 0.5]; its posterior variance
        # comes from the Hessian of the somix allele-fraction model (var_f = -1/f_hessian).
        alleliccapseg_seg_pd["f"] = somix_seg_pd["f_MAP"].copy()
        var_f = (-1.0 / somix_seg_pd["f_hessian"]).to_numpy()

        # somix reports no f quantile, so reconstruct the MINOR_ALLELE_FRACTION_POSTERIOR_90
        # that the ModelSegments MAF90 gate uses from a Gaussian (Laplace) posterior at the
        # MAP: the 90th percentile is f_MAP + z90 * sqrt(var_f), with z90 = 1.2816 (half of
        # the 2.563 = z90 - z10 spread the ModelSegments branch uses below). Snap a segment
        # to balance (f = 0.5) when this upper bound crosses --maf90_threshold, i.e. when f
        # cannot be distinguished from 0.5 -- otherwise a noisy near-balanced segment would
        # be fit as spuriously imbalanced. var_f is clipped at 0 (an ill-defined Hessian
        # gives no usable interval, so it cannot trigger the snap).
        z90 = 2.563 / 2
        f_upper90 = alleliccapseg_seg_pd["f"].to_numpy() + z90 * np.sqrt(np.clip(var_f, 0.0, None))
        alleliccapseg_seg_pd.loc[f_upper90 > args.maf90_threshold, "f"] = 0.5

    else:
        # ------------------------------------------------------------------- #
        # GATK ModelSegments input
        # ------------------------------------------------------------------- #
        model_segments_seg_pd = seg_pd

        alleliccapseg_seg_pd["Chromosome"] = model_segments_seg_pd["CONTIG"]
        alleliccapseg_seg_pd["Start.bp"] = model_segments_seg_pd["START"]
        alleliccapseg_seg_pd["End.bp"] = model_segments_seg_pd["END"]
        alleliccapseg_seg_pd["n_probes"] = model_segments_seg_pd["NUM_POINTS_COPY_RATIO"]
        alleliccapseg_seg_pd["length"] = alleliccapseg_seg_pd["End.bp"] - alleliccapseg_seg_pd["Start.bp"]
        alleliccapseg_seg_pd["n_hets"] = model_segments_seg_pd["NUM_POINTS_ALLELE_FRACTION"]

        # NOTE: ModelSegments estimates posterior credible intervals, while AllelicCapSeg
        # performs maximum a posteriori (MAP) estimation. The copy-ratio and allele-fraction
        # models fit by both also differ.

        chr_ploidy = chromosomal_ploidy(alleliccapseg_seg_pd["Chromosome"])

        alleliccapseg_seg_pd["tau"] = chr_ploidy * 2 ** model_segments_seg_pd["LOG2_COPY_RATIO_POSTERIOR_50"]
        alleliccapseg_seg_pd["f"] = model_segments_seg_pd["MINOR_ALLELE_FRACTION_POSTERIOR_50"].copy()

        # Some segments may not have HETs called due to LOH or approximate LOH.
        # alleliccapseg_seg_pd.loc[alleliccapseg_seg_pd["f"].isna(), "f"] = 1 / alleliccapseg_seg_pd["tau"]
        # alleliccapseg_seg_pd["f"] = alleliccapseg_seg_pd["f"].where(alleliccapseg_seg_pd["f"] < 0.5, 1 - alleliccapseg_seg_pd["f"]).clip(lower=0)

        # For segments with less than 10 hets, AllelicCapSeg also tries to call whether
        # a segment is "split" or not. ACS performs a simple hypothesis test on the
        # alternate-allele fractions to see if a unimodal distribution peaked at 0.5 is
        # supported over a bimodal distribution peaked at f and 1 - f. If the former is
        # supported, then AllelicCapSeg ignores the MAP estimate of f and simply sets it to 0.5.
        # For now, we replace the statistical test with a simple threshold test.
        alleliccapseg_seg_pd.loc[model_segments_seg_pd["MINOR_ALLELE_FRACTION_POSTERIOR_90"] > args.maf90_threshold, "f"] = 0.5

        # Propagate errors: If y = exp(x) and x is normally distributed, then Var(y) = exp(2*x + Var(x)) * (exp(Var(x)) - 1)
        # For a normal distribution, it is: q90 - q10 = 2.563 * sigma
        sigma_log2_tau = (model_segments_seg_pd["LOG2_COPY_RATIO_POSTERIOR_90"].to_numpy() - model_segments_seg_pd["LOG2_COPY_RATIO_POSTERIOR_10"].to_numpy()) / 2.563
        var_log_tau = np.log(2)**2 * sigma_log2_tau**2
        alleliccapseg_seg_pd["sigma.tau"] = chr_ploidy * np.exp(np.log(2) * model_segments_seg_pd["LOG2_COPY_RATIO_POSTERIOR_50"] + var_log_tau / 2) * np.sqrt(np.exp(var_log_tau) - 1)

        sigma_f = (model_segments_seg_pd["MINOR_ALLELE_FRACTION_POSTERIOR_90"].to_numpy() - model_segments_seg_pd["MINOR_ALLELE_FRACTION_POSTERIOR_10"].to_numpy()) / 2.563
        sigma_f = np.where(np.isnan(sigma_f), 1e-3, sigma_f)
        var_f = sigma_f ** 2

    # ----------------------------------------------------------------------- #
    # Shared: sex handling and allelic copy ratios (independent of input type)
    # ----------------------------------------------------------------------- #
    # Correct the diploid assumption
    if nX < 2:
        alleliccapseg_seg_pd.loc[alleliccapseg_seg_pd["Chromosome"].isin(["X", "chrX"]), "f"] = np.nan
    # Assume that Y never has heterozygous germline SNPs (technically not true, but those that arise
    # e.g. during S-phase before meiosis II are not enough to allow for stable aCR signal).
    alleliccapseg_seg_pd.loc[alleliccapseg_seg_pd["Chromosome"].isin(["Y", "chrY"]), "f"] = np.nan

    # A chromosome that is absent (ploidy 0, e.g. Y in an XX sample) has an *unobserved*,
    # not zero, total copy ratio. Set tau (and its error) to NaN so these segments are
    # dropped downstream instead of being fit as spurious homozygous deletions. Haploid
    # X/Y in a male keep their real ploidy=1 tau; only ploidy == 0 is masked. (mu/sigma
    # for the minor/major alleles inherit the NaN since they are derived from tau below.)
    alleliccapseg_seg_pd.loc[chr_ploidy == 0, "tau"] = np.nan
    alleliccapseg_seg_pd.loc[chr_ploidy == 0, "sigma.tau"] = np.nan

    alleliccapseg_seg_pd["mu.minor"] = alleliccapseg_seg_pd["f"] * alleliccapseg_seg_pd["tau"]
    alleliccapseg_seg_pd["mu.major"] = (1 - alleliccapseg_seg_pd["f"]) * alleliccapseg_seg_pd["tau"]
    # If m = f * t then Var(m) = t**2 * Var(f) + f**2 * Var(t) (+ cov)
    alleliccapseg_seg_pd["sigma.minor"] = np.sqrt(alleliccapseg_seg_pd["tau"] ** 2 * var_f + alleliccapseg_seg_pd["f"] ** 2 * alleliccapseg_seg_pd["sigma.tau"] ** 2)
    alleliccapseg_seg_pd["sigma.major"] = np.sqrt(alleliccapseg_seg_pd["tau"] ** 2 * var_f + (1 - alleliccapseg_seg_pd["f"]) ** 2 * alleliccapseg_seg_pd["sigma.tau"] ** 2)

    # AllelicCapSeg attempts to call CNLOH. It attempts to distinguish between
    # three states ("0 is flanked on both sides, 1 is one side, 2 is no cn.loh").
    alleliccapseg_seg_pd["SegLabelCNLOH"] = label_cnloh(alleliccapseg_seg_pd)

    W = alleliccapseg_seg_pd["length"] / alleliccapseg_seg_pd["length"].sum()

    good_rows = alleliccapseg_seg_pd["n_hets"] >= args.min_hets
    good_rows &= alleliccapseg_seg_pd["n_probes"] >= args.min_probes
    n = alleliccapseg_seg_pd.shape[0] - np.sum(good_rows)
    pct_genomic_drop = W.loc[~good_rows].sum() * 100
    print(f"Dropping {n}/{alleliccapseg_seg_pd.shape[0]} (-{pct_genomic_drop:.6f}% of genome) segments with min_hets < {args.min_hets} or min_probes < {args.min_probes}.")

    alleliccapseg_seg_pd = alleliccapseg_seg_pd.loc[good_rows]

    alleliccapseg_seg_pd.to_csv(output_filename, sep="\t", index=False, na_rep="NaN")

    if args.af_parameters is not None:
        # The af_parameters are only informative if allelic counts have been supplied to ModelSegments.
        # If only total copy ratios had been supplied, ModelSegments will still produce af_param output
        # but its estimates are not informed by data but by ModelSegments MCMC over the bias prior range
        # which is [0, 5], so the MEAN_BIAS will likely be close to 2.5.
        af_param = args.af_parameters
        model_segments_af_param_pd = pd.read_csv(af_param, sep="\t", comment="@")
        # ABSOLUTE requires the value of the "skew" parameter from the AllelicCapSeg
        # allele-fraction model. This parameter allows the model to account for reference bias
        # of the form f -> f / (f + (1 - f) * MEAN_BIAS).
        # We try to transform the relevant parameter in the corrected model back to a "skew",
        # but this operation is ill-defined. For WGS, the reference bias is typically negligible.
        if "PARAMETER_NAME" in model_segments_af_param_pd.columns:
            # GATK ModelSegments allele-fraction parameters: MEAN_BIAS is a named row.
            model_segments_reference_bias = model_segments_af_param_pd[model_segments_af_param_pd["PARAMETER_NAME"] == "MEAN_BIAS"]["POSTERIOR_50"]
        elif "lambda_MAP" in model_segments_af_param_pd.columns:
            # somix allele-fraction parameters (columns sample_id, s_MAP, lambda_MAP):
            # the reference bias MEAN_BIAS is the multiplier lambda_MAP. Select the row
            # for this sample; the file may hold one row per sample in a cohort.
            rows = model_segments_af_param_pd
            if "sample_id" in rows.columns and args.sample_id is not None:
                rows = rows[rows["sample_id"] == args.sample_id]
            elif rows.shape[0] > 1:
                raise Exception("Multiple rows in --af_parameters; provide --sample_id to select one.")
            if rows.empty:
                raise Exception(f"No row in --af_parameters matching --sample_id '{args.sample_id}'.")
            model_segments_reference_bias = rows["lambda_MAP"]
        else:
            raise Exception("Unknown --af_parameters input file format.")

        alleliccapseg_skew = 2. / (1. + model_segments_reference_bias)
        np.savetxt(output_skew_filename, alleliccapseg_skew)
    else:
        # tCR mode: no allele-fraction model was supplied, so there is no meaningful
        # skew. Write -1 as a 1-D array (np.savetxt rejects 0-D scalars) as a sentinel
        # for "no skew" (interpreted downstream / by ABSOLUTE as skew = 1 / disabled).
        np.savetxt(output_skew_filename, [-1.0])


if __name__ == "__main__":
    main()
