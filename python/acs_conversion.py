import argparse
import numpy as np
import os
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        prog="ModelSegmentsToAllelicCapSegConversion",
        description="""
            Convert GATK ModelSegments output to AllelicCapSeg output. 
        """,
        epilog="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.usage = "acs_conversion.py --seg <seg> --af_parameters <af_parameters> --output_dir <output_dir> [--min_hets <min_hets>] [--min_probes <min_probes>] [--maf90_threshold <maf90_threshold>] [--verbose]"
    parser.add_argument("--seg",            type=str,   required=True,  help="Path to the GATK ModelSegments modelFinal.seg output file.")
    parser.add_argument("--af_parameters",  type=str,   required=True,  help="Path to the GATK ModelSegments modelFinal.af.param output file.")
    parser.add_argument("--output_dir",     type=str,   required=True,  help="Path to the output directory.")
    parser.add_argument("--min_hets",       type=int,   default=0,     help="Minimum number of heterozygous sites for AllelicCapSeg to call a segment.")
    parser.add_argument("--min_probes",     type=int,   default=0,      help="Minimum number of target intervals for AllelicCapSeg to call a segment.")
    parser.add_argument("--maf90_threshold",type=float, default=0.485,  help="Threshold of 90% quantile for setting minor allele fraction to 0.5.")
    parser.add_argument("--verbose",        default=False,  action="store_true", help="Print information to stdout during execution.")
    return parser.parse_args()


def main():
    args = parse_args()
    print_args(args)
    convert_model_segments_to_alleliccapseg(args=args)


def print_args(args):
    if args.verbose:
        print("Calling ModelSegmentsToAllelicCapSegConversion")
        print("Arguments:")
        for key, value in vars(args).items():
            print(f"  {key}: {value}")
        print()


def convert_model_segments_to_alleliccapseg(args):
    """ Adapted from https://portal.firecloud.org/?return=terra#methods/lichtens/Model_Segments_PostProcessing/3
    """

    # get the input file names
    model_seg = args.seg
    af_param = args.af_parameters

    # get the output file names
    prefix = os.path.basename(model_seg).removesuffix(".seg").removesuffix(".modelFinal").removesuffix(".modelBegin")
    output_filename = os.path.join(args.output_dir, f"{prefix}.acs.seg")
    output_skew_filename = output_filename + ".skew"

    # read the input files
    model_segments_seg_pd = pd.read_csv(model_seg, sep='\t', comment='@', na_values='NA')
    model_segments_af_param_pd = pd.read_csv(af_param, sep='\t', comment='@')

    # define AllelicCapSeg columns
    alleliccapseg_seg_columns = [
        'Chromosome',
        'Start.bp',
        'End.bp',
        'n_probes',
        'length',
        'n_hets',
        'f',
        'tau',
        'sigma.tau',
        'mu.minor',
        'sigma.minor',
        'mu.major',
        'sigma.major',
        'SegLabelCNLOH'
    ]

    def simple_determine_allelic_fraction(model_segments_seg_pd):
        result = model_segments_seg_pd['MINOR_ALLELE_FRACTION_POSTERIOR_50'].copy()
        result.loc[model_segments_seg_pd['MINOR_ALLELE_FRACTION_POSTERIOR_90'] > args.maf90_threshold] = 0.5
        return result

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
            if data['Chromosome'][i] != data['Chromosome'][i + 1]:
                continue
            if pd.isna(data['f'][i]) or pd.isna(data['f'][i + 1]):
                continue
            if pass_t(data['tau'], data['sigma.tau'], i, i + 1) and not pass_f(data['f'], i, i + 1):
                if data['f'][i] > data['f'][i + 1]:
                    labels[i + 1] = 1
                else:
                    labels[i] = 1

        # Pass 2: Assign label 0 based on stricter CNLOH conditions for triplets
        for j in range(n_seg - 2):
            if (data['Chromosome'][j] != data['Chromosome'][j + 1] or
                data['Chromosome'][j + 1] != data['Chromosome'][j + 2]
            ):
                continue
            if pd.isna(data['f'][j]) or pd.isna(data['f'][j + 1]) or pd.isna(data['f'][j + 2]):
                continue
            if (pass_t(data['tau'], data['sigma.tau'], j, j + 2) and
                pass_t(data['tau'], data['sigma.tau'], j, j + 1) and
                pass_t(data['tau'], data['sigma.tau'], j + 1, j + 2) and
                not pass_f(data['f'], j, j + 1, 0.1) and
                not pass_f(data['f'], j + 1, j + 2, 0.1) and
                # pass_pair(data['alt'], data['ref'], data['outlier_prob'], j + 1) and
                data['f'][j + 1] < data['f'][j] and
                data['f'][j + 1] < data['f'][j + 2]
            ):
                labels[j + 1] = 0

        return labels

    def convert(model_segments_seg_pd, model_segments_af_param_pd):
        alleliccapseg_seg_pd = pd.DataFrame(columns=alleliccapseg_seg_columns)

        # The following conversions are trivial.
        alleliccapseg_seg_pd['Chromosome'] = model_segments_seg_pd['CONTIG']
        alleliccapseg_seg_pd['Start.bp'] = model_segments_seg_pd['START']
        alleliccapseg_seg_pd['End.bp'] = model_segments_seg_pd['END']
        alleliccapseg_seg_pd['n_probes'] = model_segments_seg_pd['NUM_POINTS_COPY_RATIO']
        alleliccapseg_seg_pd['length'] = alleliccapseg_seg_pd['End.bp'] - alleliccapseg_seg_pd['Start.bp']
        alleliccapseg_seg_pd['n_hets'] = model_segments_seg_pd['NUM_POINTS_ALLELE_FRACTION']

        # ModelSegments estimates posterior credible intervals, while AllelicCapSeg performs maximum a posteriori (MAP) estimation.
        # The copy-ratio and allele-fraction models fit by both also differ.
        # We will attempt a rough translation of the model fits here.

        alleliccapseg_seg_pd['f'] = simple_determine_allelic_fraction(model_segments_seg_pd)

        alleliccapseg_seg_pd['tau'] = 2. * 2 ** model_segments_seg_pd['LOG2_COPY_RATIO_POSTERIOR_50']
        alleliccapseg_seg_pd['sigma.tau'] = 2 ** model_segments_seg_pd['LOG2_COPY_RATIO_POSTERIOR_90'] - 2 ** model_segments_seg_pd['LOG2_COPY_RATIO_POSTERIOR_10']
        sigma_f = (model_segments_seg_pd['MINOR_ALLELE_FRACTION_POSTERIOR_90'].to_numpy() - model_segments_seg_pd['MINOR_ALLELE_FRACTION_POSTERIOR_10'].to_numpy()) / 2.
        sigma_mu = np.sqrt(sigma_f ** 2 + alleliccapseg_seg_pd['sigma.tau'] ** 2)  # we propagate errors in the products f * tau and (1 - f) * tau in the usual way
        alleliccapseg_seg_pd['mu.minor'] = alleliccapseg_seg_pd['f'] * alleliccapseg_seg_pd['tau']
        alleliccapseg_seg_pd['sigma.minor'] = sigma_mu
        alleliccapseg_seg_pd['mu.major'] = (1. - alleliccapseg_seg_pd['f']) * alleliccapseg_seg_pd['tau']
        alleliccapseg_seg_pd['sigma.major'] = sigma_mu

        # AllelicCapSeg attempts to call CNLOH.  It attempts to distinguish between
        # three states ("0 is flanked on both sides, 1 is one side, 2 is no cn.loh").
        alleliccapseg_seg_pd['SegLabelCNLOH'] = label_cnloh(alleliccapseg_seg_pd)

        # One important caveat: for segments with less than 10 hets, AllelicCapSeg also tries to call whether a segment is "split" or not.
        # This script will attempt to call "split" on all segments.
        # ACS performs a simple hypothesis test on the alternate-allele fractions to see if
        # a unimodal distribution peaked at 0.5 is supported over a bimodal distribution peaked at f and 1 - f.
        # If the former is supported, then AllelicCapSeg ignores the MAP estimate of f and simply sets it to be 0.5.
        # ABSOLUTE may actually be rather sensitive to this.  Again, let's ignore for now, and we can later port this
        # statistical test if necessary.

        # Finally, I believe that ABSOLUTE requires the value of the "skew" parameter from the AllelicCapSeg
        # allele-fraction model.  This parameter is supposed to allow the model to account for reference bias,
        # but the model likelihood that AllelicCapSeg uses is not valid over the entire range of the skew parameter.
        # We corrected this during the development of AllelicCNV and retain the same corrected model in ModelSegments.
        # We will try to transform the relevant parameter in the corrected model back to a "skew",
        # but this operation is ill-defined.  Luckily, for WGS, the reference bias is typically negligible.
        model_segments_reference_bias = model_segments_af_param_pd[model_segments_af_param_pd['PARAMETER_NAME'] == 'MEAN_BIAS']['POSTERIOR_50']
        alleliccapseg_skew = 2. / (1. + model_segments_reference_bias)

        # If a row has less than X (set by user) hets or number of target intervals (probes), remove:
        good_rows = alleliccapseg_seg_pd["n_hets"] >= args.min_hets
        good_rows &= alleliccapseg_seg_pd["n_probes"] >= args.min_probes
        n = alleliccapseg_seg_pd.shape[0] - np.sum(good_rows)
        pct_drop = n / alleliccapseg_seg_pd.shape[0] * 100
        print(f"Dropping {n} (-{pct_drop:.3f}%) segments with min_hets < {args.min_hets} or min_probes < {args.min_probes}.")
        alleliccapseg_seg_pd = alleliccapseg_seg_pd.loc[good_rows]

        return alleliccapseg_seg_pd, alleliccapseg_skew

    alleliccapseg_seg_pd, alleliccapseg_skew = convert(model_segments_seg_pd, model_segments_af_param_pd)

    alleliccapseg_seg_pd.to_csv(output_filename, sep='\t', index=False, na_rep='NaN')
    np.savetxt(output_skew_filename, alleliccapseg_skew)


if __name__ == "__main__":
    main()
