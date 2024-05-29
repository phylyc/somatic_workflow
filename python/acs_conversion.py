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
    parser.usage = "acs_conversion.py --seg <seg> --af_parameters <af_parameters> --output_dir <output_dir> [--min_hets <min_hets>] [--maf90_threshold <maf90_threshold>] [--verbose]"
    parser.add_argument("--seg",            type=str,   required=True,  help="Path to the GATK ModelSegments modelFinal.seg output file.")
    parser.add_argument("--af_parameters",  type=str,   required=True,  help="Path to the GATK ModelSegments modelFinal.af.param output file.")
    parser.add_argument("--output_dir",     type=str,   required=True,  help="Path to the output directory.")
    parser.add_argument("--min_hets",       type=int,   default=10,     help="Minimum number of heterozygous sites for AllelicCapSeg to call a segment.")
    parser.add_argument("--min_probes",     type=int,   default=4,      help="Minimum number of heterozygous sites for AllelicCapSeg to call a segment.")
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
        result = model_segments_seg_pd['MINOR_ALLELE_FRACTION_POSTERIOR_50']
        result[model_segments_seg_pd['MINOR_ALLELE_FRACTION_POSTERIOR_90'] > args.maf90_threshold] = 0.5
        return result

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

        # For whatever reason, AllelicCapSeg attempts to call CNLOH.  Documentation is spotty, but it seems like it attempts
        # to distinguish between three states ("0 is flanked on both sides, 1 is one side, 2 is no cn.loh").
        # Let's just set everything to 2 for now.
        # Hopefully, ABSOLUTE is robust to this ...
        alleliccapseg_seg_pd['SegLabelCNLOH'] = 2

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

        # If a row has less than X (set by user) hets or number of target intervals (probes), then assume zero
        filter_rows = alleliccapseg_seg_pd['n_hets'] < args.min_hets
        filter_rows |= alleliccapseg_seg_pd["n_probes"] < args.min_probes
        # mu.minor  sigma.minor  mu.major  sigma.major
        alleliccapseg_seg_pd.loc[filter_rows, 'n_hets'] = 0
        alleliccapseg_seg_pd.loc[filter_rows, 'f'] = np.NaN
        alleliccapseg_seg_pd.loc[filter_rows, 'mu.minor'] = np.NaN
        alleliccapseg_seg_pd.loc[filter_rows, 'sigma.minor'] = np.NaN
        alleliccapseg_seg_pd.loc[filter_rows, 'mu.major'] = np.NaN
        alleliccapseg_seg_pd.loc[filter_rows, 'sigma.major'] = np.NaN

        return alleliccapseg_seg_pd, alleliccapseg_skew

    alleliccapseg_seg_pd, alleliccapseg_skew = convert(model_segments_seg_pd, model_segments_af_param_pd)

    # ABSOLUTE doesn't like NaN segments
    alleliccapseg_seg_pd.dropna(inplace=True)

    alleliccapseg_seg_pd.to_csv(output_filename, sep='\t', index=False, na_rep='NaN')
    np.savetxt(output_skew_filename, alleliccapseg_skew)


if __name__ == "__main__":
    main()
