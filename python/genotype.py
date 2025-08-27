import argparse
from collections import Counter, defaultdict
from copy import copy
from itertools import combinations
import gzip
import multiprocessing as mp
import numpy as np
import os
import pandas as pd
import scipy.stats as st
import scipy.special as sp
import time
import warnings


def message(*args, **kwargs) -> None:
    print(f"{time.strftime('%H:%M:%S')} ", *args, **kwargs)
    return None


def parse_args():
    parser = argparse.ArgumentParser(
        prog="Genotyper",
        description="""
            Genotyper is a script that calculates genotype likelihoods from allelic counts, 
            while accounting for potential sample contamination, reference bias, and  
            segment-specific minor allele fractions. It leverages allelic pileup data 
            (from GATK GetPileupSummaries), allelic segmentation data (from GATK CalculateContamination
            or GATK ModelSegments), and contamination estimates (from GATK CalculateContamination), 
            to compute genotype likelihoods under a beta-binomial model. The script 
            computes variant genotype correlations across samples as a consistency check.

            Outputs include:
                1) A VCF with genotype information for each patient,
                2) Allelic count matrices for each allele,
                3) Sample correlation matrix, and
                4) (Optionally) sample pileups annotated with genotype likelihoods.

            Additional options allow filtering based on allele frequency, read depth, and 
            error rate thresholds, as well as output formatting and parallelization.
        """,
        epilog="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.usage = "genotype.py -O <output_dir> --patient <patient> --sample <sample> [--sample <sample> ...] -P <pileup> [-P <pileup> ...] [-C <contamination> ...] [-S <segments> ...] [-V <variant> ...] [--normal_sample <name>] [--normal_to_tumor_weight <weight>] [-D <min_read_depth>] [--min_allele_frequency <freq>] [-p <min_genotype_likelihood>] [-s <overdispersion>] [-l <ref_bias>] [--min_error_rate <rate>] [--max_error_rate <rate>] [--outlier_prior <prior>] [-F <format>] [--threads <num>] [--select_hets] [--save_sample_genotype_likelihoods] [--compress_output] [--verbose]"
    parser.add_argument("-O", "--output_dir",                   type=str,   required=True,  help="Path to the output directory.")
    parser.add_argument("--patient",                            type=str,   required=True,  help="Name of the patient.")
    parser.add_argument("--sex",                                type=str,   default="XXY",  help="Genotype sex of the patient for ploidy priors on X and Y chromosomes: {Female, Male, female, male, XX, XY, XXY, XYY, XXX, etc.}")
    parser.add_argument("--sample",                             type=str,   required=True,  action="append",    help="Assigned name of the sample. Does not have to coincide with the sample name in the pileup file header.")
    parser.add_argument("-P", "--pileup",                       type=str,   required=True,  action="append",    help="Path to the allelic pileup file (output of GATK's GetPileupSummaries; used for CalculateContamination).")
    parser.add_argument("-C", "--contamination",                type=str,   default=None,   action="append",    help="Path to the contamination estimate file (output of GATK's CalculateContamination).")
    parser.add_argument("-S", "--segments",                     type=str,   default=None,   action="append",    help="Path to the allelic copy ratio segmentation file (output of GATK's CalculateContamination OR ModelSegments).")
    parser.add_argument("-L", "--af_model_parameters",          type=str,   default=None,   action="append",    help="Path to the allelic copy ratio segmentation model parameters containing ref bias (output of GATK's ModelSegments).")
    parser.add_argument("-V", "--variant",                      type=str,   default=None,   action="append",    help="VCF file containing germline variants and population allele frequencies. (Used to get pileup summaries.)")
    parser.add_argument("--normal_sample",                      type=str,   default=None,   action="append",    help="Assigned name of the normal sample.")
    parser.add_argument("--normal_to_tumor_weight",             type=float, default=10.0,   help="Relative weight of normal samples to tumor samples when determining patient genotype.")
    parser.add_argument("-D", "--min_read_depth",               type=int,   default=10,     help="Minimum read depth per sample to consider site for genotyping.")
    parser.add_argument("--min_allele_frequency",               type=float, default=0,      help="Minimum population allele frequency to consider a site.")
    parser.add_argument("-p", "--min_genotype_likelihood",      type=float, default=0.995,  help="Probability threshold for calling and retaining genotypes.")
    parser.add_argument("-s", "--overdispersion",               type=float, default=10,     help="")
    parser.add_argument("-l", "--ref_bias",                     type=float, default=1.05,   help="")
    parser.add_argument("--min_error_rate",                     type=float, default=1e-3,   help="")
    parser.add_argument("--max_error_rate",                     type=float, default=1e-2,   help="")
    parser.add_argument("--outlier_prior",                      type=float, default=1e-3,   help="Prior probability for a variant to be an outlier.")
    parser.add_argument("--phasing_log_ratio_cap",              type=float, default=10,     help="Maximum absolute log odds ratio to avoid outlier HETs dominating the phasing signal.")
    parser.add_argument("--phasing_sample_llr_threshold",       type=float, default=0.4,    help="Sample-level log likelihood ratio threshold of 0|1 vs 1|0 to consider a HET for consensus phasing.")
    parser.add_argument("--phasing_consensus_fdr",              type=float, default=0.005,  help="Site-level FDR threshold for evidence of 0|1 vs 1|0 to consider for consensus phasing.")
    parser.add_argument("--phasing_max_num_contig_segs",        type=float, default=1000,   help="Maximum number of segments per chromosome to consider for forming phasing consensus.")
    parser.add_argument("-F", "--format",                       type=str,   default="GT",   help="VCF format field. (GT: genotype; AD: allele depth; DP: total depth; PL: phred-scaled genotype likelihoods.)")
    parser.add_argument("--threads",                            type=int,   default=1,      help="Number of threads to use for parallelization over samples.")
    parser.add_argument("--select_hets",                                    default=False,  action="store_true", help="Keep only heterozygous sites.")
    parser.add_argument("--save_sample_genotype_likelihoods",               default=False,  action="store_true", help="Save genotype likelihoods to file for each sample.")
    parser.add_argument("--compress_output",                                default=False,  action="store_true", help="Compress output files.")
    parser.add_argument("--verbose",                                        default=True,  action="store_true", help="Print information to stdout during execution.")
    return parser.parse_args()


def main():
    args = parse_args()
    print_args(args)
    genotyper = Genotyper(
        sex=args.sex,
        overdispersion=args.overdispersion,
        ref_bias=args.ref_bias,
        min_genotype_likelihood=args.min_genotype_likelihood,
        min_error_rate=args.min_error_rate,
        max_error_rate=args.max_error_rate,
        outlier_prior=args.outlier_prior,
        select_hets=args.select_hets,
        verbose=args.verbose
    )
    data = GenotypeData(
        individual_id=args.patient,
        samples=args.sample,
        variant=args.variant,
        pileup=args.pileup,
        contamination=args.contamination,
        segments=args.segments,
        af_model_parameters=args.af_model_parameters,
        min_allele_frequency=args.min_allele_frequency,
        min_read_depth=args.min_read_depth,
        default_ref_bias=args.ref_bias,
        verbose=args.verbose
    )
    phaser = Phaser(
        llr_cap=args.phasing_log_ratio_cap,
        llr_threshold=args.phasing_sample_llr_threshold,
        consensus_fdr=args.phasing_consensus_fdr,
        max_num_contig_segs=args.phasing_max_num_contig_segs
    )
    data.deduplicate_samples()
    data.subset_data()
    data.pileup_likelihoods = data.get_pileup_likelihoods_parallel(genotyper=genotyper, max_processes=args.threads)
    data.validate_pileup_likelihood_allele_frequencies()
    data.sample_correlation = data.validate_sample_correlation(genotyper=genotyper)
    data.joint_genotype_likelihood = genotyper.get_joint_genotype_likelihood(pileup_likelihoods=data.pileup_likelihoods, normal_samples=args.normal_sample, normal_to_tumor_weight=args.normal_to_tumor_weight)
    data.sort_genomic_positions()
    data.phase_hets(phaser=phaser)
    data.write_output(output_dir=args.output_dir, vcf_format=args.format, args=args)


def print_args(args):
    if args.verbose:
        message("Calling Genotyper")
        print("Arguments:")
        for key, value in vars(args).items():
            print(f"  {key}: {value}")
        print()


class Pileup(object):
    """
    Represents and processes pileup data from a specified file. (E.g. output of GATK GetPileupSummaries.)

    This class reads pileup data from a given file path, filters the data based on
    a minimum read depth, and extracts sample name information from the file header.

    Attributes:
        file_path (str): Path to the pileup file.
        columns (list of str): Column names for the DataFrame.
        df (pandas.DataFrame): DataFrame containing the pileup data.
        bam_sample_name (str): Extracted sample name from the file header.

    Args:
        file_path (str): Path to the pileup file.
        min_read_depth (int): Minimum read depth threshold for filtering data.
    """

    def __init__(self, file_path: str = None, min_read_depth: int = 0, min_allele_frequency: float = 0):
        self.file_path = file_path
        self.columns = ["contig", "position", "ref_count", "alt_count", "other_alt_count", "allele_frequency"]
        self.column_types = [str, int, int, int, int, float]
        can_read_file = True
        try:
            self.df = (
                pd.read_csv(file_path, sep="\t", comment="#", header=0, names=self.columns, low_memory=False)
                if file_path is not None
                else pd.DataFrame(columns=self.columns)
            )
        except Exception as e:
            can_read_file = False
            warnings.warn(f"Exception reading pileup file {file_path}: {e}")
            warnings.warn(f"Setting pileup to empty DataFrame.")
            self.df = pd.DataFrame(columns=self.columns)
        self.df = self.df.astype({key: dtype for key, dtype in zip(self.columns, self.column_types)})
        self.df = self.df.loc[self.df[["ref_count", "alt_count", "other_alt_count"]].sum(axis=1) >= min_read_depth]
        self.df = self.df.loc[self.df["allele_frequency"] >= min_allele_frequency]
        self.bam_sample_name = None
        if file_path is not None and can_read_file:
            open_func = gzip.open if file_path.endswith(".gz") else open
            with open_func(file_path, "rt") as pileup_file:
                pileup_header = pileup_file.readline().strip()
                if pileup_header.startswith("#<METADATA>SAMPLE="):
                    self.bam_sample_name = pileup_header.removeprefix("#<METADATA>SAMPLE=")


class PileupLikelihood(object):
    """
    Encapsulates genotype likelihoods for a given pileup.

    This class stores the genotype likelihoods along with associated sample name information.
    It also provides functionality to reindex the underlying DataFrame.

    Attributes:
        df (pandas.DataFrame): DataFrame containing the genotype likelihoods.
        assigned_sample_name (str): The name assigned to the sample by the user.
        bam_sample_name (str): The sample name obtained from the BAM file (annotated in the pileup header).

    Args:
        pileup_likelihood (pandas.DataFrame): DataFrame containing genotype likelihoods.
        assigned_sample_name (str): Assigned name of the sample.
        bam_sample_name (str): Sample name from the BAM file.
    """

    def __init__(self, pileup_likelihood: pd.DataFrame, assigned_sample_name: str, bam_sample_name: str):
        self.df = pileup_likelihood.set_index(["contig", "position"])
        self.assigned_sample_name = assigned_sample_name
        self.bam_sample_name = bam_sample_name

    def reindex(self, *args, **kwargs):
        self.df = self.df.reindex(*args, **kwargs)
        self.df["./."] = self.df["./."].fillna(1.0)  # Declare missing loci as outliers.
        self.df.fillna(0, inplace=True)  # Set all other likelihoods and read counts to 0.
        self.df = self.df.astype({"ref_count": int, "alt_count": int, "other_alt_count": int})
        return self


class Contamination(object):
    """
    Represents contamination data from a specified file. (E.g. output of GATK CalculateContamination.)

    This class reads contamination data, including contamination levels and associated
    errors, from a given file. It also extracts the sample name from the file.

    Attributes:
        file_path (str): Path to the contamination file.
        keys (list of str): Keys representing the contamination data fields.
        dict (dict): Dictionary containing contamination data.
        value (float): Contamination value.
        error (float): Error associated with the contamination value.
        bam_sample_name (str): Extracted sample name.

    Args:
        file_path (str): Path to the contamination file.
    """

    _min_error = 1e-6

    def __init__(self, file_path: str = None):
        self.file_path = file_path
        self.keys = ["sample", "contamination", "contamination_error"]
        default_dict = {key: default for key, default in zip(self.keys, [None, 0, self._min_error])}
        try:
            self.dict = (
                pd.read_csv(file_path, sep="\t", comment="#", header=0, names=self.keys, low_memory=False).loc[0].to_dict()
                if file_path is not None
                else default_dict
            )
        except Exception as e:
            warnings.warn(f"Exception reading contamination file {file_path}: {e}")
            warnings.warn(f"Setting contamination to default dictionary.")
            self.dict = default_dict
        self.value = self.dict["contamination"]
        self.error = max(self.dict["contamination_error"], self._min_error)
        self.bam_sample_name = self.dict["sample"]


class Segments(object):
    """
    Represents segmentation data. (Output of GATK CalculateContamination OR ModelSegments.)

    This class is responsible for reading segmentation data related to allele
    fractions from a file, and extracting sample name information from the file header.

    Attributes:
        file_path (str): Path to the segments file.
        columns (list of str): Column names for the DataFrame.
        df (pandas.DataFrame): DataFrame containing the segmentation data.
        bam_sample_name (str): Extracted sample name from the file header.

    Args:
        file_path (str): Path to the segments file.
    """

    def __init__(self, file_path: str = None):
        self.file_path = file_path
        self.columns = ["contig", "start", "end", "minor_allele_fraction"]
        default_df = pd.DataFrame(columns=self.columns)
        can_read_file = True
        try:
            if file_path is None:
                self.df = default_df
            elif file_path.endswith(".segments"):  # GATK CalculateContamination
                self.df = pd.read_csv(file_path, sep="\t", comment="#", header=0, names=self.columns, low_memory=False)
            elif file_path.endswith(".seg"):  # GATK ModelSegments
                self.df = pd.read_csv(file_path, sep="\t", comment="@", low_memory=False).rename(columns={"CONTIG": "contig", "START": "start", "END": "end"})
                self.df["minor_allele_fraction"] = self.df["MINOR_ALLELE_FRACTION_POSTERIOR_50"].where(
                    self.df["MINOR_ALLELE_FRACTION_POSTERIOR_50"].isna(),
                    1 / (2 * 2 ** self.df["LOG2_COPY_RATIO_POSTERIOR_50"]),
                )
                self.df["minor_allele_fraction"] = self.df["minor_allele_fraction"].where(
                    self.df["minor_allele_fraction"] < 0.5,
                    1 - self.df["minor_allele_fraction"]
                )
                self.df = self.df[self.columns]
        except Exception as e:
            can_read_file = False
            warnings.warn(f"Exception reading segments file {file_path}: {e}")
            warnings.warn(f"Setting segments to empty DataFrame.")
            self.df = default_df
        self.bam_sample_name = None
        if file_path is not None and can_read_file:
            open_func = gzip.open if file_path.endswith(".gz") else open
            with open_func(file_path, "rt") as segments_file:
                segments_header = segments_file.readline().strip()
                if segments_header.startswith("#<METADATA>SAMPLE="):
                    self.bam_sample_name = segments_header.removeprefix("#<METADATA>SAMPLE=")
                elif segments_header.startswith("@RG"):
                    self.bam_sample_name = [c.removeprefix("SM:") for c in segments_header.split("\t") if c.startswith("SM:")][0]


class AFParameters(object):
    """
    Represents allelic fraction data from a specified file. (E.g. output of GATK ModelSegments.)

    Attributes:
        file_path (str): Path to the contamination file.
        keys (list of str): Keys representing the contamination data fields.
        dict (dict): Dictionary containing contamination data.
        ref_bias (float): reference bias value lambda in f_obs = f / (f + (1 - f) * lambda)
        ref_bias_var (float): variance of reference bias posterior distribution
        outlier_probability (float): prior likelihood that a given locus is an outlier.

    Args:
        file_path (str): Path to the allelic fraction parameters file.
    """
    def __init__(self, file_path: str = None, default_ref_bias: float = 1):
        self.file_path = file_path
        self.keys = ["MEAN_BIAS", "BIAS_VARIANCE", "OUTLIER_PROBABILITY"]
        default_dict = {key: default for key, default in zip(self.keys, [default_ref_bias, 0, 0])}
        try:
            self.dict = (
                pd.read_csv(file_path, sep="\t", comment="@", low_memory=False).set_index("PARAMETER_NAME")["POSTERIOR_50"].to_dict()
                if file_path is not None
                else default_dict
            )
        except Exception as e:
            warnings.warn(f"Exception reading contamination file {file_path}: {e}")
            warnings.warn(f"Setting af model parameters to default dictionary.")
            self.dict = default_dict
        self.ref_bias = self.dict["MEAN_BIAS"]
        self.ref_bias_var = self.dict["BIAS_VARIANCE"]
        # The outlier probability is not used for this model since the outlier
        # emissions and their priors are fundamentally different. The model here
        # uses population allele frequency priors per locus.
        self.outlier_probability = self.dict["OUTLIER_PROBABILITY"]


def cross_check_sample_name(*args) -> str:
    """
    Cross-checks the sample names across multiple data sources.

    Compares the bam_sample_name attribute of each provided object. If there are
    discrepancies in sample names, it issues warnings. Returns the first unique
    sample name found.

    Args:
        *args: Variable length argument list, each argument is an object with
        bam_sample_name and file_path attributes.

    Returns:
        str: The first unique sample name found among the arguments.

    Raises:
        Warning: If there are discrepancies in sample names across files.
    """
    for a, b in combinations(args, 2):
        if a.bam_sample_name is not None and b.bam_sample_name is not None and a.bam_sample_name != b.bam_sample_name:
            warnings.warn(f"Sample name in {a.file_path} ({a.bam_sample_name}) does not match sample name in {b.file_path} ({b.bam_sample_name}).")
    sample_names = np.unique([a.bam_sample_name for a in args if a.bam_sample_name is not None])
    return str(sample_names[0]) if len(sample_names) else None


def join_pileups(pileups: list["Pileup"]):
    joint = Pileup()
    df = pd.concat([p.df for p in pileups], ignore_index=True)
    joint.df = df.groupby(["contig", "position"]).agg(
        {
            "ref_count": "sum",
            "alt_count": "sum",
            "other_alt_count": "sum",
            "allele_frequency": "max"
        }
    ).reset_index()
    joint.bam_sample_name = np.unique([p.bam_sample_name for p in pileups])[0]  # should be the same
    return joint


def join_contaminations(contaminations: list["Contamination"]):
    joint = Contamination()
    joint.value, joint.error = np.average([c.value for c in contaminations], weights=[1 / c.error ** 2 for c in contaminations], returned=True)
    joint.bam_sample_name = np.unique([c.bam_sample_name for c in contaminations])[0]  # should be the same
    joint.dict = {"sample": joint.bam_sample_name, "contamination": joint.value, "contamination_error": joint.error}
    return joint


def harmonize_segments(segments: list["Segments"]):
    joint = Segments()
    segs = pd.concat([s.df for s in segments], ignore_index=True)
    split_segs = pd.concat([non_overlapping(contig_group) for contig, contig_group in segs.groupby(by=["contig"])])
    segs = split_segs.groupby(["contig", "start", "end"]).agg({"minor_allele_fraction": "min"})
    joint.df = merge_abutting_intervals(segs).reset_index()
    joint.bam_sample_name = np.unique([s.bam_sample_name for s in segments])[0]  # should be the same
    return joint


def split_intervals(intervals: pd.DataFrame, start_col: str, end_col: str):
    starts = pd.DataFrame(data={"start": 1}, index=intervals[start_col])
    ends = pd.DataFrame(data={"end": -1}, index=intervals[end_col])
    transitions = pd.merge(starts, ends, how="outer", left_index=True, right_index=True).fillna(0)
    # This dataframe stores per position the type of transitions. Automatically sorted.
    # Now, we need to know at each position if we are at an interval start or not.
    # This can be done by counting the opening & closing parenthesis.
    transitions["is_start"] = (transitions.pop("start") + transitions.pop("end")).cumsum().astype(bool)
    # This handles overlapping intervals.
    # Create interval indices and take all intervals that have a valid start.
    new_intervals = pd.IntervalIndex.from_breaks(transitions.index, closed="both")[transitions["is_start"][:-1]]
    new_intervals = new_intervals[new_intervals.length > 0].unique()

    # Remove intervals of length 1 that are abutting with the previous and next interval.
    if new_intervals.size > 2:
        new_starts = new_intervals.left
        new_ends = new_intervals.right
        abutting = new_starts[1:] == new_ends[:-1]
        mid_abutting = np.concatenate([[False], abutting[1:] & abutting[:-1], [False]])
        length1 = new_intervals.length <= 1
        # don't remove length 1 intervals next to other length 1 intervals
        good_length1 = np.concatenate([[False], ~length1[:-2] & length1[1:-1] & ~length1[2:], [False]])
        mask = mid_abutting & good_length1
        new_intervals = new_intervals[~mask]
    return new_intervals


def non_overlapping(intervals: pd.DataFrame, start_col: str = "start", end_col: str = "end", verbose: bool = False):
    new_intervals = split_intervals(intervals, start_col, end_col)

    # Map each original interval to the new non-overlapping intervals
    non_overlapping_intervals = []
    for new_interval in new_intervals:
        contains_interval = intervals[(intervals[start_col] <= new_interval.left) & (intervals[end_col] > new_interval.left)]
        for _, interval in contains_interval.iterrows():
            new_row = copy(interval)
            new_row[start_col] = new_interval.left
            new_row[end_col] = new_interval.right
            non_overlapping_intervals.append(new_row)

    non_overlapping_df = pd.DataFrame(non_overlapping_intervals)
    return non_overlapping_df


def merge_abutting_intervals(intervals: pd.DataFrame, dist: int = 1, tol: float = 1e-3):
    new_intervals = []
    previous_interval = None
    for i, interval in intervals.iterrows():
        if previous_interval is None:
            previous_interval = interval
            new_intervals.append(interval)
            continue

        # interval.name == (contig, start, end)
        if interval.name[0] == previous_interval.name[0] and abs(interval.name[1] - previous_interval.name[2]) <= dist:
            if (abs(previous_interval - interval) < tol).all():
                new_interval = copy(previous_interval)
                new_interval.name = (new_interval.name[0], new_interval.name[1], interval.name[2])
                new_intervals[-1] = new_interval
                previous_interval = new_interval
                continue

        # Reduce end of the previous interval by 1 if it overlaps with the start of the current interval.
        # Even though we use left-closed intervals, some downstream tools expect intervals to be both left- and right-closed.
        if previous_interval.name[0] == interval.name[0] and previous_interval.name[2] == interval.name[1]:
            previous_interval.name = (previous_interval.name[0], previous_interval.name[1], previous_interval.name[2] - 1)
            new_intervals[-1] = previous_interval

        new_intervals.append(interval)
        previous_interval = interval

    merged_intervals = pd.DataFrame(new_intervals)
    merged_intervals.index.names = ["contig", "start", "end"]
    return merged_intervals


class VCF(object):
    """
    Represents and processes VCF (Variant Call Format) data.

    This class is designed to read VCF data from a file and process it by filtering
    for bi-allelic single nucleotide variants (SNVs) and reformatting it into a usable DataFrame.

    Attributes:
        file_path (str): Path to the VCF file.
        ref_dict (list): List of reference sequence information extracted from the VCF file.
        df (pandas.DataFrame): DataFrame containing processed VCF data.

    Args:
        file_path (str, optional): Path to the VCF file.
    """
    columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]

    def __init__(self, file_path: str = None, verbose: bool = False):
        self.file_path = file_path
        self.ref_dict = self.get_vcf_sequence_dict()
        self.contigs = self.get_contigs()
        if file_path is not None:
            open_func = gzip.open if file_path.endswith(".gz") or file_path.endswith(".bgz") else open
            with open_func(file_path, "rt") as f:
                self.df = pd.read_csv(f, sep="\t", comment="#", header=None, low_memory=False, names=self.columns)
        else:
            self.df = pd.DataFrame(columns=self.columns)
    
        self.df = self.df.astype({"CHROM": str, "POS": int, "ID": str, "REF": str, "ALT": str, "INFO": str})
        n_input_vars = self.df.shape[0]
        ## GetPileupSummaries only supports SNVs well.
        ## CAUTION: the input vcf may contain multi-allelic sites and other stuff!!!! (removed here)
        snv_mask = (self.df["REF"].apply(len) == 1) & (self.df["ALT"].apply(len) == 1) & ~self.df["REF"].isin(["-", "."]) & ~self.df["ALT"].isin(["-", "."])
        self.df = self.df.loc[snv_mask]
        self.df = self.df.drop_duplicates(subset=["CHROM", "POS"], keep=False)
        n_kept_vars = self.df.shape[0]
        if n_kept_vars < n_input_vars:
            print(f"    Dropping {n_input_vars - n_kept_vars} multi-allelic and non-SNP loci from Variants. {n_kept_vars} loci remaining.") if verbose else None
        self.df = self.df.set_index(["CHROM", "POS"])
        self.df.index.names = ["contig", "position"]  # align to PileupLikelihood index names

    def get_vcf_sequence_dict(self):
        """
        Extracts reference sequence information from the VCF file.

        Parses the VCF file to retrieve the sequence dictionary information
        (i.e., reference sequences) included in the VCF file header.

        Returns:
            list: A list of reference sequence information strings.
        """
        ref_dict = []
        if self.file_path is None:
            return ref_dict
        open_func = gzip.open if self.file_path.endswith(".gz") or self.file_path.endswith(".bgz") else open
        with open_func(self.file_path, "rt") as ifile:
            for line in ifile:
                if line.startswith("##contig"):
                    ref_dict.append(line.strip())
                if line.startswith("#CHROM"):
                    break
        ifile.close()
        return ref_dict

    def get_contigs(self):
        return [c.split("ID=")[-1].split(",")[0] for c in self.ref_dict]


def join_vcfs(vcfs: list["VCF"], verbose=False):
    if len(vcfs) == 0:
        return None
    elif len(vcfs) == 1:
        return vcfs[0]
    else:
        joint = copy(vcfs[0])
        df = pd.concat([vcf.df.reset_index() for vcf in vcfs])
        n_input_vars = df.shape[0]
        df.drop_duplicates(subset=["contig", "position"], keep="first")
        n_kept_vars = df.shape[0]
        if n_kept_vars < n_input_vars:
            print(f"    Dropping {n_input_vars - n_kept_vars} duplicate loci from multiple variant inputs. {n_kept_vars} loci remaining total.") if verbose else None
        joint.df = df.set_index(["contig", "position"])
        return joint


class Genotyper(object):
    """
    Performs genotype likelihood calculations using specified statistical models.

    Based on the implementation from the Broad Institute's GATK ContaminationModel.
    https://github.com/broadinstitute/gatk/blob/master/src/main/java/org/broadinstitute/hellbender/tools/walkers/contamination/ContaminationModel.java

    Attributes:
        contig_ploidies (dict): Dictionary mapping contig names to ploidy levels.
        overdispersion (float): Overdispersion parameter for the beta-binomial model.
        ref_bias (float): SNP bias of the reference allele.
        min_genotype_likelihood (float): Minimum genotype likelihood threshold for filtering.
        min_error_rate (float): Minimum error rate for a read supporting an alternate allele to be a false positive.
        max_error_rate (float): Maximum error rate for a read supporting an alternate allele to be a false positive.
        outlier_prior (float): Prior probability for a variant to be an outlier.
        select_hets (bool): Whether to select only heterozygous sites.
        verbose (bool): Whether to print information during execution

    Args:
        sex (str, optional): Genotype for sex chromosomes. (default: "XXY" for uninformed prior on X and Y)
        overdispersion (float, optional): Overdispersion parameter for the beta-binomial model.
        ref_bias (float, optional): SNP bias of the reference allele.
        min_genotype_likelihood (float, optional): Minimum genotype likelihood threshold.
        select_hets (bool, optional): Flag to select only heterozygous sites.
    """
    genotypes = ["0/0", "0/1", "1/1", "./."]

    def __init__(self, sex: str = "XXY", overdispersion: float = 100, ref_bias: float = 1.1, min_genotype_likelihood: float = 0.95, min_error_rate: float = 1e-4, max_error_rate: float = 1e-1, outlier_prior: float = 1e-4, select_hets: bool = False, verbose: bool = False):
        self.contig_ploidies = self.get_contig_ploidies(sex)
        self.overdispersion = overdispersion  # overdispersion parameter for beta-binomial model
        self.ref_bias = ref_bias  # SNP bias of the ref allele
        self.min_genotype_likelihood = min_genotype_likelihood
        self.min_error_rate = min(min_error_rate, max_error_rate)
        self.max_error_rate = max(min_error_rate, max_error_rate)
        self.outlier_prior = outlier_prior
        self.select_hets = select_hets
        self.verbose = verbose

    @staticmethod
    def get_contig_ploidies(sex):
        ploidy = {str(i): 2 for i in range(1, 23)}
        ploidy["X"] = 2
        ploidy["Y"] = 0
        ploidy["MT"] = 1
        ploidy["M"] = 1

        if sex in ["XX", "Female", "female"]:
            ploidy["X"] = 2
            ploidy["Y"] = 0
        elif sex in ["XY", "Male", "male"]:
            ploidy["X"] = 1
            ploidy["Y"] = 1

        # Handle cases with XXY, etc. We only need to count up to ploidy 2 as this
        # differentiates the presence or absence of HETs.
        if "XX" in sex:
            ploidy["X"] = 2
        elif "X" in sex:
            ploidy["X"] = 1

        if "Y" in sex:
            ploidy["Y"] = 1

        # Add "chr" prefix to contig names for other references.
        ploidy = ploidy | {f"chr{c}": p for c, p in ploidy.items()}

        return ploidy

    def get_error_prob(self, pileup: pd.DataFrame) -> float:
        """
        Estimates the error probability based on the pileup data.

        The error probability is calculated as a ratio of 'other_alt_count' to the
        total counts, clipped within specified bounds.

        Args:
            pileup (pd.DataFrame): DataFrame containing pileup data.

        Returns:
            float: Estimated error probability.
        """
        other_alt_counts = pileup["other_alt_count"].sum()
        total_counts = pileup[["ref_count", "alt_count", "other_alt_count"]].sum(axis=1).sum()
        return np.clip(3/2 * other_alt_counts / max(1, total_counts), a_min=self.min_error_rate, a_max=self.max_error_rate)

    def select_confident_calls(self, likelihoods: pd.DataFrame, genotypes: list[str] = None, normalized: bool = False) -> pd.DataFrame:
        """
        Filters genotype likelihoods based on a minimum threshold.

        Args:
            likelihoods (pd.DataFrame): DataFrame containing genotype likelihoods.
            genotypes (list of str, optional): List of genotype columns to consider.
            normalized (bool): likelihoods can be interpreted as probabilities

        Returns:
            pd.DataFrame: DataFrame containing filtered genotype likelihoods.
        """
        if genotypes is None:
            genotypes = self.genotypes
        total = likelihoods[self.genotypes].sum(axis=1) if not normalized else np.ones(likelihoods.shape[0])
        selection = likelihoods[genotypes].max(axis=1) >= self.min_genotype_likelihood * total
        return likelihoods.loc[selection]

    def calculate_genotype_likelihoods(self, pileup: pd.DataFrame, contamination: float = 0.001, segments: pd.DataFrame = None, ref_bias: float = 1) -> pd.DataFrame:
        """
        Calculates genotype likelihoods for given pileup data, and optional allelic copy ratio segmentation, and contamination estimate.

        Args:
            pileup (pd.DataFrame): DataFrame containing pileup data.
            contamination (float, optional): Contamination level to consider in calculations.
            segments (pd.DataFrame, optional): DataFrame containing segmentation data.
            ref_bias (float, optional): Bias for observing the reference allele over the alternate allele.

        Returns:
            pd.DataFrame: DataFrame with genotype likelihoods for each genomic position.
        """
        if pileup.empty:
            return pd.DataFrame(columns=["contig", "position", "ref_count", "alt_count", "other_alt_count", "allele_frequency", "0/0", "0|1", "1|0", "0/1", "1/1", "./."])

        error = self.get_error_prob(pileup=pileup)

        likelihoods = []
        in_any_seg = pd.Series(False, index=pileup.index)
        if segments is not None:
            for _, seg in segments.iterrows():
                seg_mask = (
                    (pileup["contig"] == seg["contig"])
                    & (seg["start"] <= pileup["position"])
                    & (pileup["position"] <= seg["end"])
                )
                seg_pileup = pileup.loc[seg_mask]
                if seg_pileup.empty:
                    continue
                in_any_seg |= seg_mask
                likelihoods.append(
                    self.calculate_genotype_likelihoods_per_segment(
                        pileup=seg_pileup,
                        minor_af=seg["minor_allele_fraction"],
                        contamination=contamination,
                        ref_bias=ref_bias,
                        error=error
                    )
                )
        seg_pileup = pileup.loc[~in_any_seg]
        if not seg_pileup.empty:
            likelihoods.append(
                self.calculate_genotype_likelihoods_per_segment(
                    pileup=seg_pileup,
                    contamination=contamination,
                    ref_bias=ref_bias,
                    error=error
                )
            )

        genotype_likelihoods = pd.concat(likelihoods).astype(
            {
                "contig": str,
                "position": int,
                "ref_count": int,
                "alt_count": int,
                "other_alt_count": int,
                "0/0": float,
                "0|1": float,
                "1|0": float,
                "0/1": float,
                "1/1": float,
                "./.": float,
            }
        )
        return genotype_likelihoods

    def calculate_genotype_likelihoods_per_segment(self, pileup: pd.DataFrame, minor_af: float = 0.47, contamination: float = 0.001, ref_bias: float = 1, error: float = 0.01) -> pd.DataFrame:
        """
        Calculates genotype likelihoods for a specific segment of pileup data.
        https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3167057/

        This method applies a beta-binomial model to compute the likelihoods.

        Note: Suppose X[i] ~ Binom(N, .5), and minor_af[i] = min(X[i]/N, 1 - X[i]/N),
              then median(minor_af) = 0.47 for N = 100 (sequencing depth).

        Args:
            pileup (pd.DataFrame): DataFrame containing pileup data for the segment.
            minor_af (float, optional): Minor allele fraction, default is 0.5.
            contamination (float, optional): Level of contamination to consider in the calculation.
            ref_bias (float, optional): Bias for observing the reference allele over the alternate allele.
            error (float, optional): Error probability to be used in the calculation.

        Returns:
            pd.DataFrame: DataFrame containing calculated genotype likelihoods for the segment.
        """

        contig = pileup["contig"]
        n = pileup[["ref_count", "alt_count"]].sum(axis=1)
        alt = pileup["alt_count"]
        popaf = pileup["allele_frequency"]  # population allele frequency from SNP panel

        def bias(_f):
            return _f / (_f + (1 - _f) * self.ref_bias)

        f_aa = contamination * popaf + (1 - contamination) * error / 3  # errors spread across the three wrong bases
        f_bb = contamination * popaf + (1 - contamination) * (1 - error)
        f_ab = bias(contamination * popaf + (1 - contamination) * minor_af)
        f_ba = bias(contamination * popaf + (1 - contamination) * (1 - minor_af))

        # include probability from seq error assuming error rate for alt>ref == ref>alt
        f_ab = f_ab * (1 - error) + (1 - f_ab) * error
        f_ba = f_ba * (1 - error) + (1 - f_ba) * error

        s = self.overdispersion
        # Give sites with zero read depth high outlier (no information) prior
        outlier_prior = self.outlier_prior * np.ones_like(n)
        outlier_prior[n == 0] = 1 - self.outlier_prior
        logo = np.log(outlier_prior)
        log1o = np.log1p(-outlier_prior)

        ploidy = contig.apply(lambda c: self.contig_ploidies.get(c, 2))
        # If chromosomal ploidy is 0, then only possible genotype is ./.
        log_p0_prior = np.where(ploidy == 0, -np.inf, 0)
        # If chromosomal ploidy is less than 2, then only possible genotypes are 0/0, 1/1, and ./.
        log_plt2_prior = np.where(ploidy < 2, -np.inf, 0)

        # Prior probabilities of genotypes:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=RuntimeWarning)
            log1f = np.log1p(-popaf)
            logf = np.log(popaf)

            outlier = logo + 2 * log1f + st.betabinom(n, 1, 1).logpmf(alt)  # max entropy
            aa = log_p0_prior + log1o + 2 * log1f + st.binom(n, f_aa).logpmf(alt)
            bb = log_p0_prior + log1o + 2 * logf + st.binom(n, f_bb).logpmf(alt)
            ab = log_plt2_prior + log1o + log1f + logf + st.betabinom(n, f_ab * s, (1 - f_ab) * s).logpmf(alt)
            ba = log_plt2_prior + log1o + log1f + logf + st.betabinom(n, f_ba * s, (1 - f_ba) * s).logpmf(alt)

            het = np.logaddexp(ab, ba)

        likelihoods = pileup.copy()
        likelihoods["0/0"] = np.exp(aa)
        likelihoods["0|1"] = np.exp(ab)
        likelihoods["1|0"] = np.exp(ba)
        likelihoods["0/1"] = np.exp(het)
        likelihoods["1/1"] = np.exp(bb)
        likelihoods["./."] = np.exp(outlier)

        # Output raw likelihoods since they carry information about read depth
        # or confidence for joint genotyping.

        return likelihoods

    def get_joint_genotype_likelihood(self, pileup_likelihoods: list[PileupLikelihood], normal_samples: list[str] = None, normal_to_tumor_weight: float = 1.0) -> pd.DataFrame:
        """
        Calculates joint genotype likelihoods across multiple samples from the
        same individual.

        This method aggregates genotype likelihoods from multiple PileupLikelihood
        instances, normalizes them, and applies filters based on the minimum genotype
        likelihood threshold and selection criteria for heterozygous sites.

        Args:
            pileup_likelihoods (list[PileupLikelihood]): List of PileupLikelihood objects from different samples.
            normal_samples (list[str]): List of normal sample names
            normal_to_tumor_weight (float): relative weight of normal samples to tumor samples in determining genotype

        Returns:
            pd.DataFrame: DataFrame with joint genotype likelihoods for each genomic position.
        """
        # Calculate the joint genotype likelihoods.
        message(f"Calculating joint genotype likelihoods: ", end="", flush=True) if self.verbose else None

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=RuntimeWarning)

            weights = (
                np.array([normal_to_tumor_weight if pl.assigned_sample_name in normal_samples else 1 for pl in pileup_likelihoods])
                if normal_samples is not None
                else np.ones(len(pileup_likelihoods))
            )

            def nansum(ndarray):
                # testing for NaN throws errors for some reason, so testing for 1 instead (log(p) must be < 0)
                mask = ndarray != 1
                if np.sum(mask) == 0:
                    return -np.inf
                return np.average(ndarray[mask], weights=weights[mask])

            joint_log_likelihood = {}
            for gt in self.genotypes:
                log_likelihoods = []
                for pl in pileup_likelihoods:
                    # Only use evidence in covered regions.
                    coverage_mask = pl.df[["ref_count", "alt_count"]].sum(axis=1) > 0
                    log_likelihoods.append(np.log(pl.df.loc[coverage_mask, gt]))

                joint_log_likelihood[gt] = (
                    log_likelihoods[0].fillna(1) if len(log_likelihoods) == 1
                    else pd.concat(log_likelihoods, axis=1).fillna(1).apply(nansum, axis=1)
                )

                print(gt, end=" ", flush=True) if self.verbose else None
        print("") if self.verbose else None

        # Normalize likelihoods
        total = pd.concat(joint_log_likelihood.values(), axis=1).apply(sp.logsumexp, axis=1)
        genotype_likelihood = {
            gt: np.exp(joint_log_likelihood[gt].sub(total, axis=0))
            for gt in self.genotypes
        }

        # Count total allelic read counts across all samples.
        genotype_likelihood |= {
            count: pd.concat([pl.df[count] for pl in pileup_likelihoods], axis=1).sum(axis=1)
            for count in ["ref_count", "alt_count", "other_alt_count"]
        }

        # Average allele frequencies across all samples. They should be the same.
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=RuntimeWarning)
            genotype_likelihood["allele_frequency"] = pd.concat(
                [pl.df["allele_frequency"] for pl in pileup_likelihoods],
                axis=1
            ).apply(np.nanmean, axis=1)

        df = pd.DataFrame.from_dict(genotype_likelihood)
        message(f"Processed {df.shape[0]} records.") if self.verbose else None

        df = self.select_confident_calls(likelihoods=df, normalized=True)
        message(f"Selected  {df.shape[0]} confident calls.") if self.verbose else None

        df = self.select_confident_calls(likelihoods=df, genotypes=[g for g in self.genotypes if g != "./."], normalized=True)
        message(f"Selected  {df.shape[0]} non-outliers.") if self.verbose else None

        if self.select_hets:
            df = self.select_confident_calls(likelihoods=df, genotypes=["0/1"], normalized=True)
            message(f"Selected  {df.shape[0]} heterozygous SNPs.") if self.verbose else None

        df["GT"] = df[["0/0", "0/1", "1/1", "./."]].idxmax(axis=1)

        return df


class Phaser(object):
    """
    Phases heterozygous genotype calls using per-sample phased likelihoods, with segment-aware consistency.
    """

    class _UF:
        """Union-Find data structure with path compression and union by rank."""
        def __init__(self, n):
            self.p = np.arange(n)
            self.r = np.zeros(n, dtype=int)

        def find(self, x):
            while self.p[x] != x:
                self.p[x] = self.p[self.p[x]]
                x = self.p[x]
            return x

        def union(self, a, b):
            ra, rb = self.find(a), self.find(b)
            if ra == rb: return
            if self.r[ra] < self.r[rb]:
                self.p[ra] = rb
            elif self.r[ra] > self.r[rb]:
                self.p[rb] = ra
            else:
                self.p[rb] = ra
                self.r[ra] += 1

    def __init__(self, llr_cap: float = 6, llr_threshold: float = 0.4, consensus_fdr: float = 0.01, max_num_contig_segs: int = 1000):
        self.llr_cap = llr_cap
        self.llr_threshold = llr_threshold
        self.max_num_contig_segs = max_num_contig_segs
        self.consensus_fdr = consensus_fdr
        self.consensus_threshold = sp.logit(1 - consensus_fdr)

    def phase_hets(
        self,
        joint_genotype_likelihood: pd.DataFrame,
        pileup_likelihoods: list[PileupLikelihood],
        segments: list[Segments],
    ) -> pd.DataFrame:
        """
        Replace 0/1 in self.joint_genotype_likelihood with 0|1 or 1|0 based on
        consensus of per-sample phased likelihoods, with segment-aware consistency.

        Args:
            joint_genotype_likelihood: DataFrame with joint genotype likelihoods.
            pileup_likelihoods: List of PileupLikelihood objects from different samples.
            segments: List of Segments objects corresponding to the pileup_likelihoods.

        Returns:
            DataFrame with phased genotype likelihoods.
        """
        message("Phasing HETs ...")

        if joint_genotype_likelihood.empty:
            return joint_genotype_likelihood

        jl = joint_genotype_likelihood.copy()
        contigs = jl.index.get_level_values(level="contig").drop_duplicates()

        cons_parts = []  # list[pd.Series] signed consensus per locus
        num_components_per_contig = dict()
        num_hets_per_contig = dict()
        for contig in contigs:
            is_het_in_contig = jl.index.get_level_values(level="contig") == contig
            is_het_in_contig &= jl["GT"] == "0/1"

            num_hets_per_contig[contig] = np.sum(is_het_in_contig)

            idx = jl.index[is_het_in_contig]

            if not np.any(is_het_in_contig):
                continue

            V = []
            for pl, segs in zip(pileup_likelihoods, segments):
                hets = pl.df.reindex(idx)
                # pull as numpy and sanitize
                l01 = hets["0|1"].to_numpy(dtype=float)
                l10 = hets["1|0"].to_numpy(dtype=float)

                # turn any NaN/±inf into 0 (i.e., “no evidence”)
                l01 = np.nan_to_num(l01, nan=0.0, posinf=0.0, neginf=0.0)
                l10 = np.nan_to_num(l10, nan=0.0, posinf=0.0, neginf=0.0)

                tiny = 1e-32
                llr = np.log(np.maximum(l01, tiny)) - np.log(np.maximum(l10, tiny))
                llr = np.clip(llr, -self.llr_cap, self.llr_cap)
                llr = np.nan_to_num(llr, nan=0.0, posinf=0.0, neginf=0.0)

                if segs.df.empty:
                    V.append(llr)
                    continue

                seg_c = segs.df.loc[segs.df["contig"] == contig]
                if seg_c.empty:
                    V.append(llr)
                    continue

                if seg_c.shape[0] > self.max_num_contig_segs:
                    # Don't use evidence from over-segmented samples. That data
                    # can not be trusted.
                    # V.append(llr)
                    continue

                pos = idx.get_level_values(level="position").to_numpy()
                for _, row in seg_c.iterrows():
                    het_in_seg = (row["start"] <= pos) & (pos < row["end"])
                    if np.sum(het_in_seg):
                        vseg = np.zeros_like(llr)
                        vseg[het_in_seg] = llr[het_in_seg]
                        V.append(vseg)

            if not V:
                continue

            V = [np.asarray(v, dtype=float, order="C") for v in V]

            components = self.connected_components_from_tracks(V)
            num_components_per_contig[contig] = len(components)

            for comp in components:
                Vc = np.stack([V[j] for j in comp], axis=0)  # (Jc, L)
                if not np.any(Vc):
                    continue

                cover = np.any(Vc != 0.0, axis=0)  # (L,)
                if not cover.any():
                    continue

                idx_comp = idx[cover]
                Vcov = Vc[:, cover]

                Jc, Lc = Vcov.shape

                if Jc <= 20:
                    x = self.exact_sign_max_gray(Vcov)
                else:
                    x0 = self.spectral_init(Vcov, iters=10)
                    x = self.best_improvement(Vcov, x0, max_iter=200, eps=1e-9)

                # temperature shrinkage
                if Jc == 1:
                    tau = 1
                else:
                    rho_bar = self.avg_pairwise_corr(Vcov, mask=(np.abs(Vcov).mean(axis=0) < self.llr_threshold))
                    tau = np.sqrt(1 + (Jc - 1) * rho_bar)
                    if not np.isfinite(tau) or tau < 1:
                        tau = 1

                c_comp = Vcov.T @ x / tau  # consensus LLR per locus (MAP-oriented)
                cons_parts.append(pd.Series(c_comp, index=idx_comp))

        if not cons_parts:
            message("No components with usable evidence; no HETs to phase.")
            return jl

        cons = pd.concat(cons_parts)
        accept = np.abs(cons) > self.consensus_threshold

        pick01 = (cons > 0) & accept
        pick10 = (cons < 0) & accept

        jl.loc[pick01.index[pick01], "GT"] = "0|1"
        jl.loc[pick10.index[pick10], "GT"] = "1|0"

        phased_01_hets = int(pick01.sum())
        phased_10_hets = int(pick10.sum())
        num_hets_per_contig = pd.Series(num_hets_per_contig)
        num_hets = num_hets_per_contig.sum()
        R = phased_01_hets + phased_10_hets
        pct = (R / num_hets * 100.0) if num_hets else 0.0

        message(
            f"Phased {R} / {num_hets} HETs ({pct:.1f}%, 0|1: {phased_01_hets}, 1|0: {phased_10_hets})"
            f" via |consensus LLR| > {self.consensus_threshold:.3g} (FDR < {self.consensus_fdr})"
        )

        # Per-contig breakdown (phased / total / %)
        by_contig_calls = pd.Series(accept, index=cons.index).groupby(level='contig').sum()
        contig_table = pd.DataFrame({
            "phased": by_contig_calls.astype(int),
            "n_hets": num_hets_per_contig
        }).fillna(0).astype(int)
        contig_table["pct"] = (100.0 * contig_table["phased"] / contig_table["n_hets"]).round(1)
        contig_table[["phased", "n_hets"]] = contig_table[["phased", "n_hets"]].astype(str)
        contig_table.index.name = "contig"

        print(contig_table.reindex(contigs).T.to_string())
        print("")

        # Quantiles of |consensus| for total vs phased
        def _qtab(x, name):
            qs = np.linspace(0, 1, 21)
            tab = pd.Series(np.quantile(np.abs(x), qs), index=[int(100 * q) for q in qs]).to_frame("|consensus LLR|")
            tab.index.name = f"% quantile ({name})"
            return tab.T.round(2)

        print(_qtab(cons.values, "all").to_string())
        print("")

        if R > 0:
            print(_qtab(cons.values[accept], "phased").to_string())
            print("")

        num_components_per_contig = pd.Series(num_components_per_contig).to_frame("n_comp").reindex(contigs).fillna(0).astype(int)
        num_components_per_contig.index.name = "contig"
        print(num_components_per_contig.loc[num_components_per_contig["n_comp"] > 0].astype(str).T.to_string())
        print("")

        return jl

    def connected_components_from_tracks(self, V: list[np.ndarray]) -> list[list[int]]:
        """Find connected components from tracks V (list of J arrays of length L)."""
        J = len(V)
        if J == 0:
            return []
        L = V[0].size
        pos2nodes = [[] for _ in range(L)]
        for j, v in enumerate(V):
            nz = np.flatnonzero(v)
            for i in nz:
                pos2nodes[i].append(j)
        uf = self._UF(J)
        for nodes in pos2nodes:
            if len(nodes) > 1:
                first = nodes[0]
                for k in nodes[1:]:
                    uf.union(first, k)
        groups = defaultdict(list)
        for j in range(J):
            groups[uf.find(j)].append(j)
        # drop singleton components with all-zero tracks
        comps = []
        for idxs in groups.values():
            keep = any(np.any(V[j] != 0) for j in idxs)
            if keep:
                comps.append(idxs)
        return comps

    @staticmethod
    def avg_pairwise_corr(M, mask):
        # M: (Jc, Lc) raw (unflipped) rows; mask: (Lc,) bool for near-null sites
        X = M[:, mask]
        if X.shape[1] < 20:  # too few loci → skip
            return 0.0
        # drop rows with ~zero variance
        row_var = X.var(axis=1)
        keep = row_var > 1e-8
        rho_bar = 0.0
        if np.count_nonzero(keep) >= 2:
            C = np.corrcoef(X[keep])
            iu = np.triu_indices_from(C, k=1)
            rho_bar = np.maximum(0.0, np.nanmean(C[iu]))
        return rho_bar

    @staticmethod
    def exact_sign_max_gray(Vc: np.ndarray) -> np.ndarray:
        """Exact maximize ||sum x_j v_j||^2 for small J using Gray code. Vc shape (Jc, L)."""
        Jc, _ = Vc.shape
        x = np.ones(Jc, dtype=int)
        c = Vc.sum(axis=0)  # consensus for all +1
        best_x = x.copy()
        best_F = float(np.dot(c, c))
        prev_gray = 0
        for t in range(1, 1 << Jc):
            gray = t ^ (t >> 1)
            bit = (prev_gray ^ gray).bit_length() - 1
            # flip bit
            c = c - 2 * x[bit] * Vc[bit]
            x[bit] *= -1
            F = float(np.dot(c, c))
            if F > best_F:
                best_F, best_x = F, x.copy()
            prev_gray = gray
        return best_x

    @staticmethod
    def spectral_init(Vc: np.ndarray, iters: int = 10) -> np.ndarray:
        """Top eigenvector of W=V V^T via power iteration using implicit multiplies."""
        Jc = Vc.shape[0]
        z = np.random.default_rng(42).standard_normal(Jc)
        z /= (np.linalg.norm(z) + 1e-12)
        for _ in range(iters):
            y = Vc.T @ z  # L
            z = Vc @ y  # Jc   (this equals (V V^T) z)
            nrm = np.linalg.norm(z)
            if nrm == 0:
                break
            z /= nrm
        x0 = np.sign(z)
        x0[x0 == 0] = 1
        return x0.astype(int)

    @staticmethod
    def best_improvement(Vc: np.ndarray, x: np.ndarray, max_iter: int = 200, eps: float = 1e-9) -> np.ndarray:
        """
        Coordinate ascent: always flip the j with largest positive gain until none.
        Vc: (Jc, L), x: (Jc,)
        """
        # Jc = Vc.shape[0]
        c = Vc.T @ x  # consensus
        b = np.einsum("ij,ij->i", Vc, Vc)  # ||v_j||^2
        for _ in range(max_iter):
            a = Vc @ c  # a_j = <v_j, c>
            dF = -4.0 * x * a + 4.0 * b
            j = int(np.argmax(dF))
            if dF[j] <= eps:
                break
            # flip using OLD sign
            c = c - 2.0 * x[j] * Vc[j]
            x[j] = -x[j]
        return x


class GenotypeData(object):
    """
    Manages genotype likelihood calculations and output generation.

    This class handles the workflow of processing genotype likelihoods, correlating
    sample genotypes, sub-setting data, and writing outputs for genomic analysis.

    Attributes:
        individual_id (str): Identifier for the individual.
        samples (list[str]): List of sample names.
        vcf (VCF): VCF object.
        pileups (list[Pileup]): List of Pileup objects.
        contaminations (list[Contamination]): List of Contamination objects.
        segments (list[Segments]): List of Segments objects.
        pileup_likelihoods (list[PileupLikelihood]): List of PileupLikelihood objects.
        sample_correlation (pandas.DataFrame): DataFrame containing sample genotype correlation matrix.
        joint_genotype_likelihood (pandas.DataFrame): DataFrame containing joint genotype likelihoods.

    Args:
        individual_id (str): Identifier for the individual.
        samples (list[str]): List of sample names.
        variant (list[str]): Path to the VCF files containing germline variants and population allele frequencies.
        pileup (list[str]): List of paths to the allelic pileup input files (output of GATK's GetPileupSummaries; used for CalculateContamination).
        contamination (list[str], optional): List of paths to the contamination estimate input files (output of GATK's CalculateContamination).
        segments (list[str], optional): List of paths to the allelic copy ratio segmentation input files (output of GATK's CalculateContamination).
        min_read_depth (int, optional): Minimum read depth threshold for filtering data.
    """

    def __init__(self, individual_id: str, samples: list[str], variant: list[str] = None, pileup: list[str] = None, contamination: list[str] = None, segments: list[str] = None, af_model_parameters: list[str] = None, min_read_depth: int = 0, min_allele_frequency: float = 0, default_ref_bias: float = 1, verbose: bool = False):
        self.verbose = verbose
        message("Loading data:") if verbose else None
        self.individual_id = individual_id
        self.samples = samples

        print("  Variants") if verbose else None
        self.vcf = join_vcfs([VCF(file_path=v, verbose=verbose) for v in variant], verbose=verbose)
        print("  Pileups") if verbose else None
        self.pileups = (
            [Pileup()] * len(samples)
            if pileup is None
            else [Pileup(file_path=p, min_read_depth=min_read_depth, min_allele_frequency=min_allele_frequency) for p in pileup]
        )
        print("  Contamination") if verbose else None
        self.contaminations = (
            [Contamination()] * len(samples)
            if contamination is None
            else [Contamination(file_path=c) for c in contamination]
        )
        print("  Segments") if verbose else None
        self.segments = (
            [Segments()] * len(samples)
            if segments is None
            else [Segments(file_path=s) for s in segments]
        )
        print("  AF Model Parameters") if verbose else None
        self.af_parameters = (
            [AFParameters(default_ref_bias=default_ref_bias)] * len(samples)
            if af_model_parameters is None
            else [AFParameters(file_path=a, default_ref_bias=default_ref_bias) for a in af_model_parameters]
        )
        print() if verbose else None

        self.pileup_likelihoods = None
        self.sample_correlation = None
        self.joint_genotype_likelihood = None

    def deduplicate_samples(self):
        """
        Deduplicates samples based on the sample names by merging their pileups, segments, and contaminations.

        This method removes duplicate samples from the input data based on the sample
        names. It issues a warning if duplicate samples are found.
        """
        duplicate_samples = [name for name, count in Counter(self.samples).items() if count > 1]
        if not duplicate_samples:
            return None

        if self.verbose:
            print(f"  Duplicate samples found: {duplicate_samples}.")
            print(f"  Merging pileups, segments, and contaminations for duplicate samples.\n")

        samples = []
        pileups = []
        contaminations = []
        segments = []
        for s, p, c, seg in zip(self.samples, self.pileups, self.contaminations, self.segments):
            if s in duplicate_samples:
                continue
            samples.append(s)
            pileups.append(p)
            contaminations.append(c)
            segments.append(seg)

        for s in duplicate_samples:
            s_idx = [i for i, sample in enumerate(self.samples) if sample == s]
            samples.append(s)
            pileups.append(join_pileups([self.pileups[i] for i in s_idx]))
            contaminations.append(join_contaminations([self.contaminations[i] for i in s_idx]))
            segments.append(harmonize_segments([self.segments[i] for i in s_idx]))

        self.samples = samples
        self.pileups = pileups
        self.contaminations = contaminations
        self.segments = segments

    def subset_data(self):
        """
        Subsets the data to the intersection of all pileup and variant input files.

        Aligns and trims the data to ensure it only contains information relevant
        to both the pileup and variant data.
        """
        pileup_loci = pd.MultiIndex.from_frame(df=pd.concat([p.df[["contig", "position"]] for p in self.pileups]).drop_duplicates())
        index = pileup_loci if self.vcf.df.empty else pileup_loci.join(self.vcf.df.index, how="inner")
        message(f"Subset data to {len(index)} loci (aligned to reference variant input vcf).") if self.verbose else None
        for s, p in zip(self.samples, self.pileups):
            n_loci = p.df.shape[0]
            _p = p.df.set_index(["contig", "position"])
            mask = _p.index.isin(index)
            p.df = p.df.loc[mask]
            message(f"Subset pileups for {s} to {p.df.shape[0]}/{n_loci}.")

    def get_pileup_likelihood(self, genotyper: Genotyper, assigned_sample_name: str, pileup: Pileup, contamination: Contamination, segments: Segments, af_parameters: AFParameters) -> PileupLikelihood:
        pileup_likelihood = PileupLikelihood(
            pileup_likelihood=genotyper.calculate_genotype_likelihoods(
                pileup=pileup.df,
                contamination=contamination.value,
                segments=segments.df,
                ref_bias=af_parameters.ref_bias
            ),
            assigned_sample_name=assigned_sample_name,
            bam_sample_name=cross_check_sample_name(pileup, contamination, segments)
        )
        print(".", end="", flush=True) if self.verbose else None
        return pileup_likelihood

    def get_pileup_likelihoods_parallel(self, genotyper: Genotyper, max_processes: int = 1) -> list[PileupLikelihood]:
        processes = np.min([max_processes, len(self.samples), os.cpu_count()])
        if processes > 1:
            try:
                message(f"Parallel execution with {processes} processes over samples: ", end="", flush=True) if self.verbose else None
                with mp.Pool(processes=processes) as pool:
                    pileup_likelihoods = pool.starmap(
                        func=self.get_pileup_likelihood,
                        iterable=[
                            (genotyper, assigned_sample_name, pileup, contamination, segments, af_parameters)
                            for assigned_sample_name, pileup, contamination, segments, af_parameters in zip(self.samples, self.pileups, self.contaminations, self.segments, self.af_parameters)
                        ]
                    )
                print() if self.verbose else None
                return pileup_likelihoods

            except Exception as e:
                print(f"Error in parallel execution: {e}")
        message("Serial execution over samples: ", end="", flush=True) if self.verbose else None
        pileup_likelihoods = [
            self.get_pileup_likelihood(genotyper, assigned_sample_name, pileup, contamination, segments, af_parameters)
            for assigned_sample_name, pileup, contamination, segments, af_parameters in zip(self.samples, self.pileups, self.contaminations, self.segments, self.af_parameters)
        ]
        print("") if self.verbose else None
        return pileup_likelihoods

    def validate_pileup_likelihood_allele_frequencies(self):
        """
        Validates the consistency of population allele frequencies in pileup likelihoods.

        Checks if allele frequencies are the same across all samples for each locus
        and issues a warning if inconsistencies are found. (Allele frequencies should
        be the same across samples for each locus as they are derived from the
        SNP panel used for generating the allelic pileups.)

        Raises:
            Warning: If allele frequencies are not consistent across samples.
        """
        message("Validating pileup likelihood allele frequencies ...") if self.verbose else None
        allele_frequency = pd.concat([pl.df["allele_frequency"] for pl in self.pileup_likelihoods], axis=1)
        is_same = allele_frequency.apply(lambda row: row.nunique() == 1, axis=1)
        if not is_same.all():
            warnings.warn("Allele frequencies in pileups are not the same across samples for each locus (though they should)! Please check your inputs!")

    def validate_sample_correlation(self, genotyper: Genotyper, correlation_threshold: float = 0.90) -> pd.DataFrame:
        """
        Validates the correlation of sample genotypes.

        Calculates the Kendall-tau correlation of variant genotypes between samples
        and warns if the correlation is below a specified threshold.

        Median correlation of samples between unrelated patients: 0.15 (from data)
        Median correlation of samples between 1st-degree related patients: 0.575
          0.15 = 0.575 - 0.425, so two unrelated individuals share 57.5% of their variant SNPs
          children inherit 50% of their SNPs from each parent, so two related individuals share 78.75% of their variant SNPs
          Kendall tau correlation 0.575 = 0.7875 - (1 - 0.7875)
        Median correlation of samples within same patient: 0.978 (from data)

        Args:
            genotyper (Genotyper): Genotyper object.
            correlation_threshold (float, optional): Threshold for genotype correlation to raise warning.

        Returns:
            pandas.DataFrame: DataFrame containing the sample correlation matrix.

        Raises:
            Warning: If genotype correlations between samples are below the threshold.
        """
        message("Evaluating genotype concordance ...") if self.verbose else None

        sample_gt_list = []
        for pl in self.pileup_likelihoods:
            confident_calls = genotyper.select_confident_calls(likelihoods=pl.df)
            if not confident_calls.empty:
                sample_gt_list.append(confident_calls[genotyper.genotypes].idxmax(axis=1).to_frame(pl.assigned_sample_name))
            else:
                sample_gt_list.append(pd.DataFrame(index=confident_calls.index, columns=[pl.assigned_sample_name]))

        pd.set_option("future.no_silent_downcasting", True)
        sample_genotypes = (
            pd.concat(sample_gt_list, join="outer", axis=1)
            .fillna("./.")
            .replace("./.", np.nan)
            .replace("0/0", 0)
            .replace("0/1", 1)
            .replace("1/1", 2)
        )

        if sample_genotypes.empty:
            return pd.DataFrame(columns=[pl.assigned_sample_name for pl in self.pileup_likelihoods])

        sample_names = sample_genotypes.columns
        sample_gt = [sample_genotypes[[sample]] for sample in sample_names]
        sample_n_loc = [sample_genotypes[sample].dropna().gt(0).sum() for sample in sample_names]

        sample_pairs = list(combinations(sample_names, 2))
        gt_pairs = list(combinations(sample_gt, 2))

        correlations = []
        pvalues = []
        num_loci = []
        for gt_pair in gt_pairs:
            gt = gt_pair[0].join(gt_pair[1], how="inner")
            sample_a = gt_pair[0].columns[0]
            sample_b = gt_pair[1].columns[0]
            # Subset to loci with at least one alternate allele in at least one sample,
            gt = gt.loc[(gt[sample_a] > 0) | (gt[sample_b] > 0)].dropna()
            corr = st.kendalltau(gt[sample_a].to_numpy(), gt[sample_b].to_numpy())
            correlations.append(corr.statistic)
            pvalues.append(corr.pvalue)
            num_loci.append(gt.shape[0])

        if len(sample_pairs) == 0:
            corr = pd.DataFrame([[1]], columns=sample_names, index=sample_names)
            pval = pd.DataFrame([[0]], columns=sample_names, index=sample_names)
            nloc = pd.DataFrame([[0]], columns=sample_names, index=sample_names)
        else:
            corr = pd.Series(
                {sample_pair: corr for sample_pair, corr in zip(sample_pairs, correlations)}
            ).unstack().reindex(index=sample_names, columns=sample_names)
            pval = pd.Series(
                {sample_pair: pval for sample_pair, pval in zip(sample_pairs, pvalues)}
            ).unstack().reindex(index=sample_names, columns=sample_names)
            nloc = pd.Series(
                {sample_pair: nloc for sample_pair, nloc in zip(sample_pairs, num_loci)}
            ).unstack().reindex(index=sample_names, columns=sample_names)

        sample_correlation = corr.where(~corr.isna(), corr.T)
        pval = pval.where(~pval.isna(), pval.T)
        nloc = nloc.where(~nloc.isna(), nloc.T)
        np.fill_diagonal(sample_correlation.values, 1)
        np.fill_diagonal(pval.values, 0)
        np.fill_diagonal(nloc.values, sample_n_loc)
        sample_correlation = sample_correlation.fillna(0)
        pval = pval.fillna(1)
        nloc = nloc.fillna(0).astype(int)

        if self.verbose:
            print("Kendall-tau correlation of variant genotypes between samples:")
            print(sample_correlation.round(3).to_string())
            print()
            print("Number of variant loci for evaluating concordance of samples:")
            print(nloc.to_string())
            print()
            print(f"p-values of Kendall-tau correlation of variant genotypes between samples:")
            print(pval.round(3).to_string())
            print()

        if sample_correlation[pval < 0.05].lt(correlation_threshold).any().any():
            warning_message = (
                f"\n\n"
                f"\tGenotypes between samples are not well correlated (< {correlation_threshold})!\n"
                f"\tIt is likely that not all samples are from the same patient.\n"
                f"\tPlease check your inputs!\n"
            )
            warnings.warn(warning_message)
        return sample_correlation

    def sort_genomic_positions(self) -> None:
        """
        Sorts a MultiIndex by genomic positions ('contig' and 'position' levels).

        Genomic contigs are sorted according to header in variants VCF.
        Positions within contigs are sorted numerically.

        Aligns each sample's pileup likelihood with the joint GT likelihood index.
        """
        # contig_order = self.vcf.contigs
        message("Sorting loci ... ")

        temp_df = pd.DataFrame(index=self.joint_genotype_likelihood.index).reset_index()
        contigs = self.vcf.contigs
        # Add all remaining contigs that are not in vcf reference
        contigs += list(sorted(set(temp_df["contig"]) - set(self.vcf.contigs)))
        temp_df["contig"] = pd.Categorical(temp_df["contig"], categories=contigs, ordered=True)
        temp_df.sort_values(by=["contig", "position"], inplace=True)
        sorted_index = pd.MultiIndex.from_frame(temp_df[["contig", "position"]])

        self.joint_genotype_likelihood = self.joint_genotype_likelihood.loc[sorted_index]
        self.vcf.df = self.vcf.df.loc[sorted_index]
        for pl in self.pileup_likelihoods:
            pl.df = pl.df.reindex(sorted_index)

    def phase_hets(self, phaser: Phaser):
        self.joint_genotype_likelihood = phaser.phase_hets(
            joint_genotype_likelihood=self.joint_genotype_likelihood,
            pileup_likelihoods=self.pileup_likelihoods,
            segments=self.segments
        )

    def write_output(self, output_dir: str, vcf_format: str = "GT", args: argparse.Namespace = None):
        """
        Writes the processed data to output files.

        Generates and saves various output files based on the processed data,
        including (optional) sample genotype likelihoods and VCF files.

        Args:
            output_dir (str): Directory where the output files will be saved.
            vcf_format (str, optional): Format of the VCF FORMAT output.
            args (argparse.Namespace, optional): Arguments passed to the script, used for additional information in the output.
        """
        message(f"Writing output to {output_dir} ...") if args.verbose else None
        os.makedirs(output_dir, exist_ok=True)
        self.write_sample_correlation(output_dir=output_dir, args=args)
        self.write_sample_likelihoods(output_dir=output_dir, args=args)
        self.write_counts_tables(output_dir=output_dir, args=args)
        self.write_genotype_vcf(output_dir=output_dir, vcf_format=vcf_format, args=args)

    def write_sample_correlation(self, output_dir: str, args: argparse.Namespace = None):
        file = f"{output_dir}/{self.individual_id}.sample_correlation.tsv" + (".gz" if args.compress_output else "")
        self.sample_correlation.to_csv(file, sep="\t")
        print(f"  {file}") if args.verbose else None
        file = f"{output_dir}/{self.individual_id}.sample_correlation.min.txt"
        np.savetxt(file, [np.min(self.sample_correlation.values)])
        print(f"  {file}") if args.verbose else None

    def write_sample_likelihoods(self, output_dir: str, args: argparse.Namespace = None):
        if args.save_sample_genotype_likelihoods:
            for pl in self.pileup_likelihoods:
                sample_file = f"{output_dir}/{pl.assigned_sample_name}.likelihoods.pileup" + (".gz" if args.compress_output else "")
                open_func = gzip.open if args.compress_output else open
                with open_func(sample_file, "wt") as output_file:
                    output_file.write(f"#<METADATA>SAMPLE={pl.bam_sample_name}\n")
                    pl.df.fillna(0).astype({"ref_count": int, "alt_count": int, "other_alt_count": int}).to_csv(output_file, sep="\t")
                print(f"  {sample_file}") if args.verbose else None

    def write_counts_tables(self, output_dir: str, args: argparse.Namespace = None):
        for count in ["ref_count", "alt_count", "other_alt_count"]:
            table = pd.concat(
                [pl.df[count].to_frame(pl.assigned_sample_name).fillna(0) for pl in self.pileup_likelihoods],
                axis=1
            ).fillna(0).astype(int)
            file = f"{output_dir}/{self.individual_id}.germline.{count}.tsv" + (".gz" if args.compress_output else "")
            table.to_csv(file, sep="\t", index=False)
            print(f"  {file}") if args.verbose else None

    def write_genotype_vcf(self, output_dir: str, vcf_format: str = "GT", args: argparse.Namespace = None):
        """
        Converts the genotype data to VCF format and writes to a file.

        Args:
            output_dir (str): Directory where the output files will be saved.
            vcf_format (str, optional): The format for the VCF FORMAT output.
            args (argparse.Namespace, optional): Arguments namespace containing additional information for VCF header.
        """
        file = f"{output_dir}/{self.individual_id}.germline.vcf" + (".gz" if args.compress_output else "")

        df = self.joint_genotype_likelihood

        if not df.empty:
            df["id"] = "."
            df["ref"] = self.vcf.df["REF"] if not self.vcf.df.empty else "N"
            df["alt"] = self.vcf.df["ALT"] if not self.vcf.df.empty else "N"
            df["qual"] = "."
            df["filter"] = "."
            df["info"] = df["allele_frequency"].apply(lambda af: f"AF={af:.6f}")
            df["format"] = vcf_format

            # create fields for "SAMPLE" column:
            # df["GT"] = df[['0/0', '0/1', '1/1', './.']].astype(float).idxmax(axis=1)
            if "AD" in vcf_format:
                df["AD"] = df[['ref_count', 'alt_count']].astype(int).apply(lambda row: f"{row['ref_count']},{row['alt_count']}", axis=1)
            if "DP" in vcf_format:
                df["DP"] = df[['ref_count', 'alt_count', 'other_alt_count']].sum(axis=1).astype(int).astype(str)
            if "PL" in vcf_format:
                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore", category=RuntimeWarning)
                    phred = np.clip(-10 * np.log(df[['0/0', '0/1', '1/1', './.']]), 0, 99).astype(int)
                pl_min = phred.min(axis=1)
                df["PL"] = phred.sub(pl_min, axis=0).apply(lambda row: ",".join([str(p) for p in row]), axis=1)
            df["sample"] = df[vcf_format.split(":")].apply(lambda row: ":".join([str(f) for f in row]), axis=1)
            df = df.dropna()

        open_func = gzip.open if args.compress_output else open
        with open_func(file, "wt") as vcf:
            vcf.write("##fileformat=VCFv4.2\n")
            vcf.write("##FILTER=<ID=PASS,Description=\"All filters passed\">\n")
            vcf.write("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency, for each ALT allele, in the same order as listed\">\n")
            if "AD" in vcf_format:
                vcf.write("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n")
            if "DP" in vcf_format:
                vcf.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth; some reads may have been filtered\">\n")
            if "GQ" in vcf_format:
                vcf.write("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n")
            vcf.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
            if "PL" in vcf_format:
                vcf.write("##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">\n")
            for ref_line in self.vcf.ref_dict:
                vcf.write(f"{ref_line}\n")
            source_comment = "##source=genotype.py"
            if args is not None:
                for key, value in vars(args).items():
                    source_comment += f" --{key} {value}"
            vcf.write(f"{source_comment}\n")
            vcf.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{self.individual_id}\n")
            if not df.empty:
                df.reset_index()[["contig", "position", "id", "ref", "alt", "qual", "filter", "info", "format", "sample"]].to_csv(vcf, sep="\t", index=False, header=False, na_rep=".")

        print(f"  {file}") if args.verbose else None


if __name__ == "__main__":
    main()
