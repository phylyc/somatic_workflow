import argparse
from functools import reduce
from itertools import combinations
import gzip
import multiprocessing as mp
import numpy as np
import os
import pandas as pd
import scipy.stats as st
import warnings


def parse_args():
    parser = argparse.ArgumentParser(
        prog="Genotyper",
        description="""
            This script calculates genotype likelihoods from allelic counts, taking into account potential sample contamination and segment-specific minor allele frequency. 
            It uses allelic pileup data (GATK GetPileupSummaries), segments data (GATK CalculateContamination), and contamination estimates (GATK CalculateContamination), 
            to compute genotype likelihoods under specified models (binomial or beta-binomial). 
            Variant genotype correlation between all samples is computed. 
            Outputs 1) a VCF with patient genotype information, 2) allelic count matrices for each allele, 3) sample correlation matrix, and (optionally) 4) sample pileups with genotype likelihood annotations.
        """,
        epilog="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.usage = "genotype.py -O <output_dir> -V <variant> --patient <patient> --sample <sample> [--sample <sample> ...]-P <pileup> [-P <pileup> ...] [-S <segments> ...] [-C <contamination> ...] [-M <model>] [-D <min_read_depth>] [-p <min_genotype_likelihood>] [-F <format>] [--select_hets] [--save_sample_genotype_likelihoods] [--verbose]"
    parser.add_argument("-O", "--output_dir",                   type=str,   required=True,  help="Path to the output directory.")
    parser.add_argument("--patient",                            type=str,   required=True,  help="Name of the patient.")
    parser.add_argument("--sample",                             type=str,   required=True,  action="append",    help="Assigned name of the sample. Does not have to coincide with the sample name in the pileup file.")
    parser.add_argument("-P", "--pileup",                       type=str,   required=True,  action="append",    help="Path to the allelic pileup input file (output of GATK's GetPileupSummaries; used for CalculateContamination).")
    parser.add_argument("-S", "--segments",                     type=str,   default=None,   action="append",    help="Path to the allelic copy ratio segmentation input file (output of GATK's CalculateContamination).")
    parser.add_argument("-C", "--contamination",                type=str,   default=None,   action="append",    help="Path to the contamination estimate input file (output of GATK's CalculateContamination).")
    parser.add_argument("-V", "--variant",                      type=str,   default=None,   help="A VCF file containing common germline variants and population allele frequencies. (Used to get pileup summaries.)")
    parser.add_argument("-M", "--model",                        type=str,   default="betabinom", choices=["binom", "betabinom"], help="Genotype likelihood model.")
    parser.add_argument("-D", "--min_read_depth",               type=int,   default=10,     help="Minimum read depth per sample to consider site for genotyping.")
    parser.add_argument("-p", "--min_genotype_likelihood",      type=float, default=0.95,   help="Probability threshold for calling and retaining genotypes.")
    parser.add_argument("-F", "--format",                       type=str,   default="GT",   help="VCF format field. (GT: genotype; AD: allele depth; DP: total depth; PL: phred-scaled genotype likelihoods.)")
    parser.add_argument("--threads",                            type=int,   default=1,      help="Number of threads to use for parallelization.")
    parser.add_argument("--select_hets",                                    default=False,  action="store_true", help="Keep only heterozygous sites.")
    parser.add_argument("--save_sample_genotype_likelihoods",               default=False,  action="store_true", help="Save genotype likelihoods to file for each sample.")
    parser.add_argument("--verbose",                                        default=False,  action="store_true", help="Print information to stdout during execution.")
    return parser.parse_args()


def main():
    args = parse_args()
    print_args(args)
    genotyper = Genotyper(
        model=args.model,
        min_genotype_likelihood=args.min_genotype_likelihood,
        select_hets=args.select_hets
    )
    data = GenotypeData(
        individual_id=args.patient,
        samples=args.sample,
        variant=args.variant,
        pileup=args.pileup,
        segments=args.segments,
        contamination=args.contamination,
        min_read_depth=args.min_read_depth,
    )
    data.pileup_likelihoods = data.get_pileup_likelihoods_parallel(genotyper=genotyper, max_processes=args.threads)
    data.validate_pileup_likelihood_allele_frequencies()
    data.sample_correlation = data.validate_sample_correlation(verbose=args.verbose)
    data.joint_genotype_likelihood = genotyper.get_joint_genotype_likelihood(pileup_likelihoods=data.pileup_likelihoods)
    data.subset_data()
    data.write_output(output_dir=args.output_dir, vcf_format=args.format, args=args)


def print_args(args):
    if args.verbose:
        print("Calling Genotyper")
        print("Arguments:")
        for key, value in vars(args).items():
            print(f"  {key}: {value}")


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

    def __init__(self, file_path: str = None, min_read_depth: int = 0):
        self.file_path = file_path
        self.columns = ["contig", "position", "ref_count", "alt_count", "other_alt_count", "allele_frequency"]
        self.df = (
            pd.read_csv(file_path, sep="\t", comment="#", header=0, names=self.columns, low_memory=False)
            if file_path is not None
            else pd.DataFrame(columns=self.columns)
        )
        self.df = self.df.loc[self.df[["ref_count", "alt_count", "other_alt_count"]].sum(axis=1) >= min_read_depth]
        self.bam_sample_name = None
        if file_path is not None:
            with open(file_path, "r") as pileup_file:
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
        self.df["./."].fillna(1.0, inplace=True)  # Declare missing loci as outliers.
        self.df.fillna(0, inplace=True)  # Set all other likelihoods and read counts to 0.
        self.df = self.df.astype({"ref_count": int, "alt_count": int, "other_alt_count": int})
        return self


class Segments(object):
    """
    Represents segmentation data from a specified file. (E.g. output of GATK CalculateContamination.)

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
        self.columns = ["contig", "start_pos", "end_pos", "minor_allele_fraction"]
        self.df = (
            pd.read_csv(file_path, sep="\t", comment="#", header=0, names=self.columns, low_memory=False)
            if file_path is not None
            else pd.DataFrame(columns=self.columns)
        )
        self.bam_sample_name = None
        if file_path is not None:
            with open(file_path, "r") as segments_file:
                segments_header = segments_file.readline().strip()
                if segments_header.startswith("#<METADATA>SAMPLE="):
                    self.bam_sample_name = segments_header.removeprefix("#<METADATA>SAMPLE=")


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

    def __init__(self, file_path: str = None):
        self.file_path = file_path
        self.keys = ["sample", "contamination", "contamination_error"]
        self.dict = (
            pd.read_csv(file_path, sep="\t", comment="#", header=0, names=self.keys, low_memory=False).loc[0].to_dict()
            if file_path is not None
            else {key: default for key, default in zip(self.keys, [None, 0, 0])}
        )
        self.value = self.dict["contamination"]
        self.error = self.dict["contamination_error"]
        self.bam_sample_name = self.dict["sample"]


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

    def __init__(self, file_path: str = None):
        self.file_path = file_path
        self.ref_dict = self.get_vcf_sequence_dict()
        self.df = (
            pd.read_csv(file_path, sep="\t", comment="#", header=None, low_memory=False, names=self.columns)
            if file_path is not None
            else pd.DataFrame(columns=self.columns)
        )
        self.df = self.df.astype({"CHROM": str, "POS": int, "ID": str, "REF": str, "ALT": str, "INFO": str})
        # GetPileupSummaries only supports SNVs well.
        # CAUTION: the input vcf may contain multi-allelic sites and other stuff!!!! (removed here)
        snv_mask = (self.df["REF"].apply(len) == 1) & (self.df["ALT"].apply(len) == 1) & ~self.df["REF"].isin(["-", "."]) & ~self.df["ALT"].isin(["-", "."])
        self.df = self.df.loc[snv_mask].drop_duplicates(subset=["CHROM", "POS"], keep=False)
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
        open_func = gzip.open if self.file_path.endswith(".gz") else open
        with open_func(self.file_path, "rt") as ifile:
            for line in ifile:
                if line.startswith("##contig"):
                    ref_dict.append(line.strip())
                if line.startswith("#CHROM"):
                    break
        ifile.close()
        return ref_dict


class Genotyper(object):
    """
    Performs genotype likelihood calculations using specified statistical models.

    Based on the implementation from the Broad Institute's GATK ContaminationModel.
    https://github.com/broadinstitute/gatk/blob/master/src/main/java/org/broadinstitute/hellbender/tools/walkers/contamination/ContaminationModel.java

    Attributes:
        model (str): The model to use for calculating genotype likelihoods ('binom' or 'betabinom').
        min_genotype_likelihood (float): Minimum genotype likelihood threshold for filtering.
        select_hets (bool): Whether to select only heterozygous sites.

    Args:
        model (str): The model for calculating genotype likelihoods.
        min_genotype_likelihood (float, optional): Minimum genotype likelihood threshold.
        select_hets (bool, optional): Flag to select only heterozygous sites.
    """

    def __init__(self, model: str, min_genotype_likelihood: float = 0.95, select_hets: bool = True):
        self.model = model
        self.min_genotype_likelihood = min_genotype_likelihood
        self.select_hets = select_hets

    @staticmethod
    def get_error_prob(pileup: pd.DataFrame, min_error: float = 1e-6, max_error: float = 1e-1) -> float:
        """
        Estimates the error probability based on the pileup data.

        The error probability is calculated as a ratio of 'other_alt_count' to the
        total counts, clipped within specified bounds.

        Args:
            pileup (pd.DataFrame): DataFrame containing pileup data.
            min_error (float, optional): Minimum bound for the error probability.
            max_error (float, optional): Maximum bound for the error probability.

        Returns:
            float: Estimated error probability.
        """
        other_alt_counts = pileup["other_alt_count"].sum()
        total_counts = pileup[["ref_count", "alt_count", "other_alt_count"]].sum(axis=1).sum()
        return np.clip(1.5 * other_alt_counts / min(1, total_counts), a_min=min_error, a_max=max_error)

    def calculate_genotype_likelihoods(self, pileup: pd.DataFrame, segments: pd.DataFrame = None, contamination: float = 0) -> pd.DataFrame:
        """
        Calculates genotype likelihoods for given pileup data, optionally segmented by genomic regions.

        Args:
            pileup (pd.DataFrame): DataFrame containing pileup data.
            segments (pd.DataFrame, optional): DataFrame containing segmentation data.
            contamination (float, optional): Contamination level to consider in calculations.

        Returns:
            pd.DataFrame: DataFrame with genotype likelihoods for each genomic position.
        """
        if pileup.empty:
            return pd.DataFrame(columns=["contig", "position", "ref_count", "alt_count", "other_alt_count", "allele_frequency", "0/0", "0/1", "1/1", "./."])

        error = self.get_error_prob(pileup=pileup)

        likelihoods = []
        in_any_seg = False
        for _, seg in segments.iterrows():
            seg_mask = (pileup["contig"] == seg["contig"]) & (seg["start_pos"] <= pileup["position"]) & (pileup["position"] <= seg["end_pos"])
            in_any_seg |= seg_mask
            seg_pileup = pileup.loc[seg_mask]
            if seg_pileup.empty:
                continue
            likelihoods.append(
                self.calculate_genotype_likelihoods_per_segment(
                    pileup=seg_pileup,
                    minor_af=seg["minor_allele_fraction"],
                    contamination=contamination,
                    error=error
                )
            )
        seg_pileup = pileup.loc[~in_any_seg]
        if not seg_pileup.empty:
            likelihoods.append(
                self.calculate_genotype_likelihoods_per_segment(
                    pileup=seg_pileup,
                    contamination=contamination,
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
                "0/1": float,
                "1/1": float,
                "./.": float,
            }
        )
        return genotype_likelihoods

    def calculate_genotype_likelihoods_per_segment(self, pileup: pd.DataFrame, minor_af: float = 0.4434, contamination: float = 0, error: float = 0.01) -> pd.DataFrame:
        """
        Calculates genotype likelihoods for a specific segment of pileup data.

        This method applies either a binomial or beta-binomial model to compute the
        likelihoods, based on the model specified during the initialization of the
        Genotyper.

        Note: Median minor allele frequency as inferred by the GATK's CalculateContamination
              is 0.4434 across a thousand normal samples.

        Args:
            pileup (pd.DataFrame): DataFrame containing pileup data for the segment.
            minor_af (float, optional): Minor allele frequency, default is 0.5.
            contamination (float, optional): Level of contamination to consider in the calculation.
            error (float, optional): Error probability to be used in the calculation.

        Returns:
            pd.DataFrame: DataFrame containing calculated genotype likelihoods for the segment.
        """

        n = pileup[["ref_count", "alt_count", "other_alt_count"]].sum(axis=1)
        alt = pileup["alt_count"]
        f = pileup["allele_frequency"]  # population allele frequency from SNP panel

        f_aa = contamination * f + (1 - contamination) * error / 3  # / 3 is in GATK code
        f_ab = contamination * f + (1 - contamination) * minor_af
        f_ba = contamination * f + (1 - contamination) * (1 - minor_af)
        f_bb = contamination * f + (1 - contamination) * (1 - error)
        s = 100  # overdispersion
        lam = 1.1  # SNP bias of the ref allele

        def bias(_f):
            return _f / (_f + (1 - _f) * lam)

        # Calculate log-likelihoods:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=RuntimeWarning)
            log1f = np.log1p(-f)
            logf = np.log(f)

        if self.model == "binom":
            aa = 2 * log1f + st.binom(n, f_aa).logpmf(alt)
            ab = log1f + logf + st.binom(n, f_ab).logpmf(alt)
            ba = log1f + logf + st.binom(n, f_ba).logpmf(alt)
            bb = 2 * logf + st.binom(n, f_bb).logpmf(alt)
            outlier = -np.inf
        elif self.model == "betabinom":
            aa = 2 * log1f + st.binom(n, bias(f_aa)).logpmf(alt)
            ab = log1f + logf + st.betabinom(n, bias(f_ab) * s, (1 - bias(f_ab)) * s).logpmf(alt)
            ba = log1f + logf + st.betabinom(n, bias(f_ba) * s, (1 - bias(f_ba)) * s).logpmf(alt)
            bb = 2 * logf + st.binom(n, bias(f_bb)).logpmf(alt)
            # outlier = st.betabinom(N, 1, 1).logpmf(alt)  # max entropy; todo: needs prior / weight
            outlier = -np.inf
        else:
            raise ValueError(f"model is {self.model} but has to be one of binom, betabinom.")

        het = np.logaddexp(ab, ba)
        total = reduce(np.logaddexp, [aa, het, bb, outlier])

        likelihoods = pileup.copy()
        likelihoods["0/0"] = np.exp(aa - total)
        likelihoods["0/1"] = np.exp(het - total)
        likelihoods["1/1"] = np.exp(bb - total)
        likelihoods["./."] = np.exp(outlier - total)

        return likelihoods

    def get_joint_genotype_likelihood(self, pileup_likelihoods: list[PileupLikelihood]) -> pd.DataFrame:
        """
        Calculates joint genotype likelihoods across multiple samples from the
        same individual.

        This method aggregates genotype likelihoods from multiple PileupLikelihood
        instances, normalizes them, and applies filters based on the minimum genotype
        likelihood threshold and selection criteria for heterozygous sites.

        Args:
            pileup_likelihoods (list[PileupLikelihood]): List of PileupLikelihood objects from different samples.

        Returns:
            pd.DataFrame: DataFrame with joint genotype likelihoods for each genomic position.
        """
        # Calculate the joint genotype likelihoods.
        genotypes = ["0/0", "0/1", "1/1", "./."]
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=RuntimeWarning)
            joint_log_likelihood = {
                gt: pd.concat([np.log(pl.df[gt]) for pl in pileup_likelihoods], axis=1).apply(np.nansum, axis=1)
                for gt in genotypes
            }
        # Normalize likelihoods
        total = pd.concat(joint_log_likelihood.values(), axis=1).apply(lambda row: reduce(np.logaddexp, row), axis=1)
        genotype_likelihood = {
            gt: np.exp(joint_log_likelihood[gt].sub(total, axis=0))
            for gt in genotypes
        }

        # Count total allelic read counts across all samples.
        genotype_likelihood |= {
            count: pd.concat([pl.df[count] for pl in pileup_likelihoods], axis=1).sum(axis=1)
            for count in ["ref_count", "alt_count", "other_alt_count"]
        }

        # Average allele frequencies across all samples. They should be the same.
        genotype_likelihood["allele_frequency"] = pd.concat(
            [pl.df["allele_frequency"] for pl in pileup_likelihoods],
            axis=1
        ).apply(np.nanmean, axis=1)

        df = pd.DataFrame.from_dict(genotype_likelihood)

        # subset to confident calls
        df = df.loc[(df[genotypes] >= self.min_genotype_likelihood).any(axis=1)]
        if self.select_hets:
            df = df.loc[df["0/1"] >= self.min_genotype_likelihood]

        return df


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
        segments (list[Segments]): List of Segments objects.
        contaminations (list[Contamination]): List of Contamination objects.
        pileup_likelihoods (list[PileupLikelihood]): List of PileupLikelihood objects.
        sample_correlation (pandas.DataFrame): DataFrame containing sample genotype correlation matrix.
        joint_genotype_likelihood (pandas.DataFrame): DataFrame containing joint genotype likelihoods.

    Args:
        individual_id (str): Identifier for the individual.
        samples (list[str]): List of sample names.
        variant (str): Path to the VCF file containing common germline variants and population allele frequencies.
        pileup (list[str]): List of paths to the allelic pileup input files (output of GATK's GetPileupSummaries; used for CalculateContamination).
        segments (list[str], optional): List of paths to the allelic copy ratio segmentation input files (output of GATK's CalculateContamination).
        contamination (list[str], optional): List of paths to the contamination estimate input files (output of GATK's CalculateContamination).
        min_read_depth (int, optional): Minimum read depth threshold for filtering data.
    """

    def __init__(self, individual_id: str, samples: list[str], variant: str, pileup: list[str], segments: list[str] = None, contamination: list[str] = None, min_read_depth: int = 10):
        self.individual_id = individual_id
        self.samples = samples
        self.vcf = VCF(file_path=variant)
        self.pileups = [Pileup(file_path=p, min_read_depth=min_read_depth) for p in pileup]
        self.segments = [Segments(file_path=s) for s in segments] if segments is not None else [Segments() for _ in pileup]
        self.contaminations = [Contamination(file_path=c) for c in contamination] if contamination is not None else [Contamination() for _ in pileup]
        self.pileup_likelihoods = None
        self.sample_correlation = None
        self.joint_genotype_likelihood = None

    @staticmethod
    def get_pileup_likelihood(genotyper: Genotyper, assigned_sample_name: str, pileup: Pileup, segments: Segments, contamination: Contamination) -> PileupLikelihood:
        pileup_likelihood = PileupLikelihood(
            pileup_likelihood=genotyper.calculate_genotype_likelihoods(
                pileup=pileup.df,
                segments=segments.df,
                contamination=contamination.value,
            ),
            assigned_sample_name=assigned_sample_name,
            bam_sample_name=cross_check_sample_name(pileup, segments, contamination)
        )
        print(".", end="", flush=True)
        return pileup_likelihood

    def get_pileup_likelihoods_parallel(self, genotyper: Genotyper, max_processes: int = 1) -> list[PileupLikelihood]:
        processes = np.min([max_processes, len(self.samples), os.cpu_count()])
        try:
            print(f"Parallel execution with {processes} processes: ", end="")
            with mp.Pool(processes=processes) as pool:
                pileup_likelihoods = pool.starmap(
                    func=self.get_pileup_likelihood,
                    iterable=[
                        (genotyper, assigned_sample_name, pileup, segments, contamination)
                        for assigned_sample_name, pileup, segments, contamination in zip(self.samples, self.pileups, self.segments, self.contaminations)
                    ]
                )
            print()
            return pileup_likelihoods
        except Exception as e:
            print(f"Error in parallel execution: {e}")
            print("Falling back to serial execution.")
            return [
                self.get_pileup_likelihood(genotyper, assigned_sample_name, pileup, segments, contamination)
                for assigned_sample_name, pileup, segments, contamination in zip(self.samples, self.pileups, self.segments, self.contaminations)
            ]

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
        allele_frequency = pd.concat([pl.df["allele_frequency"] for pl in self.pileup_likelihoods], axis=1)
        is_same = allele_frequency.apply(lambda row: row.nunique() == 1, axis=1)
        if not is_same.all():
            warnings.warn("Allele frequencies in pileups are not the same across samples for each locus (though they should)! Please check your inputs!")

    def validate_sample_correlation(self, correlation_threshold: float = 0.95, verbose=False) -> pd.DataFrame:
        """
        Validates the correlation of sample genotypes.

        Calculates the Kendall-tau correlation of variant genotypes between samples
        and warns if the correlation is below a specified threshold.

        Args:
            correlation_threshold (float, optional): Threshold for genotype correlation.
            verbose (bool, optional): If True, prints detailed correlation information.

        Returns:
            pandas.DataFrame: DataFrame containing the sample correlation matrix.

        Raises:
            Warning: If genotype correlations between samples are below the threshold.
        """
        sample_genotypes = pd.concat(
            [
                pl.df[["0/0", "0/1", "1/1", "./."]].idxmax(axis=1).to_frame(pl.assigned_sample_name)
                for pl in self.pileup_likelihoods
            ],
            axis=1
        ).fillna("./.")

        if sample_genotypes.empty:
            return pd.DataFrame(columns=[pl.assigned_sample_name for pl in self.pileup_likelihoods])

        # Subset to loci with at least one alternate allele in at least one sample.
        sample_genotypes = sample_genotypes.loc[sample_genotypes.apply(lambda row: any(["1" in gt for gt in row.values]), axis=1)]
        sample_genotypes = sample_genotypes.replace("0/0", 0).replace("0/1", 1).replace("1/1", 2).replace("./.", np.nan)
        sample_correlation = sample_genotypes.corr(method="kendall")
        # sample_correlation.sort_values(by=sample_correlation.columns[0], axis=0, ascending=False, inplace=True)
        # sample_correlation = sample_correlation.reindex(sample_correlation.index, axis=1)
        if verbose:
            print(f"Concordance of samples is evaluated at {sample_genotypes.dropna().shape[0]} to {sample_genotypes.shape[0]} variant loci.")
            print("Kendall-tau correlation of variant genotypes between samples:")
            print(sample_correlation.round(3))
        if sample_correlation.lt(correlation_threshold).any().any():
            message = (
                f"\n\n"
                f"\tGenotypes between samples are not well correlated (< {correlation_threshold})!\n"
                f"\tIt is likely that not all samples are from the same individual.\n"
                f"\tPlease check your inputs!\n"
            )
            warnings.warn(message)
        return sample_correlation

    @staticmethod
    def sort_genomic_positions(index: pd.MultiIndex) -> pd.MultiIndex:
        """
        Sorts a MultiIndex by genomic positions ('contig' and 'position' levels).

        Genomic contigs are sorted numerically and then by 'X', 'Y', 'MT'. Positions
        within contigs are sorted numerically.

        Args:
            index (pd.MultiIndex): MultiIndex to be sorted, expected to have 'contig' and 'position' levels.

        Returns:
            pd.MultiIndex: Sorted MultiIndex.
        """
        contig_order = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
        temp_df = pd.DataFrame(index=index).reset_index()
        temp_df["contig"] = pd.Categorical(temp_df["contig"], categories=contig_order, ordered=True)
        temp_df.sort_values(by=["contig", "position"], inplace=True)
        return pd.MultiIndex.from_frame(temp_df[["contig", "position"]])

    def subset_data(self):
        """
        Subsets the data to the intersection of all pileup and variant input files.

        Aligns and trims the data to ensure it only contains information relevant
        to both the pileup and variant data.
        """
        index = (
            self.joint_genotype_likelihood.index.join(self.vcf.df.index, how="inner")
            if not self.vcf.df.empty
            else self.joint_genotype_likelihood.index
        )
        index = self.sort_genomic_positions(index=index)
        self.pileup_likelihoods = [pl.reindex(index) for pl in self.pileup_likelihoods]
        # If pileup summaries of one sample do not cover a locus in another pileup,
        # the allele frequency is set to 0. This is rescuing those allele frequencies:
        for pl in self.pileup_likelihoods:
            pl.df["allele_frequency"] = self.joint_genotype_likelihood["allele_frequency"]
        self.joint_genotype_likelihood = self.joint_genotype_likelihood.loc[index]
        if not self.vcf.df.empty:
            self.vcf.df = self.vcf.df.reindex(index)

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
        print(f"Writing output to {output_dir} ...") if args.verbose else None
        os.makedirs(output_dir, exist_ok=True)
        if args.save_sample_genotype_likelihoods:
            for pl in self.pileup_likelihoods:
                sample_file = f"{output_dir}/{pl.assigned_sample_name}.likelihoods.pileup"
                with open(sample_file, "w") as output_file:
                    output_file.write(f"#<METADATA>SAMPLE={pl.bam_sample_name}\n")
                    pl.df.to_csv(output_file, sep="\t")
                print(f"  {sample_file}") if args.verbose else None

        self.sample_correlation.to_csv(f"{output_dir}/{self.individual_id}.sample_correlation.tsv", sep="\t")
        print(f"  {output_dir}/{self.individual_id}.sample_correlation.tsv") if args.verbose else None

        for count in ["ref_count", "alt_count", "other_alt_count"]:
            pd.concat(
                [pl.df[count].to_frame(pl.assigned_sample_name) for pl in self.pileup_likelihoods],
                axis=1
            ).fillna(0).astype(int).to_csv(f"{output_dir}/{self.individual_id}.hets.{count}.tsv", sep="\t", index=False)
            print(f"  {output_dir}/{self.individual_id}.hets.{count}.tsv") if args.verbose else None

        self.to_vcf(output=f"{output_dir}/{self.individual_id}.hets.vcf", vcf_format=vcf_format, args=args)
        print(f"  {output_dir}/{self.individual_id}.hets.vcf") if args.verbose else None

    def to_vcf(self, output: str, vcf_format: str = "GT", args: argparse.Namespace = None):
        """
        Converts the genotype data to VCF format and writes to a file.

        Args:
            output (str): Path to the output VCF file.
            vcf_format (str, optional): The format for the VCF FORMAT output.
            args (argparse.Namespace, optional): Arguments namespace containing additional information for VCF header.
        """
        df = self.joint_genotype_likelihood

        df["id"] = "."
        df["ref"] = self.vcf.df["REF"] if not self.vcf.df.empty else "N"
        df["alt"] = self.vcf.df["ALT"] if not self.vcf.df.empty else "N"
        df["qual"] = "."
        df["filter"] = "."
        df["info"] = df["allele_frequency"].apply(lambda af: f"AF={af:.6f}")
        # df["info"] = self.variant.df["INFO"]  # may contain more than AF (e.g. AN, AC, AF, etc.)
        df["format"] = vcf_format

        # create fields for "SAMPLE" column:
        df["GT"] = df[['0/0', '0/1', '1/1', './.']].astype(float).idxmax(axis=1)
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

        with open(output, "w") as vcf:
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
            df.reset_index()[["contig", "position", "id", "ref", "alt", "qual", "filter", "info", "format", "sample"]].to_csv(vcf, sep="\t", index=False, header=False, na_rep=".")


if __name__ == "__main__":
    main()
