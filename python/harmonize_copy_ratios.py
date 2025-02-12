import argparse
from collections import defaultdict
from copy import copy
import gzip
import multiprocessing as mp
import numpy as np
import os
import pandas as pd
import time
from typing import Union
import warnings


def parse_args():
    parser = argparse.ArgumentParser(
        prog="HarmonizeCopyRatios",
        description="""
            
        """,
        epilog="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.usage = "harmonize_copy_ratios.py -I <copy_ratio> [-I <copy_ratio> ...] [-S <sample> ...] [-P <patient>] [-O <output_dir>] [-D <ref_dict>] [--normal_sample <name>] [--suffix <suffix>] [--threads <threads>] [--compress_output] [--verbose]"
    parser.add_argument("-I", "--copy_ratio",       type=str,   required=True,      action="append",        help="Path to the (denoised) copy ratio file (e.g. output of GATK's DenoiseReadCounts).")
    parser.add_argument("-S", "--sample",           type=str,                       action="append",        help="Assigned name of the sample paired to the input file. Does not have to coincide with the sample name in the file header.")
    parser.add_argument("-P", "--patient",          type=str,   default="Patient",                          help="Assigned name of the patient.")
    parser.add_argument("-O", "--output_dir",       type=str,   default=".",                                help="Path to the output directory.")
    parser.add_argument("-D", "--ref_dict",         type=str,                                               help="Path to the reference dictionary. This will be used as header instead of the copy_ratio input headers.")
    parser.add_argument("--normal_sample",          type=str,                                               help="Assigned name of the normal sample.")
    parser.add_argument("--suffix",                 type=str,   default=".harmonized.CR.tsv",               help="Suffix to append to each sample name to create the output file name.")
    parser.add_argument("--threads",                type=int,   default=1,                                  help="Number of threads to use for parallelization over chromosomes.")
    parser.add_argument("--min_target_length",      type=int,   default=100,                                help="Minimum length of a harmonized target to be included in output.")
    parser.add_argument("--merge_abutting_distance", type=int,  default=1,                                  help="Maximum distance of abutting intervals to merge.")
    parser.add_argument("--merge_abutting_tolerance", type=float,   default=1e-3,                           help="Maximum copy ratio difference between abutting intervals to merge.")
    parser.add_argument("--interval_set_rule",      type=str,   default="INTERSECTION",                     help="Output the INTERSECTION / UNION of all harmonized intervals across all samples.")
    parser.add_argument("--compress_output",                    default=False,      action="store_true",    help="Compress output files.")
    parser.add_argument("--column_names",           type=str,               nargs="*",                      help="Column names for input files")
    parser.add_argument("--column_types",           type=str,               nargs="*",                      help="Column types for input files. Must be the same number of arguments as for --column_names.")
    parser.add_argument("--contig_col",             type=str,   default="CONTIG",                           help="Name of the column containing the contig information.")
    parser.add_argument("--start_col",              type=str,   default="START",                            help="Name of the column containing the start position of the interval.")
    parser.add_argument("--end_col",                type=str,   default="END",                              help="Name of the column containing the end position of the interval.")
    parser.add_argument("--copy_ratio_col",         type=str,   default="LOG2_COPY_RATIO",                  help="Name of the column containing the log2 copy ratio.")
    parser.add_argument("--sample_col",             type=str,   default="SAMPLE",                           help="Name of the column containing the sample name.")
    parser.add_argument("--comment_char",           type=str,   default="@",                                help="Character that starts a comment line in the input files.")
    parser.add_argument("--agg_col",                type=str,   default=None,       action="append",        help="Name of column to aggregate over. Can be listed more than once.")
    parser.add_argument("--agg_func",               type=str,   default=None,       action="append",        help="Function to aggregate over column. Needs to be listed in same order as --agg_col.")
    parser.add_argument("--filter_germline_calls",              default=False,      action="store_true",    help="If --normal_sample is given, filter all CALL segments with + or - assignment.")
    parser.add_argument("--verbose",                            default=False,      action="store_true",    help="Print information to stdout during execution.")
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    print_args(args)
    h = Harmonizer.from_args(args)
    h.harmonize()


def print_args(args):
    if args.verbose:
        print("Calling HarmonizeCopyRatios")
        print("Arguments:")
        for key, value in vars(args).items():
            print(f"  {key}: {value}")
        print()


class Harmonizer(object):
    @staticmethod
    def from_args(args):
        return Harmonizer(
            copy_ratio_files=args.copy_ratio,
            sample_names=args.sample,
            normal_sample=args.normal_sample,
            patient=args.patient,
            output_dir=args.output_dir,
            suffix=args.suffix,
            threads=args.threads,
            min_target_length=args.min_target_length,
            compress_output=args.compress_output,
            column_names=args.column_names,
            column_types=args.column_types,
            contig_col=args.contig_col,
            start_col=args.start_col,
            end_col=args.end_col,
            sample_col=args.sample_col,
            comment_char=args.comment_char,
            ref_dict_file=args.ref_dict,
            agg_col=args.agg_col,
            agg_func=args.agg_func,
            interval_set_rule=args.interval_set_rule,
            merge_abutting_distance=args.merge_abutting_distance,
            merge_abutting_tolerance=args.merge_abutting_tolerance,
            filter_germline_calls=args.filter_germline_calls,
            verbose=args.verbose,
        )

    def __init__(
        self,
        copy_ratio_files: Union[tuple[str], list[str]],
        sample_names: Union[tuple[str], list[str]] = None,
        normal_sample: str = None,
        patient: str = "Patient",
        output_dir: str = ".",
        suffix: str = ".harmonized.CR",
        threads: int = 1,
        min_target_length: int = 100,
        compress_output: bool = False,
        column_names: Union[tuple[str], list[str]] = ("CONTIG", "START", "END", "LOG2_COPY_RATIO"),
        column_types: Union[tuple[str], list[str]] = (str, int, int, float),
        contig_col: str = "CONTIG",
        start_col: str = "START",
        end_col: str = "END",
        sample_col: str = "SAMPLE",
        comment_char: str = "@",
        ref_dict_file: str = None,
        agg_col: Union[tuple[str], list[str]] = None,
        agg_func: Union[tuple[str], list[str]] = None,
        interval_set_rule: str = "INTERSECTION",
        merge_abutting_distance: int = 1,
        merge_abutting_tolerance: float = 1e-3,
        filter_germline_calls: bool = False,
        verbose: bool = False,
    ):
        self.copy_ratio_files = copy_ratio_files
        self.headers = defaultdict(list)
        self.ref_dict_file = ref_dict_file
        self.ref_dict = None
        self.contigs = []
        self.copy_ratios = []
        self.harmonized_copy_ratios = None
        self.sample_names = sample_names
        self.normal_sample = normal_sample
        self.patient = patient
        self.output_dir = output_dir
        self.suffix = suffix
        self.threads = threads
        self.min_target_length = min_target_length
        self.compress_output = compress_output
        self.column_names = column_names
        self.column_types = {c: t for c, t in zip(column_names, column_types) if t is not None}
        self.contig_col = contig_col
        self.start_col = start_col
        self.end_col = end_col
        self.sample_col = sample_col
        self.comment_char = comment_char
        self.agg_col = agg_col if agg_col is not None else []
        self.agg_func = agg_func if agg_func is not None else []
        self.interval_set_rule = interval_set_rule
        self.merge_abutting_distance = merge_abutting_distance
        self.merge_abutting_tolerance = merge_abutting_tolerance
        self.filter_germline_calls = filter_germline_calls
        self.verbose = verbose

    def harmonize(self):
        self.load_data()
        self.harmonized_copy_ratios = self._harmonize()
        self.write_output_files()

    def load_data(self):
        def load_file(file_path, sample_name=None):
            header, df = self.get_header_and_df(file_path=file_path)
            if sample_name is not None:
                df[self.sample_col] = sample_name
            self.copy_ratios.append(df)
            print(f"Number of loci in {cr_file_path}: {df.shape[0]}") if self.verbose else None
            return header, df

        if self.sample_names is None:
            for cr_file_path in self.copy_ratio_files:
                header, df = load_file(file_path=cr_file_path)
                self.headers[self.patient].append(header)
        else:
            for sample_name, cr_file_path in zip(self.sample_names, self.copy_ratio_files):
                header, df = load_file(file_path=cr_file_path, sample_name=sample_name)
                self.headers[sample_name].append(header)
        if self.ref_dict_file is not None:
            self.ref_dict = self.get_header(self.ref_dict_file)
        print() if self.verbose else None

        self.contigs = self.get_contigs()

    def get_contigs(self):
        contigs = []
        if self.ref_dict is not None:
            contigs = [
                h.split("\t")[1].removeprefix("SN:")
                for h in self.ref_dict.split("\n")
                if h.startswith(f"{self.comment_char}SQ")
            ]
        for _, header in self.headers.items():
            contigs = [
                h.split("\t")[1].removeprefix("SN:")
                for h in header[0].split("\n")
                if h.startswith(f"{self.comment_char}SQ")
            ]
            if len(contigs):
                break
        if not len(contigs):
            # Add default contigs for hg19 and hg38 in case no headers were found:
            contigs = (
                [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
                + [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
            )
        return contigs

    def get_header(self, file_path: str):
        open_func = gzip.open if file_path.endswith(".gz") else open
        with open_func(file_path, "rt") as file:
            # save all lines starting with comment_char as header
            header = "".join([line for line in file if line.startswith(self.comment_char)])
        return header

    def get_header_and_df(self, file_path: str):
        try:
            header = self.get_header(file_path)
            # TODO: support hdf5
            df = pd.read_csv(file_path, sep="\t", comment=self.comment_char, header=0, names=self.column_names, low_memory=False)
            if self.column_types:
                df = df.astype({key: self.column_types[key] for key in self.column_names if key in self.column_types})
        except Exception as e:
            warnings.warn(f"Exception reading file {file_path}: {e}")
            warnings.warn(f"Setting header to None and df to empty DataFrame.")
            header = None
            df = pd.DataFrame(columns=self.column_names)
        return header, df

    def write_output_files(self):
        os.makedirs(self.output_dir, exist_ok=True)

        def get_header(name):
            # propagate the metadata from copy ratios as they contain bam names which may differ from assigned sample names
            metadata = (
                "\n".join([line for line in self.headers[name][0].split("\n") if not line.startswith("@HD") and not line.startswith("@SQ")])
                if len(self.headers[name]) > 0
                else ""
            )
            return (
                self.ref_dict + metadata
                if self.ref_dict is not None
                else (self.headers[name][0] if name in self.headers and len(self.headers[name]) > 0 else None)
            )

        if self.sample_names is None:
            if self.harmonized_copy_ratios.empty:
                df = self.harmonized_copy_ratios
            else:
                with warnings.catch_warnings():
                    # .stack has new implementation in pandas >= 2.1.0 that raises FutureWarning
                    warnings.simplefilter(action="ignore", category=FutureWarning)
                    df = self.harmonized_copy_ratios.stack(level=self.sample_col)
                sort_by = [self.sample_col, self.contig_col, self.start_col, self.end_col]
                df = df.reindex(self.sort_genomic_positions(index=df.index, by=sort_by)).reset_index()
            self.write_header_and_df(
                header=get_header(self.patient),
                df=df,
                file_path=f"{self.output_dir}/{self.patient}{self.suffix}" + (".gz" if self.compress_output else ""),
            )
        else:
            for sample_name in np.unique(self.sample_names):
                df = (
                    self.harmonized_copy_ratios[sample_name].reset_index()
                    if not self.harmonized_copy_ratios.empty
                    else pd.DataFrame(columns=self.column_names)
                )
                self.write_header_and_df(
                    header=get_header(sample_name),
                    df=df,
                    file_path=f"{self.output_dir}/{sample_name}{self.suffix}" + (".gz" if self.compress_output else ""),
                )

    def write_header_and_df(self, header: str, df: pd.DataFrame, file_path: str):
        print(f"Writing output file {file_path}") if self.verbose else None
        # TODO: support hdf5
        open_func = gzip.open if file_path.endswith(".gz") else open
        with open_func(file_path, "wt") as file:
            if header is not None:
                file.write(header)
            df = df.astype({c: t for c, t in self.column_types.items() if c in df.columns})
            df.to_csv(file, sep="\t", index=False)

    def split_intervals(self, intervals: pd.DataFrame):
        starts = pd.DataFrame(data={"start": 1}, index=intervals[self.start_col])
        ends = pd.DataFrame(data={"end": -1}, index=intervals[self.end_col])
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

    def non_overlapping(self, intervals: pd.DataFrame):
        new_intervals = self.split_intervals(intervals)

        # Map each original interval to the new non-overlapping intervals
        non_overlapping_intervals = []
        for new_interval in new_intervals:
            contains_interval = intervals[(intervals[self.start_col] <= new_interval.left) & (intervals[self.end_col] > new_interval.left)]
            for _, interval in contains_interval.iterrows():
                new_row = copy(interval)
                new_row[self.start_col] = new_interval.left
                new_row[self.end_col] = new_interval.right
                non_overlapping_intervals.append(new_row)

        non_overlapping_df = pd.DataFrame(non_overlapping_intervals)
        print(".", end="", flush=True) if self.verbose else None
        return non_overlapping_df

    def sort_genomic_positions(self, index: pd.MultiIndex, by=None) -> pd.MultiIndex:
        by = [self.contig_col, self.start_col, self.end_col] if by is None else by
        contig_order = self.contigs
        contig_order += list(sorted(set(index.get_level_values(self.contig_col)) - set(self.contigs)))
        temp_df = pd.DataFrame(index=index).reset_index().astype({c: t for c, t in self.column_types.items() if c in index.names})
        temp_df[self.contig_col] = pd.Categorical(temp_df[self.contig_col], categories=contig_order, ordered=True)
        temp_df.sort_values(by=by, inplace=True)
        temp_df = temp_df.astype({c: t for c, t in self.column_types.items() if c in temp_df.columns})
        return pd.MultiIndex.from_frame(temp_df[index.names])

    def merge_abutting_intervals(self, intervals: pd.DataFrame):
        new_intervals = []
        previous_interval = None
        for i, interval in intervals.iterrows():
            if previous_interval is None:
                previous_interval = interval
                new_intervals.append(interval)
                continue

            # interval.name == (contig, start, end)
            if interval.name[0] == previous_interval.name[0] and abs(interval.name[1] - previous_interval.name[2]) <= self.merge_abutting_distance:
                equal = previous_interval.shape[0] > 0
                for _pi, _i in zip(previous_interval, interval):
                    try:
                        # if interval is just a number
                        equal &= abs(_pi - _i) < self.merge_abutting_tolerance
                    except:
                        try:
                            # if interval is a vector of numbers
                            equal &= _pi == _i
                        except:
                            pass
                if equal:
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
        merged_intervals.index.names = [self.contig_col, self.start_col, self.end_col]
        return merged_intervals

    def select_intervals(self, intervals: pd.DataFrame, min_length: int = 1):
        lengths = intervals.index.get_level_values(self.end_col) - intervals.index.get_level_values(self.start_col) + 1
        return intervals[lengths >= min_length]

    def harmonize_loci_parallel(self, copy_ratios: pd.DataFrame):
        n_unique_contigs = copy_ratios[self.contig_col].nunique()
        processes = np.min([self.threads, n_unique_contigs, os.cpu_count()])
        success = False
        if processes > 1:
            try:
                print(f"Parallel execution with {processes} processes over contigs: ", end="", flush=True) if self.verbose else None
                with mp.Pool(processes=processes) as pool:
                    harmonized_loci_by_contig = pool.starmap(
                        func=self.non_overlapping,
                        iterable=[[contig_group] for contig, contig_group in copy_ratios.groupby(by=[self.contig_col])]
                    )
                success = True
            except Exception as e:
                print(f"Error in parallel execution: {e}")
        if not success:
            print("Serial execution over contigs: ", end="", flush=True) if self.verbose else None
            harmonized_loci_by_contig = [
                self.non_overlapping(contig_group)
                for contig, contig_group in copy_ratios.groupby(by=[self.contig_col])
            ]
        print() if self.verbose else None
        harmonized_loci = pd.concat(harmonized_loci_by_contig)
        harmonized_loci = harmonized_loci.astype({c: t for c, t in self.column_types.items() if c in harmonized_loci.columns})
        return harmonized_loci

    @staticmethod
    def agg_func_map(func):
        if func == "strongest_signal":
            return lambda x: x.values[np.argmax(np.abs(x.values))]
        elif func == "error_mean":
            return lambda x: np.sqrt(np.mean(x ** 2))
        elif func == "first":
            return lambda x: x.values[0]
        elif func == "last":
            return lambda x: x.values[-1]
        return func

    def _harmonize(self):
        print("Harmonizing copy ratio intervals:") if self.verbose else None
        cr = pd.concat(self.copy_ratios, ignore_index=True)
        if cr.empty:
            print("  No loci found to harmonize in any of the input files.") if self.verbose else None
            return pd.DataFrame(columns=[self.contig_col, self.start_col, self.end_col, self.sample_col] + self.agg_col)

        n_intervals_per_file = [
            cr_f[[self.contig_col, self.start_col, self.end_col]].drop_duplicates().shape[0]
            for cr_f in self.copy_ratios
        ]
        n_intervals = cr[[self.contig_col, self.start_col, self.end_col]].drop_duplicates().shape[0]
        print(f"  Number of loci per file before harmonization: {n_intervals_per_file}") if self.verbose else None

        t = time.perf_counter()
        harmonized_loci = self.harmonize_loci_parallel(copy_ratios=cr)
        print(f"Harmonizing copy ratio intervals took {time.perf_counter() - t:.3f} seconds.") if self.verbose else None

        # Bring into table format with samples as columns
        aggregate = (
            harmonized_loci
            .groupby(by=[self.contig_col, self.start_col, self.end_col, self.sample_col])
            .agg(
                **{
                    col: pd.NamedAgg(
                        column=col,
                        aggfunc=self.agg_func_map(func)
                    )
                    for col, func in zip(self.agg_col, self.agg_func)
                }
            )
            .unstack(level=self.sample_col)
            .swaplevel(axis=1)
        )
        aggregate = aggregate.reindex(self.sort_genomic_positions(index=aggregate.index))

        n_harmonized = aggregate.shape[0]
        pct_increase = 100 * (n_harmonized / np.max(n_intervals) - 1)
        print(f"  Number of loci after harmonization: {n_harmonized} (+{pct_increase:.3f}%)") if self.verbose else None

        # Drop intervals that are not present in all samples:
        if self.interval_set_rule == "INTERSECTION":
            consensus_intervals = aggregate.dropna()
            n_consensus = consensus_intervals.shape[0]
            pct_drop = 100 * (1 - n_consensus / n_harmonized)
            print(f"  Number of loci after dropping loci not present in all samples: {n_consensus} (-{pct_drop:.3f}%)") if self.verbose else None
            if n_consensus == 0:
                return consensus_intervals
        else:
            consensus_intervals = aggregate
            n_consensus = n_harmonized

        # Drop intervals that are called as amp or del in the matched normal sample:
        if self.filter_germline_calls and self.normal_sample is not None and "CALL" in self.column_names:
            germline_mask = consensus_intervals[(self.normal_sample, "CALL")].isin(["-", "+"])
            somatic_intervals = consensus_intervals.loc[~germline_mask]
            n_somatic = somatic_intervals.shape[0]
            pct_drop = 100 * (1 - n_somatic / n_consensus)
            print(f"  Number of loci after dropping loci called as CNV in germline sample: {n_somatic} (-{pct_drop:.3f}%)") if self.verbose else None
            if n_somatic == 0:
                return somatic_intervals
        else:
            somatic_intervals = consensus_intervals
            n_somatic = n_consensus

        # merge abutting intervals with same copy ratio
        merged_consensus_intervals = self.merge_abutting_intervals(somatic_intervals)
        n_merged = merged_consensus_intervals.shape[0]
        pct_drop = 100 * (1 - n_merged / n_somatic)
        print(f"  Number of loci after merging abutting intervals (dist <= {self.merge_abutting_distance}) with the same copy ratio (tol < {self.merge_abutting_tolerance}): {n_merged} (-{pct_drop:.3f}%)") if self.verbose else None

        # Drop intervals that are too short:
        selected_merged_consensus_intervals = self.select_intervals(merged_consensus_intervals, min_length=self.min_target_length)
        n_selected = selected_merged_consensus_intervals.shape[0]
        pct_drop = 100 * (1 - n_selected / n_merged)
        print(f"  Number of loci after dropping loci shorter than {self.min_target_length} bp: {n_selected} (-{pct_drop:.3f}%)") if self.verbose else None

        # print() if self.verbose else None
        return selected_merged_consensus_intervals


if __name__ == "__main__":
    main()
