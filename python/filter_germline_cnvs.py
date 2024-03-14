import argparse
from collections import defaultdict
from copy import copy
import gzip
import multiprocessing as mp
import numpy as np
import os
import pandas as pd
import time
import warnings


def parse_args():
    parser = argparse.ArgumentParser(
        prog="FilterGermlineCNVs",
        description="""

        """,
        epilog="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.usage = "python filter_germline_cnvs.py --tumor_called_copy_ratio_segmentation <tumor_called_copy_ratio_segmentation> --normal_called_copy_ratio_segmentation <normal_called_copy_ratio_segmentation> [-O <output_dir>] [--verbose]"
    parser.add_argument("--tumor_called_copy_ratio_segmentation",   type=str,   required=True,  help="Path to the called segmented copy ratio file of the tumor sample (output of GATK's ModelSegments + CallCopyRatioSegments).")
    parser.add_argument("--normal_called_copy_ratio_segmentation",  type=str,   required=True,  help="Path to the called segmented copy ratio file of the normal sample (output of GATK's ModelSegments + CallCopyRatioSegments).")
    parser.add_argument("-O", "--output",                           type=str,   required=True,  help="Path to the output file.")
    parser.add_argument("--threads",                                type=int,   default=1,      help="Number of threads to use for parallelization over chromosomes.")
    parser.add_argument("--min_segment_length",                     type=int,   default=100,    help="Minimum length of a segment to be included in output.")
    parser.add_argument("--compress_output",                                    default=False,  action="store_true", help="Compress output files.")
    parser.add_argument("--verbose",                                            default=False,  action="store_true", help="Print information to stdout during execution.")
    return parser.parse_args()


def main():
    args = parse_args()
    print_args(args)
    harmonize_intervals(args=args)


def print_args(args):
    if args.verbose:
        print("Calling HarmonizeCopyRatios")
        print("Arguments:")
        for key, value in vars(args).items():
            print(f"  {key}: {value}")
        print()


def get_header_and_df(file_path: str, columns: list[str] = None, column_types: dict[str, type] = None, comment_char: str = "@"):
    open_func = gzip.open if file_path.endswith(".gz") else open
    try:
        with open_func(file_path, "rt") as file:
            # save all lines starting with comment_char as header
            header = "".join([line for line in file if line.startswith(comment_char)])
        # TODO: support hdf5
        df = pd.read_csv(file_path, sep="\t", comment=comment_char, header=0, names=columns, low_memory=False)
        if column_types is not None:
            df = df.astype({key: column_types[key] for key in columns if key in column_types})
    except Exception as e:
        warnings.warn(f"Exception reading file {file_path}: {e}")
        warnings.warn(f"Setting header to None and df to empty DataFrame.")
        header = None
        df = pd.DataFrame(columns=columns)
    return header, df


def write_header_and_df(header: str, df: pd.DataFrame, file_path: str, verbose: bool = False):
    print(f"Writing output file {file_path}") if verbose else None
    # TODO: support hdf5
    open_func = gzip.open if file_path.endswith(".gz") else open
    with open_func(file_path, "wt") as file:
        file.write(header) if header is not None else None
        df.to_csv(file, sep="\t", index=False)


def split_intervals(intervals: pd.DataFrame, start_col: str, end_col: str):
    # Create a sorted list of all unique interval boundaries
    boundaries = sorted(set(intervals[start_col]).union(set(intervals[end_col])))
    # Use the boundaries to form non-overlapping intervals
    new_intervals = pd.IntervalIndex.from_breaks(boundaries, closed='both')
    return new_intervals


def non_overlapping(intervals: pd.DataFrame, start_col: str = "START", end_col: str = "END", verbose: bool = False):
    new_intervals = split_intervals(intervals, start_col, end_col)

    # Map each original interval to the new non-overlapping intervals
    interval_mapping = defaultdict(list)
    for interval in new_intervals:
        contained_by_interval = intervals[(intervals[start_col] <= interval.left) & (intervals[end_col] > interval.left)]
        for _, row in contained_by_interval.iterrows():
            interval_mapping[interval].append(row)

    # Flatten the mapping into a list of non-overlapping intervals with their corresponding data
    non_overlapping_intervals = []
    for interval, rows in interval_mapping.items():
        for row in rows:
            new_row = copy(row)
            new_row[start_col] = interval.left
            new_row[end_col] = interval.right
            non_overlapping_intervals.append(new_row)

    non_overlapping_df = pd.DataFrame(non_overlapping_intervals)
    print(".", end="", flush=True) if verbose else None
    return non_overlapping_df


def sort_genomic_positions(index: pd.MultiIndex) -> pd.MultiIndex:
    contig_order = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
    temp_df = pd.DataFrame(index=index).reset_index()
    temp_df["CONTIG"] = pd.Categorical(temp_df["CONTIG"], categories=contig_order, ordered=True)
    temp_df.sort_values(by=["CONTIG", "START", "END"], inplace=True)
    return pd.MultiIndex.from_frame(temp_df[["CONTIG", "START", "END"]])


def merge_abutting_intervals(intervals: pd.DataFrame, cols: list[str], dist: int = 1, tol: float = 1e-3):
    new_intervals = []
    previous_interval = None
    for i, interval in intervals.iterrows():
        if previous_interval is None:
            previous_interval = interval
            new_intervals.append(interval)
            continue

        # interval.name == (contig, start, end)
        if interval.name[0] == previous_interval.name[0] and abs(interval.name[1] - previous_interval.name[2]) <= dist:
            if (abs(previous_interval[cols] - interval[cols]) < tol).all() & (previous_interval["CALL"] == interval["CALL"]):
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
    merged_intervals.index.names = ["CONTIG", "START", "END"]
    return merged_intervals


def select_intervals(intervals: pd.DataFrame, min_length: int = 1):
    lengths = intervals.index.get_level_values("END") - intervals.index.get_level_values("START") + 1
    return intervals[lengths >= min_length]


def harmonize_loci_parallel(copy_ratios: pd.DataFrame, verbose: bool = False, max_processes: int = 1):
    n_unique_contigs = copy_ratios["CONTIG"].nunique()
    processes = np.min([max_processes, n_unique_contigs, os.cpu_count()])
    if processes > 1:
        try:
            print(f"Parallel execution with {processes} processes over contigs: ", end="", flush=True) if verbose else None
            with mp.Pool(processes=processes) as pool:
                harmonized_loci_by_contig = pool.starmap(
                    func=non_overlapping,
                    iterable=[
                        (contig_group, "START", "END", verbose)
                        for contig, contig_group in copy_ratios.groupby(by=["CONTIG"])
                    ]
                )
            print() if verbose else None
            return pd.concat(harmonized_loci_by_contig)

        except Exception as e:
            print(f"Error in parallel execution: {e}")
    print("Serial execution over contigs: ", end="", flush=True) if verbose else None
    harmonized_loci_by_contig = [
        non_overlapping(contig_group, "START", "END", verbose)
        for contig, contig_group in copy_ratios.groupby(by=["CONTIG"])
    ]
    print() if verbose else None
    return pd.concat(harmonized_loci_by_contig)


def harmonize_and_filter(copy_ratios: list[pd.DataFrame], sample_name: str, normal_sample_name: str, merge_on_cols: list[str], min_target_length: int = 1, max_processes: int = 1, verbose: bool = False):
    cr = pd.concat(copy_ratios, ignore_index=True)
    if cr.empty:
        print("  No loci found to harmonize in any of the input files.") if verbose else None
        return cr

    t = time.perf_counter()
    harmonized_loci = harmonize_loci_parallel(copy_ratios=cr, verbose=verbose, max_processes=max_processes)
    print(f"Harmonizing intervals took {time.perf_counter() - t:.3f} seconds.") if verbose else None

    aggregate = harmonized_loci.groupby(by=["CONTIG", "START", "END", "SAMPLE"]).agg(
        {key: "sum" for key in harmonized_loci.columns if key not in ["CONTIG", "START", "END", "SAMPLE"]}
    ).unstack(level="SAMPLE").swaplevel(axis=1)
    aggregate = aggregate.reindex(sort_genomic_positions(index=aggregate.index))

    max_n_intervals = max([cr.shape[0] for cr in copy_ratios])
    n_harmonized = aggregate.shape[0]
    pct_increase = 100 * (n_harmonized / max_n_intervals - 1)
    print(f"  Number of loci after harmonization: {n_harmonized} (+{pct_increase:.3f}%)") if verbose else None

    # Drop intervals that are called as amp or del in the normal sample:
    germline_mask = aggregate[(normal_sample_name, "CALL")].isin(["-", "+"])
    somatic_intervals = aggregate.loc[~germline_mask, sample_name].dropna(how="all")
    n_somatic = somatic_intervals.shape[0]
    pct_drop = 100 * (1 - n_somatic / n_harmonized)
    print(f"  Number of loci after dropping loci called as CNV in germline sample: {n_somatic} (-{pct_drop:.3f}%)") if verbose else None
    if n_somatic == 0:
        return somatic_intervals

    # merge abutting intervals with same copy ratio
    merged_consensus_intervals = merge_abutting_intervals(somatic_intervals, cols=merge_on_cols)
    n_merged = merged_consensus_intervals.shape[0]
    pct_drop = 100 * (1 - n_merged / n_somatic)
    print(f"  Number of loci after merging abutting intervals with the same copy ratio: {n_merged} (-{pct_drop:.3f}%)") if verbose else None

    # Drop intervals that are too short:
    selected_merged_consensus_intervals = select_intervals(merged_consensus_intervals, min_length=min_target_length)
    n_selected = selected_merged_consensus_intervals.shape[0]
    pct_drop = 100 * (1 - n_selected / n_merged)
    print(f"  Number of loci after dropping loci shorter than {min_target_length} bp: {n_selected} (-{pct_drop:.3f}%)") if verbose else None

    return selected_merged_consensus_intervals


def harmonize_intervals(args):
    columns = [
        "CONTIG", "START", "END",
        "NUM_POINTS_COPY_RATIO", "NUM_POINTS_ALLELE_FRACTION",
        "LOG2_COPY_RATIO_POSTERIOR_10", "LOG2_COPY_RATIO_POSTERIOR_50", "LOG2_COPY_RATIO_POSTERIOR_90",
        "MINOR_ALLELE_FRACTION_POSTERIOR_10", "MINOR_ALLELE_FRACTION_POSTERIOR_50", "MINOR_ALLELE_FRACTION_POSTERIOR_90",
        "CALL"
    ]
    column_types = {
        "CONTIG": str,
        "START": int,
        "END": int,
        "NUM_POINTS_COPY_RATIO": int,
        "NUM_POINTS_ALLELE_FRACTION": int,
        "LOG2_COPY_RATIO_POSTERIOR_10": float,
        "LOG2_COPY_RATIO_POSTERIOR_50": float,
        "LOG2_COPY_RATIO_POSTERIOR_90": float,
        "MINOR_ALLELE_FRACTION_POSTERIOR_10": float,
        "MINOR_ALLELE_FRACTION_POSTERIOR_50": float,
        "MINOR_ALLELE_FRACTION_POSTERIOR_90": float,
        "CALL": str
    }
    cr_header, cr_df = get_header_and_df(file_path=args.tumor_called_copy_ratio_segmentation, columns=columns, column_types=column_types)
    _, normal_cr_df = get_header_and_df(file_path=args.normal_called_copy_ratio_segmentation, columns=columns, column_types=column_types)

    sample_name = os.path.basename(args.tumor_called_copy_ratio_segmentation)
    normal_sample_name = os.path.basename(args.normal_called_copy_ratio_segmentation)

    cr_df["SAMPLE"] = sample_name
    normal_cr_df["SAMPLE"] = normal_sample_name

    print(f"Number of loci in {sample_name}: {cr_df.shape[0]}") if args.verbose else None
    print(f"Number of loci in {normal_sample_name}: {normal_cr_df.shape[0]}") if args.verbose else None
    print() if args.verbose else None

    print("Harmonizing and filter copy ratio intervals:") if args.verbose else None
    if not cr_df.empty and not normal_cr_df.empty:
        filtered_cr = harmonize_and_filter(
            copy_ratios=[cr_df, normal_cr_df],
            sample_name=sample_name,
            normal_sample_name=normal_sample_name,
            merge_on_cols=["LOG2_COPY_RATIO_POSTERIOR_10", "LOG2_COPY_RATIO_POSTERIOR_50", "LOG2_COPY_RATIO_POSTERIOR_90", "MINOR_ALLELE_FRACTION_POSTERIOR_10", "MINOR_ALLELE_FRACTION_POSTERIOR_50", "MINOR_ALLELE_FRACTION_POSTERIOR_90"],
            min_target_length=args.min_segment_length,
            max_processes=args.threads,
            verbose=args.verbose
        )
        print() if args.verbose else None
    else:
        filtered_cr = cr_df

    write_header_and_df(
        header=cr_header,
        df=filtered_cr.reset_index()[columns] if not filtered_cr.empty else pd.DataFrame(columns=columns),
        file_path=args.output,
        verbose=args.verbose
    )


if __name__ == "__main__":
    main()
