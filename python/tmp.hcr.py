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
        prog="HarmonizeCopyRatios",
        description="""

        """,
        epilog="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.usage = "harmonize_copy_ratios.py --sample <sample> [--sample <sample> ...] -I <copy_ratio> [-I <copy_ratio> ...] [-O <output_dir>] [--compress_output] [--verbose]"
    parser.add_argument("--sample", type=str, required=True, action="append", help="Assigned name of the sample. Does not have to coincide with the sample name in the pileup file header.")
    parser.add_argument("-I", "--copy_ratio", type=str, required=True, action="append", help="Path to the (denoised) copy ratio file (output of GATK's DenoiseReadCounts).")
    parser.add_argument("-O", "--output_dir", type=str, default=".", help="Path to the output directory.")
    parser.add_argument("--suffix", type=str, default=".harmonized.CR", help="Suffix to append to each sample name to create the output file name.")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads to use for parallelization over samples.")
    parser.add_argument("--compress_output", default=False, action="store_true", help="Compress output files.")
    parser.add_argument("--verbose", default=False, action="store_true", help="Print information to stdout during execution.")
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


def merge_abutting_intervals(intervals: pd.DataFrame, dist: int = 1, cr_tol: float = 1e-3):
    new_intervals = []
    previous_interval = None
    for i, interval in intervals.iterrows():
        if previous_interval is None:
            previous_interval = interval
            new_intervals.append(interval)
            continue

        # interval.name == (contig, start, end)
        if interval.name[0] == previous_interval.name[0] and abs(interval.name[1] - previous_interval.name[2]) <= dist:
            if (abs(previous_interval - interval) < cr_tol).all():
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


def harmonize(copy_ratios: list[pd.DataFrame], max_processes: int = 1, verbose: bool = False):
    cr = pd.concat(copy_ratios, ignore_index=True)
    if cr.empty:
        print("  No loci found to harmonize in any of the input files.") if verbose else None
        return cr

    t = time.perf_counter()
    harmonized_loci = harmonize_loci_parallel(copy_ratios=cr, verbose=verbose, max_processes=max_processes)
    print(f"Harmonizing copy ratio intervals took {time.perf_counter() - t:.3f} seconds.") if verbose else None

    aggregate = harmonized_loci.groupby(by=["CONTIG", "START", "END", "SAMPLE"]).agg(
        {
            "LOG2_COPY_RATIO": "mean"
        }
    ).unstack(level="SAMPLE").swaplevel(axis=1)
    aggregate = aggregate.reindex(sort_genomic_positions(index=aggregate.index))

    max_n_intervals = max([cr.shape[0] for cr in copy_ratios])
    n_harmonized = aggregate.shape[0]
    pct_increase = 100 * (n_harmonized / max_n_intervals - 1)
    print(f"  Number of loci after harmonization: {n_harmonized} (+{pct_increase:.3f}%)") if verbose else None

    # Drop intervals that are not present in all samples:
    consensus_intervals = aggregate.dropna()
    n_consensus = consensus_intervals.shape[0]
    pct_drop = 100 * (1 - n_consensus / n_harmonized)
    print(f"  Number of loci after dropping loci not present in all samples: {n_consensus} (-{pct_drop:.3f}%)") if verbose else None
    if n_consensus == 0:
        return consensus_intervals

    # merge abutting intervals with same copy ratio
    merged_consensus_intervals = merge_abutting_intervals(consensus_intervals)
    n_merged = merged_consensus_intervals.shape[0]
    pct_drop = 100 * (1 - n_merged / n_consensus)
    print(f"  Number of loci after merging abutting intervals with the same copy ratio: {n_merged} (-{pct_drop:.3f}%)") if verbose else None

    return merged_consensus_intervals


def harmonize_intervals(args):
    columns = ["CONTIG", "START", "END", "LOG2_COPY_RATIO"]
    column_types = {
        "CONTIG": str,
        "START": int,
        "END": int,
        "LOG2_COPY_RATIO": float,
    }
    cr_headers = defaultdict(list)
    cr = []

    # load data (per sequencing run):
    for sample_name, cr_file_path in zip(args.sample, args.copy_ratio):
        cr_header, cr_df = get_header_and_df(file_path=cr_file_path, columns=columns, column_types=column_types)
        cr_headers[sample_name].append(cr_header)
        cr_df["SAMPLE"] = sample_name
        cr.append(cr_df)
        print(f"Number of loci in {cr_file_path}: {cr_df.shape[0]}") if args.verbose else None
    print() if args.verbose else None

    # harmonize intervals:
    print("Harmonizing copy ratio intervals:") if args.verbose else None
    harmonized_cr = harmonize(copy_ratios=cr, max_processes=args.threads, verbose=args.verbose)
    print() if args.verbose else None

    # write output files:
    os.makedirs(args.output_dir, exist_ok=True)
    for sample_name in np.unique(args.sample):
        cr_header = cr_headers[sample_name][0] if len(cr_headers[sample_name]) > 0 else None
        cr_df = harmonized_cr[sample_name].reset_index() if not harmonized_cr.empty else pd.DataFrame(columns=columns)
        write_header_and_df(
            header=cr_header,
            df=cr_df,
            file_path=f"{args.output_dir}/{sample_name}{args.suffix}.tsv" + (".gz" if args.compress_output else ""),
            verbose=args.verbose
        )


if __name__ == "__main__":
    main()
