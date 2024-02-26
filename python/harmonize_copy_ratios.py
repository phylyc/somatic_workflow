import argparse
from collections import defaultdict
from copy import copy
import gzip
import pandas as pd
import warnings


def parse_args():
    parser = argparse.ArgumentParser(
        prog="HarmonizeCopyRatios",
        description="""
            
        """,
        epilog="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.usage = "genotype.py --sample <sample> [--sample <sample> ...] -D <denoised_cr> [-D <denoised_cr> ...] -S <standardized_cr> [-S <standardized_cr>] [-O <output_dir>] [--compress_output] [--verbose]"
    parser.add_argument("--sample",                 type=str,   required=True,  action="append",    help="Assigned name of the sample. Does not have to coincide with the sample name in the pileup file header.")
    parser.add_argument("-D", "--denoised_cr",      type=str,   required=True,  action="append",    help="Path to the denoised copy ratio file (output of GATK's DenoiseReadCounts).")
    parser.add_argument("-S", "--standardized_cr",  type=str,   required=True,  action="append",    help="Path to the standardized copy ratio file (output of GATK's DenoiseReadCounts).")
    parser.add_argument("-O", "--output_dir",       type=str,   default=".",    help="Path to the output directory.")
    parser.add_argument("--compress_output",                    default=False,  action="store_true", help="Compress output files.")
    parser.add_argument("--verbose",                            default=False,  action="store_true", help="Print information to stdout during execution.")
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


def split_intervals(intervals: pd.DataFrame, start_col: str = "START", end_col: str = "END"):
    starts = pd.DataFrame(data={start_col: 1}, index=intervals[start_col])
    ends = pd.DataFrame(data={end_col: -1}, index=intervals[end_col])
    transitions = pd.merge(starts, ends, how="outer", left_index=True, right_index=True).fillna(0)
    # This dataframe stores per position the type of transitions. Automatically sorted.
    # Now, we need to know at each position if we are at an interval start or not.
    # This can be done by counting the opening & closing parenthesis.
    transitions["is_start"] = (transitions.pop(start_col) + transitions.pop(end_col)).cumsum().astype(bool)
    # This handles overlapping intervals.
    # Create interval indices and take all intervals that have a valid start.
    new_intervals = pd.IntervalIndex.from_breaks(transitions.index, closed="left")[transitions["is_start"][:-1]]
    return new_intervals[~new_intervals.is_empty]


def non_overlapping(intervals: pd.DataFrame, start_col: str = "START", end_col: str = "END"):
    new_intervals = split_intervals(intervals, start_col=start_col, end_col=end_col)
    non_overlapping_interval_list = []
    for i, interval in intervals.iterrows():
        pd_interval = pd.Interval(interval[start_col], interval[end_col], closed="left")
        overlapping_new_intervals = [n_i for n_i in new_intervals if pd_interval.overlaps(n_i)]
        # Split interval and content into new intervals:
        for new_interval in overlapping_new_intervals:
            new_split_interval = copy(interval)
            new_split_interval[start_col] = new_interval.left
            new_split_interval[end_col] = new_interval.right
            non_overlapping_interval_list.append(new_split_interval)
    return pd.DataFrame(non_overlapping_interval_list)


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

        new_intervals.append(interval)
        previous_interval = interval

    merged_intervals = pd.DataFrame(new_intervals)
    merged_intervals.index.names = ["CONTIG", "START", "END"]
    return merged_intervals


def harmonize(copy_ratios: list[pd.DataFrame], verbose: bool = False):
    cr = pd.concat(copy_ratios, ignore_index=True)
    if cr.empty:
        print("  No loci found to harmonize in any of the input files.") if verbose else None
        return cr

    harmonized_loci = pd.concat(
        [
            non_overlapping(contig_group, start_col="START", end_col="END")
            for contig, contig_group in cr.groupby(by=["CONTIG"])
        ]
    )
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
    dcr_headers = defaultdict(list)
    scr_headers = defaultdict(list)
    dcr = []
    scr = []

    # load data (per sequencing run):
    for sample_name, dcr_file_path, scr_file_path in zip(args.sample, args.denoised_cr, args.standardized_cr):
        dcr_header, dcr_df = get_header_and_df(file_path=dcr_file_path, columns=columns, column_types=column_types)
        scr_header, scr_df = get_header_and_df(file_path=scr_file_path, columns=columns, column_types=column_types)
        dcr_headers[sample_name].append(dcr_header)
        scr_headers[sample_name].append(scr_header)
        dcr_df["SAMPLE"] = sample_name
        scr_df["SAMPLE"] = sample_name
        dcr.append(dcr_df)
        scr.append(scr_df)
        print(f"Number of loci in {dcr_file_path}: {dcr_df.shape[0]}") if args.verbose else None
    print() if args.verbose else None

    # harmonize intervals:
    print("Harmonizing denoised copy ratio intervals:") if args.verbose else None
    harmonized_dcr = harmonize(copy_ratios=dcr, verbose=args.verbose)
    print() if args.verbose else None
    print("Harmonizing standardized copy ratio intervals:") if args.verbose else None
    harmonized_scr = harmonize(copy_ratios=scr, verbose=args.verbose)
    print() if args.verbose else None

    # write output files:
    for sample_name in args.sample:
        dcr_header = dcr_headers[sample_name][0] if len(dcr_headers[sample_name]) > 0 else None
        scr_header = scr_headers[sample_name][0] if len(scr_headers[sample_name]) > 0 else None
        dcr_df = harmonized_dcr[sample_name].reset_index() if not harmonized_dcr.empty else pd.DataFrame(columns=columns)
        scr_df = harmonized_scr[sample_name].reset_index() if not harmonized_scr.empty else pd.DataFrame(columns=columns)
        write_header_and_df(
            header=dcr_header,
            df=dcr_df,
            file_path=f"{args.output_dir}/{sample_name}.denoised_CR.tsv" + (".gz" if args.compress_output else ""),
            verbose=args.verbose
        )
        write_header_and_df(
            header=scr_header,
            df=scr_df,
            file_path=f"{args.output_dir}/{sample_name}.standardized_CR.tsv" + (".gz" if args.compress_output else ""),
            verbose=args.verbose
        )


if __name__ == "__main__":
    main()
