import argparse
from collections import defaultdict
import gzip
import pandas as pd
import warnings


def parse_args():
    parser = argparse.ArgumentParser(
        prog="MergePileups",
        description="""
            Merge multiple allelic pileup files into one file per sample. Pileups with the same sample name (e.g. from different sequencing runs) 
            are merged by summing ref_count, alt_count, other_alt_count and taking the maximum allele_frequency.
            The input files are expected to be in the format of GATK's GetPileupSummaries output.
            The output files will be in the same format.
        """,
        epilog="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.usage = "merge_pileups.py --sample <sample> [--sample <sample> ...] -P <pileup> [-P <pileup> ...] [-D <ref_dict>] [-O <output_dir>] [--compress_output] [--verbose]"
    parser.add_argument("--sample",             type=str,   required=True,  action="append",    help="Assigned name of the sample. Does not have to coincide with the sample name in the pileup file header.")
    parser.add_argument("-P", "--pileup",       type=str,   required=True,  action="append",    help="Path to the allelic pileup file (output of GATK's GetPileupSummaries; used for CalculateContamination).")
    parser.add_argument("-D", "--ref_dict",     type=str,                   help="Path to the reference dictionary to sort rows.")
    parser.add_argument("-O", "--output_dir",   type=str,   default=".",    help="Path to the output directory.")
    parser.add_argument("--compress_output",                default=False,  action="store_true", help="Compress output files.")
    parser.add_argument("--verbose",                        default=False,  action="store_true", help="Print information to stdout during execution.")
    return parser.parse_args()


def main():
    args = parse_args()
    print_args(args)
    merge_pileups(args=args)


def print_args(args):
    if args.verbose:
        print("Calling MergePileups")
        print("Arguments:")
        for key, value in vars(args).items():
            print(f"  {key}: {value}")
        print()


def get_header_and_df(file_path: str, columns: list[str] = None, column_types:  dict[str, type] = None, comment_char: str = "#"):
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
    open_func = gzip.open if file_path.endswith(".gz") else open
    with open_func(file_path, "wt") as file:
        file.write(header) if header is not None else None
        df.to_csv(file, sep="\t", index=False)


def sort_genomic_positions(index: pd.MultiIndex, contig_order: list[str]) -> pd.MultiIndex:
    contig_order += list(set(index.get_level_values("contig")) - set(contig_order))
    temp_df = pd.DataFrame(index=index).reset_index()
    temp_df["contig"] = pd.Categorical(temp_df["contig"], categories=contig_order, ordered=True)
    temp_df.sort_values(by=["contig", "position"], inplace=True)
    return pd.MultiIndex.from_frame(temp_df[["contig", "position"]])


def get_contigs(ref_dict: str) -> list[str]:
    if ref_dict is None:
        return (
            [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
            + [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
        )
    with open(ref_dict, "rt") as file:
        contig_header = [line for line in file if line.startswith("@SQ")]
    return [h.split("\t")[1].removeprefix("SN:") for h in contig_header]


def merge_pileups(args):
    columns = ["contig", "position", "ref_count", "alt_count", "other_alt_count", "allele_frequency"]
    column_types = {
        "contig": str,
        "position": int,
        "ref_count": int,
        "alt_count": int,
        "other_alt_count": int,
        "allele_frequency": float
    }
    headers = defaultdict(list)
    pileups = defaultdict(list)

    contig_order = get_contigs(args.ref_dict)

    # load data (per sequencing run):
    for sample_name, file_path in zip(args.sample, args.pileup):
        header, df = get_header_and_df(file_path=file_path, columns=columns, column_types=column_types)
        headers[sample_name].append(header)
        pileups[sample_name].append(df)
        print(f"Number of loci in {file_path}: {df.shape[0]}") if args.verbose else None
    print() if args.verbose else None

    # group by sample:
    merged_headers = {}
    merged_pileups = {}
    for sample_name, dfs in pileups.items():
        aggregate = pd.concat(
            dfs, ignore_index=True
        ).groupby(by=["contig", "position"]).agg(
            {
                "ref_count": "sum",
                "alt_count": "sum",
                "other_alt_count": "sum",
                "allele_frequency": "max"
            }
        )
        merged_pileup = aggregate.reindex(sort_genomic_positions(index=aggregate.index, contig_order=contig_order)).reset_index()
        print(f"Number of loci in merged pileup for {sample_name}: {merged_pileup.shape[0]}") if args.verbose else None
        merged_headers[sample_name] = headers[sample_name][0] if len(headers[sample_name]) > 0 else None
        merged_pileups[sample_name] = merged_pileup.astype(column_types)
    print() if args.verbose else None

    # write output files:
    for sample_name, merged_pileup in merged_pileups.items():
        write_header_and_df(
            header=merged_headers[sample_name],
            df=merged_pileup,
            file_path=f"{args.output_dir}/{sample_name}.pileup" + (".gz" if args.compress_output else ""),
            verbose=args.verbose
        )


if __name__ == "__main__":
    main()
