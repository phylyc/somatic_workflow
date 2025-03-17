import argparse
import numpy as np
import pandas as pd
import time
import warnings


def message(*args, **kwargs) -> None:
    print(f"{time.strftime('%H:%M:%S')} ", *args, **kwargs)
    return None


def parse_args():
    parser = argparse.ArgumentParser(
        prog="PileupToAllelicCounts",
        description="""
        """,
        epilog="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.usage = "pileup_to_allelic_counts.py --pileup <pileup> --gvcf <gvcf> --output <output> [--intervals <intervals>] [--select_hets] [--verbose]"
    parser.add_argument("--pileup",         type=str,   required=True,  help="Path to a GATK GetPileupSummaries output-like file.")
    parser.add_argument("--gvcf",           type=str,   required=True,  help="Path to a genotyped germline vcf that contains ref and alt allele information for each pileup locus and contains the GT field.")
    parser.add_argument("-D", "--ref_dict", type=str,   help="Path to the reference dictionary to sort rows.")
    parser.add_argument("--intervals",      type=str,   default=None,   help="Path to the (denoised total copy ratio) intervals to map nearby pileups to intervals. Required columns: CONTIG, START, END")
    parser.add_argument("--het_to_interval_mapping_max_distance", type=int, default=0, help="If a pileup location does not map to an interval provided in the intervals, it will be mapped to the closest interval, which is at most this many base pairs away.")
    parser.add_argument("--aggregate_hets", default=False, action="store_true", help="Aggregate pileups in intervals.")
    parser.add_argument("--min_read_depth", type=int,   default=0, help="Minimum read depth.")
    parser.add_argument("--output",         type=str,   required=True,  help="Path to the output file.")
    parser.add_argument("--error_output",   type=str,                   help="Path to the error rate output file.")
    parser.add_argument("--select_hets",    default=False, action="store_true", help="Keep only heterozygous sites.")
    parser.add_argument("--verbose",        default=False, action="store_true", help="Print information to stdout during execution.")
    return parser.parse_args()


def main():
    args = parse_args()
    print_args(args)
    convert_pileup_to_allelic_counts(args=args)


def print_args(args):
    if args.verbose:
        message("Calling PileupToAllelicCounts")
        print("Arguments:")
        for key, value in vars(args).items():
            print(f"  {key}: {value}")
        print()


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


def convert_pileup_to_allelic_counts(args):
    pileup = pd.read_csv(
        f"{args.pileup}", sep="\t", comment="#", low_memory=False
    ).astype({"contig": str, "position": int, "ref_count": int, "alt_count": int, "other_alt_count": int, "allele_frequency": float})

    gvcf = pd.read_csv(f"{args.gvcf}", sep="\t", comment="#", header=None, low_memory=False)
    gvcf.columns = ["contig", "position", "id", "ref", "alt", "qual", "filter", "info", "format", "genotype"][:gvcf.shape[1]]
    for col, dtype in {"contig": str, "position": int, "id": str, "ref": str, "alt": str, "info": str, "genotype": str}.items():
        if col in gvcf.columns:
            gvcf = gvcf.astype({col: dtype})

    intervals = pd.read_csv(
        f"{args.intervals}", sep="\t", comment="@", low_memory=False
    ).astype({"CONTIG": str, "START": int, "END": int}) if args.intervals is not None else None

    contig_order = get_contigs(args.ref_dict)

    if pileup.empty:
        print("No pileups found.")
        return None

    df = pd.merge(pileup, gvcf, how="inner", on=["contig", "position"])
    df = df.loc[df[["ref_count", "alt_count"]].sum(axis=1) >= args.min_read_depth]
    df.drop_duplicates(subset=["contig", "position"], inplace=True)

    if args.select_hets:
        if "genotype" in df.columns:
            het_mask = df["genotype"].str.contains(r"0[/|]1|1[/|]0")
            print(f"Selecting {het_mask.sum()} / {het_mask.shape[0]} HETs")
            df = df.loc[het_mask]
        else:
            args.select_hets = False
            print(f"Skipping HET selection since no genotype data is available.")

    if intervals is not None:
        # 1) Map each pileup locus to an interval, choosing the closest interval
        #    up to args.padding distance away from an interval start/end.
        # 2) Aggregate the allelic read counts for all pileups mapped to each interval

        dfs = []
        for contig in df["contig"].unique():
            _df = df.loc[df["contig"] == contig]
            i_contig = intervals.loc[intervals["CONTIG"] == contig]
            if i_contig.empty:
                continue

            # Find positions in intervals:
            _df = pd.merge_asof(_df, i_contig[["START", "END"]], left_on="position", right_on="START", direction="backward")
            _df = pd.merge_asof(_df, i_contig[["START", "END"]], left_on="position", right_on="START", direction="forward", suffixes=("_pre", "_post"))

            # Compute distance to both intervals
            _df["distance_pre"] = (_df["position"] - _df["END_pre"]).clip(lower=0)
            _df["distance_post"] = _df["START_post"] - _df["position"]

            # Get smallest distance and drop pileups too far away.
            _df["distance"] = _df[["distance_pre", "distance_post"]].min(axis=1)
            _df = _df.loc[_df["distance"] <= args.het_to_interval_mapping_max_distance]

            # Select the closest interval
            use_post = _df["distance_post"] < _df["distance_pre"]
            _df["START"] = np.where(use_post, _df["START_post"], _df["START_pre"])
            _df["END"] = np.where(use_post, _df["END_post"], _df["END_pre"])

            mapped_hets = _df["distance"] > 0
            # Move pileup position to center of interval so that ModelSegments doesn't throw it away
            _df.loc[mapped_hets, "position"] = (_df["START"] + _df["END"]) // 2

            # Aggregate read counts per interval
            if args.aggregate_hets:
                if not args.select_hets:
                    warnings.warn("Aggregating pileups in intervals, but not explicitly pre-selected for HETs!")

                # TODO: infer relative haplotype phase of HETs within each interval
                # Account for haplotype phase if available
                if "genotype" in _df.columns:
                    phased_ref_count = np.where(_df["genotype"] == "1|0", _df["alt_count"], _df["ref_count"])
                    _df["alt_count"] = np.where(_df["genotype"] == "1|0", _df["ref_count"], _df["alt_count"])
                    _df["ref_count"] = phased_ref_count

                _df = _df.groupby(["contig", "START", "END"]).agg(
                    ref_count=("ref_count", "sum"),
                    alt_count=("alt_count", "sum"),
                    ref=("ref", "first"),
                    alt=("alt", "first")
                ).reset_index()

                _df["position"] = (_df["START"] + _df["END"]) // 2

            dfs.append(_df)

        df = pd.concat(dfs).astype({"contig": str, "position": int, "ref_count": int, "alt_count": int})

        print(f"Remaining HETs after aggregating and mapping to intervals: {df.shape[0]}")

    df["depth"] = df[["ref_count", "alt_count"]].sum(axis=1)
    df = df.sort_values(by=["depth"], ascending=False).drop_duplicates(subset=["contig", "position"]).set_index(["contig", "position"], drop=True)
    df = df.reindex(sort_genomic_positions(index=df.index, contig_order=contig_order)).reset_index()
    df[["contig", "position", "ref_count", "alt_count", "ref", "alt"]].to_csv(f"{args.output}", sep="\t", index=False, header=False, mode="a")

    if args.error_output is not None:
        other_alt_counts = pileup["other_alt_count"].sum()
        total_counts = pileup[["ref_count", "alt_count", "other_alt_count"]].sum(axis=1).sum()
        error_probability = np.clip(1.5 * other_alt_counts / max(1, total_counts), a_min=0.001, a_max=0.05)
        np.savetxt(f"{args.error_output}", [error_probability])


if __name__ == "__main__":
    main()
