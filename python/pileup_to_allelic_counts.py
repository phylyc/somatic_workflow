import argparse
import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        prog="PileupToAllelicCounts",
        description="""
        """,
        epilog="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.usage = "pileup_to_allelic_counts.py --pileup <pileup> --gvcf <gvcf> --output <output> [--select_hets]"
    parser.add_argument("--pileup",         type=str,   required=True,  help="Path to a GATK GetPileupSummaries output-like file.")
    parser.add_argument("--gvcf",           type=str,   required=True,  help="Path to a genotyped germline vcf that contains ref and alt allele information for each pileup locus and contains the GT field.")
    parser.add_argument("--intervals",      type=str,   default=None,   help="Path to the (denoised total copy ratio) intervals to aggregate pileups into one allelic count. Required columns: CONTIG, START, END")
    parser.add_argument("--het_to_interval_mapping_max_distance", type=int, default=0, help="If a pileup location does not map to an interval provided in the intervals, it will be mapped to the closest interval, which is at most this many base pairs away.")
    parser.add_argument("--min_read_depth", type=int,   default=0, help="Minimum read depth.")
    parser.add_argument("--output",         type=str,   required=True,  help="Path to the output file.")
    parser.add_argument("--error_output",   type=str,                   help="Path to the error rate output file.")
    parser.add_argument("--select_hets",    default=False, action="store_true", help="Keep only heterozygous sites.")
    return parser.parse_args()


def main():
    args = parse_args()
    convert_pileup_to_allelic_counts(args=args)


def convert_pileup_to_allelic_counts(args):
    vcf_columns = ["contig", "position", "id", "ref", "alt", "qual", "filter", "info", "format", "genotype"]

    pileup = pd.read_csv(
        f"{args.pileup}", sep="\t", comment="#", low_memory=False
    ).astype({"contig": str, "position": int, "ref_count": int, "alt_count": int, "other_alt_count": int, "allele_frequency": float})
    gvcf = pd.read_csv(
        f"{args.gvcf}", sep="\t", comment="#", header=None, low_memory=False, names=vcf_columns
    ).astype({"contig": str, "position": int, "id": str, "ref": str, "alt": str, "info": str, "genotype": str})
    intervals = pd.read_csv(
        f"{args.intervals}", sep="\t", comment="@", low_memory=False
    ).astype({"CONTIG": str, "START": int, "END": int}) if args.intervals is not None else None

    if pileup.empty:
        print("No pileups found.")
        return None

    df = pd.merge(pileup, gvcf, how="inner", on=["contig", "position"])
    df = df.loc[df[["ref_count", "alt_count"]].sum(axis=1) >= args.min_read_depth]

    if args.select_hets:
        het_mask = df["genotype"] == "0/1"
        print(f"Selecting {het_mask.sum()} / {het_mask.shape[0]} HETs")
        df = df.loc[het_mask]

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

            # Aggregate read counts per interval
            _df = _df.groupby(["contig", "START", "END"]).agg(
                ref_count=("ref_count", "sum"),
                alt_count=("alt_count", "sum"),
                ref=("ref", "first"),
                alt=("alt", "first")
            ).reset_index()

            # Move pileups to center of interval
            _df["position"] = (_df["START"] + _df["END"]) // 2

            dfs.append(_df)

        df = pd.concat(dfs).astype({"contig": str, "position": int, "ref_count": int, "alt_count": int})

        print(f"Remaining HETs after aggregating and mapping to intervals: {df.shape[0]}")

    df[["contig", "position", "ref_count", "alt_count", "ref", "alt"]].to_csv(f"{args.output}", sep="\t", index=False, header=False, mode="a")

    if args.error_output is not None:
        other_alt_counts = pileup["other_alt_count"].sum()
        total_counts = pileup[["ref_count", "alt_count", "other_alt_count"]].sum(axis=1).sum()
        error_probability = np.clip(1.5 * other_alt_counts / max(1, total_counts), a_min=0.001, a_max=0.05)
        np.savetxt(f"{args.error_output}", [error_probability])


if __name__ == "__main__":
    main()
