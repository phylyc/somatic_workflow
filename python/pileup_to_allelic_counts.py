import argparse
import numpy as np
import pandas as pd
import time
import warnings
import gzip


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
    parser.add_argument("--phase_q_tolerance", type=float, default=0.2, help="")
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


def read_gvcf(filename: str):
    if filename.endswith(".gz") or filename.endswith(".bgz"):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")


def convert_pileup_to_allelic_counts(args):
    pileup = pd.read_csv(
        f"{args.pileup}", sep="\t", comment="#", low_memory=False
    ).astype({"contig": str, "position": int, "ref_count": int, "alt_count": int, "other_alt_count": int, "allele_frequency": float})

    with read_gvcf(args.gvcf) as gvcf_file:
        gvcf = pd.read_csv(gvcf_file, sep="\t", comment="#", header=None, low_memory=False)
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
            print(f"Selecting {het_mask.sum()} HETs / {het_mask.shape[0]} sites.")
            df = df.loc[het_mask]
        else:
            args.select_hets = False
            print(f"Skipping HET selection since no genotype data is available.")

    if intervals is not None:
        # Map each pileup locus to an interval, choosing the closest interval up to
        # args.het_to_interval_mapping_max_distance distance away from an interval start/end.

        if args.aggregate_hets:
            if not args.select_hets:
                warnings.warn("Aggregating pileups in intervals, but not explicitly pre-selected for HETs!")

        message("Mapping nearby pileups to intervals:", end="", flush=True)
        dfs = []
        for contig in df["contig"].unique():
            _df = df.loc[df["contig"] == contig]
            i_contig = intervals.loc[intervals["CONTIG"] == contig]
            if _df.empty or i_contig.empty:
                continue

            _df = _df.sort_values("position")
            i_contig = i_contig.sort_values("START")

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

            # Move pileup into unoccupied position in interval so that ModelSegments doesn't throw it away
            pileups = []
            for (start, end), group in _df.groupby(["START", "END"]):
                mapped_hets = group["distance"] > 0
                pileups_within_interval = group.loc[~mapped_hets]
                pileups.append(pileups_within_interval)
                if not any(mapped_hets):
                    continue

                occupied_positions = pileups_within_interval["position"].to_numpy()
                # Map to the last free positions since ModelSegments uses the
                # first appearing HET, and HETs within the interval are a bit
                # more informative than HETs mapped to the interval from nearby.
                free_positions = [p for p in np.arange(end - 1, start, -1) if p not in occupied_positions]
                n_mapped = min(len(free_positions), np.sum(mapped_hets))
                if n_mapped > 0:
                    mapped_positions = free_positions[:n_mapped]
                    for new_pos, (_, p) in zip(mapped_positions, group.loc[mapped_hets].iterrows()):
                        p["position"] = new_pos
                        pileups.append(p.to_frame(0).T)

            if len(pileups):
                _df = pd.concat(pileups, ignore_index=True)

                # Aggregate read counts per interval
                if args.aggregate_hets:
                    # Account for haplotype phase if available
                    if "0|1" in _df.columns and "1|0" in _df.columns:
                        p01 = _df["0|1"].to_numpy(dtype=float)
                        w   = _df["0/1"].to_numpy(dtype=float)

                        q = np.full_like(w, 0.5)  # default to 0.5 when w==0
                        nz = w > 0
                        q[nz] = p01[nz] / w[nz]
                        q[np.abs(q - 0.5) < args.phase_q_tolerance] = 0.5

                        a = _df["alt_count"].astype(float).to_numpy()
                        r = _df["ref_count"].astype(float).to_numpy()

                        ref_oriented = q * r + (1 - q) * a
                        alt_oriented = q * a + (1 - q) * r
                        _df = _df.assign(ref_count=ref_oriented, alt_count=alt_oriented)

                    elif "genotype" in _df.columns:
                        phased_ref_count = np.where(_df["genotype"] == "1|0", _df["alt_count"], _df["ref_count"])
                        phased_alt_count = np.where(_df["genotype"] == "1|0", _df["ref_count"], _df["alt_count"])
                        _df = _df.assign(ref_count=phased_ref_count, alt_count=phased_alt_count)

                    else:
                        warnings.warn("Aggregating pileups not accounting for haplotypes!")

                    _df = _df.groupby(by=["contig", "START", "END"], sort=False).agg(
                        ref_count=("ref_count", "sum"),
                        alt_count=("alt_count", "sum"),
                        ref=("ref", "first"),
                        alt=("alt", "first")
                    ).reset_index()

                    _df["position"] = (_df["START"] + _df["END"]) // 2
                    _df["ref_count"] = _df["ref_count"].astype(int)
                    _df["alt_count"] = _df["alt_count"].astype(int)

                    print("+", end="", flush=True) if args.verbose else None

                if not _df.empty:
                    dfs.append(_df)

            print(".", end="", flush=True) if args.verbose else None

        print("")

        if len(dfs):
            df = pd.concat(dfs).astype({"contig": str, "position": int, "ref_count": int, "alt_count": int})
        else:
            df = pd.DataFrame(columns=["contig", "position", "ref_count", "alt_count", "ref", "alt"])

        print(f"Remaining pileups after mapping to intervals: {df.shape[0]}")

    if not df.empty:
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
