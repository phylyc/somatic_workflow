import argparse
import numpy as np
import pandas as pd
import os
import scipy.stats as st


def parse_args():
    parser = argparse.ArgumentParser(
        prog="MapToAbsoluteCopyNumber",
        description="""
            Adds back segments 
        """,
        epilog="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.usage = "map_to_absolute_copy_number.py [--sample <sample>] --absolute_seg <absolute_seg> --cr_seg <cr_seg> --gvcf <gvcf> --purity <purity> --ploidy <ploidy> --output <output>"
    parser.add_argument("--sample",         type=str,   required=False, help="Sample name.")
    parser.add_argument("--sex",            type=str,   required=False, help="Patient's sex genotype.")
    parser.add_argument("--absolute_seg",   type=str,   required=True,  help="Path to a ABSOLUTE segtab output file.")
    parser.add_argument("--cr_seg",         type=str,   required=True,  help="Path to a ACS segmentation output file.")
    parser.add_argument("--gvcf",           type=str,   required=True,  help="Path to a genotyped germline vcf that contains ref and alt allele information for each pileup locus and contains the GT field.")
    parser.add_argument("--purity",         type=float, required=True,  help="Tumor purity as inferred by ABSOLUTE")
    parser.add_argument("--ploidy",         type=float, required=True,  help="Tumor ploidy as inferred by ABSOLUTE")
    parser.add_argument("--outdir",         type=str,   required=True,  help="Path to the output directory to write the extended segmentations.")
    return parser.parse_args()


def main():
    args = parse_args()
    map_to_cn(args=args)


def map_to_cn(args):
    abs_dtypes = {
        "sample": str,
        "Chromosome": str, "Start.bp": int, "End.bp": int,
        "n_probes": int, "length": int, "seg_sigma": float, "W": float,
        "total_copy_ratio": float, "modal_total_cn": int, "expected_total_cn": float, "total_HZ": int, "total_amp": int, "corrected_total_cn": float, "rescaled_total_cn": float,
        "bi.allelic": int, "copy.ratio": float,
        "hscr.a1": float, "hscr.a2": float,
        "modal.a1": int, "modal.a2": int,
        "expected.a1": float, "expected.a2": float,
        "subclonal.a1": int, "subclonal.a2": int,
        "cancer.cell.frac.a1": float, "ccf.ci95.low.a1": float, "ccf.ci95.high.a1": float,
        "cancer.cell.frac.a2": float, "ccf.ci95.low.a2": float, "ccf.ci95.high.a2": float,
        "LOH": int, "HZ": int, "SC_HZ": int,
        "amp.a1": int, "amp.a2": int,
        "rescaled.cn.a1": float, "rescaled.cn.a2": float
    }
    acs_dtypes = {
        "Chromosome": str, "Start.bp": int, "End.bp": int,
        "n_probes": int, "length": int, "n_hets": int,
        "f": float, "tau": float, "sigma.tau": float,
        "mu.minor": float, "sigma.minor": float, "mu.major": float, "sigma.major": float,
        "SegLabelCNLOH": int
    }
    vcf_columns = ["contig", "position", "id", "ref", "alt", "qual", "filter", "info", "format", "genotype"]

    abs_seg = pd.read_csv(
        f"{args.absolute_seg}", sep="\t", comment="#", low_memory=False
    ).astype(abs_dtypes)
    abs_seg_cols = abs_seg.columns
    abs_seg = abs_seg.set_index(["Chromosome", "Start.bp", "End.bp"])
    cr_seg = pd.read_csv(
        f"{args.cr_seg}", sep="\t", comment="@", low_memory=False
    ).astype(acs_dtypes).set_index(["Chromosome", "Start.bp", "End.bp"])
    gvcf = pd.read_csv(
        f"{args.gvcf}", sep="\t", comment="#", header=None, low_memory=False, names=vcf_columns
    ).astype({"contig": str, "position": int, "id": str, "ref": str, "alt": str, "info": str, "genotype": str})

    seg = pd.concat([abs_seg.drop(columns=["n_probes", "length"]), cr_seg], axis=1).sort_index()

    if args.sample is None:
        args.sample = seg["Sample"].dropna().unique()[0]

    D = args.purity * args.ploidy + (1 - args.purity) * 2
    b = (1 - args.purity) * 2 / D
    delta = args.purity / D

    # rescale copy number
    seg["CN"] = (seg["tau"] - b).clip(lower=0) / delta / 2
    # correct offset
    diff = seg["CN"].median() - seg["rescaled_total_cn"].median()
    seg["CN"] -= diff
    seg["CN"] = seg["CN"].clip(lower=0)
    seg["CN.sigma"] = seg["sigma.tau"] / delta / 2

    r_corr = seg[["CN", "rescaled_total_cn"]].corr().loc["CN", "rescaled_total_cn"]
    print(f"Correlation between rescaled total copy number and ABSOLUTE output: {r_corr}")

    # map to cluster
    cluster_values = seg["corrected_total_cn"].round(1).dropna().unique()
    is_integer = np.modf(cluster_values)[0] == 0

    def map_to_cluster(row, p_threshold=0.05, log_odds_ratio=-1):
        log_p_threshold = np.log(p_threshold)
        norm = st.norm(loc=row["CN"], scale=row["CN.sigma"])
        logcdf = norm.logcdf(cluster_values)
        logsf = norm.logsf(cluster_values)

        is_valid_cluster = (logcdf > log_p_threshold) & (logsf > log_p_threshold)
        valid_clusters = cluster_values[is_valid_cluster]
        if not len(valid_clusters):
            return row["CN"]

        valid_cluster_is_integer = is_integer[is_valid_cluster]
        logpdf = norm.logpdf(valid_clusters)
        if 0 < np.sum(valid_cluster_is_integer) < len(valid_cluster_is_integer):
            logpdf_int_idx = np.argmax(logpdf[valid_cluster_is_integer])
            logpdf_frac_idx = np.argmax(logpdf[~valid_cluster_is_integer])
            if logpdf[valid_cluster_is_integer][logpdf_int_idx] - logpdf[~valid_cluster_is_integer][logpdf_frac_idx] > log_odds_ratio:
                return valid_clusters[valid_cluster_is_integer][logpdf_int_idx]
            else:
                return valid_clusters[~valid_cluster_is_integer][logpdf_frac_idx]
        else:
            logpdf_idx = np.argmax(logpdf)
            return valid_clusters[logpdf_idx]

    seg["corrected_CN"] = seg.apply(map_to_cluster, axis=1)

    c_corr = seg[["corrected_CN", "corrected_total_cn"]].corr().loc["corrected_CN", "corrected_total_cn"]
    print(f"Correlation between corrected total copy number and ABSOLUTE output: {c_corr}")

    new_segs = seg["corrected_total_cn"].isna()
    seg.loc[new_segs, "sample"] = args.sample
    seg.loc[new_segs, "total_copy_ratio"] = seg.loc[new_segs, "tau"] / 2
    seg.loc[new_segs, "copy.ratio"] = seg.loc[new_segs, "tau"] / 2
    # seg.loc[new_segs, "hscr.a1"] = seg.loc[new_segs, "tau"] * seg.loc[new_segs, "f"]
    # seg.loc[new_segs, "hscr.a2"] = seg.loc[new_segs, "tau"] * (1 - seg.loc[new_segs, "f"])
    seg.loc[new_segs, "rescaled_total_cn"] = seg.loc[new_segs, "CN"]
    seg.loc[new_segs, "corrected_total_cn"] = seg.loc[new_segs, "corrected_CN"]
    seg.loc[new_segs, "modal_total_cn"] = seg.loc[new_segs, "corrected_CN"].round()
    seg.loc[new_segs, "LOH"] = 1
    seg.loc[new_segs, "rescaled.cn.a1"] = 0
    seg.loc[new_segs, "rescaled.cn.a2"] = seg.loc[new_segs, "CN"]

    print(f"Number of rescued segments: {new_segs.sum()}")

    def sort_genomic_positions(index: pd.MultiIndex) -> pd.MultiIndex:
        by = ["Chromosome", "Start.bp", "End.bp"]
        contig_order = (
            [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
            + [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
        )
        contig_order += list(sorted(set(index.get_level_values("Chromosome")) - set(contig_order)))
        temp_df = pd.DataFrame(index=index).reset_index().astype({c: t for c, t in abs_dtypes.items() if c in index.names})
        temp_df["Chromosome"] = pd.Categorical(temp_df["Chromosome"], categories=contig_order, ordered=True)
        temp_df.sort_values(by=by, inplace=True)
        temp_df = temp_df.astype({c: t for c, t in abs_dtypes.items() if c in temp_df.columns})
        return pd.MultiIndex.from_frame(temp_df[index.names])

    seg = seg.reindex(sort_genomic_positions(index=seg.index))

    seg["Segment_Mean"] = np.log2(seg["rescaled_total_cn"] + 1e-4) - np.log2(args.ploidy)
    if args.sex in ["XY", "Male"]:
        seg.loc[seg["Chromosome"].isin(["X", "Y"]), "Segment_Mean"] += 1

    seg.loc[new_segs].reset_index()[["Chromosome", "Start.bp", "End.bp"]].to_csv(f"{args.outdir}/{args.sample}.rescued_intervals.txt", sep="\t", index=False)

    seg = seg.reset_index()

    seg[abs_seg_cols].to_csv(f"{args.outdir}/{args.sample}.segtab.completed.txt", sep="\t", index=False)
    seg[["sample", "Chromosome", "Start.bp", "End.bp", "Segment_Mean", "rescaled_total_cn"]].to_csv(f"{args.outdir}/{args.sample}.IGV.seg.completed.txt", sep="\t", index=False)


if __name__ == "__main__":
    main()
