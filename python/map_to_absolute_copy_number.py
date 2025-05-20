import argparse
import numpy as np
import pandas as pd
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
    parser.add_argument("--sex",            type=str,   default="XXY",  help="Patient's sex genotype.")
    parser.add_argument("--absolute_seg",   type=str,   required=True,  help="Path to a ABSOLUTE segtab output file.")
    parser.add_argument("--cr_seg",         type=str,   required=True,  help="Path to a ACS segmentation output file.")
    parser.add_argument("--gvcf",           type=str,   required=False, help="Path to a genotyped germline vcf that contains ref and alt allele information for each pileup locus and contains the GT field.")
    parser.add_argument("--purity",         type=float, required=True,  help="Tumor purity as inferred by ABSOLUTE")
    parser.add_argument("--ploidy",         type=float, required=True,  help="Tumor ploidy as inferred by ABSOLUTE")
    parser.add_argument("--outdir",         type=str,   required=True,  help="Path to the output directory to write the extended segmentations.")
    return parser.parse_args()


def main():
    args = parse_args()
    map_to_cn(args=args)


def map_to_cn(args):
    if args.sex in ["Female", "female"]:
        args.sex = "XX"
    if args.sex in ["Male", "male"]:
        args.sex = "XY"

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

    abs_seg = pd.read_csv(f"{args.absolute_seg}", sep="\t", comment="#", low_memory=False)
    for col, dtype in abs_dtypes.items():
        try:
            abs_seg = abs_seg.astype({col: dtype})
        except:
            continue
    abs_seg_cols = abs_seg.columns
    abs_seg = abs_seg.set_index(["Chromosome", "Start.bp", "End.bp"])

    cr_seg = pd.read_csv(f"{args.cr_seg}", sep="\t", comment="@", low_memory=False)
    for col, dtype in acs_dtypes.items():
        try:
            cr_seg = cr_seg.astype({col: dtype})
        except:
            continue
    cr_seg = cr_seg.set_index(["Chromosome", "Start.bp", "End.bp"])

    gvcf = None
    if args.gvcf is not None:
        gvcf = pd.read_csv(
            f"{args.gvcf}", sep="\t", comment="#", header=None, low_memory=False, names=vcf_columns
        ).astype({"contig": str, "position": int, "id": str, "ref": str, "alt": str, "info": str, "genotype": str})

    seg = pd.concat([abs_seg.drop(columns=["n_probes", "length"]), cr_seg], axis=1).sort_index()

    nX = args.sex.count("X")
    nY = args.sex.count("Y")

    def correct_diploid_assumtion(col):
        seg.loc[seg["Chromosome"].isin(["X", "chrX"]), col] *= 2 / nX
        seg.loc[seg["Chromosome"].isin(["Y", "chrY"]), col] *= 2 / nY if nY > 0 else 0

    if args.sample is None:
        args.sample = seg["Sample"].dropna().unique()[0]

    D = args.purity * args.ploidy + (1 - args.purity) * 2
    b = (1 - args.purity) * 2 / D
    delta = args.purity / D

    # rescale copy number, assuming diploidy
    seg["CN"] = (seg["tau"] - b).clip(lower=0) / delta / 2
    correct_diploid_assumtion(col="CN")
    # correct offset
    diff = seg["CN"].median() - seg["rescaled_total_cn"].median()
    seg["CN"] -= diff
    seg["CN"] = seg["CN"].clip(lower=0)
    seg["CN.sigma"] = seg["sigma.tau"] / delta / 2
    correct_diploid_assumtion(col="CN.sigma")

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

    new_segs = seg[["corrected_total_cn", "rescaled_total_cn"]].isna().any(axis=1)
    seg.loc[new_segs, "sample"] = args.sample
    seg.loc[new_segs, "total_copy_ratio"] = seg.loc[new_segs, "tau"] / 2
    seg.loc[new_segs, "copy.ratio"] = seg.loc[new_segs, "tau"] / 2
    correct_diploid_assumtion(col="total_copy_ratio")
    correct_diploid_assumtion(col="copy.ratio")

    # seg.loc[new_segs, "hscr.a1"] = seg.loc[new_segs, "tau"] * seg.loc[new_segs, "f"]
    # seg.loc[new_segs, "hscr.a2"] = seg.loc[new_segs, "tau"] * (1 - seg.loc[new_segs, "f"])
    seg.loc[new_segs, "rescaled_total_cn"] = seg.loc[new_segs, "CN"]
    seg.loc[new_segs, "expected_total_cn"] = seg.loc[new_segs, "CN"]
    seg.loc[new_segs, "corrected_total_cn"] = seg.loc[new_segs, "corrected_CN"]
    seg.loc[new_segs, "modal_total_cn"] = seg.loc[new_segs, "corrected_CN"].round()
    seg.loc[new_segs, "total_HZ"] = (seg.loc[new_segs, "rescaled_total_cn"] == 0).astype(int)
    # extracted from default settings in the ABSOLUTE package
    seg.loc[new_segs, "total_amp"] = (seg.loc[new_segs, "rescaled_total_cn"] >= 7).astype(int)
    seg.loc[new_segs, "rescaled.cn.a1"] = 0
    seg.loc[new_segs, "rescaled.cn.a2"] = seg.loc[new_segs, "CN"]
    seg.loc[new_segs, "LOH"] = ((seg.loc[new_segs, "rescaled.cn.a1"] == 0) | (seg.loc[new_segs, "rescaled.cn.a1"] == 0)).astype(int)
    seg.loc[new_segs, "HZ"] = ((seg.loc[new_segs, "rescaled.cn.a1"] == 0) & (seg.loc[new_segs, "rescaled.cn.a1"] == 0)).astype(int)

    seg["W"] = seg["length"] / seg["length"].sum()


    # sample	Chromosome	Start.bp	End.bp	n_probes	length	seg_sigma	W	total_copy_ratio	modal_total_cn	expected_total_cn	total_HZ	total_amp	corrected_total_cn	rescaled_total_cn	bi.allelic	copy.ratio	hscr.a1	hscr.a2	modal.a1	modal.a2	expected.a1	expected.a2	subclonal.a1	subclonal.a2	cancer.cell.frac.a1	ccf.ci95.low.a1	ccf.ci95.high.a1	cancer.cell.frac.a2	ccf.ci95.low.a2	ccf.ci95.high.a2	LOH	HZ	SC_HZ	amp.a1	amp.a2	rescaled.cn.a1	rescaled.cn.a2

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

    seg.loc[new_segs].reset_index()[["Chromosome", "Start.bp", "End.bp"]].to_csv(f"{args.outdir}/{args.sample}.rescued_intervals.txt", sep="\t", index=False)

    seg = seg.reset_index()
    seg["Segment_Mean"] = np.log2(seg["rescaled_total_cn"].clip(lower=1e-2)) - np.log2(args.ploidy)
    seg[abs_seg_cols].to_csv(f"{args.outdir}/{args.sample}.segtab.completed.txt", sep="\t", index=False)
    seg[["sample", "Chromosome", "Start.bp", "End.bp", "Segment_Mean", "rescaled_total_cn"]].to_csv(f"{args.outdir}/{args.sample}.IGV.seg.completed.txt", sep="\t", index=False)


if __name__ == "__main__":
    main()
