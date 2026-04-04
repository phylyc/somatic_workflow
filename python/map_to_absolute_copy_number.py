import argparse
import numpy as np
import pandas as pd
import os
import scipy.special as sp
import scipy.stats as st
import time
import warnings

warnings.filterwarnings('ignore')


def message(*args, **kwargs) -> None:
    print(f"{time.strftime('%H:%M:%S')} ", *args, **kwargs)
    return None


def parse_args():
    parser = argparse.ArgumentParser(
        prog="MapToAbsoluteCopyNumber",
        description="""
            Adds back segments 
        """,
        epilog="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.usage = "map_to_absolute_copy_number.py --outdir <outdir> --purity <purity> --ploidy <ploidy> [--acs_cr_seg <acs_cr_seg> / --somix_cr_seg <somix_cr_seg>] [--absolute_segtab <absolute_segtab>] [--sample <sample>] [--sex <sex>] [--normal_ploidy <normal_ploidy>] [--min_hets <min_hets>] [--min_probes <min_probes>]"
    parser.add_argument("--outdir",             type=str,   required=True,  help="Path to the output directory to write the extended segmentations.")
    parser.add_argument("--purity",             type=float, required=True,  help="Tumor purity as inferred by ABSOLUTE")
    parser.add_argument("--ploidy",             type=float, required=True,  help="Tumor ploidy as inferred by ABSOLUTE")
    parser.add_argument("--sample",             type=str,   required=False, help="Sample name.")
    parser.add_argument("--sex",                type=str,   default="XXY",  help="Patient's sex genotype.")
    parser.add_argument("--acs_cr_seg",         type=str,   required=False, help="Path to a ACS segmentation output file.")
    parser.add_argument("--somix_cr_seg",       type=str,   required=False, help="Path to a somix segmentation output file.")
    parser.add_argument("--absolute_segtab",    type=str,   required=False, help="Path to a ABSOLUTE segtab output file.")
    parser.add_argument("--normal_ploidy",      type=int,   required=False, default=2, help="Normal/germline ploidy of that organism.")
    parser.add_argument("--min_hets",           type=int,   default=0,      help="Minimum number of heterozygous sites for AllelicCapSeg to call a segment.")
    parser.add_argument("--min_probes",         type=int,   default=0,      help="Minimum number of target intervals for AllelicCapSeg to call a segment.")
    return parser.parse_args()


def main():
    args = parse_args()
    print_args(args)
    map_to_cn(args=args)


def print_args(args):
    message("Calling MapToAbsoluteCopyNumber")
    print("Arguments:")
    for key, value in vars(args).items():
        print(f"  {key}: {value}")
    print()


def map_to_cn(args):
    s = args.sex.upper()
    if s in ["FEMALE"]:
        s = "XX"
    elif s in ["MALE"]:
        s = "XY"
    elif s in ["UNKNOWN"]:
        s = "XXY"

    if args.purity <= 0:
        args.purity = 1
    if args.ploidy <= 0:
        args.ploidy = args.normal_ploidy

    nX = s.count("X")
    nY = s.count("Y")

    def get_chromosomal_ploidy(chr):
        return nX if chr in ["X", "chrX"] else nY if chr in ["Y", "chrY"] else args.normal_ploidy

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

    extra_columns = [
        "ccf_mean",
        "ccf_median",
        "ccf_mode",
        "ccf_ci95_low",
        "ccf_ci95_high",
        "ccf_ci95_width",
        "ccf_entropy",
        "ccf_max_p",
        "ccf_support_50",
        "ccf_support_90",
        "current_pair_entropy",
        "current_pair_effn",
        "tail_mass_high_cn",
        "high_copy_flag",
        "ambiguous_ccf_flag",
        "baseline_q1",
        "baseline_q2",
    ]

    ###########################################################################
    ### LOADING DATA
    ###########################################################################

    try:
        abs_seg = pd.read_csv(f"{args.absolute_segtab}", sep="\t", comment="#", low_memory=False)
    except Exception as e:
        message(e)
        message("Using empty dataframe instead.")
        abs_seg = pd.DataFrame(None, columns=list(abs_dtypes.keys()))

    for col, dtype in abs_dtypes.items():
        if col in abs_seg.columns:
            abs_seg = abs_seg.astype({col: dtype}, errors="ignore")
    abs_seg = abs_seg.set_index(["Chromosome", "Start.bp", "End.bp"])

    if args.acs_cr_seg is not None:
        cr_seg = pd.read_csv(f"{args.acs_cr_seg}", sep="\t", comment="@", low_memory=False)

    elif args.somix_cr_seg is not None:
        cr_seg = pd.read_csv(f"{args.somix_cr_seg}", sep="\t", low_memory=False)
        cr_seg = cr_seg.rename(columns={"contig": "Chromosome", "n_markers": "n_probes", "n_snps": "n_hets", "f_MAP": "f"})
        if args.sample is not None:
            cr_seg = cr_seg.loc[cr_seg["sample_id"] == args.sample]
        chr_ploidy = cr_seg["Chromosome"].map(get_chromosomal_ploidy)
        cr_seg["tau"] = np.exp(cr_seg["log_tCR"]) * chr_ploidy
        cr_seg["sigma.tau"] = chr_ploidy * np.exp(cr_seg["log_tCR"] + cr_seg["sem_log_tCR"]**2 / 2) * np.sqrt(np.exp(cr_seg["sem_log_tCR"]**2) - 1)
        if "f" in cr_seg.columns:
            cr_seg["mu.minor"] = cr_seg["f"] * cr_seg["tau"]
            cr_seg["mu.major"] = (1 - cr_seg["f"]) * cr_seg["tau"]
            var_f = - 1 / cr_seg["f_hessian"]
            cr_seg["sigma.minor"] = np.sqrt(cr_seg["tau"]**2 * var_f + cr_seg["f"]**2 * cr_seg["sigma.tau"]**2)
            cr_seg["sigma.major"] = np.sqrt(cr_seg["tau"]**2 * var_f + (1 - cr_seg["f"])**2 * cr_seg["sigma.tau"]**2)
        cr_seg = cr_seg.reindex(columns=acs_dtypes.keys())

    else:
        message("One of acs_cr_seg or somix_cr_seg input arguments must be defined!")
        return None

    for col, dtype in acs_dtypes.items():
        cr_seg = cr_seg.astype({col: dtype}, errors="ignore")
    cr_seg = cr_seg.set_index(["Chromosome", "Start.bp", "End.bp"])

    seg = pd.concat([abs_seg.drop(columns=["n_probes", "length"], errors="ignore"), cr_seg], axis=1).sort_index().reset_index()

    if args.sample is None:
        args.sample = seg["sample"].dropna().unique()[0]

    if seg.empty:
        message("No segments to map.")
        return None

    ###########################################################################
    ### DEFINING UTILITY FUNCTIONS
    ###########################################################################

    # map to cluster
    cluster_values = seg["corrected_total_cn"].round(1).dropna().unique()
    is_integer = np.modf(cluster_values)[0] == 0

    def map_to_cluster(cn, sigma, p_threshold=0.05, log_odds_ratio=-1):
        log_p_threshold = np.log(p_threshold)
        norm = st.norm(loc=cn, scale=sigma)
        logcdf = norm.logcdf(cluster_values)
        logsf = norm.logsf(cluster_values)

        is_valid_cluster = (logcdf > log_p_threshold) & (logsf > log_p_threshold)
        valid_clusters = cluster_values[is_valid_cluster]
        if not len(valid_clusters):
            return cn

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

    seg["chr_ploidy"] = chr_ploidy = seg["Chromosome"].map(get_chromosomal_ploidy)
    seg["is_parental_haploid"] = (chr_ploidy == 1) & ((nX == 1) | (nY > 0))
    seg["sigma.a1"] = seg["sigma.minor"]
    seg["sigma.a2"] = seg["sigma.major"]
    seg["mu.a1"] = seg["mu.minor"]
    seg["mu.a2"] = seg["mu.major"]

    args.purity = args.purity if args.purity > 0 else 1
    D = (1 - args.purity) * chr_ploidy + args.purity * args.ploidy * chr_ploidy / args.normal_ploidy
    b = (1 - args.purity) * chr_ploidy / D
    delta = args.purity / D

    seg["W"] = seg["length"] / seg["length"].sum()

    # correct offset: Find alpha for
    # CN = (tau - b - alpha) / delta / c
    # sum(w * CN) = ploidy
    num = np.where(chr_ploidy > 0, seg["W"] * (seg["tau"] - b) / delta / chr_ploidy, 0)
    den = np.where(chr_ploidy > 0, seg["W"] / delta / chr_ploidy, 0)
    alpha = (np.sum(num) - args.ploidy) / np.sum(den)

    message(f"Shift total copy number (tau) by {-alpha} to fit onto comb.")
    seg["tau"] -= alpha
    seg["tau"] = seg["tau"].clip(lower=0)

    # rescale copy number
    seg["CN"] = np.where(chr_ploidy > 0, (seg["tau"] - b).clip(lower=0) / delta / chr_ploidy, 0)
    seg["CN.sigma"] = np.where(chr_ploidy > 0, seg["sigma.tau"] / delta / chr_ploidy, 0)
    seg["corrected_CN"] = seg.apply(lambda row: map_to_cluster(row["CN"], row["CN.sigma"]), axis=1)

    r_corr = seg[["CN", "rescaled_total_cn"]].corr().loc["CN", "rescaled_total_cn"]
    diff = np.abs(seg['CN'] - seg['rescaled_total_cn'])
    message(f"Correlation between rescaled total copy number and ABSOLUTE output:       {r_corr:.3f}, "
            f"median |diff| = {np.nanmedian(diff):.3f} ± {np.nanstd(diff):.3f}")

    c_corr = seg[["corrected_CN", "corrected_total_cn"]].corr().loc["corrected_CN", "corrected_total_cn"]
    diff = np.abs(seg['corrected_CN'] - seg['corrected_total_cn'])
    message(f"Correlation between corrected total copy number and ABSOLUTE output:      {c_corr:.3f}, "
            f"median |diff| = {np.nanmedian(diff):.3f} ± {np.nanstd(diff):.3f}")

    # RESCUE SEGMENTS
    # NOTE: This may (very likely) also "rescue" artifactual homozygous deletions.
    # If you have a cohort of samples, check for recurrent segment boundaries and
    # un-rescue segments within boundaries that are shared for >10% of the cohort.

    new_segs = seg[["corrected_total_cn", "rescaled_total_cn"]].isna().any(axis=1)
    seg["is_new"] = new_segs
    seg.loc[new_segs, "sample"] = args.sample
    seg.loc[new_segs, "total_copy_ratio"] = np.where(chr_ploidy.loc[new_segs] > 0, seg.loc[new_segs, "tau"].div(chr_ploidy.loc[new_segs]).fillna(0), 0)
    seg.loc[new_segs, "copy.ratio"] = np.where(chr_ploidy.loc[new_segs] > 0, seg.loc[new_segs, "tau"].div(chr_ploidy.loc[new_segs]).fillna(0), 0)
    seg.loc[new_segs, "seg_sigma"] = np.where(chr_ploidy.loc[new_segs] > 0, seg.loc[new_segs, "sigma.tau"].div(chr_ploidy.loc[new_segs]).fillna(0), 0)

    seg.loc[new_segs, "rescaled_total_cn"] = seg.loc[new_segs, "CN"]
    seg.loc[new_segs, "expected_total_cn"] = seg.loc[new_segs, "CN"]
    seg.loc[new_segs, "corrected_total_cn"] = seg.loc[new_segs, "corrected_CN"]
    seg.loc[new_segs, "modal_total_cn"] = seg.loc[new_segs, "corrected_CN"].round()
    seg.loc[new_segs, "total_HZ"] = (seg.loc[new_segs, "rescaled_total_cn"] == 0).astype(int)
    # extracted from default settings in the ABSOLUTE package
    seg.loc[new_segs, "total_amp"] = (seg.loc[new_segs, "rescaled_total_cn"] >= 7).astype(int)

    def wmode(values, weights):
        val = np.rint(values).astype(int)
        counts = np.bincount(val, weights, minlength=val.max() + 1)
        candidates = np.flatnonzero(counts == counts.max())
        return candidates[np.argmin(np.abs(candidates - 2))]

    C0_by_chr = seg.groupby("Chromosome").apply(lambda g: wmode(g["modal_total_cn"].values, g["W"].values), include_groups=False).to_dict()

    def split_alleles(row) -> tuple[float, float]:
        CN = row["rescaled_total_cn"]
        c0 = C0_by_chr.get(row["Chromosome"], 2)
        if CN < 1 or row["is_parental_haploid"] or not row["SegLabelCNLOH"]:
            return 0, CN
        elif (c0 % 2 == 0) and (abs(CN - c0) <= 0.25):  # balanced plateau
            return CN / 2, CN / 2
        elif pd.notna(row["mu.minor.abs"]):
            m = np.clip(row["mu.minor.abs"], 0, CN / 2)
            return min(m, CN - m), max(m, CN - m)
        else:
            return min(1, CN - 1), max(1, CN - 1)

    seg["mu.minor.abs"] = seg["f"] * seg["rescaled_total_cn"]
    seg["mu.major.abs"] = (1 - seg["f"]) * seg["rescaled_total_cn"]

    mn, mj = zip(*seg.apply(split_alleles, axis=1))
    mn, mj = np.array(mn), np.array(mj)

    # allele haploid-specific copy ratio (HSCR) from CN using per-chrom ploidy mixture
    seg["mn"], seg["mj"] = mn, mj
    seg["_hscr.a1"] = seg["mn"] * delta * chr_ploidy + b
    seg["_hscr.a2"] = seg["mj"] * delta * chr_ploidy + b

    c_corr_cn_a1 = seg[["rescaled.cn.a1", "mn"]].corr().fillna(1).loc["rescaled.cn.a1", "mn"]
    diff = np.abs(seg['rescaled.cn.a1'] - seg['mn'])
    message(f"Correlation between rescaled allelic copy number a1 and ABSOLUTE output:  {c_corr_cn_a1:.3f}, "
            f"median |diff| = {np.nanmedian(diff):.3f} ± {np.nanstd(diff):.3f}")
    c_corr_cn_a2 = seg[["rescaled.cn.a2", "mj"]].corr().fillna(1).loc["rescaled.cn.a2", "mj"]
    diff = np.abs(seg['rescaled.cn.a2'] - seg['mj'])
    message(f"Correlation between rescaled allelic copy number a2 and ABSOLUTE output:  {c_corr_cn_a2:.3f}, "
            f"median |diff| = {np.nanmedian(diff):.3f} ± {np.nanstd(diff):.3f}")
    c_corr_hscr_a1 = seg[["hscr.a1", "_hscr.a1"]].corr().fillna(1).loc["hscr.a1", "_hscr.a1"]
    diff = np.abs(seg['hscr.a1'] - seg['_hscr.a1'])
    message(f"Correlation between haplotype-specific copy ratio a1 and ABSOLUTE output: {c_corr_hscr_a1:.3f}, "
            f"median |diff| = {np.nanmedian(diff):.3f} ± {np.nanstd(diff):.3f}")
    c_corr_hscr_a2 = seg[["hscr.a2", "_hscr.a2"]].corr().fillna(1).loc["hscr.a2", "_hscr.a2"]
    diff = np.abs(seg['hscr.a2'] - seg['_hscr.a2'])
    message(f"Correlation between haplotype-specific copy ratio a2 and ABSOLUTE output: {c_corr_hscr_a2:.3f}, "
            f"median |diff| = {np.nanmedian(diff):.3f} ± {np.nanstd(diff):.3f}")

    new_segs_idx = new_segs.loc[new_segs].index.values

    seg.loc[new_segs, "rescaled.cn.a1"] = mn[new_segs_idx]
    seg.loc[new_segs, "rescaled.cn.a2"] = mj[new_segs_idx]
    seg.loc[new_segs, "hscr.a1"] = seg.loc[new_segs, "_hscr.a1"].values
    seg.loc[new_segs, "hscr.a2"] = seg.loc[new_segs, "_hscr.a2"].values
    seg.loc[new_segs, "expected.a1"] = seg.loc[new_segs, "rescaled.cn.a1"]
    seg.loc[new_segs, "expected.a2"] = seg.loc[new_segs, "rescaled.cn.a2"]
    corrected_a1 = seg.apply(lambda row: map_to_cluster(row["rescaled.cn.a1"], row["sigma.a1"]), axis=1)
    corrected_a2 = seg.apply(lambda row: map_to_cluster(row["rescaled.cn.a2"], row["sigma.a2"]), axis=1)
    seg.loc[new_segs, "modal.a1"] = corrected_a1.loc[new_segs].round()
    seg.loc[new_segs, "modal.a2"] = corrected_a2.loc[new_segs].round()
    seg.loc[new_segs, "amp.a1"] = (seg.loc[new_segs, "rescaled.cn.a1"] >= 4).astype(int)
    seg.loc[new_segs, "amp.a2"] = (seg.loc[new_segs, "rescaled.cn.a2"] >= 4).astype(int)
    seg.loc[new_segs, "LOH"] = (((seg.loc[new_segs, "rescaled.cn.a1"] == 0) | (seg.loc[new_segs, "rescaled.cn.a2"] == 0)) & ~seg["is_parental_haploid"]).astype(int)
    seg.loc[new_segs, "HZ"] = ((seg.loc[new_segs, "rescaled.cn.a1"] == 0) & (seg.loc[new_segs, "rescaled.cn.a2"] == 0)).astype(int)

    seg["bi.allelic"] = ((seg["modal.a1"] > 0) & (seg["modal.a2"] > 0)).astype(int)


    def get_ccf(allele, nu=10, max_cn=20):
        hscr = seg[f"hscr.{allele}"].to_numpy()  # haplotype-specific copy number
        seg_sigma = seg[f"sigma.{allele}"].to_numpy()
        scale = seg_sigma * np.sqrt((nu - 2) / nu)

        cn_grid = np.arange(max_cn)
        # allelic copy number comb
        comb = cn_grid[None, :] * delta.to_numpy()[:, None] * chr_ploidy.to_numpy()[:, None] + b.to_numpy()[:, None]

        # Get global modal CN, averaged across both alleles, to determine deletion vs amplification state
        _mu = seg["total_copy_ratio"].to_numpy()  # number of haploid copies
        _seg_sigma = seg[f"seg_sigma"].to_numpy()
        _scale = _seg_sigma * np.sqrt((nu - 2) / nu)

        comb_max = comb.max(axis=1)
        use_out = np.where(comb_max >= 1, 1, comb_max)  # assume haploid neutral
        mu_neutral = np.where(chr_ploidy > 0, (use_out - b) / delta / chr_ploidy, 0)
        Wq0 = st.norm.pdf(mu_neutral[:, None], loc=cn_grid[None, :], scale=1000) + 1e-10
        Wq0 /= Wq0.sum(axis=1, keepdims=True)
        log_prior = np.log(Wq0)
        ll = st.t.logpdf(comb, df=nu, loc=_mu[:, None], scale=_scale[:, None])
        log_mat = ll + log_prior
        log_mat -= sp.logsumexp(log_mat, axis=1, keepdims=True)
        seg_Q_post = np.exp(log_mat)
        col_sums = np.nansum(seg_Q_post, axis=0)
        modal_cn = 1 + int(np.nanargmax(col_sums[1:]))

        del_ix = hscr < comb[:, modal_cn]
        idx = np.arange(seg.shape[0])

        #TODO: Some amplification states 2->3 are not in concordance with how
        # ABSOLUTE assigns them as 1->2, leading to much lower CCF estimate for
        # those segments (~0 as opposed to ~1). It is not clear to me where the
        # discrepancy lies. PH

        qc = seg[f"modal.{allele}"].astype(int).to_numpy()
        qs = qc.copy()

        rows = np.where(del_ix)[0]
        for i in rows:
            # we want largest v with comb[v] < mu, v in [0, modal_cn - 1]
            j = np.searchsorted(comb[i, :], hscr[i], side="left") - 1
            v = np.clip(j, 0, modal_cn - 1)
            qc[i] = v + 1
            qs[i] = v

        rows = np.where(~del_ix)[0]
        for i in rows:
            # we want smallest v with mu < comb[v], v in [modal_cn + 1, max_cn - 1]
            j = np.searchsorted(comb[i, :], hscr[i], side="right")
            v = np.clip(j, modal_cn + 1, max_cn - 1)
            qc[i] = v - 1
            qs[i] = v

        ccf_grid = np.linspace(0, 1, 101)
        # For a subclonal CNA that moves the integer CN from qc (background)
        # to qs (altered) in a fraction f of cancer cells (the CCF), ABSOLUTE models
        # the expected copy-ratio as a linear function of f (ccf_grid):
        dd = comb[idx, qs] - comb[idx, qc]
        cr_grid = dd[:, None] * ccf_grid[None, :] + comb[idx, qc][:, None]
        cr_dens = st.t.logpdf(hscr[:, None], df=nu, loc=cr_grid, scale=scale[:, None])
        p = np.exp(cr_dens - sp.logsumexp(cr_dens, axis=1, keepdims=True))
        ecdf = np.cumsum(p, axis=1)

        # The CCF is the most likely change from background to altered state:
        hat = ccf_grid[np.argmax(p, axis=1)]  # mode
        # hat = (p * ccf_grid[None, :]).sum(axis=1)  # mean
        low = np.full(hat.shape, np.nan)
        high = np.full(hat.shape, np.nan)
        for i in range(hat.size):
            low[i] = np.interp(0.025, ecdf[i], ccf_grid, left=0)
            high[i] = np.interp(0.975, ecdf[i], ccf_grid, right=1)
        low = np.where(low < hat, low, hat)
        high = np.where(hat < high, high, hat)

        return hat, low, high

    def get_segment_ccf(
        nu=10,
        max_cn=200,
        pair_post_thresh=0.02,
        pair_window=2,
        high_cn_threshold=12,
        high_cn_sigma_threshold=1.5,
        high_cn_tail_mass_threshold=0.10,
        high_cn_pair_post_thresh=1e-3,
        baseline_total_cn_cap=8,
        ambiguous_width_threshold=0.50,
        ambiguous_max_p_threshold=0.05,
        ambiguous_pair_effn_threshold=4.0,
    ):
        """
        Joint segment-level CCF estimation using both alleles simultaneously.

        ---------------------------------------------
        1) Reports robust summaries of the CCF posterior:
           - mean, median, mode, equal-tail 95% interval
           The median is the recommended working point estimate for downstream use.

        2) Uses a chromosome-specific baseline *allele pair* when possible,
           rather than only a chromosome-specific modal total CN split.
           This is still approximate, but behaves better for gains / LOH backgrounds.

        3) Emits diagnostics for downstream use in CLaDE:
           - posterior width / entropy / max posterior mass
           - effective number of current-state pairs
           - high-copy flag
           - ambiguity flag

        Returns
        -------
        pd.DataFrame
            Indexed like `seg`, with columns:
                ccf_mean
                ccf_median
                ccf_mode
                ccf_ci95_low
                ccf_ci95_high
                ccf_ci95_width
                ccf_entropy
                ccf_max_p
                ccf_support_50
                ccf_support_90
                current_pair_entropy
                current_pair_effn
                tail_mass_high_cn
                high_copy_flag
                ambiguous_ccf_flag
                baseline_q1
                baseline_q2
        """
        eps = 1e-12

        hscr1 = seg["hscr.a1"].to_numpy(dtype=float)
        hscr2 = seg["hscr.a2"].to_numpy(dtype=float)

        sigma1 = seg["sigma.a1"].to_numpy(dtype=float)
        sigma2 = seg["sigma.a2"].to_numpy(dtype=float)

        total_cn = seg["CN"].to_numpy(dtype=float)
        total_cn_sigma = seg["CN.sigma"].to_numpy(dtype=float)

        cp = chr_ploidy.to_numpy(dtype=float)
        d = delta.to_numpy(dtype=float)
        bb = b.to_numpy(dtype=float)

        chrom = seg["Chromosome"].astype(str).to_numpy()
        is_parental_haploid = seg["is_parental_haploid"].to_numpy(dtype=bool)

        # current-state seeds
        q1_cur0 = pd.to_numeric(seg["modal.a1"], errors="coerce").to_numpy()
        q2_cur0 = pd.to_numeric(seg["modal.a2"], errors="coerce").to_numpy()

        q1_fallback = np.rint(pd.to_numeric(seg["rescaled.cn.a1"], errors="coerce").to_numpy())
        q2_fallback = np.rint(pd.to_numeric(seg["rescaled.cn.a2"], errors="coerce").to_numpy())

        q1_cur0 = np.where(np.isfinite(q1_cur0), q1_cur0, q1_fallback)
        q2_cur0 = np.where(np.isfinite(q2_cur0), q2_cur0, q2_fallback)

        q1_cur0 = np.where(np.isfinite(q1_cur0), q1_cur0, 0.0)
        q2_cur0 = np.where(np.isfinite(q2_cur0), q2_cur0, 0.0)

        q1_cur0 = np.clip(np.rint(q1_cur0), 0, max_cn - 1).astype(int)
        q2_cur0 = np.clip(np.rint(q2_cur0), 0, max_cn - 1).astype(int)

        cn_grid = np.arange(max_cn, dtype=float)
        ccf_grid = np.linspace(0.0, 1.0, 101)

        nseg = seg.shape[0]

        ccf_mean = np.full(nseg, np.nan, dtype=float)
        ccf_median = np.full(nseg, np.nan, dtype=float)
        ccf_mode = np.full(nseg, np.nan, dtype=float)
        ccf_low = np.full(nseg, np.nan, dtype=float)
        ccf_high = np.full(nseg, np.nan, dtype=float)
        ccf_width = np.full(nseg, np.nan, dtype=float)
        ccf_entropy = np.full(nseg, np.nan, dtype=float)
        ccf_max_p = np.full(nseg, np.nan, dtype=float)
        ccf_support_50 = np.full(nseg, np.nan, dtype=float)
        ccf_support_90 = np.full(nseg, np.nan, dtype=float)

        pair_entropy = np.full(nseg, np.nan, dtype=float)
        pair_effn = np.full(nseg, np.nan, dtype=float)

        tail_mass_high_cn = np.full(nseg, np.nan, dtype=float)
        high_copy_flag = np.zeros(nseg, dtype=int)
        ambiguous_ccf_flag = np.zeros(nseg, dtype=int)

        baseline_q1 = np.full(nseg, np.nan, dtype=float)
        baseline_q2 = np.full(nseg, np.nan, dtype=float)

        # ------------------------------------------------------------------
        # chromosome-specific baseline allele-pair estimate
        # ------------------------------------------------------------------
        def fallback_baseline_pair(ch):
            c0 = int(round(C0_by_chr.get(ch, 2)))
            chr_p = get_chromosomal_ploidy(ch)
            if chr_p == 1 and ((nX == 1) or (nY > 0)):
                return (0, max(c0, 1))
            if c0 <= 0:
                return (0, 0)
            if c0 % 2 == 0:
                return (c0 // 2, c0 // 2)
            q1 = c0 // 2
            q2 = c0 - q1
            return (q1, q2)

        q1_base_seed = np.where(np.isfinite(q1_cur0), q1_cur0, 0).astype(int)
        q2_base_seed = np.where(np.isfinite(q2_cur0), q2_cur0, 0).astype(int)

        # Use reasonably low/moderate CN states to estimate the broad/background allele pair.
        usable_for_baseline = (
                np.isfinite(q1_base_seed)
                & np.isfinite(q2_base_seed)
                & (q1_base_seed >= 0)
                & (q2_base_seed >= 0)
                & ((q1_base_seed + q2_base_seed) <= baseline_total_cn_cap)
        )

        chrom_baseline_pair = {}
        weights_all = seg["W"].to_numpy(dtype=float)

        for ch in pd.unique(chrom):
            idx = np.where((chrom == ch) & usable_for_baseline)[0]

            if idx.size == 0:
                chrom_baseline_pair[ch] = fallback_baseline_pair(ch)
                continue

            pair_weight = {}
            for j in idx:
                a = int(q1_base_seed[j])
                b_ = int(q2_base_seed[j])
                key = (a, b_)
                pair_weight[key] = pair_weight.get(key, 0.0) + float(weights_all[j])

            max_w = max(pair_weight.values())
            candidates = [k for k, v in pair_weight.items() if np.isclose(v, max_w)]

            c0 = int(round(C0_by_chr.get(ch, 2)))

            # Tie-breaker:
            #   1) total CN closest to chromosome baseline total CN
            #   2) more balanced split
            #   3) smaller total CN
            def tie_score(pair):
                a, b_ = pair
                return (abs((a + b_) - c0), abs(a - b_), a + b_)

            chrom_baseline_pair[ch] = min(candidates, key=tie_score)

        # ------------------------------------------------------------------
        # main loop
        # ------------------------------------------------------------------
        for i in range(nseg):
            if not (np.isfinite(hscr1[i]) and np.isfinite(hscr2[i])):
                continue
            if not (np.isfinite(cp[i]) and cp[i] > 0 and np.isfinite(d[i]) and d[i] > 0):
                continue

            sc1 = sigma1[i] * np.sqrt((nu - 2) / nu) if np.isfinite(sigma1[i]) and sigma1[i] > 0 else 1e-3
            sc2 = sigma2[i] * np.sqrt((nu - 2) / nu) if np.isfinite(sigma2[i]) and sigma2[i] > 0 else 1e-3

            comb_i = cn_grid * d[i] * cp[i] + bb[i]

            # total-CN posterior
            if np.isfinite(total_cn[i]) and np.isfinite(total_cn_sigma[i]) and total_cn_sigma[i] > 0:
                log_t = st.norm.logpdf(cn_grid, loc=total_cn[i], scale=max(total_cn_sigma[i], 1e-3))
            elif np.isfinite(total_cn[i]):
                log_t = st.norm.logpdf(cn_grid, loc=total_cn[i], scale=0.25)
            else:
                log_t = np.zeros_like(cn_grid)

            log_t -= sp.logsumexp(log_t)
            t_post = np.exp(log_t)

            # baseline allele pair
            if is_parental_haploid[i]:
                qc10, qc20 = fallback_baseline_pair(chrom[i])
            else:
                qc10, qc20 = chrom_baseline_pair.get(chrom[i], fallback_baseline_pair(chrom[i]))

            qc10 = int(np.clip(qc10, 0, max_cn - 1))
            qc20 = int(np.clip(qc20, 0, max_cn - 1))

            baseline_q1[i] = qc10
            baseline_q2[i] = qc20

            # posterior over current allele states from observed HSCR and total-CN feasibility
            feasible = np.flip(np.cumsum(np.flip(t_post)))
            feasible = np.clip(feasible, eps, 1.0)

            log_q1 = st.t.logpdf(hscr1[i], df=nu, loc=comb_i, scale=sc1) + np.log(feasible)
            log_q2 = st.t.logpdf(hscr2[i], df=nu, loc=comb_i, scale=sc2) + np.log(feasible)

            log_q1 -= sp.logsumexp(log_q1)
            log_q2 -= sp.logsumexp(log_q2)

            q1_post = np.exp(log_q1)
            q2_post = np.exp(log_q2)

            # high-copy / ecDNA-like regime
            tail_mass = t_post[int(high_cn_threshold):].sum() if high_cn_threshold < len(t_post) else 0.0
            tail_mass_high_cn[i] = tail_mass

            is_high_copy = (
                    (np.isfinite(total_cn[i]) and total_cn[i] >= high_cn_threshold)
                    or (
                            np.isfinite(total_cn_sigma[i])
                            and total_cn_sigma[i] >= high_cn_sigma_threshold
                            and tail_mass >= high_cn_tail_mass_threshold
                    )
                    or (tail_mass >= 0.5)
            )
            high_copy_flag[i] = int(is_high_copy)

            if not is_high_copy:
                # ordinary regime: stay close to anchored current calls
                cand1 = {q1_cur0[i], int(np.argmax(q1_post))}
                cand2 = {q2_cur0[i], int(np.argmax(q2_post))}

                for q in range(max(0, q1_cur0[i] - pair_window), min(max_cn, q1_cur0[i] + pair_window + 1)):
                    if q1_post[q] >= pair_post_thresh:
                        cand1.add(int(q))
                for q in range(max(0, q2_cur0[i] - pair_window), min(max_cn, q2_cur0[i] + pair_window + 1)):
                    if q2_post[q] >= pair_post_thresh:
                        cand2.add(int(q))

                cand1 = sorted(cand1)
                cand2 = sorted(cand2)
            else:
                # high-copy regime: allow a broad set of current states
                cand1 = np.where(q1_post >= high_cn_pair_post_thresh)[0].tolist()
                cand2 = np.where(q2_post >= high_cn_pair_post_thresh)[0].tolist()

                if not cand1:
                    cand1 = [int(np.argmax(q1_post))]
                if not cand2:
                    cand2 = [int(np.argmax(q2_post))]

                # retain anchored state if itself high-CN
                if q1_cur0[i] >= high_cn_threshold:
                    cand1.append(int(q1_cur0[i]))
                if q2_cur0[i] >= high_cn_threshold:
                    cand2.append(int(q2_cur0[i]))

                cand1 = sorted(set(cand1))
                cand2 = sorted(set(cand2))

            # current-state pair posterior
            pair_states = []
            pair_weights = []

            for qs1 in cand1:
                for qs2 in cand2:
                    tcur = qs1 + qs2
                    if tcur >= max_cn:
                        continue

                    w = q1_post[qs1] * q2_post[qs2] * t_post[tcur]
                    if w <= 0:
                        continue

                    pair_states.append((qs1, qs2))
                    pair_weights.append(w)

            if not pair_states:
                pair_states = [(q1_cur0[i], q2_cur0[i])]
                pair_weights = [1.0]

            pair_weights = np.asarray(pair_weights, dtype=float)
            pair_weights /= pair_weights.sum()

            pe = -np.sum(pair_weights * np.log(np.clip(pair_weights, eps, 1.0)))
            pair_entropy[i] = pe
            pair_effn[i] = np.exp(pe)

            # mixture posterior over segment-level CCF
            mix = np.zeros_like(ccf_grid, dtype=float)

            for (qs1, qs2), w in zip(pair_states, pair_weights):
                cr1_grid = ((qs1 - qc10) * ccf_grid + qc10) * d[i] * cp[i] + bb[i]
                cr2_grid = ((qs2 - qc20) * ccf_grid + qc20) * d[i] * cp[i] + bb[i]

                logp1 = st.t.logpdf(hscr1[i], df=nu, loc=cr1_grid, scale=sc1)
                logp2 = st.t.logpdf(hscr2[i], df=nu, loc=cr2_grid, scale=sc2)

                logp = logp1 + logp2
                logp -= sp.logsumexp(logp)
                p = np.exp(logp)

                mix += w * p

            if not np.isfinite(mix).all() or mix.sum() <= 0:
                qs1, qs2 = q1_cur0[i], q2_cur0[i]
                cr1_grid = ((qs1 - qc10) * ccf_grid + qc10) * d[i] * cp[i] + bb[i]
                cr2_grid = ((qs2 - qc20) * ccf_grid + qc20) * d[i] * cp[i] + bb[i]
                logp1 = st.t.logpdf(hscr1[i], df=nu, loc=cr1_grid, scale=sc1)
                logp2 = st.t.logpdf(hscr2[i], df=nu, loc=cr2_grid, scale=sc2)
                logp = logp1 + logp2
                logp -= sp.logsumexp(logp)
                mix = np.exp(logp)

            mix = np.clip(mix, 0, np.inf)
            mix /= mix.sum()

            ecdf = np.cumsum(mix)

            # robust summaries
            ccf_mean[i] = np.sum(mix * ccf_grid)
            ccf_median[i] = np.interp(0.5, ecdf, ccf_grid, left=0.0, right=1.0)
            ccf_mode[i] = ccf_grid[np.argmax(mix)]

            ccf_low[i] = np.interp(0.025, ecdf, ccf_grid, left=0.0, right=1.0)
            ccf_high[i] = np.interp(0.975, ecdf, ccf_grid, left=0.0, right=1.0)
            ccf_width[i] = ccf_high[i] - ccf_low[i]

            ce = -np.sum(mix * np.log(np.clip(mix, eps, 1.0)))
            ccf_entropy[i] = ce
            ccf_max_p[i] = mix.max()
            ccf_support_50[i] = np.sum(mix >= 0.5 * mix.max())
            ccf_support_90[i] = np.sum(mix >= 0.1 * mix.max())

            ambiguous_ccf_flag[i] = int(
                is_high_copy
                or (ccf_width[i] >= ambiguous_width_threshold)
                or (ccf_max_p[i] <= ambiguous_max_p_threshold)
                or (pair_effn[i] >= ambiguous_pair_effn_threshold)
            )

        precision = 3
        return pd.DataFrame(
            {
                "ccf_mean": ccf_mean.round(precision),
                "ccf_median": ccf_median.round(precision),
                "ccf_mode": ccf_mode.round(precision),
                "ccf_ci95_low": ccf_low.round(precision),
                "ccf_ci95_high": ccf_high.round(precision),
                "ccf_ci95_width": ccf_width.round(precision),
                "ccf_entropy": ccf_entropy.round(precision),
                "ccf_max_p": ccf_max_p.round(precision),
                "ccf_support_50": ccf_support_50,
                "ccf_support_90": ccf_support_90,
                "current_pair_entropy": pair_entropy.round(precision),
                "current_pair_effn": pair_effn.round(precision),
                "tail_mass_high_cn": tail_mass_high_cn.round(precision),
                "high_copy_flag": high_copy_flag,
                "ambiguous_ccf_flag": ambiguous_ccf_flag,
                "baseline_q1": baseline_q1,
                "baseline_q2": baseline_q2,
            },
            index=seg.index,
        )

    for allele in ["a1", "a2"]:
        ccf_hat, ccf_low, ccf_hi = get_ccf(allele)
        ccf_hat = pd.Series(ccf_hat, index=seg.index).round(2)
        ccf_low = pd.Series(ccf_low, index=seg.index).round(2)
        ccf_hi = pd.Series(ccf_hi, index=seg.index).round(2)
        if allele == "a1":
            ccf_hat.loc[seg["is_parental_haploid"]] = 1
            ccf_low.loc[seg["is_parental_haploid"]] = 1
            ccf_hi.loc[seg["is_parental_haploid"]] = 1
        is_subclonal = ((seg[f"rescaled.cn.{allele}"] != 1) & (0 < ccf_hat) & (ccf_hat < 1)).astype(int)

        _df = pd.concat([ccf_hat.to_frame("ccf"), seg[f"cancer.cell.frac.{allele}"]], axis=1)
        ccf_corr = _df.corr().loc["ccf", f"cancer.cell.frac.{allele}"]
        diff = np.abs(_df["ccf"] - _df[f"cancer.cell.frac.{allele}"])
        message(f"Correlation between inferred CCF for {allele} and ABSOLUTE output: {ccf_corr:.3f}, "
                f"median |diff| = {np.nanmedian(diff):.3f} ± {np.nanstd(diff):.3f}")

        for col, series in zip([f"cancer.cell.frac.{allele}", f"ccf.ci95.low.{allele}", f"ccf.ci95.high.{allele}", f"subclonal.{allele}"], [ccf_hat, ccf_low, ccf_hi, is_subclonal]):
            seg.loc[seg[col].isna(), col] = series.loc[seg[col].isna()]

    seg_ccf = get_segment_ccf()

    seg = pd.concat([seg, seg_ccf], axis=1)

    schz_na = seg["SC_HZ"].isna()
    seg.loc[schz_na, "SC_HZ"] = (
        seg.loc[schz_na, "HZ"].astype(bool)
        & ((seg.loc[schz_na, "cancer.cell.frac.a1"] < 1)
         | (seg.loc[schz_na, "cancer.cell.frac.a2"] < 1))
    ).astype(int)

    message(f"Number of rescued segments: {new_segs.sum()}")

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

    seg = seg.set_index(["Chromosome", "Start.bp", "End.bp"])
    seg = seg.reindex(sort_genomic_positions(index=seg.index))
    seg = seg.reset_index()
    seg["Segment_Mean"] = np.log2(seg["rescaled_total_cn"].clip(lower=1e-2)) - np.log2(np.where(chr_ploidy > 0, args.ploidy * chr_ploidy / args.normal_ploidy, 1))

    good_rows = (seg["n_hets"] >= args.min_hets) | seg["n_hets"].isna()
    good_rows &= (seg["n_probes"] >= args.min_probes) | seg["n_probes"].isna()
    n = seg.shape[0] - np.sum(good_rows)
    pct_genomic_drop = seg.loc[~good_rows, "W"].sum() * 100
    print(f"Dropping {n}/{seg.shape[0]} (-{pct_genomic_drop:.6f}% of genome) segments with min_hets < {args.min_hets} or min_probes < {args.min_probes}.")

    seg = seg.loc[good_rows]
    seg["W"] = seg["length"] / seg["length"].sum()

    os.makedirs(args.outdir, exist_ok=True)

    seg.loc[seg["is_new"]].reset_index()[["Chromosome", "Start.bp", "End.bp"]].to_csv(f"{args.outdir}/{args.sample}.rescued_intervals.txt", sep="\t", index=False)
    abs_seg_cols = [c for c in list(abs_dtypes.keys()) + extra_columns if c in seg.columns]
    seg[abs_seg_cols].to_csv(f"{args.outdir}/{args.sample}.segtab.completed.txt", sep="\t", index=False)
    seg[["sample", "Chromosome", "Start.bp", "End.bp", "Segment_Mean", "rescaled_total_cn"]].to_csv(f"{args.outdir}/{args.sample}.IGV.seg.completed.txt", sep="\t", index=False)


if __name__ == "__main__":
    main()
