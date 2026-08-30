import argparse
import csv
import gzip
import os
import re
import time
import warnings

import numpy as np
import pandas as pd
import scipy.special as sp
import scipy.stats as st

warnings.filterwarnings("ignore")


ABSOLUTE_EXTRA_COLUMNS = """
observed_alt	normal_allele_count	A1.ix	A2.ix	T.seg.ix	Number of times codon change is in COSMIC	q_hat	HS_q_hat_1	HS_q_hat_2	total_q_hat	SC.ix	C.ix	total_qc	total_qs	SC_Aq_d	SC_Aq_a	C_Aq	Pr_SC_SC.ix	Pr_SC_C.ix	clonal_scna_mut_ix	Pr_somatic_clonal	Pr_germline	Pr_subclonal	Pr_subclonal_wt0	Pr_wt0	Pr_ge2	Pr_GL_som_HZ_alt	Pr_GL_som_HZ_ref	Pr_cryptic_SCNA	modal_q_s	LL	subclonal.ix	subclonal_wt0.ix	clonal.ix	wt0.ix	clonal_het.ix	ge2.ix	homozygous.ix	ccf_hat	ccf_CI95_low	ccf_CI95_high	H1	H2	H3	H4	detection_power	detection_power_for_single_read	purity	SSNV_skew	0	0.01	0.02	0.03	0.04	0.05	0.06	0.07	0.08	0.09	0.1	0.11	0.12	0.13	0.14	0.15	0.16	0.17	0.18	0.19	0.2	0.21	0.22	0.23	0.24	0.25	0.26	0.27	0.28	0.29	0.3	0.31	0.32	0.33	0.34	0.35	0.36	0.37	0.38	0.39	0.4	0.41	0.42	0.43	0.44	0.45	0.46	0.47	0.48	0.49	0.5	0.51	0.52	0.53	0.54	0.55	0.56	0.57	0.58	0.59	0.6	0.61	0.62	0.63	0.64	0.65	0.66	0.67	0.68	0.69	0.7	0.71	0.72	0.73	0.74	0.75	0.76	0.77	0.78	0.79	0.8	0.81	0.82	0.83	0.84	0.85	0.86	0.87	0.88	0.89	0.9	0.91	0.92	0.93	0.94	0.95	0.96	0.97	0.98	0.99	1	local_cn_a1	local_cn_a2
""".strip().split("\t")

CCF_GRID = np.round(np.linspace(0, 1, 101), 2)
CCF_COLUMNS = [f"{x:.2f}".rstrip("0").rstrip(".") for x in CCF_GRID]
MAF_BASE_COLUMNS = [
    "Tumor_Sample_Barcode",
    "Chromosome",
    "Start_position",
    "End_position",
    "Reference_Allele",
    "Tumor_Seq_Allele2",
]


def message(*args, **kwargs) -> None:
    print(f"{time.strftime('%H:%M:%S')} ", *args, **kwargs)


def parse_args():
    parser = argparse.ArgumentParser(
        prog="CalculateCancerCellFraction",
        description="""
            Assign each somatic SNV / INDEL a cancer-cell-fraction (CCF) posterior on the
            0.01 grid, given the sample's absolute copy-number segtab (from
            map_to_absolute_copy_number) and --purity/--ploidy. Each variant is mapped onto
            its segment to obtain the local copy number, then its observed VAF is converted
            to a CCF posterior; outputs include the point estimate ccf_hat, the 95% interval
            ccf_CI95_low/high, the full per-grid posterior, and detection-power columns.
            Multiplicity (alt copies per carrying cell) is folded into the VAF->CCF mapping,
            so a clonal multi-copy mutation still yields ccf ~ 1 (never > 1).

            Variants are EXCLUDED if they cannot be placed on a segment, fall on MT, or have
            no tumor coverage. Variants already present in --absolute_maf keep their ABSOLUTE
            CCF; caller-only variants (--snv_maf / --indel_maf) are "rescued" (inferred here).
            Unlike map_to_absolute's per-segment CCF, variant CCF is identifiable because VAF
            directly encodes purity * ccf * multiplicity / local_copy_number.
        """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.usage = (
        "calculate_cancer_cell_fraction.py --outdir <outdir> --absolute_segtab <segtab> "
        "--purity <purity> --ploidy <ploidy> --ssnv_skew <ssnv_skew> "
        "[--absolute_maf <absolute_maf>] [--snv_maf <snv_maf>] [--indel_maf <indel_maf>] "
        "[--sample <sample>] [--sex <sex>] [--normal_ploidy <normal_ploidy>]"
    )
    parser.add_argument("--outdir", type=str, required=True, help="Output directory.")
    parser.add_argument("--purity", type=float, required=True, help="Tumor purity as inferred by ABSOLUTE.")
    parser.add_argument("--ploidy", type=float, required=True, help="Tumor ploidy as inferred by ABSOLUTE.")
    parser.add_argument("--ssnv_skew", type=float, required=False, default=0.9883274, help="SSNV skew parameter used by ABSOLUTE.")
    parser.add_argument("--sample", type=str, required=False, help="Sample name.")
    parser.add_argument("--sex", type=str, default="XXY", help="Patient's sex genotype.")
    parser.add_argument("--absolute_maf", type=str, required=False, help="Path to ABSOLUTE ABS_MAF output.")
    parser.add_argument("--absolute_segtab", type=str, required=True, help="Path to completed segtab from map_to_absolute_copy_number.py.")
    parser.add_argument("--snv_maf", type=str, required=False, help="Path to caller SNV MAF.")
    parser.add_argument("--indel_maf", type=str, required=False, help="Path to caller INDEL MAF.")
    parser.add_argument("--normal_ploidy", type=int, default=2, help="Normal/germline ploidy of that organism.")
    parser.add_argument("--copy_num_type", type=str, default="allelic", help="Type of copy number mode {allelic, total}.")
    return parser.parse_args()


def print_args(args):
    message("Calling CalculateCancerCellFraction")
    print("Arguments:")
    for key, value in vars(args).items():
        print(f"  {key}: {value}")
    print()


def standardize_sex(sex: str) -> str:
    s = str(sex).upper()
    if s == "FEMALE":
        return "XX"
    if s == "MALE":
        return "XY"
    if s in ["UNKNOWN", "NAN", "NONE"]:
        return "XXY"
    return s


def normalize_chromosome(chrom) -> str:
    c = str(chrom).strip()
    if c.lower().startswith("chr"):
        c = c[3:]
    return "MT" if c == "M" else c


def get_chromosomal_ploidy(chrom, sex: str, normal_ploidy: int = 2) -> int:
    s = standardize_sex(sex)
    nX = s.count("X")
    nY = s.count("Y")
    chrom = normalize_chromosome(chrom)
    if chrom == "X":
        return nX
    if chrom == "Y":
        return nY
    return normal_ploidy


def sort_genomic_df(df: pd.DataFrame, chr_col="Chromosome", start_col="Start_position", end_col="End_position") -> pd.DataFrame:
    missing = [col for col in [chr_col, start_col, end_col] if col not in df.columns]
    if missing:
        if df.empty:
            return df.copy()
        raise KeyError(missing[0])
    contig_order = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
    observed = set(df[chr_col].astype(str))
    contig_order += list(sorted(observed - set(contig_order)))
    out = df.copy()
    out["_chrom_sort"] = pd.Categorical(out[chr_col].astype(str), categories=contig_order, ordered=True)
    out = out.sort_values(by=["_chrom_sort", start_col, end_col], kind="stable").drop(columns="_chrom_sort")
    return out


def normalize_maf_columns(df: pd.DataFrame) -> pd.DataFrame:
    return df.rename(columns={c: c.replace("(", ".").replace(")", ".") for c in df.columns})


def count_leading_comment_lines(path: str, comment: str = "#") -> int:
    """Number of leading comment lines, counted at line start only.

    ``pd.read_csv(comment=...)`` truncates a line at the first occurrence of the
    comment character *anywhere* in that line, not just at the start. MAF files
    use ``#`` as a leading header-block marker but also carry it inside free-text
    annotation values (gene descriptions, ClinVar terms), so passing ``comment``
    to the parser silently amputates those rows. Count the leading block here and
    pass ``skiprows`` instead.
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    n = 0
    with opener(path, "rt", newline="") as fh:
        for line in fh:
            if not line.startswith(comment):
                break
            n += 1
    return n


#: A field as ``csv.QUOTE_MINIMAL`` would have written it: wrapped in quotes, with
#: every interior quote doubled. Deliberately strict -- see :func:`unquote_fields`.
QUOTE_MINIMAL_FIELD = re.compile(r'"(?:[^"]|"")*"')


def unquote_fields(df: pd.DataFrame) -> pd.DataFrame:
    """Undo CSV-style field quoting *after* the row has already been split on tabs.

    Reading with ``QUOTE_NONE`` keeps every tab as a delimiter, which is what stops
    records merging, but it also leaves in place the quote characters a writer added
    around a value that contains a literal ``"`` (``\"\"\"X\"`` for the value ``"X``).
    Undoing that here is safe in a way that doing it in the parser is not: the field
    boundaries are already fixed, so no amount of quote imbalance can swallow a tab
    or a newline.

    The match has to be *exact*, not merely "starts and ends with a quote". MAF has
    no quoting convention, so a quote is ordinarily literal data, and the values that
    matter most here look deceptively like quoted fields: Funcotator writes HGNC
    alias/previous names as quoted comma-separated lists,
    ``"plexin 2", "plexin-A2"``. That starts and ends with a quote but is not a
    quoted field, and stripping its outer characters yields the lopsided
    ``plexin 2", "plexin-A2``. So a value is only unquoted when it is a complete
    ``QUOTE_MINIMAL`` encoding -- wrapped, with every interior quote doubled -- which
    the list form is not. Anything else is passed through untouched.
    """
    out = df.copy()
    for col in out.columns:
        if out[col].dtype != object:
            continue
        values = out[col]
        is_str = values.map(lambda v: isinstance(v, str))
        if not bool(is_str.any()):
            continue
        target = values.where(~is_str, values.astype(str))
        encoded = is_str & target.str.fullmatch(QUOTE_MINIMAL_FIELD).fillna(False)
        if bool(encoded.any()):
            target = target.where(
                ~encoded,
                target.str.slice(1, -1).str.replace('""', '"', regex=False),
            )
        out[col] = values.where(~is_str, target)
    return out


def read_table(path: str, comment="#") -> pd.DataFrame:
    """Read a MAF/segtab TSV without letting quote or comment handling corrupt it.

    Two parser settings matter here and both default to the wrong thing for MAF:

    * ``quoting=csv.QUOTE_NONE`` -- MAF is a plain tab-delimited format with no
      quoting convention, so ``"`` is a literal data character. Funcotator writes
      HGNC alias/previous-name values that contain quoted, comma-separated lists
      (``... (yeast)", "origin recognition complex, subunit 6"``). Read with
      pandas' default ``quotechar='"'``, an odd number of quotes on a line opens a
      quoted region that swallows every following tab *and newline* until the next
      quote. That merges whole records into one giant field: the surviving row
      loses field boundaries (its values shift left, tail padded with NaN) and the
      absorbed rows vanish. Observed in production output as ABS_MAF rows whose
      CCF-grid block lands 72 columns early and whose ``individual_id`` holds a
      posterior density value instead of a sample name.
    * ``skiprows`` rather than ``comment`` -- see
      :func:`count_leading_comment_lines`.

    Genuinely absent or empty inputs still degrade to an empty frame because
    several inputs are optional, but a *parse* failure is re-raised: silently
    substituting an empty frame for a malformed file is how this class of bug
    stays invisible.
    """
    try:
        skiprows = count_leading_comment_lines(path, comment=comment) if comment else 0
        return unquote_fields(
            pd.read_csv(
                path,
                sep="\t",
                skiprows=skiprows,
                quoting=csv.QUOTE_NONE,
                low_memory=False,
            )
        )
    except (FileNotFoundError, pd.errors.EmptyDataError) as e:
        message("Caught exception:", e)
        message("Using empty dataframe instead.")
        return pd.DataFrame()


def sanitize_delimiters(df: pd.DataFrame) -> pd.DataFrame:
    """Replace embedded tabs/newlines in string fields with single spaces.

    Guarantees the written table is a rectangular TSV that any tab-splitter can
    read, so no downstream consumer needs quote-aware parsing. With the read path
    fixed this should be a no-op; if it is not, the offending columns are named
    so the upstream annotation source can be corrected.
    """
    out = df.copy()
    offenders: dict[str, int] = {}
    for col in out.columns:
        if out[col].dtype != object:
            continue
        as_str = out[col].astype(str)
        hit = as_str.str.contains(r"[\t\r\n]", regex=True, na=False)
        n = int(hit.sum())
        if not n:
            continue
        offenders[str(col)] = n
        out.loc[hit, col] = as_str.loc[hit].str.replace(r"[\t\r\n]+", " ", regex=True)
    if offenders:
        message(
            f"WARNING: replaced embedded tab/newline characters in {sum(offenders.values())} "
            f"field(s) before writing: {offenders}"
        )
    return out


def write_table(df: pd.DataFrame, path: str) -> None:
    """Write a plain, quote-free, rectangular TSV."""
    sanitize_delimiters(df).to_csv(
        path, sep="\t", index=False, quoting=csv.QUOTE_NONE, escapechar=None
    )


def read_optional_maf(path: str | None) -> pd.DataFrame:
    if path is None:
        return pd.DataFrame()
    return normalize_maf_columns(read_table(path, comment="#"))


def to_numeric(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    out = df.copy()
    for col in columns:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")
    return out


def choose_output_columns(caller_df: pd.DataFrame, abs_maf: pd.DataFrame) -> list[str]:
    if not caller_df.empty:
        caller_cols = [c for c in caller_df.columns if not c.startswith("_")]
    elif len(abs_maf.columns):
        caller_cols = [c for c in abs_maf.columns if not c.startswith("_")]
    else:
        caller_cols = MAF_BASE_COLUMNS.copy()
    extra_cols = [c for c in ABSOLUTE_EXTRA_COLUMNS if c not in caller_cols]
    return caller_cols + extra_cols


def ensure_columns(df: pd.DataFrame, columns: list[str], fill_value=np.nan) -> pd.DataFrame:
    out = df.copy()
    for col in columns:
        if col not in out.columns:
            out[col] = fill_value
    return out


def infer_sample_name(args, seg: pd.DataFrame, caller_df: pd.DataFrame, abs_maf: pd.DataFrame) -> str:
    if args.sample is not None:
        return args.sample
    for df in [caller_df, abs_maf]:
        if not df.empty and "Tumor_Sample_Barcode" in df.columns:
            vals = df["Tumor_Sample_Barcode"].dropna().astype(str).unique()
            if len(vals):
                return vals[0]
    vals = seg["sample"].dropna().astype(str).unique()
    return vals[0] if len(vals) else "sample"


def make_variant_key(df: pd.DataFrame) -> pd.Series:
    return (
        df["_chrom"].astype(str)
        + ":" + df["Start_position"].astype(int).astype(str)
        + ":" + df["End_position"].astype(int).astype(str)
        + ":" + df["Reference_Allele"].astype(str)
        + ":" + df["Tumor_Seq_Allele2"].astype(str)
    )


def prepare_caller_maf(df: pd.DataFrame, label: str) -> pd.DataFrame:
    if df.empty:
        return df
    out = df.copy()
    out = to_numeric(out, ["Start_position", "End_position", "t_alt_count", "t_ref_count", "n_alt_count", "n_ref_count", "tumor_f"])
    out["_source"] = label
    out["_chrom"] = out["Chromosome"].map(normalize_chromosome)
    out["_variant_key"] = make_variant_key(out)
    return out


def prepare_abs_maf(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    out = df.copy()
    out = to_numeric(
        out,
        [
            "Start_position", "End_position", "t_alt_count", "t_ref_count", "n_alt_count", "n_ref_count",
            "ccf_hat", "ccf_CI95_low", "ccf_CI95_high", "local_cn_a1", "local_cn_a2",
            "detection_power", "detection_power_for_single_read",
        ],
    )
    out["_chrom"] = out["Chromosome"].map(normalize_chromosome)
    out["_variant_key"] = make_variant_key(out)
    return out


def prepare_segtab(df: pd.DataFrame, sex: str, normal_ploidy: int) -> pd.DataFrame:
    out = df.copy()
    # A failed ABSOLUTE extraction may leave a zero-byte segtab, especially when
    # an older cached Postprocess task invokes this script after copy-number mapping
    # has gracefully skipped. Treat it as a valid empty segmentation instead of
    # indexing absent coordinate columns.
    out = ensure_columns(out, ["sample", "Chromosome", "Start.bp", "End.bp"])
    # A total-copy-ratio segtab names the total-CN columns without the "total" qualifier
    # and carries no allelic calls; normalize the names so total-CN lookups (modal_total_cn)
    # resolve. Guard against clobbering a canonical column an allelic segtab already has.
    tcr_to_canonical = {"modal_cn": "modal_total_cn", "expected_cn": "expected_total_cn"}
    out = out.rename(columns={s: d for s, d in tcr_to_canonical.items()
                              if s in out.columns and d not in out.columns})
    out = to_numeric(
        out,
        ["Start.bp", "End.bp", "rescaled.cn.a1", "rescaled.cn.a2", "rescaled_total_cn", "modal.a1", "modal.a2", "modal_total_cn"],
    )
    out["_chrom"] = out["Chromosome"].map(normalize_chromosome)
    out["_chr_ploidy"] = out["_chrom"].map(lambda c: get_chromosomal_ploidy(c, sex=sex, normal_ploidy=normal_ploidy))
    out = sort_genomic_df(out, chr_col="_chrom", start_col="Start.bp", end_col="End.bp").reset_index(drop=True)
    out["_seg_idx1"] = np.arange(1, out.shape[0] + 1, dtype=int)
    return out


def concat_callers(snv_maf: pd.DataFrame, indel_maf: pd.DataFrame) -> pd.DataFrame:
    dfs = [df for df in [snv_maf, indel_maf] if not df.empty]
    if not dfs:
        return pd.DataFrame()
    out = pd.concat(dfs, axis=0, ignore_index=True, sort=False)
    before = out.shape[0]
    out = out.drop_duplicates(subset=["_variant_key"], keep="first").reset_index(drop=True)
    dropped = before - out.shape[0]
    if dropped:
        message(f"Dropped {dropped} duplicate caller rows by variant key.")
    return out


def map_variants_to_segments(variants: pd.DataFrame, seg: pd.DataFrame) -> pd.DataFrame:
    out = variants.copy()
    out["T.seg.ix"] = np.nan
    out["A1.ix"] = np.nan
    out["A2.ix"] = np.nan
    out["local_cn_a1"] = np.nan
    out["local_cn_a2"] = np.nan
    out["q_hat"] = np.nan
    out["HS_q_hat_1"] = np.nan
    out["HS_q_hat_2"] = np.nan
    out["total_q_hat"] = np.nan
    out["normal_allele_count"] = np.nan

    for chrom, v_idx in out.groupby("_chrom").groups.items():
        s = seg.loc[seg["_chrom"] == chrom]
        if s.empty:
            continue

        starts = s["Start.bp"].to_numpy()
        ends = s["End.bp"].to_numpy()
        seg_rows = s.index.to_numpy()

        vv = out.loc[list(v_idx)]
        pos = vv["Start_position"].to_numpy()
        pos_end = vv["End_position"].to_numpy()

        idx = np.searchsorted(starts, pos, side="right") - 1
        valid = (idx >= 0) & (pos_end <= ends[np.clip(idx, 0, len(ends) - 1)])
        variant_rows = vv.index.to_numpy()[valid]
        matched_seg_rows = seg_rows[idx[valid]]

        # A total-copy-ratio segtab has no allelic call columns; fall back to NaN for
        # those so total mode degrades to total-CN CCF instead of raising. The CCF
        # posterior only needs the total copy number (q_hat) and chromosomal ploidy.
        def seg_value(name):
            if name in seg.columns:
                return seg.loc[matched_seg_rows, name].to_numpy()
            return np.full(matched_seg_rows.shape[0], np.nan)

        out.loc[variant_rows, "T.seg.ix"] = seg.loc[matched_seg_rows, "_seg_idx1"].to_numpy()
        out.loc[variant_rows, "A1.ix"] = seg.loc[matched_seg_rows, "_seg_idx1"].to_numpy()
        out.loc[variant_rows, "A2.ix"] = seg.loc[matched_seg_rows, "_seg_idx1"].to_numpy()
        out.loc[variant_rows, "local_cn_a1"] = seg_value("rescaled.cn.a1")
        out.loc[variant_rows, "local_cn_a2"] = seg_value("rescaled.cn.a2")
        out.loc[variant_rows, "HS_q_hat_1"] = seg_value("modal.a1")
        out.loc[variant_rows, "HS_q_hat_2"] = seg_value("modal.a2")
        out.loc[variant_rows, "total_q_hat"] = seg_value("modal_total_cn")
        out.loc[variant_rows, "q_hat"] = seg_value("modal_total_cn")
        out.loc[variant_rows, "normal_allele_count"] = seg.loc[matched_seg_rows, "_chr_ploidy"].to_numpy()

    return out


def d_beta_binom(k, a, b, n, log=False):
    log_d = sp.betaln(k + a, n - k + b) + sp.gammaln(n + 1) - sp.gammaln(k + 1) - sp.gammaln(n - k + 1) - sp.betaln(a, b)
    return log_d if log else np.exp(log_d)


def log_density_alt_cond_coverage(alt, cov, f, ssnv_skew, rho, epsilon):
    f = np.asarray(f, dtype=float)
    skewed_f = np.clip(ssnv_skew * f, 0.0, 1.0)
    tiny = np.finfo(float).tiny
    A = np.clip(skewed_f * rho, tiny, None)
    B = np.clip((1 - skewed_f) * rho, tiny, None)

    l_dens = d_beta_binom(alt, A, B, cov, log=True)
    l_dens = np.asarray(l_dens, dtype=float)
    l_dens[f == 0.0] = st.binom.logpmf(alt, cov, epsilon / 3)
    l_dens[f == 1.0] = st.binom.logpmf(alt, cov, 1 - epsilon / 3)
    l_dens[np.isnan(l_dens)] = -np.inf
    return l_dens


def calc_ccf_posterior_ll_grid(alt, ref, alpha, q, normal_allele_count, ssnv_skew, rho, epsilon):
    f = (alpha * CCF_GRID) / (normal_allele_count * (1 - alpha) + alpha * q)
    cov = alt + ref
    return log_density_alt_cond_coverage(alt, cov, f, ssnv_skew=ssnv_skew, rho=rho, epsilon=epsilon)


def calc_ccf_dens(mut_df: pd.DataFrame, alpha: float, ssnv_skew: float, rho: float, epsilon: float) -> np.ndarray:
    ccf_ll = np.full((mut_df.shape[0], CCF_GRID.size), np.nan)
    for i, row in enumerate(mut_df.itertuples(index=False)):
        if pd.isna(row.q_hat):
            continue
        # An observed somatic variant implies at least one tumor copy at the locus, so a
        # modal total copy number of 0 (e.g. local_cn_a1 = 0 and local_cn_a2 < 0.5, which
        # rounds q_hat down to 0) is a copy-number rounding artifact, not a homozygous
        # deletion of the mutated allele. Floor the copy number at 1 so the variant still
        # gets a CCF posterior instead of being dropped with an empty grid.
        q = max(float(row.q_hat), 1.0)
        ccf_ll[i, :] = calc_ccf_posterior_ll_grid(
            alt=int(row.t_alt_count),
            ref=int(row.t_ref_count),
            alpha=alpha,
            q=q,
            normal_allele_count=float(row.normal_allele_count),
            ssnv_skew=ssnv_skew,
            rho=rho,
            epsilon=epsilon,
        )

    ccf_dens = np.full_like(ccf_ll, np.nan, dtype=float)
    for i in range(ccf_ll.shape[0]):
        row_ll = ccf_ll[i, :]
        if np.all(np.isnan(row_ll)):
            continue
        finite = np.isfinite(row_ll)
        if not np.any(finite):
            continue
        log_z = sp.logsumexp(row_ll[finite])
        ccf_dens[i, finite] = np.exp(row_ll[finite] - log_z)
    return ccf_dens


def calc_ccf_95ci(ccf_dens: np.ndarray) -> pd.DataFrame:
    ccf_ci95 = np.full((ccf_dens.shape[0], 2), np.nan)
    ccf_hat = np.full(ccf_dens.shape[0], np.nan)

    for i in range(ccf_dens.shape[0]):
        if np.any(np.isnan(ccf_dens[i, :])):
            continue
        max_ix = np.argmax(ccf_dens[i, :])
        ccf_hat[i] = CCF_GRID[max_ix]
        ecdf = np.cumsum(ccf_dens[i, :])
        ccf_ci95[i, :] = np.interp([0.025, 0.975], ecdf, CCF_GRID)

    # Only rows with an actual posterior (finite ccf_hat) get boundary defaults; rows
    # with no CCF posterior (e.g. missing copy number) stay all-NaN so ccf_hat,
    # ccf_CI95_low, and ccf_CI95_high are reported consistently as NaN rather than as a
    # spurious [0, 1] interval paired with an NaN point estimate.
    computed = ~np.isnan(ccf_hat)
    nix1 = computed & np.isnan(ccf_ci95[:, 0])
    ccf_ci95[nix1, 0] = CCF_GRID.min()
    nix2 = computed & np.isnan(ccf_ci95[:, 1])
    ccf_ci95[nix2, 1] = CCF_GRID.max()
    ix = computed & (ccf_ci95[:, 1] > CCF_GRID[-2])
    ccf_ci95[ix, 1] = 1.0

    return pd.DataFrame({
        "ccf_hat": np.round(ccf_hat, 2),
        "ccf_CI95_low": np.round(ccf_ci95[:, 0], 2),
        "ccf_CI95_high": np.round(ccf_ci95[:, 1], 2),
    })


def power_calc(n: int, eps: float, fdr: float, delta: float, ssnv_skew: float, rho: float) -> float:
    if pd.isna(n) or pd.isna(delta) or n < 0:
        return np.nan
    n = int(n)
    k = np.arange(n + 1)

    p = st.binom.pmf(k, n=n, p=eps)
    pv = 1 - np.cumsum(np.concatenate(([0.0], p)))
    ks_candidates = np.flatnonzero(pv <= fdr)
    if len(ks_candidates) == 0:
        return 0.0
    ks = int(ks_candidates[0])
    if ks <= 1:
        return 1.0

    denom = pv[ks - 1] - pv[ks]
    cval = 0.0 if denom == 0 else (fdr - pv[ks]) / denom

    skewed_delta = np.clip(ssnv_skew * delta, 0.0, 1.0)
    tiny = np.finfo(float).tiny
    A = max(skewed_delta * rho, tiny)
    B = max((1 - skewed_delta) * rho, tiny)

    p1 = d_beta_binom(k, A, B, n)
    power = 1 - np.sum(p1[: ks - 1]) + cval * p1[ks - 2]
    return float(np.clip(power, 0.0, 1.0))


def power_calc_for_single_read(n: int, delta: float, ssnv_skew: float, rho: float) -> float:
    if pd.isna(n) or pd.isna(delta) or n < 0:
        return np.nan
    n = int(n)
    skewed_delta = np.clip(ssnv_skew * delta, 0.0, 1.0)
    tiny = np.finfo(float).tiny
    A = max(skewed_delta * rho, tiny)
    B = max((1 - skewed_delta) * rho, tiny)
    p0 = d_beta_binom(0, A, B, n)
    return float(np.clip(1 - p0, 0.0, 1.0))


def mode_ssnv_power_calc(mut_df: pd.DataFrame, alpha: float, ssnv_skew: float, rho: float, single_read: bool = False) -> np.ndarray:
    eps = 1e-3 / 3
    fdr = 0.5e-7

    cov = (mut_df["t_alt_count"] + mut_df["t_ref_count"]).to_numpy(dtype=float)
    qt = mut_df["q_hat"].to_numpy(dtype=float)
    total_q = mut_df["total_q_hat"].to_numpy(dtype=float)
    qt[np.isnan(qt)] = total_q[np.isnan(qt)]
    # Floor a rounded-to-zero copy number at 1, matching calc_ccf_dens: an observed
    # variant implies >=1 tumor copy, so the locus is not truly deleted.
    qt[np.isfinite(qt) & (qt < 1.0)] = 1.0
    nc = mut_df["normal_allele_count"].to_numpy(dtype=float)
    delta = alpha / (nc * (1 - alpha) + alpha * qt)

    power = np.full(mut_df.shape[0], np.nan, dtype=float)
    for i in range(mut_df.shape[0]):
        if single_read:
            power[i] = power_calc_for_single_read(cov[i], delta[i], ssnv_skew=ssnv_skew, rho=rho)
        else:
            power[i] = power_calc(cov[i], eps=eps, fdr=fdr, delta=delta[i], ssnv_skew=ssnv_skew, rho=rho)
    return power


def fill_immediate_columns(df: pd.DataFrame, purity: float, ssnv_skew: float):
    out = df.copy()
    out["observed_alt"] = out["t_alt_count"]
    out["purity"] = purity
    out["SSNV_skew"] = ssnv_skew
    out["clonal_scna_mut_ix"] = 1
    out["LL"] = out["_ll"]
    return out


def compare_to_absolute(inferred: pd.DataFrame, abs_maf: pd.DataFrame) -> None:
    if inferred.empty or abs_maf.empty:
        return

    # A total-mode ABS_MAF carries no allelic columns; only compare what it has.
    compare_cols = [c for c in [
        "ccf_hat", "ccf_CI95_low", "ccf_CI95_high", "local_cn_a1", "local_cn_a2",
        "detection_power", "detection_power_for_single_read",
    ] if c in abs_maf.columns]
    merged = inferred.merge(
        abs_maf[["_variant_key"] + compare_cols],
        on="_variant_key",
        suffixes=("_infer", "_abs"),
        how="inner",
    )
    if merged.empty:
        message("No overlap between inferred caller variants and ABS_MAF for concordance check.")
        return

    def summarize(col):
        if f"{col}_infer" not in merged.columns or f"{col}_abs" not in merged.columns:
            return
        x = pd.to_numeric(merged[f"{col}_infer"], errors="coerce")
        y = pd.to_numeric(merged[f"{col}_abs"], errors="coerce")
        keep = x.notna() & y.notna()
        if keep.sum() < 2:
            message(f"Not enough overlapping non-missing values for {col}.")
            return
        x = x.loc[keep]
        y = y.loc[keep]
        corr = pd.concat([x, y], axis=1).corr().iloc[0, 1]
        diff = np.abs(x - y)
        message(
            f"Correlation between inferred {col} and ABSOLUTE output: {corr:.3f}, "
            f"median |diff| = {np.nanmedian(diff):.3f} ± {np.nanstd(diff):.3f}"
        )

    message(f"Concordance check on {merged.shape[0]} overlapping variants.")
    for col in compare_cols:
        summarize(col)


def calculate_ccf(args):
    sex = standardize_sex(args.sex)
    if args.purity <= 0:
        args.purity = 1
    if args.ploidy <= 0:
        args.ploidy = args.normal_ploidy

    epsilon = 1e-3
    b = -33.67232
    m = 82.77464
    rho = float(np.exp(m * (args.ssnv_skew / 2) + b))

    abs_maf = prepare_abs_maf(read_optional_maf(args.absolute_maf))
    snv_maf = prepare_caller_maf(read_optional_maf(args.snv_maf), "snv_maf") if args.snv_maf else pd.DataFrame()
    indel_maf = prepare_caller_maf(read_optional_maf(args.indel_maf), "indel_maf") if args.indel_maf else pd.DataFrame()
    caller_df = concat_callers(snv_maf, indel_maf)
    seg = prepare_segtab(read_table(args.absolute_segtab, comment="#"), sex=sex, normal_ploidy=args.normal_ploidy)

    # Update local allelic CN in ABS_MAF from the (possibly corrected) segtab.
    # The segtab may have been updated by --allelic_split_focals in
    # map_to_absolute_copy_number, making the original ABSOLUTE local_cn_a1/a2 stale.
    # Only applies in allelic mode, where both segtab and ABS_MAF carry allelic columns.
    if not abs_maf.empty and "rescaled.cn.a1" in seg.columns and "local_cn_a1" in abs_maf.columns:
        abs_maf_remapped = map_variants_to_segments(abs_maf, seg)
        remapped_mask = abs_maf_remapped["T.seg.ix"].notna()
        stale_a1 = remapped_mask & (abs_maf_remapped["local_cn_a1"] != abs_maf["local_cn_a1"])
        stale_a2 = remapped_mask & (abs_maf_remapped["local_cn_a2"] != abs_maf["local_cn_a2"])
        n_updated = int((stale_a1 | stale_a2).sum())
        if n_updated:
            message(f"Updating local_cn_a1/a2 from segtab for {n_updated} ABS_MAF variants.")
            abs_maf.loc[remapped_mask, "local_cn_a1"] = abs_maf_remapped.loc[remapped_mask, "local_cn_a1"]
            abs_maf.loc[remapped_mask, "local_cn_a2"] = abs_maf_remapped.loc[remapped_mask, "local_cn_a2"]

    sample_name = infer_sample_name(args=args, seg=seg, caller_df=caller_df, abs_maf=abs_maf)
    output_maf = os.path.join(args.outdir, f"{sample_name}.ABS_MAF.{args.copy_num_type}.completed.txt")
    os.makedirs(args.outdir, exist_ok=True)

    if caller_df.empty:
        output_columns = choose_output_columns(caller_df=caller_df, abs_maf=abs_maf)
        out = ensure_columns(abs_maf, output_columns)
        out = sort_genomic_df(out, chr_col="Chromosome", start_col="Start_position", end_col="End_position")
        write_table(out, output_maf)
        message(f"Wrote {out.shape[0]} variants to {output_maf}")
        return None

    message(f"Loaded {caller_df.shape[0]} caller variants.")
    if abs_maf.empty:
        message(
            "WARNING: no ABSOLUTE ABS_MAF provided; inferring CCF de novo for all caller "
            "variants from the supplied total copy number (no allelic multiplicity "
            "resolution). If the segtab was produced without ABSOLUTE, its copy number "
            "assumes the segment-length-weighted mean copy number equals ploidy "
            "(sum(W * CN) = ploidy), which can miss whole-genome doubling."
        )
    else:
        message(f"Loaded {abs_maf.shape[0]} ABS_MAF variants.")

    mapped = map_variants_to_segments(caller_df, seg)
    n_mapped = int(mapped["T.seg.ix"].notna().sum())
    message(f"Mapped {n_mapped}/{mapped.shape[0]} caller variants to segments.")

    mt_mask = mapped["_chrom"].eq("MT")
    if mt_mask.any():
        message(f"Skipping {int(mt_mask.sum())} mitochondrial caller variants; MT is out of scope for allelic CCF inference.")

    unmapped = mapped["T.seg.ix"].isna()
    if unmapped.any():
        message(f"Skipping {int(unmapped.sum())} caller variants that do not overlap any segment.")

    invalid_counts = mapped[["t_alt_count", "t_ref_count"]].isna().any(axis=1)
    if invalid_counts.any():
        message(f"Skipping {int(invalid_counts.sum())} caller variants with missing tumor read counts.")

    coverage = mapped["t_alt_count"] + mapped["t_ref_count"]
    nonpositive_coverage = coverage <= 0
    if nonpositive_coverage.any():
        message(f"Skipping {int(nonpositive_coverage.sum())} caller variants with non-positive tumor coverage.")

    skip = mt_mask | unmapped | invalid_counts | nonpositive_coverage
    mapped = mapped.loc[~skip].copy().reset_index(drop=True)

    output_columns = choose_output_columns(caller_df=caller_df, abs_maf=abs_maf)
    if mapped.empty:
        if abs_maf.empty:
            out = ensure_columns(pd.DataFrame(), output_columns)
        else:
            out = ensure_columns(abs_maf, output_columns)
            out = sort_genomic_df(out, chr_col="Chromosome", start_col="Start_position", end_col="End_position")
        write_table(out, output_maf)
        message(f"Wrote {out.shape[0]} variants to {output_maf}")
        return None

    ccf_dens = calc_ccf_dens(mapped, alpha=args.purity, ssnv_skew=args.ssnv_skew, rho=rho, epsilon=epsilon)
    ccf_ci = calc_ccf_95ci(ccf_dens)

    for i, col in enumerate(CCF_COLUMNS):
        mapped[col] = ccf_dens[:, i]
    mapped[["ccf_hat", "ccf_CI95_low", "ccf_CI95_high"]] = ccf_ci

    detection_power = mode_ssnv_power_calc(mapped, alpha=args.purity, ssnv_skew=args.ssnv_skew, rho=rho, single_read=False)
    detection_power_single = mode_ssnv_power_calc(mapped, alpha=args.purity, ssnv_skew=args.ssnv_skew, rho=rho, single_read=True)
    mapped["detection_power"] = detection_power
    mapped["detection_power_for_single_read"] = detection_power_single

    log_dens = np.full(mapped.shape[0], np.nan)
    finite_rows = np.all(np.isfinite(ccf_dens), axis=1)
    log_dens[finite_rows] = np.nanmax(np.log(ccf_dens[finite_rows]), axis=1)
    mapped["_ll"] = log_dens

    if "tumor_f" not in mapped.columns:
        mapped["tumor_f"] = np.nan
    missing_tumor_f = mapped["tumor_f"].isna()
    if missing_tumor_f.any():
        denom = mapped["t_alt_count"] + mapped["t_ref_count"]
        mapped.loc[missing_tumor_f, "tumor_f"] = np.where(
            denom.loc[missing_tumor_f] > 0,
            mapped.loc[missing_tumor_f, "t_alt_count"] / denom.loc[missing_tumor_f],
            np.nan,
        )

    inferred = fill_immediate_columns(mapped, purity=args.purity, ssnv_skew=args.ssnv_skew)
    compare_to_absolute(inferred=inferred, abs_maf=abs_maf)
    inferred = ensure_columns(inferred, output_columns)

    if abs_maf.empty:
        out = inferred[output_columns].copy()
        message(f"ABS_MAF absent; writing inferred values for all {out.shape[0]} caller variants.")
    else:
        abs_keys = set(abs_maf["_variant_key"].astype(str))
        rescued = inferred.loc[~inferred["_variant_key"].isin(abs_keys)].copy()
        out = pd.concat([ensure_columns(abs_maf, output_columns), rescued[output_columns]], axis=0, ignore_index=True, sort=False)
        message(f"Number of rescued variants: {rescued.shape[0]}")

    out["Chromosome"] = out["Chromosome"].map(normalize_chromosome)
    out = to_numeric(out, ["Start_position", "End_position"])
    out = sort_genomic_df(out, chr_col="Chromosome", start_col="Start_position", end_col="End_position")
    write_table(out.drop(columns=[c for c in out.columns if c.startswith("_")]), output_maf)
    message(f"Wrote {out.shape[0]} variants to {output_maf}")
    return None


def main():
    args = parse_args()
    print_args(args)
    calculate_ccf(args=args)


if __name__ == "__main__":
    main()
