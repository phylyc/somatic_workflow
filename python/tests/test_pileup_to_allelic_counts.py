"""Unit / IO / stat tests for pileup_to_allelic_counts.py."""

import gzip

import numpy as np
import pandas as pd
import pytest

import pileup_to_allelic_counts as p2a
from helpers import args_ns, pileup_df, write_pileup, run_script, synth_wes_pileup


def write_gvcf(path, records, gz=False):
    """records: (contig, pos, ref, alt, gt) -> 10-column gvcf body (no header)."""
    lines = ["\t".join([str(c), str(p), ".", ref, alt, ".", "PASS", ".", "GT", gt])
             for (c, p, ref, alt, gt) in records]
    text = "\n".join(lines) + "\n"
    if gz:
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        path.write_text(text)
    return str(path)


def write_intervals(path, rows):
    df = pd.DataFrame(rows, columns=["CONTIG", "START", "END"])
    with open(path, "w") as fh:
        fh.write("@HD\tVN:1.6\n")
        df.to_csv(fh, sep="\t", index=False)
    return str(path)


def base_args(pileup, gvcf, output, **kw):
    defaults = dict(
        pileup=pileup, gvcf=gvcf, output=output, intervals=None, ref_dict=None,
        het_to_interval_mapping_max_distance=0, min_read_depth=0,
        error_output=None, select_hets=False, aggregate_hets=False, verbose=False,
    )
    defaults.update(kw)
    return args_ns(**defaults)


def read_output(path):
    return pd.read_csv(path, sep="\t", header=None,
                       names=["contig", "position", "ref_count", "alt_count",
                              "ref", "alt"])


# --------------------------------------------------------------------------- #
# pure helpers
# --------------------------------------------------------------------------- #
@pytest.mark.unit
def test_get_contigs_default_and_dict(ref_dict):
    assert p2a.get_contigs(None)[:2] == ["1", "2"]
    assert p2a.get_contigs(ref_dict) == ["1", "2", "X"]


@pytest.mark.unit
def test_sort_genomic_positions():
    idx = pd.MultiIndex.from_tuples([("2", 1), ("1", 5), ("1", 2)],
                                    names=["contig", "position"])
    out = p2a.sort_genomic_positions(idx, ["1", "2"])
    assert list(out) == [("1", 2), ("1", 5), ("2", 1)]


@pytest.mark.io
def test_read_gvcf_plain_and_gzip(tmp_path):
    plain = write_gvcf(tmp_path / "g.vcf", [("1", 1, "A", "G", "0/1")])
    with p2a.read_gvcf(plain) as fh:
        assert "A\tG" in fh.read()
    gz = write_gvcf(tmp_path / "g.vcf.gz", [("1", 1, "A", "G", "0/1")], gz=True)
    with p2a.read_gvcf(gz) as fh:
        assert "A\tG" in fh.read()


# --------------------------------------------------------------------------- #
# convert core
# --------------------------------------------------------------------------- #
@pytest.mark.stat
def test_basic_join_and_min_depth(tmp_path):
    pu = write_pileup(tmp_path / "p.pileup",
                      pileup_df([("1", 100, 10, 10, 0, 0.2),
                                 ("1", 200, 1, 0, 0, 0.1)]))
    gv = write_gvcf(tmp_path / "g.vcf",
                    [("1", 100, "A", "G", "0/1"), ("1", 200, "C", "T", "1/1")])
    out = tmp_path / "out.tsv"
    p2a.convert_pileup_to_allelic_counts(
        base_args(pu, gv, str(out), min_read_depth=2))
    df = read_output(out)
    assert list(df["position"]) == [100]
    assert df.iloc[0]["ref"] == "A" and df.iloc[0]["alt"] == "G"


@pytest.mark.stat
def test_select_hets_keeps_only_hets(tmp_path):
    pu = write_pileup(tmp_path / "p.pileup",
                      pileup_df([("1", 100, 10, 10, 0, 0.2),
                                 ("1", 200, 0, 20, 0, 0.9)]))
    gv = write_gvcf(tmp_path / "g.vcf",
                    [("1", 100, "A", "G", "0/1"), ("1", 200, "C", "T", "1/1")])
    out = tmp_path / "out.tsv"
    p2a.convert_pileup_to_allelic_counts(
        base_args(pu, gv, str(out), select_hets=True))
    df = read_output(out)
    assert list(df["position"]) == [100]


@pytest.mark.stat
def test_interval_remapping_to_free_position(tmp_path):
    # pos 50 inside [1,100]; pos 150 is 50bp past the end -> remapped within max_dist
    pu = write_pileup(tmp_path / "p.pileup",
                      pileup_df([("1", 50, 5, 5, 0, 0.2),
                                 ("1", 150, 5, 5, 0, 0.2)]))
    gv = write_gvcf(tmp_path / "g.vcf",
                    [("1", 50, "A", "G", "0/1"), ("1", 150, "C", "T", "0/1")])
    iv = write_intervals(tmp_path / "iv.tsv", [("1", 1, 100)])
    out = tmp_path / "out.tsv"
    p2a.convert_pileup_to_allelic_counts(
        base_args(pu, gv, str(out), intervals=iv,
                  het_to_interval_mapping_max_distance=60))
    df = read_output(out)
    # original 50 kept, 150 remapped to the last free position end-1 == 99
    assert set(df["position"]) == {50, 99}


@pytest.mark.stat
def test_interval_drops_far_loci(tmp_path):
    pu = write_pileup(tmp_path / "p.pileup",
                      pileup_df([("1", 50, 5, 5, 0, 0.2),
                                 ("1", 5000, 5, 5, 0, 0.2)]))
    gv = write_gvcf(tmp_path / "g.vcf",
                    [("1", 50, "A", "G", "0/1"), ("1", 5000, "C", "T", "0/1")])
    iv = write_intervals(tmp_path / "iv.tsv", [("1", 1, 100)])
    out = tmp_path / "out.tsv"
    p2a.convert_pileup_to_allelic_counts(
        base_args(pu, gv, str(out), intervals=iv,
                  het_to_interval_mapping_max_distance=10))
    df = read_output(out)
    assert set(df["position"]) == {50}


@pytest.mark.stat
def test_aggregate_hets_oriented_sums(tmp_path):
    # phased hets: 0|1 (alt on B), 1|0 (alt on A) -> oriented aggregation
    pu = write_pileup(tmp_path / "p.pileup",
                      pileup_df([("1", 10, 8, 2, 0, 0.2),
                                 ("1", 20, 3, 7, 0, 0.7)]))
    gv = write_gvcf(tmp_path / "g.vcf",
                    [("1", 10, "A", "G", "0|1"), ("1", 20, "C", "T", "1|0")])
    iv = write_intervals(tmp_path / "iv.tsv", [("1", 1, 100)])
    out = tmp_path / "out.tsv"
    p2a.convert_pileup_to_allelic_counts(
        base_args(pu, gv, str(out), intervals=iv, select_hets=True,
                  aggregate_hets=True))
    df = read_output(out)
    # alt_sum = alt[0|1] + ref[1|0] = 2 + 3 = 5; ref_sum = ref[0|1] + alt[1|0] = 8 + 7 = 15
    row = df.loc[df["position"] == 1].iloc[0]
    assert row["ref_count"] == 15 and row["alt_count"] == 5


@pytest.mark.stat
def test_error_output_clip(tmp_path):
    pu = write_pileup(tmp_path / "p.pileup",
                      pileup_df([("1", 100, 10, 10, 10, 0.2)]))
    gv = write_gvcf(tmp_path / "g.vcf", [("1", 100, "A", "G", "0/1")])
    out = tmp_path / "out.tsv"
    err = tmp_path / "err.txt"
    p2a.convert_pileup_to_allelic_counts(
        base_args(pu, gv, str(out), error_output=str(err)))
    # 1.5 * 10 / 30 = 0.5 -> clipped to a_max 0.05
    assert float(np.loadtxt(err)) == pytest.approx(0.05)


@pytest.mark.stat
def test_empty_pileup_no_output(tmp_path):
    pu = write_pileup(tmp_path / "p.pileup", pileup_df([]))
    gv = write_gvcf(tmp_path / "g.vcf", [("1", 100, "A", "G", "0/1")])
    out = tmp_path / "out.tsv"
    result = p2a.convert_pileup_to_allelic_counts(base_args(pu, gv, str(out)))
    assert result is None
    assert not out.exists()


# --------------------------------------------------------------------------- #
# Black-box expectations (synthetic cohort, contract-derived).
# --------------------------------------------------------------------------- #
@pytest.mark.regression
def test_p2a_alleles_from_gvcf(synth, tmp_path):
    # output ref/alt allele letters come from the gvcf REF/ALT (the pileup has none)
    pileup = synth_wes_pileup(synth, "PXX", "PXX-T1")
    gvcf = synth.files["patients"]["PXX"]["gvcf"]
    out = tmp_path / "ac.tsv"
    run_script("pileup_to_allelic_counts.py",
               ["--pileup", pileup, "--gvcf", gvcf, "--output", str(out)])
    ac = pd.read_csv(out, sep="\t", header=None,
                     names=["contig", "position", "ref_count", "alt_count", "ref", "alt"])
    gv = pd.read_csv(gvcf, sep="\t", comment="#", header=None,
                     names=["CHROM", "POS", "ID", "REF", "ALT", "Q", "F", "I", "FMT", "S"])
    gmap = {(str(r.CHROM), int(r.POS)): (r.REF, r.ALT) for r in gv.itertuples()}
    assert len(ac) > 0
    for r in ac.itertuples():
        assert (r.ref, r.alt) == gmap[(str(r.contig), int(r.position))], (r.contig, r.position)


@pytest.mark.regression
def test_p2a_includes_homozygous_sites_without_select_hets(synth, tmp_path):
    pileup = synth_wes_pileup(synth, "PXX", "PXX-T1")
    gvcf = synth.files["patients"]["PXX"]["gvcf"]
    out = tmp_path / "ac.tsv"
    run_script("pileup_to_allelic_counts.py",
               ["--pileup", pileup, "--gvcf", gvcf, "--output", str(out)])
    ac = pd.read_csv(out, sep="\t", header=None,
                     names=["contig", "position", "ref_count", "alt_count", "ref", "alt"])
    pos = set(zip(ac["contig"].astype(str), ac["position"].astype(int)))
    # 1:300000 is 1/1 and 3:200000 is 0/0 — both present without --select_hets
    assert ("1", 300000) in pos and ("3", 200000) in pos


@pytest.mark.regression
def test_p2a_min_read_depth_filters_output(synth, tmp_path):
    pileup = synth_wes_pileup(synth, "PXX", "PXX-T1")
    gvcf = synth.files["patients"]["PXX"]["gvcf"]
    pu = pd.read_csv(pileup, sep="\t", comment="#")
    pu["depth"] = pu["ref_count"] + pu["alt_count"]
    thr = int(pu["depth"].median())
    out = tmp_path / "ac.tsv"
    run_script("pileup_to_allelic_counts.py",
               ["--pileup", pileup, "--gvcf", gvcf, "--output", str(out),
                "--min_read_depth", str(thr)])
    ac = pd.read_csv(out, sep="\t", header=None,
                     names=["contig", "position", "ref_count", "alt_count", "ref", "alt"])
    kept = set(zip(ac["contig"].astype(str), ac["position"].astype(int)))
    depth = {(str(r.contig), int(r.position)): r.depth for r in pu.itertuples()}
    assert all(depth[k] >= thr for k in kept)
    assert any(d < thr for d in depth.values()) and len(kept) < len(depth)


@pytest.mark.regression
def test_p2a_deterministic_and_filters_only_shrink(synth, tmp_path):
    # rerunning is byte-identical; adding --select_hets can only shrink the row set
    pileup = synth_wes_pileup(synth, "PXX", "PXX-T1")
    gvcf = synth.files["patients"]["PXX"]["gvcf"]
    a, b = tmp_path / "a.tsv", tmp_path / "b.tsv"
    for o in (a, b):
        run_script("pileup_to_allelic_counts.py",
                   ["--pileup", pileup, "--gvcf", gvcf, "--output", str(o)])
    assert a.read_bytes() == b.read_bytes()
    hets = tmp_path / "h.tsv"
    run_script("pileup_to_allelic_counts.py",
               ["--pileup", pileup, "--gvcf", gvcf, "--output", str(hets), "--select_hets"])
    full_pos = {tuple(x.split("\t")[:2]) for x in a.read_text().splitlines()}
    het_pos = {tuple(x.split("\t")[:2]) for x in hets.read_text().splitlines()}
    assert het_pos <= full_pos and len(het_pos) < len(full_pos)


@pytest.mark.regression
def test_pileup_to_allelic_counts_hets(synth, tmp_path):
    # het allelic counts from pileup x gvcf: only hets emitted, counts from pileup
    pileup = synth_wes_pileup(synth, "PXX", "PXX-T1")
    gvcf = synth.files["patients"]["PXX"]["gvcf"]
    out = tmp_path / "ac.tsv"
    err = tmp_path / "err.txt"
    run_script("pileup_to_allelic_counts.py", [
        "--pileup", pileup, "--gvcf", gvcf, "--output", str(out),
        "--error_output", str(err), "--select_hets", "--min_read_depth", "10",
    ])
    ac = pd.read_csv(out, sep="\t", header=None,
                     names=["contig", "position", "ref_count", "alt_count", "ref", "alt"])
    pu = pd.read_csv(pileup, sep="\t", comment="#")

    truth = synth.truth["germline_sites"]
    hets = truth[(truth.patient_id == "PXX") & (truth.genotype == "0/1") & truth.is_snv]
    out_pos = set(zip(ac["contig"].astype(str), ac["position"].astype(int)))
    for _, g in hets.iterrows():
        assert (str(g.contig), int(g.position)) in out_pos, (g.contig, g.position)
    # only hets (no hom sites) selected
    hom_pos = truth[(truth.patient_id == "PXX") & (truth.genotype != "0/1")]
    for _, g in hom_pos.iterrows():
        assert (str(g.contig), int(g.position)) not in out_pos

    # allelic counts copied from the pileup
    merged = ac.merge(pu, on=["contig", "position"], suffixes=("", "_pu"))
    assert (merged["ref_count"] == merged["ref_count_pu"]).all()
    assert 0.001 <= float(np.loadtxt(err)) <= 0.05


@pytest.mark.regression
def test_p2a_sex_chromosome_hets_like_autosomal(synth, tmp_path):
    # chrX hets in a female sample are emitted exactly like autosomal hets — no
    # special-casing of sex contigs in the allelic-count conversion.
    pileup = synth_wes_pileup(synth, "PXX", "PXX-T1")
    gvcf = synth.files["patients"]["PXX"]["gvcf"]
    out = tmp_path / "ac.tsv"
    run_script("pileup_to_allelic_counts.py",
               ["--pileup", pileup, "--gvcf", gvcf, "--output", str(out), "--select_hets"])
    ac = pd.read_csv(out, sep="\t", header=None,
                     names=["contig", "position", "ref_count", "alt_count", "ref", "alt"])
    pos = set(zip(ac["contig"].astype(str), ac["position"].astype(int)))
    # PXX (female, chrX diploid) has 0/1 germline sites at X:400000 and X:500000
    truth = synth.truth["germline_sites"]
    x_hets = truth[(truth.patient_id == "PXX") & (truth.contig.astype(str) == "X") &
                   (truth.genotype == "0/1") & truth.is_snv]
    assert len(x_hets) >= 1
    for _, g in x_hets.iterrows():
        assert ("X", int(g.position)) in pos, g.position
    # the chrX het rows carry the same pileup-derived allelic counts as any other row
    pu = pd.read_csv(pileup, sep="\t", comment="#")
    merged = ac.merge(pu, on=["contig", "position"], suffixes=("", "_pu"))
    assert (merged["ref_count"] == merged["ref_count_pu"]).all()
