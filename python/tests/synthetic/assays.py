"""Assay profiles and the genomic target/bin layout.

Minimal, not biologically realistic: just enough structure that each assay's
characteristic data shape (WES = CR+allelic, ULP-WGS = CR-only, deep-panel =
allelic-only) and the harmonize target-overlap cases fall out.
"""

from __future__ import annotations

from dataclasses import dataclass

from .model import CONTIGS, Assay

BIN_SIZE = 50_000  # copy-ratio bin width on the synthetic reference


@dataclass(frozen=True)
class AssayProfile:
    name: str
    has_cr: bool         # provides copy-ratio bins (tCR)
    has_allelic: bool    # provides usable allelic depth (aCR)
    het_depth: int       # mean read depth at SNP sites (0 => ~no allelic data)
    cr_depth: int        # mean read depth backing a CR bin (affects log2CR noise)


PROFILES: dict[Assay, AssayProfile] = {
    Assay.WES: AssayProfile("WES", has_cr=True, has_allelic=True, het_depth=120, cr_depth=150),
    Assay.ULP_WGS: AssayProfile("ULP_WGS", has_cr=True, has_allelic=False, het_depth=0, cr_depth=40),
    Assay.DEEP_PANEL: AssayProfile("DEEP_PANEL", has_cr=False, has_allelic=True, het_depth=800, cr_depth=0),
}


def all_cr_bins() -> list[tuple[str, int, int, int]]:
    """Genome-wide CR bins (contig, start, end, per-contig index). MT excluded."""
    bins = []
    for contig, length in CONTIGS:
        if contig == "MT":
            continue
        idx, start = 0, 1
        while start + BIN_SIZE - 1 <= length:
            bins.append((contig, start, start + BIN_SIZE - 1, idx))
            start += BIN_SIZE
            idx += 1
    return bins


def target_bins(target_set: str | None) -> list[tuple[str, int, int, int]]:
    """CR bins covered by a target set.

    ULP-WGS (None) is genome-wide. WES panels are subsets chosen so that:
      WES_A (even idx) and WES_C (odd idx) are DISJOINT  -> empty intersection,
      WES_A (even idx) and WES_B (idx % 3 != 0) PARTIALLY overlap.
    PANEL is a handful of bins (used for allelic targets, not segmentation).
    """
    bins = all_cr_bins()
    if target_set is None:
        return bins
    if target_set == "WES_A":
        return [b for b in bins if b[3] % 2 == 0]
    if target_set == "WES_B":
        return [b for b in bins if b[3] % 3 != 0]
    if target_set == "WES_C":
        return [b for b in bins if b[3] % 2 == 1]
    if target_set == "PANEL":
        return [b for b in bins if b[3] % 7 == 0][:6]
    raise ValueError(f"unknown target set {target_set!r}")
