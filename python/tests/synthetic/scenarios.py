"""Derive degenerate-input scenarios from the synthetic cohort.

These reproduce, under the one synthetic-data design, the edge cases the legacy
test_*.sh scripts exercised: empty inputs, single sample, duplicate sample names,
no normal, empty / partial copy-ratio interval intersections, and chr-prefixed
contigs. Most scenarios are just particular *invocations* over the generated files
(handled in the tests); the few that need new files are produced here.
"""

from __future__ import annotations

import pathlib

from . import emit
from .assays import PROFILES
from .model import CONTIGS, Assay, Patient, Sample, chr_name


def write_empty_pileup(path, sample_name: str = "EMPTY") -> str:
    path = pathlib.Path(path)
    with open(path, "w") as fh:
        fh.write(f"#<METADATA>SAMPLE={sample_name}\n")
        fh.write(emit.PILEUP_HEADER)
    return str(path)


def write_empty_cr(path, chr_prefixed: bool = False) -> str:
    path = pathlib.Path(path)
    with open(path, "w") as fh:
        for contig, length in CONTIGS:
            fh.write(f"@SQ\tSN:{chr_name(contig, chr_prefixed)}\tLN:{length}\n")
        fh.write("CONTIG\tSTART\tEND\tLOG2_COPY_RATIO\n")
    return str(path)


def write_cr_for_target_set(path, patient: Patient, sample: Sample, target_set,
                            rng, chr_prefixed: bool = False) -> str:
    """A denoised-CR file for the sample restricted to a specific target set,
    used to build disjoint (WES_A vs WES_C) and partial (WES_A vs WES_B) overlaps."""
    emit.write_denoised_cr(path, patient, sample, target_set, rng, chr_prefixed,
                           PROFILES[Assay.WES].cr_depth)
    return str(path)
