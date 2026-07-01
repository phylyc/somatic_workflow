"""CLI entry point for the synthetic golden-cohort generator.

From one ground-truth tumour model it builds the cohort, emits the ground-truth
answer tables + manifest, and forward-simulates every upstream tool-output file
(pileups, denoised CR, ModelSegments, VCFs, ABSOLUTE/SOMIX segtabs, MAFs). The
test fixtures regenerate this into a temp dir each run; see python/TESTING.md.

Usage:
    python tests/synthetic/make_synthetic_testing_data.py --outdir tests/synthetic/data [--seed 0] [--chr]
"""

from __future__ import annotations

import argparse
import sys

# Support both `python -m tests.synthetic.make_...` and direct invocation.
if __package__ in (None, ""):
    import pathlib
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent.parent))
    from tests.synthetic.cohort import build_cohort
    from tests.synthetic.truth import write_truth
    from tests.synthetic.generate import generate_cohort
else:
    from .cohort import build_cohort
    from .truth import write_truth
    from .generate import generate_cohort


def parse_args(argv=None):
    p = argparse.ArgumentParser(prog="make_synthetic_testing_data",
                                description="Generate the synthetic golden cohort.")
    p.add_argument("--outdir", required=True, help="Output directory.")
    p.add_argument("--seed", type=int, default=0, help="RNG seed for the simulation (tests use 0).")
    p.add_argument("--chr", action="store_true", help="Emit chr-prefixed contig naming.")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    cohort = build_cohort()
    written = write_truth(cohort, args.outdir)
    index = generate_cohort(cohort, args.outdir, seed=args.seed, chr_prefixed=args.chr)
    n_files = sum(len(s["runs"]) * 2 + 4 for p in index["patients"].values()
                  for s in p["samples"].values())
    print(f"Wrote truth + manifest + ~{n_files} tool-output files for "
          f"{len(cohort.patients)} patients to {args.outdir}")
    for name, path in written.items():
        print(f"  {name}: {path}")


if __name__ == "__main__":
    main()
