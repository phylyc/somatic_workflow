import argparse
import pandas as pd
import os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--patient_id", required=True)
    parser.add_argument("--absolute_maf", nargs="+", required=True)
    parser.add_argument("--absolute_segtab", nargs="*", required=False)
    parser.add_argument("--absolute_purity", nargs="+", required=True, type=float)
    parser.add_argument("--timepoints", nargs='*', type=int, required=False, help="Optional explicit timepoint ordering (one per sample)")
    parser.add_argument("--outfile", required=True)
    args = parser.parse_args()

    # Assume sample_id is just the MAF filename minus ".ABS_MAF.txt"
    sample_ids = [os.path.basename(x).replace(".ABS_MAF.txt", "") for x in args.absolute_maf]
    n = len(sample_ids)

    maf_fns = args.absolute_maf
    seg_fns = args.absolute_segtab if args.absolute_segtab else [""] * n
    purities = args.absolute_purity
    timepoints = args.timepoints if args.timepoints else list(range(1, n+1)) # Simple sequential timepoint assignment (can be changed)

    if len(maf_fns) != n:
        raise ValueError(f"Number of maf_fns ({len(maf_fns)}) does not match number of samples ({n}).")
    if len(seg_fns) != n:
        raise ValueError(f"Number of seg_fns ({len(seg_fns)}) does not match number of samples ({n}).")
    if len(purities) != n:
        raise ValueError(f"Number of purities ({len(purities)}) does not match number of samples ({n}).")
    if len(timepoints) != n:
        raise ValueError(f"Number of timepoints ({len(timepoints)}) does not match number of samples ({n}).")

    df = pd.DataFrame({
        "sample_id": sample_ids,
        "maf_fn": maf_fns,
        "seg_fn": seg_fns,
        "purity": purities,
        "timepoint": timepoints,
    })

    df.to_csv(args.outfile, sep="\t", index=False)

if __name__ == "__main__":
    main()