# CNV calling

## Goal

This stage infers denoised total copy ratios, allelic copy ratios, segmentation, and sample-level copy-number events.

## Inputs

Main inputs include:

- target intervals
- total read counts
- allelic pileups at common germline sites
- optional rare germline sites from the SNV workflow
- CNV panel of normals, if available
- contamination estimates
- reference annotations for target intervals, where applicable

## Outputs

Key outputs include:

- contamination tables
- segmentation tables
- copy-ratio plots
- modeled segments
- sample-level copy-ratio event calls
- SNP correlation QC outputs

## User-facing method summary

CNV calling is built around:

1. collecting total and allelic read counts
2. harmonizing targets across samples/platforms
3. genotyping germline SNPs across all samples
4. segmenting total and allelic copy-ratio signal jointly across samples
5. fitting per-sample copy-ratio profiles on the shared segmentation
6. generating event calls and plots

The multi-sample design improves segmentation stability and helps recover signal in low-purity samples.

## Panels of normals

A CNV PoN is preferred but not strictly mandatory.

Recommended practice:

- use a platform-specific PoN
- use a sex-specific PoN as well if accurate X/Y copy-number interpretation matters

If no suitable PoN exists for a run, an empty file can be supplied so the workflow falls back to denoising based on target-annotation information.

## Important caveat for mixed platforms

For CNV analysis, signal is restricted to the intersection of targets across the included runs/samples. If the interval sets differ too much, CNV output can become poor or unusable.

Practical guidance:

- avoid combining WES and narrow targeted panels for CNV analysis
- in mixed-modality settings, consider excluding less informative runs from `use_sample_for_tCR`
- treat mixed-platform CNV results conservatively unless interval overlap is substantial

## Key QC output: SNP sample correlation

The genotype-correlation outputs are a critical identity check.

Users should inspect:

- `snp_sample_correlation`
- `snp_sample_correlation_min`

Suggested investigation thresholds:

- WES: investigate if `snp_sample_correlation_min < 0.9`
- WGS: investigate if `snp_sample_correlation_min < 0.99`

Values below these levels should prompt a review of whether all samples truly come from the same individual.

## Assay-specific side notes

- small targeted panels are generally too sparse for robust CNV analysis
- ULP WGS can support CNV inference, but purity/ploidy estimation is not considered supported in that setting
