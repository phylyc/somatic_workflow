# CNV calling

CNV calling is built around:

1. collecting total and allelic read counts
2. harmonizing targets across samples/platforms
3. genotyping germline SNPs across all samples
4. segmenting total and allelic copy-ratio signal jointly across samples
5. fitting per-sample copy-ratio profiles on the shared segmentation
6. generating event calls and plots

The multi-sample design improves segmentation stability and helps recover signal in low-purity samples.

## Goal

From BAM-level read counts and germline SNP pileups, produce per-sample segmented copy-ratio profiles (total and allelic).

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

## Panels of normals

A CNV PoN is **highly recommended** but not mandatory.

Recommended practice:

- use samples sequenced on the same platform to create the PoN
- Ideally, the PoN is created using at least 40 normal samples
- use a sex-specific PoN as well if accurate X/Y copy-number interpretation matters
- After generating the `cnv_panel_of_normals`, load the HDF5 file and create an elbow plot to visualize the singular values of the eigensamples. Determine the appropriate number of eigensamples by identifying where the elbow plot begins to plateau. If this number is less than the total number of normal BAMs, rerun the workflow with an updated `number_of_eigensamples`.

If no suitable PoN exists for a run, an empty file can be supplied so the workflow falls back to denoising based on target-annotation information.

## Assay-specific side notes

- small targeted panels are generally too sparse for robust CNV analysis
- ULP WGS can support CNV inference, but purity/ploidy estimation is not considered supported in that setting

## Important caveat for mixed platforms

For CNV analysis, signal is restricted to the intersection of targets across the included runs/samples. If the interval sets differ too much, CNV output can become poor or unusable.

Practical guidance:

- avoid combining WES and narrow targeted panels for CNV analysis
- in mixed-modality settings, consider excluding less informative runs using `use_sample_for_tCR`
- treat mixed-platform CNV results conservatively unless interval overlap is substantial
