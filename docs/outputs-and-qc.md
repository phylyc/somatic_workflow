# Outputs and QC

This page focuses on the files that most users should look at first.

## Open these first

For most analyses, start with:

- `absolute_acr_plot` (PDF)
- `cr_plot` (PDF)
- `phylogic_report` (HTML), if present

These are the most directly interpretable summary outputs.

## Primary deliverables

The main deliverables are:

- `absolute_segtab`
- `absolute_maf`
- `contamination_table`
- `gvcf`
- `snp_sample_correlation`
- PhylogicNDT outputs, especially `phylogic_report`

Depending on the project, somatic and germline VCF outputs may also be primary deliverables.

## Important QC checks

### 1. Sample identity / consistency

Inspect:

- `snp_sample_correlation`
- `snp_sample_correlation_min`

Suggested thresholds for further investigation:

- WES: investigate below `0.9`
- WGS: investigate below `0.99`

### 2. Contamination

Inspect the contamination tables.

A contamination estimate above roughly `0.05` should prompt closer review.

### 3. Segmentation quality

Inspect `cr_plot` and related segmentation outputs.

Potential warning signs:

- extremely fragmented segmentation
- more than roughly `1000` segments
- noisy profiles with weak plateau structure

### 4. ABSOLUTE plausibility

Inspect `absolute_acr_plot`.

Potential warning signs:

- poor agreement between proposed copy states and visible segmental plateaus
- poor agreement between purity expectations and SNV VAF modes
- inconsistent whole-genome doubling interpretation across samples

### 5. Normal sample behavior

If a designated normal sample shows visible aneuploidy in the copy-ratio plots, that does not automatically invalidate the run. In this workflow, a "normal" sample is only weakly privileged as germline normal; if the data support substantial tumor content, it can effectively behave like a tumor sample.

## Mostly intermediate or cache-style outputs

Many other workflow outputs are primarily intermediate data products or cacheable artifacts. They remain useful for debugging and reruns, but are not the first files most users should inspect.
