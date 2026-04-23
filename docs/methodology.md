# Methodology

This page gives a scientific overview of the workflow. It is adapted from publication-style methods text and intentionally keeps the user-facing documentation pages lighter.

## High-level scientific model

Genomic variation is characterized with a multi-sample analysis pipeline that jointly analyzes all samples from the same patient. Joint somatic calling leverages evidence across samples, increasing sensitivity in low-purity settings and for subclonal variants. The workflow combines somatic short-variant calling, copy-number inference, purity/ploidy estimation, and clonal reconstruction.

## Stage summary

### SNV calling

The workflow calls somatic SNVs and indels using a combination of Mutect1-derived candidates and multi-sample Mutect2 calling. Filtering incorporates:

- read orientation bias
- contamination estimates
- allelic copy-ratio segmentations
- germline resources
- sequencing-artifact resources
- hard filters on depth and quality metrics
- realignment-artifact filtering

Passing variants are functionally annotated. The workflow is designed to support both matched-normal and tumor-only analysis.

### CNV calling

Copy-number analysis combines:

- total read counts in target intervals
- allelic counts at common germline sites
- harmonization across sequencing platforms
- contamination estimation
- multi-sample SNP genotyping
- multi-sample segmentation
- per-sample total and allelic copy-ratio inference

The multi-sample genotyping step is central to robust allelic copy-ratio modeling across all samples from the same patient.

### ABSOLUTE

Purity, ploidy, and absolute copy number are inferred per sample from copy-ratio segmentation and somatic short variants. The workflow then supports rescue of segments and variants that fall outside the raw ABSOLUTE export but can be completed using the chosen purity/ploidy solution.

### PhylogicNDT

Somatic mutations are clustered across samples from the same patient using CCF information. Phylogenetic trees are then reconstructed, and relative molecular timing is estimated. In this implementation, focal CNAs can also be incorporated as timing events.

## Intended level of detail

This page is a narrative scientific overview rather than a mathematical derivation. For practical usage, see the stage-specific pages:

- [SNV calling](snv-calling.md)
- [CNV calling](cnv-calling.md)
- [ABSOLUTE review and rerun](absolute-review-and-rerun.md)
- [PhylogicNDT](phylogicndt.md)

## Future additions

The following are expected future extensions or documentation additions:

- TMB calculation details
- mutational-signature estimation details
- a publication link once available
