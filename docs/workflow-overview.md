# Workflow overview

This workflow turns bulk DNA sequencing data into somatic and germline calls, copy-number profiles, purity/ploidy estimates, and optional clonal reconstructions.

## End-to-end flow

```text
BAM/BAI + intervals + references/resources
    ↓
coverage and allelic counts
    ↓
multi-sample SNV/indel discovery and filtering
    ↓
multi-sample CNV modeling
    ↓
ABSOLUTE purity/ploidy analysis
    ↓
manual review of ABSOLUTE plots
    ↓
ABSOLUTE extraction and rescue of dropped segments/variants
    ↓
optional PhylogicNDT clustering, tree building, and timing
```

## Stage 1: coverage and allelic counts

The workflow first collects:

- callable loci
- total read counts in target intervals
- allelic pileups at common germline sites
- contamination estimates
- first-pass segmentation used downstream

These outputs support both variant filtering and CNV modeling.

## Stage 2: SNV/indel calling

Somatic short variants are called jointly across samples. The workflow combines a Mutect1-based force-call strategy with multi-sample Mutect2 calling, followed by statistical filtering, hard filtering, optional realignment-artifact filtering, and functional annotation.

## Stage 3: CNV calling

Copy-number inference uses denoised total copy ratios plus allelic counts. It includes harmonization across platforms, multi-sample SNP genotyping, multi-sample segmentation, and per-sample copy-ratio/event calling.

## Stage 4: ABSOLUTE

ABSOLUTE is run per sample to infer purity, ploidy, and absolute copy number.

This stage is intentionally split into two parts:

1. generate candidate ABSOLUTE solutions and plots
2. manually choose one solution per sample and rerun the workflow to extract final ABSOLUTE-based outputs

## Stage 5: PhylogicNDT

When appropriate inputs exist, the workflow can run PhylogicNDT to cluster somatic events across samples, build phylogenetic trees, and estimate molecular timing. In this implementation, focal CNAs can also be incorporated as timing events.

## Why multi-sample analysis?

Multi-sample analysis improves sensitivity by borrowing evidence across samples from the same patient. This is especially useful in low-purity settings and for subclonal variants that may be weak in one sample but clearer in another.
