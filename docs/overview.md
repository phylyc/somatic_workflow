# Documentation overview

This documentation describes the `somatic_workflow` multi-sample bulk DNA sequencing workflow for somatic variant discovery and clonal analysis.

## What this workflow is for

`somatic_workflow` is an all-in-one best-practice pipeline for genomic characterization of bulk DNA next-generation sequencing data. Starting from aligned sequencing data (`BAM`/`BAI`) plus a reference genome and auxiliary resources, it produces somatic and germline variant calls, copy-number profiles, contamination estimates, purity/ploidy solutions, and optional clonal-phylogenetic outputs.

At a high level, the workflow performs:

1. coverage and allelic-count collection
2. multi-sample somatic SNV/indel calling and filtering
3. multi-sample CNV modeling
4. per-sample ABSOLUTE purity/ploidy analysis
5. optional PhylogicNDT clonal reconstruction and timing

## Intended audience

The primary audience is analysts and bioinformatics users working with somatic bulk DNA sequencing data. The documentation assumes familiarity with basic concepts such as:

- reference genomes
- aligned BAM/BAI files
- SNVs, indels, and CNVs
- matched-normal versus tumor-only analysis

## Documentation map

- [Getting started](getting-started.md)
- [Supported use cases](supported-use-cases.md)
- [Input model](input-model.md)
- [Workflow overview](workflow-overview.md)
- [SNV calling](snv-calling.md)
- [CNV calling](cnv-calling.md)
- [ABSOLUTE review and rerun](absolute-review-and-rerun.md)
- [PhylogicNDT](phylogicndt.md)
- [Outputs and QC](outputs-and-qc.md)
- [Reruns and failure recovery](reruns-and-failure-recovery.md)
- [Reference resources](resources.md)
- [Methodology](methodology.md)

## Recommended reading order

For a new user:

1. read [Getting started](getting-started.md)
2. skim [Supported use cases](supported-use-cases.md)
3. read [Workflow overview](workflow-overview.md)
4. read [Outputs and QC](outputs-and-qc.md)
5. use [ABSOLUTE review and rerun](absolute-review-and-rerun.md) for the second-pass run

## Status and release model

These docs are written in release-style language and should describe the workflow as the canonical user-facing version. For scientific details beyond the operational behavior documented here, see [Methodology](methodology.md).

## Placeholders

A few sections still contain placeholders such as:

- `<TODO_PUBLIC_DATASET_URL>`
- `<TODO_TERRA_WORKFLOW_LINK>`
- `<TODO_JSON_TEMPLATE>`

These are intentional and can be filled in later.
