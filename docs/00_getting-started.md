# Getting started

This page is your starting point for understanding the workflow. If you are a new user, please review the documentation in the following order:

- [README.md](https://github.com/phylyc/somatic_workflow/blob/v2.0.0/README.md)
- [Getting started](00_getting-started.md)
- [Supported use cases](01_supported-use-cases.md)
- [Inputs](02_inputs.md)
- [Outputs](03_outputs.md)
- [CNV calling](04_cnv-calling.md)
- [SNV calling](05_snv-calling.md)
- [ABSOLUTE review and rerun](06_absolute-review-and-rerun.md)
- [PhylogicNDT](07_phylogicndt.md)
- [Ancestry calling](08_ancestry-calling.md)
- [Failure recovery](09_failure-recovery.md)
- [Reference resources](10_resources.md)
- [References / citations](11_references.md)

## Intended audience

The primary audience is analysts and bioinformaticians working with bulk DNA sequencing data. The documentation assumes familiarity with basic concepts such as:

- reference genomes
- aligned BAM/BAI files
- SNVs, indels, and CNVs
- matched-normal versus tumor-only analysis

## Minimal input model

A minimal run generally includes:

- one `patient_id`
- one or more tumor `BAM`/`BAI` pairs
- per-run target intervals
- a reference FASTA, FASTA index, and sequence dictionary
- germline-resource files
- common germline alleles
- an SNV panel of normals (see docs/05_snv-calling.md for more details)
- an optional CNV panel of normals per run (see docs/04_cnv-calling.md)
- a realignment index image if the realignment artifact filter is enabled

## Primary supported modes

The workflow is designed primarily for:
- multi-sample (usually including matched normal sample) whole-exome / whole-genome sequencing data
- matched tumor-normal analysis for whole-exome / whole-genome sequencing data
- multi-sample tumor-only whole-exome / whole genome-sequencing
- single sample tumor only whole-exome / whole-genome sequencing

Matched normals improve results, but they are not required. Tumor-only analysis is a first-class mode and is expected to yield sensible results.

## Common first-run suggestions

- Prefer a platform-appropriate panel of normals where available.
- Treat the first full run as the setup for the ABSOLUTE review step. Rerun with ABSOLUTE Cache inputs filled out (review docs/06_absolute-review-and-rerun.md for more details)
- Use Terra/Cromwell call caching where helpful, but if inputs or core settings change substantially, rerunning is safer than assuming cache reuse is valid.
