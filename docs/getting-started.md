# Getting started

This page is a practical starting point for running the workflow.

## Recommended execution path

The main public execution path is:

1. import the workflow into Terra from Dockstore  
2. provide workflow inputs and reference resources  
3. run the workflow once from BAMs  
4. inspect the ABSOLUTE plots and choose one solution per sample  
5. rerun the workflow with cached outputs plus `Cache.absolute_solution` to complete ABSOLUTE extraction and downstream clonal analysis

`<TODO_TERRA_WORKFLOW_LINK>`

## Minimal input model

A minimal sensible run generally includes:

- one `patient_id`
- one or more tumor `BAM`/`BAI` pairs
- per-run target intervals
- a reference FASTA, FASTA index, and sequence dictionary
- germline-resource files
- common germline alleles
- an SNV panel of normals
- an optional CNV panel of normals per run
- a realignment index image if the realignment artifact filter is enabled

Matched normals improve results, but tumor-only analysis is supported and still yields sensible output.

## Minimal example JSON

This example shows a matched tumor-normal WES-style structure. Replace all file paths with real locations.

```json
{
  "MultiSampleSomaticWorkflow.patient_id": "Patient",
  "MultiSampleSomaticWorkflow.sex": "XX",
  "MultiSampleSomaticWorkflow.bams": [
    "file1.bam",
    "file2.bam"
  ],
  "MultiSampleSomaticWorkflow.bais": [
    "file1.bai",
    "file2.bai"
  ],
  "MultiSampleSomaticWorkflow.sample_names": [
    "Sample1",
    "Sample2"
  ],
  "MultiSampleSomaticWorkflow.normal_sample_names": [
    "Sample2"
  ],
  "MultiSampleSomaticWorkflow.target_intervals": [
    "WES_target_intervals.interval_list",
    "WES_target_intervals.interval_list"
  ],
  "MultiSampleSomaticWorkflow.cnv_panel_of_normals": [
    "WES_pon.hdf5",
    "WES_pon.hdf5"
  ],
  "MultiSampleSomaticWorkflow.Files.ref_dict": "gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.dict",
  "MultiSampleSomaticWorkflow.Files.ref_fasta": "gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta",
  "MultiSampleSomaticWorkflow.Files.ref_fasta_index": "gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai",
  "MultiSampleSomaticWorkflow.Files.germline_resource": "gnomad_sites.af_only.vcf.gz",
  "MultiSampleSomaticWorkflow.Files.germline_resource_idx": "gnomad_sites.af_only.vcf.gz.tbi",
  "MultiSampleSomaticWorkflow.Files.common_germline_alleles": "gnomad_sites.AFgt0.05.vcf.gz",
  "MultiSampleSomaticWorkflow.Files.common_germline_alleles_idx": "gnomad_sites.AFgt0.05.vcf.gz.tbi",
  "MultiSampleSomaticWorkflow.Files.snv_panel_of_normals": "snv_pon.vcf.gz",
  "MultiSampleSomaticWorkflow.Files.snv_panel_of_normals_idx": "snv_pon.vcf.gz.tbi",
  "MultiSampleSomaticWorkflow.Files.realignment_bwa_mem_index_image": "reference.fasta.img"
}
```

For tumor-only analysis, leave `normal_sample_names` empty or omit it.

## What to open first after the first run

The most important first-pass review files are:

- ABSOLUTE PDF plots
- copy-ratio PDF plots
- PhylogicNDT HTML report, if produced

See [Outputs and QC](outputs-and-qc.md) for details.

## Common first-run advice

- Start with a matched tumor-normal WES case if possible.
- Use tumor-only mode if no matched normal exists; the workflow is designed to still produce useful output.
- Prefer a platform-appropriate panel of normals where available.
- Treat the first full run as the setup for the ABSOLUTE review step.
- Use Terra/Cromwell call caching where helpful, but if inputs or core settings change substantially, rerunning is safer than assuming cache reuse is valid.

## Public test dataset

`<TODO_PUBLIC_DATASET_URL>`

A small public smoke-test dataset is still to be selected.
