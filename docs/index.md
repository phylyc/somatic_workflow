# Multi-Sample Somatic Workflow

The **MultiSampleSomaticWorkflow (MSSW)** is a pipeline for somatic analysis of tissue
or liquid biopsy samples from cancer patients (WES, WGS, and panel DNA sequencing).
It is a WDL workflow built on GATK4 covering preprocessing, SNV calling, CNV calling,
and clonal analysis.

!!! tip "New here?"
    Start with **[Getting started](00_getting-started.md)**, then read the
    **[Inputs](02_inputs.md)** and **[Outputs](03_outputs.md)** pages.

## Documentation

| Page | What it covers |
| --- | --- |
| [Getting started](00_getting-started.md) | First-run walkthrough and orientation |
| [Supported use cases](01_supported-use-cases.md) | What data and scenarios are supported |
| [Inputs](02_inputs.md) | Full input reference and minimal run template |
| [Outputs](03_outputs.md) | Reviewing and understanding outputs |
| [CNV calling](04_cnv-calling.md) | Copy number variation calling details |
| [SNV calling](05_snv-calling.md) | Single nucleotide variant calling details |
| [ABSOLUTE review and rerun](06_absolute-review-and-rerun.md) | Manual solution selection and re-run |
| [PhylogicNDT](07_phylogicndt.md) | Phylogenetic and clonal dynamics analysis |
| [Ancestry calling](08_ancestry-calling.md) | Genetic ancestry inference |
| [Failure recovery](09_failure-recovery.md) | Recovering from failed shards/tasks |
| [Reference resources](10_resources.md) | Reference and resource bundles |
| [References](11_references.md) | Citations for underlying tools |

## Useful links

- Source repository: [phylyc/somatic_workflow](https://github.com/phylyc/somatic_workflow)
- Import via [Dockstore](https://dockstore.org/workflows/github.com/phylyc/somatic_workflow/MultiSampleSomaticWorkflow:master)

