# Workflow for Somatic Variant Discovery with GATK4+

This repository contains a WDL (Workflow Description Language) workflow for performing multi-sample somatic variant calling using GATK4. The workflow covers preprocessing, single nucleotide variant (SNV) calling, copy number variation (CNV) calling, and clonal analysis. It is available on [Dockstore](https://dockstore.org/workflows/github.com/phylyc/somatic_workflow/MultiSampleSomaticWorkflow:master) to import into e.g. a [Terra](https://app.terra.bio/) workspace. Alternatively, clone this repository and run the workflow locally using [Cromwell](https://cromwell.readthedocs.io/en/stable/):
```bash
java -jar cromwell.jar run multi-sample_somatic_workflow.wdl --inputs multi-sample_somatic_workflow.inputs.json
```

---
## Concepts

- Patient → Sample → SequencingRuns: a patient has ≥1 sample(s), and a sample may have multiple sequencing runs. Runs are grouped by `sample_names` which are inferred from bam read group names if not provided.
- tCR / aCR: denoised **total** copy ratio and **allelic** copy ratio estimated from harmonized targets and SNP read counts.


## Workflow Structure

The workflow is organized into the following main tasks:

### 1. Preprocessing
- **1.1 Preprocess and Split Calling Intervals** Prepare the intervals for scatter-gather tasks in shards for SNV calling.
- **1.2 Define Patient/Sample/SequencingRun sets**: Define a patient as a set of samples, and each sample as a set of sequencing runs. Sequencing runs are grouped by sample name.

### 2. Coverage Aggregation
- **2.1 Callable loci per sample**: Collect genomic regions per sample with sufficient coverage for SNV detection
- **2.2 Total read counts**: Collect total read counts in target intervals and perform denoising of total copy ratios via tangent-normalization.
- **2.3 Allelic read count**: Collect allelic read counts at common germline sites (SNP panel: `common_germline_alleles`).
- **2.4 Harmonization of target intervals and allelic counts**: Harmonize target intervals across samples, subset to intersection; merge read count data from multiple sequencing runs per sample.
- **2.5 Contamination estimation**: Estimate *out-of-patient* contamination in each sample.
- **2.6 First-pass copy ratio segmentation**: Get a prior single-sample allelic copy ratio segmentation for filtering SNVs (3.3) and genotyping of germline sites (4.1).
  - **2.6.1 Filter total read counts**: Filters total read counts based on first-pass segmentation (turned off by default).

### 3. SNV Calling
- **3.1 Mutect1 single-sample calling**: Use Mutect1 for single-sample mutation calling in tumor-normal mode if a matched normal sample is available, otherwise in tumor-only mode.
- **3.2 Mutect2 multi-sample calling**: Use Mutect2 for multi-sample mutation calling. Force-call alleles that were called via Mutect1.
- **3.3 Filter Variant Calls**: Annotate and select somatic vs germline vs artifactual variants based on various filters.
  - **3.3a Filter**: Apply statistical filters for sequencing artifacts, germline variants, read orientation bias, and contamination, among others.
  - **3.3b Hard Filter**: Apply hard filters based on base quality, mapping quality, fragment length, read depth, read orientation quality, position on the read, and population allele frequency.
  - **3.3c Realignment Filter** Filter based on realignment success (to hg38 or whichever reference is given by `realignment_bwa_mem_index_image`).
- **3.4 Annotate SNVs**: Annotate SNVs with functional information.
- **3.5 Tumor mutational burden (TMB) estimation**: (coming soon)

### 4. CNV Calling
- **4.1 SNP genotyping**: Genotype allelic count data at common (from 2.3) and rare (from 3.3) germline sites using evidence across all samples; harmonize loci across samples; phase HETs using allelic imbalance.
- **4.2 Multi-sample segmentation**: Segment denoised total copy ratios and allelic copy ratio across multiple samples
- **4.3 Per-sample total and allelic copy ratio inference**: Infer total and allelic copy ratios for each sample.
- **4.4 Per-sample event calling**: Call amplifications and deletions for each sample. 
- **4.5 Per-sample segmentation plotting**: Plot the segmented denoised copy ratios and allelic copy ratios for each sample.

### 5. Clonal Analysis
- **5.1 ABSOLUTE** Perform per-sample clonal analysis of the identified somatic variants and estimate tumor purity, ploidy, and absolute copy number.
- **5.2 ABSOLUTE extraction**: Extract results for one chosen solution (needs manual input).
  - **5.2a**: Rescue dropped segments.
- **5.3 PhylogicNDT**: Build phylogenetic trees from all samples, perform growth kinetics, and timing analysis.

Please remember to always review the intermediate results to ensure that the final results are as expected. Inappropriate filtering or parameter settings can lead to misleading output.

---
## Expected Inputs

The input arguments are grouped into 
```wdl
MultiSampleSomaticWorkflow # (mandatory + optional)
Cache # (workflow output)
Files # (reference data resources)
Parameters # (workflow parameters)
RuntimeParameters
```
### A) Core input (`MultiSampleSomaticWorkflow`)

| Name                                        | Type             | Required | Notes                                                                                     |
| ------------------------------------------- |------------------| :------: |-------------------------------------------------------------------------------------------|
| `patient_id`                                | String           |     ✓    | Label used across outputs                                                                 |
| `sex`                                       | String?          |          | “XX” or “XY”                                                                              |
| `sample_names`                              | Array\[String]?  |          | Groups runs into samples; if omitted, inferred from BAMs                                  |
| `timepoints`                                | Array\[Int]?     |          | Ordinal or days for phylogenetic timing                                                   |
| `bams` / `bais`                             | Array\[File]+    |     ✓    | BAMs and BAI indices (per run)                                                            |
| `target_intervals`                          | Array\[File]+    |     ✓    | WES padded targets **or** WGS binned genomic intervals (prefer blacklist-filtered)        |
| `annotated_target_intervals`                | Array\[File]?    |          | Targets with GC, mappability, optional segdups                                            |
| `cnv_panel_of_normals`                      | Array\[File]?    |          | If unavailable for a platform, point to **0-byte** file to fall back to annotated targets |
| `is_paired_end`                             | Array\[Boolean]? |          | Avoids double-counting for short cfDNA inserts                                            |
| `use_sample_for_tCR` / `use_sample_for_aCR` | Array\[Boolean]? |          | Choose samples contributing to tCR/aCR estimation                                         |
| `normal_sample_names`                       | Array\[String]?  |          | First entry used as the matched normal (choose highest coverage if known)                 |

### B) Reference data (`Files`)

| Name                                                      | Type    | Required | Notes                                                                                                       |
| --------------------------------------------------------- | ------- | :------: |-------------------------------------------------------------------------------------------------------------|
| `interval_list` / `interval_lists` / `interval_blacklist` | File(s) |          | Short-variant calling intervals (apply blacklist); if not provided, SNV are called across the whole genome. |
| `ref_fasta` / `_index` / `_dict`                          | File    |     ✓    | Reference genome                                                                                            |
| `force_call_alleles` / `_idx`                             | File?   |          | COSMIC/driver sites to force-call                                                                           |
| `snv_panel_of_normals` / `_idx` (+ v4.1 variants)         | File?   |          | Platform-specific artifact PoN                                                                              |
| `germline_resource` / `_idx` (+ v4.1 variants)            | File?   |          | gnomAD with AF for filtering                                                                                |
| `common_germline_alleles` / `_idx`                        | File?   |          | Common biallelic SNPs for aCR + contamination                                                               |
| `realignment_bwa_mem_index_image`                         | File?   |          | For realignment filter (e.g., hg38)                                                                         |
| `funcotator_transcript_list`                              | File?   |          | Transcript whitelist for annotation                                                                         |
| `funcotator_data_sources_tar_gz`                          | File?   |          | If omitted, downloaded automatically (much slower)                                                          |

### C) Workflow Parameters (`Parameters`)

| Name                                                      | Type  | Required | Notes                                                                                                                         |
| --------------------------------------------------------- |-------| :------: |-------------------------------------------------------------------------------------------------------------------------------|
| `total_mean_read_depth` | Int   |          | Approximate total mean read depth across all samples (used to dynamically scale scattter-gather shard number for SNV calling) |

---
## Outputs

### Per sequencing **run**

| Name                                | Type                 | Meaning                          |
| ----------------------------------- |----------------------| -------------------------------- |
| `callable_loci`                     | Array\[Array\[File]] | Regions passing coverage filters |
| `total_read_counts`                 | Array\[Array\[File]] | Counts in targets (pre-denoise)  |
| `denoised_total_copy_ratios`        | Array\[Array\[File]] | tCR after tangent-normalization  |
| `snppanel_allelic_pileup_summaries` | Array\[Array\[File]] | Allelic counts summaries         |

### Per **sample**

| Name                                           | Type               | Meaning                                                                          |
| ---------------------------------------------- |--------------------|----------------------------------------------------------------------------------|
| `sample_names_ordered`                         | Array\[String]     | Sample order used downstream                                                     |
| `harmonized_callable_loci`                     | Array\[File]       | Union of callable loci per seq run                                               |
| `harmonized_denoised_total_copy_ratios`        | Array\[File]       | Harmonized tCR per sample                                                        |
| `harmonized_snppanel_allelic_pileup_summaries` | Array\[File]       | Harmonized allelic summaries                                                     |
| `contamination_table`                          | Array\[File]       | Per-sample contamination estimates                                               |
| `af_segmentation_table`                        | Array\[File]       | Allelic copy-ratio segmentation (intervals annotated with minor allele fraction) |
| `allelic_pileup_summaries`                     | Array\[File]       | Per-sample allelic pileups                                                       |
| `aggregated_allelic_read_counts`               | Array\[File]       | Per target interval aggregated allelic counts                                    |
| `genotype_error_probabilities`                 | Array\[Float]      | Genotyping error metrics                                                         |
| `af_model_parameters` / `cr_model_parameters`  | Array\[File]       | Model params (AF & copy-ratio)                                                   |
| `called_copy_ratio_segmentation`               | Array\[File]       | Segmentation calls (amp/del)                                                     |
| `cr_plot`                                      | Array\[File]       | Segmentation plot                                                                |
| `acs_copy_ratio_segmentation`                  | Array\[File]       | Allelic copy-ratio segmentation                                                  |
| `acs_copy_ratio_skew`                          | Array\[Float]      | Allelic skew per sample                                                          |
| `annotated_somatic_variants` / `_idx`          | Array\[File]       | Annotated somatic VCF(+ index)                                                   |
| `absolute_*` (incl. `purity`, `ploidy`)        | Array\[File/Float] | ABSOLUTE outputs/metrics                                                         |

### Per **interval shard**

| Name                                                   | Type         | Meaning                                |
| ------------------------------------------------------ | ------------ |----------------------------------------|
| `raw_calls_mutect2_vcf_scattered` / `_idx_scattered`   | Array\[File] | Shard VCF + index                      |
| `raw_mutect2_stats_scattered`                          | Array\[File] | Mutect2 stats per shard                |
| `raw_mutect2_bam_out_scattered` / `_bai_out_scattered` | Array\[File] | BAMs with evidence of called mutations |
| `raw_mutect2_artifact_priors_scattered`                | Array\[File] | Artifact priors                        |

### **Patient** level

| Name                                                         | Type | Meaning                                        |
| ------------------------------------------------------------ | ---- |------------------------------------------------|
| `raw_snv_calls_vcf` / `_idx`                                 | File | Merged raw SNV calls                           |
| `mutect2_stats`                                              | File | Combined stats                                 |
| `orientation_bias`                                           | File | Read-orientation artifacts                     |
| `filtered_vcf` / `_idx`                                      | File | After statistical/hard/realignment filters     |
| `filtering_stats`                                            | File | Filtering summaries                            |
| `somatic_vcf` / `_idx`                                       | File | Final selected somatic calls                   |
| `germline_vcf` / `_idx`                                      | File | Germline calls retained                        |
| `rare_germline_alleles` / `_idx`                             | File | Rare germline sites used for CNV               |
| `somatic_calls_bam` / `_bai`                                 | File | Evidence BAM for selected somatic calls        |
| `gvcf` / `_idx`                                              | File | Genotyping VCF                                 |
| `snp_ref_counts` / `snp_alt_counts` / `snp_other_alt_counts` | File | Allele-specific count tables                   |
| `snp_sample_correlation`                                     | File | Cross-sample variant SNP genotype correlations |
| `modeled_segments`                                           | File | CNV modeled segments                           |
| `phylogic_*` (report/timing/graph/growth)                    | File | PhylogicNDT outputs                            |


All outputs can be provided as `Cache` input in which case the corresponding task to generate that output is skipped if possible.

If `Array[File] Cache.absolute_acr_rdata` is supplied as input, the workflow skips straight to the clonal decomposition step. Use that together with `Array[Int] Cache.absolute_solution` to extract the chosen purity/ploidy mode for each sample and run phylogenetic inference. (Other required cached inputs for this to run successfully are `Array[File] Cache.acs_copy_ratio_segmentation` and `Array[Float] Cache.acs_copy_ratio_skew`.)


---
## Troubleshooting

**A few Mutect2 shards failed while most succeeded**
- Retry problematic shards with more memory: set `Cache.high_mem_shards` to the list of shard IDs and re-run (call-caching will reuse prior successes if enabled).
- If needed, increase `Parameters.mutect2_high_mem_factor` (e.g., 3).
- As a last resort, skip shards with `Cache.skip_shards` (acknowledging lost calls in those regions).


---
## Important Notes:
- Mutect2 `use_linked_de_bruijn_graph`, while increasing sensitivity, has trouble calling variants in complex regions. It is strongly recommended (necessary) to use it together with `recover_all_dangling_branches` (both turned on by default). Pre-calling with Mutect1 also helps.
- Runtime parameters are optimized for implementations on Google Cloud Platform (GCP) and HPC cluster with SLURM scheduler and for applications to whole exome sequencing data.
- For assistance running workflows on GCP or locally, refer to the [GATK tutorial](https://gatk.broadinstitute.org/hc/en-us/articles/360035530952).
- Access necessary reference and resources bundles via the [GATK Resource Bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360036212652).


---
## FAQ

**Q: Can I run tumor-only (no matched normal)?**
A: Yes, Mutect1/2 and CNV tumor-only calling are supported. Expect more rare germline calls in the somatic SNV call set and consider stronger filtering (e.g. GERMQ >= 40) for SNV burden tests. 

**Q: No CNV panel of normals (PoN) for my sequencing platform?**
A: Provide a **0-byte** file for that entry. The workflow will fall back to `annotated_target_intervals` (using gc-content) for denoising. This is still provides decent results for WGS samples, but WES samples tend to have more batch structure. 

**Q: My targeted panels differ across samples—OK?**
A: Harmonization intersects targets across samples and mixed panels reduce effective target space. Prefer consistent panels for multi-sample CNV detection.


---
## Software version requirements :
- **GATK**: Version 4.6.2.0. 
- **Cromwell**: Tested successfully on version 86.
