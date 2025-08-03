# Workflow for Somatic Variant Discovery with GATK4+

This repository contains a WDL (Workflow Description Language) workflow for performing multi-sample somatic variant calling using GATK4. The workflow covers preprocessing, single nucleotide variant (SNV) calling, copy number variation (CNV) calling, and clonal analysis. It is available on [Dockstore](https://dockstore.org/workflows/github.com/phylyc/somatic_workflow/MultiSampleSomaticWorkflow:master) to import into e.g. a [Terra](https://app.terra.bio/) workspace. Alternatively, clone this repository and run the workflow locally using [Cromwell](https://cromwell.readthedocs.io/en/stable/):
```bash
java -jar cromwell.jar run mult-sample_somatic_workflow.wdl --inputs mult-sample_somatic_workflow.inputs.json
```

## Workflow Structure

The workflow is organized into the following main tasks:

### 1. Preprocessing
- **1.1 Preprocess and Split Calling Intervals**: Prepare the intervals for scatter-gather tasks in shards for variant calling.
- **1.2 Define Patients**: Define a patient as a set of samples, and each sample as a set of sequencing runs. Sequencing runs are grouped by sample name.

### 2. Coverage Aggregation
- **2.1 Callable Loci Collection**: Collect genomic regions per sample with sufficient coverage for SNV discovery.
- **2.2 Total Read Count Collection**: Collect total read counts in target intervals and perform denoising of total copy ratios via tangent-normalization.
- **2.3 Allelic Read Count Collection**: Collect allelic read counts at common germline sites (SNP panel: `common_germline_alleles`).
- **2.4 Harmonization of Target Intervals and Allelic Counts**: Harmonize target intervals across samples, subset to intersection; merge read count data from multiple sequencing runs per sample.
- **2.5 Contamination Estimation**: Estimate out-of-patient contamination in each sample.
- **2.6 First-pass copy ratio segmentation**: get a prior single-sample allelic copy ratio segmentation which is used for filtering SNVs (3.3) and for genotyping germline sites (4.1).

### 3. SNV Calling
Tasks involved in the detection and analysis of short nucleotide variations:
- **3.1 Mutect1 Single-sample Calling**: Use Mutect1 for single-sample mutation calling in tumor-normal mode if a matched normal sample is available, otherwise in tumor-only mode.
- **3.2 Mutect2 Multi-sample Calling**: Use Mutect2 for multi-sample mutation calling. Force-call alleles that were called via Mutect1.
- **3.3 Filter Variant Calls**: Annotate and select somatic vs germline vs artifactual variants based on various filters.
  - **3.3a Filter**: Apply statistical filters for sequencing artifacts, germline variants, read orientation bias, and contamination, among others.
  - **3.3b Hard Filter**: Apply hard filters based on base quality, mappability, fragment length, read depth, read orientation quality, position on the read, and population allele frequency.
  - **3.3c Realignment Filter**: Filter based on realignment success (to hg38 or whichever reference is given by `realignment_bwa_mem_index_image`).
- **3.4 Annotate SNVs**: Annotate short nucleotide variants with functional information.
- **3.5 Tumor mutational burden (TMB) estimation**: (coming soon)

### 4. CNV Calling
Tasks involved in the detection and analysis of copy number variations:
- **4.1 SNP Genotyping**: Genotype allelic count data at common (from 2.3) and rare (from 3.3) germline sites using evidence across all samples; harmonize loci across samples.
- **4.2 Multi-sample Segmentation**: Segment denoised total copy ratios and allelic copy ratio across multiple samples.
- **4.3 Per-sample Copy Ratio Inference**: Infer total and allelic copy ratios for each sample.
- **4.4 Per-sample Event Calling**: Call amplifications and deletions for each sample. 
- **4.5 Per-sample Segmentation Plotting**: Plot the segmented denoised copy ratios and allelic copy ratios for each sample.

### 5. Clonal Analysis
- **5.1 ABSOLUTE**: Perform per-sample clonal analysis of the identified somatic variants and estimate tumor purity and ploidy.
- **5.2 ABSOLUTE extraction**: Extract results for one chosen solution (needs manual input).
  - **5.2a**: Rescue dropped segments.
- **5.3 PhylogicNDT**: (Coming soon)

Please remember to always review the intermediate results to ensure that the final results are as expected. Inappropriate filtering or parameter settings can lead to misleading output.


## Expected Inputs

The input arguments are grouped into 
```wdl
MultiSampleSomaticWorkflow # (mandatory + optional)
Cache # (workflow output)
Files # (reference data resources)
Parameters # (workflow parameters)
RuntimeParameters
```
Detailed descriptions of expected inputs (MultiSampleSomaticWorkflow):

```wdl
# This string is used to label the outputs of the workflow.
String patient_id
# Sex genotype: XX or XY.
String? sex
# If defined, all arrays must have the same length. Each entry corresponds to a
# sequencing run, and the same index in each array corresponds to the same run.
# Several sequencing runs from the same physical sample are allowed and will be
# grouped based on the sample name. If the sample name is not provided, it will
# be inferred from the bam.
Array[String]? sample_names
Array[File]+ bams
Array[File]+ bais
# For targeted sequencing, the (possibly padded and ideally mappability blacklist-removed)
# target intervals must be supplied. For whole genome sequencing, the intervals
# are just the chromosomal intervals (ideally mappability blacklist-removed).
Array[File]+ target_intervals
# The target_intervals annotated with gc content, mappability, and segmental duplications. 
Array[File]? annotated_target_intervals
# If a panel of normals is not available for the sequencing platform of a sample,
# its corresponding path must point to an empty file (of size 0B). The 
# annotated_target_intervals will instead be used for denoising of total read counts.
Array[File]? cnv_panel_of_normals
# Setting this avoids double counting evidence from paired-end reads. This is
# particularly important for cell-free DNA samples, where the majority of
# templates is shorter than twice the read length.
Array[Boolean]? is_paired_end
# Whether to use the sample for denoised copy ratio (dCR) and allelic copy ratio 
# (aCR) estimation. If not provided, all samples will be used.
Array[Boolean]? use_sample_for_tCR  # Boolean inputs, ensure correct format!
Array[Boolean]? use_sample_for_aCR  # Boolean inputs, ensure correct format!
    
# A list of normal sample names. If not provided, all samples will be treated as
# tumor samples. The first sample in this list will be used as the "matched" normal
# sample. If known, it is recommended to choose the one with the highest coverage.
Array[String]? normal_sample_names
```
Reference data (Files):
```wdl
# Intervals for short variant calling. The interval_list and interval_lists will
# be combined and interval_blacklist will be subtracted from the result.
File? interval_list
File? interval_blacklist
Array[File]? interval_lists

# Reference fasta, index, and dictionary files.
File ref_fasta
File ref_fasta_index
File ref_dict

# VCF file of sites to force mutation calling at, e.g. COSMIC or known driver mutations.
File? force_call_alleles
File? force_call_alleles_idx
# VCF file of common SNV sequencing artifacts to filter out. This should match the 
# sequencing platform(s) of the samples.
File? snv_panel_of_normals
File? snv_panel_of_normals_idx
# Same VCF files in vcf v4.1 format for Mutect1
File? snv_panel_of_normals_v4_1
File? snv_panel_of_normals_v4_1_idx
# VCF file with AF field for annotating / filtering germline alleles (gnomAD).
File? germline_resource
File? germline_resource_idx
# Same VCF files in vcf v4.1 format for Mutect1
File? germline_resource_v4_1
File? germline_resource_v4_1_idx
# VCF file of common biallelic germline alleles (e.g. population allele frequency > 5%) to
# collect allelic counts at for allelic copy ratio (aCR) and contamination estimation.
File? common_germline_alleles
File? common_germline_alleles_idx
# BWA index image file for realignment task (to hg38).
File? realignment_bwa_mem_index_image
# List of transcript names to use for the Funcotator annotation.
File? funcotator_transcript_list
# Tarball of data sources for Funcotator. If not provided, the tarball will automatically
# be downloaded from the GATK resource bundle, which is much slower.
File? funcotator_data_sources_tar_gz
```

## Expected Final Outputs

These outputs can be used as cached workflow inputs with the same name to skip the corresponding tasks.

```wdl
# for each sequencing run:
# CACHE (as returned by the workflow)
Array[Array[File]]? callable_loci = Output.callable_loci
Array[Array[File]]? total_read_counts = Output.total_read_counts
Array[Array[File]]? denoised_total_copy_ratios = Output.denoised_total_copy_ratios
Array[Array[File]]? snppanel_allelic_pileup_summaries = Output.snppanel_allelic_pileup_summaries

# for each sample:
# CACHE (as returned by the workflow)
Array[String]? sample_names_ordered = Output.sample_names_ordered
Array[File]? harmonized_callable_loci = Output.harmonized_callable_loci
Array[File]? harmonized_denoised_total_copy_ratios = Output.harmonized_denoised_total_copy_ratios
Array[File]? harmonized_snppanel_allelic_pileup_summaries = Output.harmonized_snppanel_allelic_pileup_summaries
Array[File]? contamination_table = Output.contamination_table
Array[File]? af_segmentation_table = Output.af_segmentation_table
Array[File]? allelic_pileup_summaries = Output.allelic_pileup_summaries
Array[File]? aggregated_allelic_read_counts = Output.aggregated_allelic_read_counts
Array[Float]? genotype_error_probabilities = Output.genotype_error_probabilities
Array[File]? af_model_parameters = Output.af_model_parameters
Array[File]? cr_model_parameters = Output.cr_model_parameters
Array[File]? called_copy_ratio_segmentation = Output.called_copy_ratio_segmentation
Array[File]? cr_plot = Output.cr_plot
Array[File]? acs_copy_ratio_segmentation = Output.acs_copy_ratio_segmentation
Array[Float]? acs_copy_ratio_skew = Output.acs_copy_ratio_skew
Array[File]? annotated_somatic_variants = Output.annotated_somatic_variants
Array[File?]? annotated_somatic_variants_idx = Output.annotated_somatic_variants_idx
Array[File]? absolute_acr_rdata = Output.absolute_acr_rdata
Array[File]? absolute_acr_plot = Output.absolute_acr_plot
Array[File]? absolute_snv_maf = Output.absolute_snv_maf
Array[File]? absolute_indel_maf = Output.absolute_indel_maf
Array[Int]? absolute_solution = Output.absolute_solution
Array[File]? absolute_maf = Output.absolute_maf
Array[File]? absolute_segtab = Output.absolute_segtab
Array[File]? absolute_table = Output.absolute_table
Array[Float]? purity = Output.purity
Array[Float]? ploidy = Output.ploidy
Array[Int]? timepoint = Output.timepoint

Array[File]? first_pass_cr_segmentations = CoverageWorkflow.first_pass_cr_segmentations
Array[File]? first_pass_cr_plots = CoverageWorkflow.first_pass_cr_plots
Array[File]? first_pass_af_model_parameters = CoverageWorkflow.first_pass_af_model_parameters
Array[File]? first_pass_cr_model_parameters = CoverageWorkflow.first_pass_cr_model_parameters

# for each interval shard:
# CACHE (as returned by the workflow)
Array[File]? raw_calls_mutect2_vcf_scattered = Output.raw_calls_mutect2_vcf_scattered
Array[File]? raw_calls_mutect2_vcf_idx_scattered = Output.raw_calls_mutect2_vcf_idx_scattered
Array[File]? raw_mutect2_stats_scattered = Output.raw_mutect2_stats_scattered
Array[File]? raw_mutect2_bam_out_scattered = Output.raw_mutect2_bam_out_scattered
Array[File]? raw_mutect2_bai_out_scattered = Output.raw_mutect2_bai_out_scattered
Array[File]? raw_mutect2_artifact_priors_scattered = Output.raw_mutect2_artifact_priors_scattered

# for patient
File? raw_snv_calls_vcf = out_patient.raw_snv_calls_vcf
File? raw_snv_calls_vcf_idx = out_patient.raw_snv_calls_vcf_idx
File? mutect2_stats = out_patient.mutect2_stats
File? orientation_bias = out_patient.orientation_bias
File? filtered_vcf = out_patient.filtered_vcf
File? filtered_vcf_idx = out_patient.filtered_vcf_idx
File? filtering_stats = out_patient.filtering_stats
File? somatic_vcf = out_patient.somatic_vcf
File? somatic_vcf_idx = out_patient.somatic_vcf_idx
File? germline_vcf = out_patient.germline_vcf
File? germline_vcf_idx = out_patient.germline_vcf_idx
File? rare_germline_alleles = out_patient.rare_germline_alleles
File? rare_germline_alleles_idx = out_patient.rare_germline_alleles_idx
File? somatic_calls_bam = out_patient.somatic_calls_bam
File? somatic_calls_bai = out_patient.somatic_calls_bai
File? gvcf = out_patient.gvcf
File? gvcf_idx = out_patient.gvcf_idx
File? snp_ref_counts = out_patient.snp_ref_counts
File? snp_alt_counts = out_patient.snp_alt_counts
File? snp_other_alt_counts = out_patient.snp_other_alt_counts
File? snp_sample_correlation = out_patient.snp_sample_correlation
File? modeled_segments = out_patient.modeled_segments
```


## Important Notes:
- Mutect2 `use_linked_de_bruijn_graph`, while increasing sensitivity, has trouble calling variants in complex regions. It is strongly recommended (necessary) to use it together with `recover_all_dangling_branches` (both turned on by default).
- Runtime parameters are optimized for implementations on Google Cloud Platform (GCP) and HPC cluster with SLURM scheduler and for applications to whole exome sequencing data.
- For assistance running workflows on GCP or locally, refer to the [GATK tutorial](https://gatk.broadinstitute.org/hc/en-us/articles/360035530952).
- Access necessary reference and resources bundles via the [GATK Resource Bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360036212652).

## Software version requirements :
- **GATK**: Version 4.6.2.0. 
- **Cromwell**: Tested successfully on version 86.
