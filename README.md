# Workflow for Somatic Variant Discovery with GATK4+

This repository contains a WDL (Workflow Description Language) workflow for performing multi-sample somatic variant calling using GATK4. The workflow covers preprocessing, single nucleotide variant (SNV) calling, copy number variation (CNV) calling, and clonal analysis. It is available on [Dockstore](https://dockstore.org/workflows/github.com/phylyc/somatic_workflow/MultiSampleSomaticWorkflow:master) to import into e.g. a [Terra](https://app.terra.bio/) workspace. Alternatively, clone this repository and run the workflow locally using [Cromwell](https://cromwell.readthedocs.io/en/stable/):
```bash
java -jar cromwell.jar run mult-sample_somatic_workflow.wdl --inputs mult-sample_somatic_workflow.inputs.json
```

## Workflow Structure

The workflow is organized into the following main tasks:

### 1. Preprocessing
- **1.1 Define Patients**: Define a patient as a set of samples, and each sample as a set of sequencing runs. Sequencing runs are grouped by sample name.
- **1.2 Preprocess and Split Calling Intervals**: Prepare the intervals for scatter-gather tasks for variant calling, realignment, and annotation.

### 2. Coverage Aggregation
- **2.1 Total Read Count Collection**: Collect total read counts in target intervals and perform denoising of total copy ratios via tangent-normalization.
- **2.2 Allelic Read Count Collection**: Collect allelic read counts at common germline sites (SNP panel: `common_germline_alleles`).
- **2.3 Harmonization of Target Intervals and Allelic Counts**: Harmonize target intervals across samples, subset to intersection; merge read count data from multiple sequencing runs per sample.
- **2.4 Contamination Estimation**: Estimate out-of-patient contamination in each sample.
- **2.5 First-pass copy ratio segmentation**: Same as steps 4.2-4.6 to get a prior allelic copy ratio segmentation which is used for filtering SNVs (3.2) and for genotyping germline sites (4.1).

### 3. SNV Calling
Tasks involved in the detection and analysis of short nucleotide variations:
- **3.1 Mutect2 Multi-sample Calling**: Use Mutect2 for multi-sample mutation calling.
- **3.2 Filter Variant Calls**: Apply statistical filters for sequencing artifacts, germline variants, read orientation bias, and contamination, among others.
- **3.3 Hard Filtering**: Apply hard filters based on base quality, mappability, fragment length, read depth, read orientation quality, position on the read, and population allele frequency.
- **3.4 Realignment Filter**: Filter based on realignment success (to hg38 or whichever reference is given by `realignment_bwa_mem_index_image`).
- **3.5 Annotate SNVs**: Annotate short nucleotide variants with functional information.

### 4. CNV Calling
Tasks involved in the detection and analysis of copy number variations:
- **4.1 SNP Genotyping**: Genotype allelic count data at common (from 2.3) and rare (from 3.2) germline sites using evidence across all samples; harmonize loci across samples.
- **4.2 Multi-sample Segmentation**: Segment denoised total copy ratios and allelic copy ratio across multiple samples.
- **4.3 Per-sample Copy Ratio Inference**: Infer copy ratios for each sample.
- **4.4 Per-sample Event Calling**: Call amplifications and deletions for each sample. 
- **4.5 Per-sample Germline Filtering**: Filter segments called as amp/del in the matched normal sample from the tumor segmentations.
- **4.6 Per-sample Segmentation Plotting**: Plot the segmented denoised copy ratios and allelic copy ratios for each sample.

### 5. Clonal Analysis
- **5.1 ABSOLUTE**: Perform per-sample clonal analysis of the identified variants and estimate tumor purity and ploidy.

Please remember to always review the intermediate results to ensure that the final results are as expected. Inappropriate filtering or parameter settings can lead to misleading output.

## Expected Inputs

Detailed descriptions of expected inputs:

```wdl
# This string is used to label the patient-specific outputs of the workflow.
String individual_id
# If defined, all arrays must have the same length. Each entry corresponds to a 
# sequencing run, and the same index in each array corresponds to the same run.
# Several sequencing runs from the same physical sample are allowed and will be 
# grouped based on the sample name. If the sample name is not provided, it will
# be inferred from the bam.
Array[String]? sample_names
Array[File]+ bams
Array[File]+ bais
# For targeted sequencing, the (possibly padded and ideally blacklist-removed) 
# target intervals must be supplied. For whole genome sequencing, the intervals 
# are just the binned chromosomal intervals (ideally blacklist-removed). 
Array[File]+ target_intervals
# The target_intervals annotated with gc content, mappability, and segmental duplications. 
Array[File]? annotated_target_intervals
# If a panel of normals is not available for the sequencing platform of a sample,
# its corresponding path must point to an empty file (of size 0B). The 
# annotated_target_intervals will instead be used for denoising.
Array[File]? cnv_panel_of_normals
# Paired-end reads that overlap at some sites of interest lead to double counting. 
# This process in general is more of an issue in cell-free DNA where the vast majority 
# of templates are ~166bp long, which is shorter than twice the read length. The result 
# is that many bases on a given template are reported twice, once from each paired-end 
# read. This is mitigated by setting the respective `is_paired_end` entry to `true`.
Array[Boolean]? is_paired_end       # Boolean inputs, ensure correct format!
# Whether to use the sample for denoised copy ratio (dCR) and allelic copy ratio 
# (aCR) estimation. If not provided, all samples will be used.
Array[Boolean]? use_sample_for_dCR  # Boolean inputs, ensure correct format!
Array[Boolean]? use_sample_for_aCR  # Boolean inputs, ensure correct format!
    
# A list of normal sample names. If not provided, all samples will be treated as
# tumor samples. The first sample in this list will be used as the "matched" normal
# sample. If known, it is recommended to choose the one with the highest coverage.
Array[String]? normal_sample_names
```
Reference data:
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
# VCF file of common sequencing artifacts to filter out. This should match the 
# sequencing platform(s) of the samples.
File? snv_panel_of_normals
File? snv_panel_of_normals_idx
# VCF file of annotating / filtering germline alleles (gnomAD).
File? germline_resource
File? germline_resource_idx
# VCF file of common biallelic germline alleles (e.g. population allele frequency > 5%) to
# collect allelic counts at for first-pass allelic copy ratio (aCR) and contamination estimation.
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
```wdl
### COVERAGE WORKFLOW:
# Read count summarization of bam files.
# Total read counts in supplied target intervals for each sequencing run (!)
Array[File?]? target_read_counts
# Harmonized denoised total copy ratios for each sample
Array[File?]? denoised_copy_ratios
# Allelic read count pileup tables at sites of common and rare germline sites
Array[File]? germline_pileups
# Interval list of bases covered by at least 4 reads
Array[File?]? covered_regions_interval_list
# Contamination estimates for each sample
Array[File]? contamination_tables
    
# FIRST-PASS copy ratio segmentations; see below for final segmentations
# Multi-sample segmentation intervals
File? first_pass_modeled_segments
# Denoised copy ratio and allelic copy ratio segmentations per sample
Array[File]? first_pass_cr_segmentations
# Denoised copy ratio and allelic copy ratio segmentation model plots
Array[File]? first_pass_cr_plots
# Model parameters from GATK's ModelSegments 
Array[File]? first_pass_af_model_parameters
Array[File]? first_pass_cr_model_parameters

### SNV WORKFLOW
# Raw Mutect2 output of called variants
File? unfiltered_vcf
File? unfiltered_vcf_idx
File? mutect_stats
# Learned read orientation bias for filtering FFPE or OxOG artifacts
File? orientation_bias
# unfiltered_vcf with filter annotations for several artifacts, germline, and PASS (somatic)
File? filtered_vcf
File? filtered_vcf_idx
File? filtering_stats
# Final, selected multi-sample somatic vcf
File? somatic_vcf
File? somatic_vcf_idx
# Locally realigned multi-sample bam containing only reads covering the filtered, selected somatic calls
File? somatic_calls_bam
File? somatic_calls_bai
# MAF with functional annotations from Funcotator 
Array[File]? annotated_variants
Array[File?]? annotated_variants_idx

### CNV WORKFLOW
# gVCF file with genotype information for sites from the SNP panel (common germline alleles) 
# and rare germline alleles that were called by the SNV workflow
File? germline_vcf
File? germline_vcf_idx
# Allelic read count matrices at sites of the germline_vcf
File? germline_ref_counts
File? germline_alt_counts
File? germline_other_alt_counts
# Sample cross-correlation table of variant genotypes. Low correlation indicates patient-mismatch or sample swaps. 
File? germline_sample_correlation
# Allelic read count pileup tables at sites of the germline_vcf with genotype likelihoods
Array[File]? germline_sample_genotype_likelihoods
# Germline pileups in allelic count format
Array[File]? germline_allelic_counts

# FINAL copy ratio segmentations
# Multi-sample segmentation intervals
File? modeled_segments
# Denoised copy ratio and allelic copy ratio segmentations per sample
Array[File]? cr_segmentations
# Denoised copy ratio and allelic copy ratio segmentation model plots
Array[File]? cr_plots
# Model parameters from GATK's ModelSegments 
Array[File]? af_model_parameters
Array[File]? cr_model_parameters
# Allelic read count tables as used by ModelSegments to infer allelic copy ratio.
# The genomic loci are moved to match the target intervals, and read counts are 
# aggregated across ach interval (haplotype-aware).
Array[File]? segmentation_allelic_counts

# ABSOLUTE output files per sample. Manual review necessary for downstream analysis!
Array[File]? absolute_acr_plots
Array[File]? absolute_acr_rdata
```

## Important Notes:
- Up until GATK v4.5.0.0, force-calling alleles led to mis-classification of filtered variants in the same way as [the flag `--genotype-germline-sites` does](https://github.com/broadinstitute/gatk/issues/7391). Use GATK v4.6.0.0 or higher.
- `use_linked_de_bruijn_graph`, while increasing sensitivity, has trouble calling variants in complex regions. It is strongly recommended (necessary) to use it together with `recover_all_dangling_branches` (both turned on by default).
- Runtime parameters are optimized for implementations on Google Cloud Platform (GCP) and HPC cluster with SLURM scheduler.
- For assistance running workflows on GCP or locally, refer to the [GATK tutorial](https://gatk.broadinstitute.org/hc/en-us/articles/360035530952).
- Access necessary reference and resources bundles via the [GATK Resource Bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360036212652).

## Software version requirements :
- **GATK**: Version 4.6.1.0. 
- **Cromwell**: Tested successfully on version 86.
