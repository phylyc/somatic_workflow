# Somatic Workflow for Variant Discovery with GATK4

This repository contains a WDL (Workflow Description Language) workflow for performing multi-sample somatic variant calling using GATK4. The workflow covers preprocessing, copy number variation (CNV) calling, single nucleotide variant (SNV) calling, and clonal analysis. It is available on [Dockstore](https://dockstore.org/workflows/github.com/phylyc/somatic_workflow/MultiSampleSomaticWorkflow:master) to import into e.g. a [Terra](https://app.terra.bio/) workspace. Alternatively, clone this repository and run the workflow locally using [Cromwell](https://cromwell.readthedocs.io/en/stable/):
```bash
java -jar cromwell.jar run mult-sample_somatic_workflow.wdl --inputs mult-sample_somatic_workflow.inputs.json
```

## Workflow Structure

The workflow is organized into the following main tasks:

### 1. Preprocessing
- **1.1 Define Patients**: Define a patient as a set of samples, and each sample as a set of sequencing runs. Sequencing runs are grouped by sample name.
- **1.2 Preprocess and Split Calling Intervals**: Prepare the intervals for scatter-gather tasks for variant calling, realignment, and annotation.

### 2. CNV Calling
Tasks involved in the detection and analysis of copy number variations:
- **2.1 Read Count Collection**: Collect read counts in target intervals and perform denoising of total copy ratios per target interval.
- **2.2 Allelic Count Collection**: Collect allelic counts at common germline sites (SNP panel).
- **2.3 Harmonization of Target Intervals**: Ensure consistency of target intervals across samples, removing intervals not present in all samples; merge data from multiple sequencing runs per sample.
- **2.4 Genotyping and Contamination Estimation**: Genotype allelic count data at common germline sites, harmonize allelic loci across samples, and estimate out-of-patient contamination in each sample.
- **2.5 Multi-sample Segmentation**: Segment denoised copy ratios and allelic counts across multiple samples.
- **2.6 Per-sample Copy Ratio Inference**: Infer copy ratios for each sample.
- **2.7 Per-sample Event Calling**: Call amplifications and deletions for each sample. 
- **2.8 Per-sample Germline Filtering**: Filter segments called as amp/del in the matched normal sample from the tumor segmentations.
- **2.9 Per-sample Segmentation Plotting**: Plot the segmented denoised copy ratios and allelic copy ratios for each sample.

### 3. SNV Calling
Tasks involved in the detection and analysis of short nucleotide variations:
- **3.1 Mutect2 Multi-sample Calling**: Use Mutect2 for multi-sample mutation calling.
- **3.2 Filter Mutect2 Calls**: Apply filters for germline variants, read orientation bias, and contamination (from CNV workflow), among others.
- **3.3 Hard Filtering**: Apply hard filters based on base quality, mappability, fragment length, read depth, read orientation quality, and position on the read.
- **3.4 Realignment Filter**: Filter based on realignment success (to hg38 or whichever reference is given by `realignment_bwa_mem_index_image`).
- **3.5 Annotate SNVs**: Annotate short nucleotide variants with functional information.

### 4. Clonal Analysis
- **4.1 ABSOLUTE**: Perform per-sample clonal analysis of the identified variants.

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
# VCF file of germline alleles to filter out (gnomAD).
File? germline_resource
File? germline_resource_idx
# VCF file of common biallelic germline alleles (population allele frequency > 5%) to
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

## Expected Outputs
```wdl
# VCF file with genotype information for sites of the SNP panel (common germline alleles)
File? genotyped_snparray_vcf
File? genotyped_snparray_vcf_idx
# Count matrices of allelic counts at sites of the SNP panel (can be used for AllelicCapSeg)
File? snparray_ref_counts
File? snparray_alt_counts
File? snparray_other_alt_counts
# Sample cross-correlation table of variant genotypes. Low correlation indicates
# high contamination or sample swaps.
File? sample_snp_correlation
# Pileup tables at sites of the SNP panel with genotype likelihoods
Array[File]? sample_snparray_genotype_likelihoods
# Pileup tables at sites of the SNP panel
Array[File]? snparray_pileups
# Allelic count format of pileup tables at sites of the SNP panel
Array[File]? snparray_allelic_counts
# Contamination estimates for each sample
Array[File]? contamination_tables
# Coarse allelic copy ratio segmentations for each sample
Array[File]? segmentation_tables
# Read counts in supplied target intervals for each sequencing run (!)
Array[File?]? target_read_counts
# Harmonized denoised copy ratios for each sample
Array[File?]? denoised_copy_ratios

# Multi-sample segmentation intervals
File? modeled_segments
# Denoised copy ratio and allelic copy ratio segmentations per sample
Array[File]? cr_segments
# -0+ called segments
Array[File]? called_cr_segments
# Denoised copy ratio and allelic copy ratio segmentation model plots
Array[File]? cr_plots
# Model parameters from GATK's ModelSegments 
Array[File]? af_model_parameters
Array[File]? cr_model_parameters

# Mutect2 output of called variants
File? unfiltered_vcf
File? unfiltered_vcf_idx
File? mutect_stats
File? locally_realigned_bam
File? locally_realigned_bai
# Learned read orientation bias for filtering FFPE or OxOG artifacts
File? orientation_bias
# filtered output of called variants (+realignment & hard filters)
File? filtered_vcf
File? filtered_vcf_idx
# Filtered somatic variants
File? somatic_vcf
File? somatic_vcf_idx
# Filtered germline variants (if keep_germline=true)
File? germline_vcf
File? germline_vcf_idx
File? filtering_stats
# Pileup tables of allelic counts at called somatic sites
Array[File?]? called_somatic_allelic_counts
# Pileup tables of allelic counts at called germline sites
Array[File?]? called_germline_allelic_counts
# Functionally annotated variants (Funcotator) in MAF or VCF format per sample
Array[File]? annotated_variants
# Its index if in VCF format.
Array[File?]? annotated_variants_idx

# ABSOLUTE output files per sample. Since ABSOLUTE provides many solutions, further 
# review of this output is necessary.
Array[File]? absolute_plots
Array[File]? absolute_rdata
```

## Important Notes :
- It is recommended to move all output files to a separate bucket or directory after workflow completion as intermediate files can be large and unnecessarily consume storage ($$$). 
- Allelic counts are collected via GetPileupSummaries. This tool ignores certain MNVs and INDELs and should eventually be replaced by [CollectByBaseCounts](https://github.com/broadinstitute/gatk/pull/6545). The GenotypeSNPArray workflow corrects some of its shortcomings.
- As of GATK v4.3.0.0 force-calling alleles leads to mis-classification of filtered variants in the same way as [the flag `--genotype-germline-sites` does](https://github.com/broadinstitute/gatk/issues/7391). It is therefore recommended to use `keep_germline=false` (default) until this has been fixed. Ideally, run Mutect2 also only the force-calling sites and merge the results with the main call set.
- `use_linked_de_bruijn_graph`, while increasing sensitivity, has trouble calling variants in complex regions. Strongly recommended to use with `recover_all_dangling_branches` (both turned on by default). Ideally, run with and without these options and use the joint call set.
- Runtime parameters are optimized for implementations on Google Cloud Platform (GCP) and HPC cluster with SLURM scheduler.
- For assistance running workflows on GCP or locally, refer to the [GATK tutorial](https://gatk.broadinstitute.org/hc/en-us/articles/360035530952).
- Access necessary reference and resources bundles via the [GATK Resource Bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360036212652).

## Software version requirements :
- **GATK**: Version 4.3.0.0. 
- **Cromwell**: Tested successfully on version 86.
