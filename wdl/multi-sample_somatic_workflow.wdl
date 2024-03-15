version development

import "patient.wdl" as p
import "patient.define.wdl" as p_def
import "workflow_arguments.wdl" as wfargs
import "workflow_resources.wdl" as wfres
import "runtime_collection.wdl" as rtc
import "cnv_workflow.wdl" as cnv
import "snv_workflow.wdl" as snv
import "clonal_analysis_workflow.wdl" as clone


workflow MultiSampleSomaticWorkflow {
    input {
        # todo: add options for inputting cached output

        Patient patient = GetPatient.patient

        # This string is used to label the outputs of the workflow.
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
        # are just the chromosomal intervals (ideally blacklist-removed).
        Array[File]+ target_intervals
        Array[File]? annotated_target_intervals
        # If a panel of normals is not available for the sequencing platform of a sample,
        # it corresponding path must point to an empty file. The annotated_target_intervals
        # will instead be used for denoising.
        Array[File]? cnv_panel_of_normals
        Array[Boolean]? is_paired_end       # Boolean inputs, ensure correct format!
        # Whether to use the sample for denoised copy ratio (dCR) and allelic copy ratio
        # (aCR) estimation. If not provided, all samples will be used.
        Array[Boolean]? use_sample_for_dCR  # Boolean inputs, ensure correct format!
        Array[Boolean]? use_sample_for_aCR  # Boolean inputs, ensure correct format!

        # A list of normal sample names. If not provided, all samples will be treated as
        # tumor samples.
        Array[String]? normal_sample_names

        Int scatter_count = 10

        WorkflowArguments args = Parameters.arguments
        WorkflowResources resources = Files.resources
        RuntimeCollection runtime_collection = RuntimeParameters.rtc
    }

    call rtc.DefineRuntimeCollection as RuntimeParameters {
        input:
            num_bams = length(bams),
            scatter_count = scatter_count,
    }

    call wfres.DefineWorkflowResources as Files

    call wfargs.DefineWorkflowArguments as Parameters {
        input:
            scatter_count = scatter_count,
            resources = resources,
            runtime_collection = runtime_collection,
    }

    call p_def.DefinePatient as GetPatient {
        input:
            individual_id = individual_id,

            sample_names = sample_names,
            bams = bams,
            bais = bais,
            target_intervals = target_intervals,
            annotated_target_intervals = annotated_target_intervals,
            cnv_panel_of_normals = cnv_panel_of_normals,
            is_paired_end = is_paired_end,
            use_for_dCR = use_sample_for_dCR,
            use_for_aCR = use_sample_for_aCR,
            normal_sample_names = normal_sample_names,

            runtime_collection = runtime_collection,
    }

    call cnv.CNVWorkflow {
        input:
            args = args,
            patient = patient,
            runtime_collection = runtime_collection,
    }

    call snv.SNVWorkflow {
        input:
            args = args,
            patient = CNVWorkflow.updated_patient,
            runtime_collection = runtime_collection,
    }

    call clone.ClonalAnalysisWorkflow {
        input:
            args = args,
            patient = SNVWorkflow.updated_patient,
            runtime_collection = runtime_collection,
    }

    output {
        File? genotyped_snparray_vcf = CNVWorkflow.genotyped_snparray_vcf
        File? genotyped_snparray_vcf_idx = CNVWorkflow.genotyped_snparray_vcf_idx
        File? snparray_ref_counts = CNVWorkflow.snparray_ref_counts
        File? snparray_alt_counts = CNVWorkflow.snparray_alt_counts
        File? snparray_other_alt_counts = CNVWorkflow.snparray_other_alt_counts
        File? sample_snp_correlation = CNVWorkflow.sample_snp_correlation
        Array[File]? sample_snparray_genotype_likelihoods = CNVWorkflow.sample_snparray_genotype_likelihoods
        Array[File]? snparray_pileups = CNVWorkflow.snparray_pileups
        Array[File]? snparray_allelic_counts = CNVWorkflow.snparray_allelic_counts
        Array[File]? contamination_tables = CNVWorkflow.contamination_tables
        Array[File]? segmentation_tables = CNVWorkflow.segmentation_tables
        Array[File?]? target_read_counts = CNVWorkflow.target_read_counts
        Array[File?]? denoised_copy_ratios = CNVWorkflow.denoised_copy_ratios

        File? modeled_segments = CNVWorkflow.modeled_segments
        Array[File]? filtered_called_copy_ratio_segmentations = CNVWorkflow.filtered_called_copy_ratio_segmentations
        Array[File]? cr_plots = CNVWorkflow.cr_plots
        Array[File]? af_model_parameters = CNVWorkflow.af_model_parameters
        Array[File]? cr_model_parameters = CNVWorkflow.cr_model_parameters

        File? unfiltered_vcf = SNVWorkflow.unfiltered_vcf
        File? unfiltered_vcf_idx = SNVWorkflow.unfiltered_vcf_idx
        File? mutect_stats = SNVWorkflow.mutect_stats
        File? locally_realigned_bam = SNVWorkflow.locally_realigned_bam
        File? locally_realigned_bai = SNVWorkflow.locally_realigned_bai
        File? orientation_bias = SNVWorkflow.orientation_bias
        File? filtered_vcf = SNVWorkflow.filtered_vcf
        File? filtered_vcf_idx = SNVWorkflow.filtered_vcf_idx
        File? somatic_vcf = SNVWorkflow.somatic_vcf
        File? somatic_vcf_idx = SNVWorkflow.somatic_vcf_idx
        File? germline_vcf = SNVWorkflow.germline_vcf
        File? germline_vcf_idx = SNVWorkflow.germline_vcf_idx
        File? filtering_stats = SNVWorkflow.filtering_stats
        Array[File?]? called_somatic_allelic_counts = SNVWorkflow.called_somatic_allelic_counts
        Array[File?]? called_germline_allelic_counts = SNVWorkflow.called_germline_allelic_counts
        Array[File]? annotated_variants = SNVWorkflow.annotated_variants
        Array[File?]? annotated_variants_idx = SNVWorkflow.annotated_variants_idx

#        Array[File]? covered_regions_bed = SNVWorkflow.covered_regions_bed
#        Array[File?]? covered_regions_bam = SNVWorkflow.covered_regions_bam
#        Array[File?]? covered_regions_bai = SNVWorkflow.covered_regions_bai
#        Array[File?]? covered_regions_interval_list = SNVWorkflow.covered_regions_interval_list

        Array[File]? absolute_plots = ClonalAnalysisWorkflow.absolute_plots
        Array[File]? absolute_rdata = ClonalAnalysisWorkflow.absolute_rdata
    }
}
