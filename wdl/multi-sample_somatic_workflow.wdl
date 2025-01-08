version development

import "patient.wdl" as p
import "patient.define.wdl" as p_def
import "workflow_arguments.wdl" as wfargs
import "workflow_resources.wdl" as wfres
import "runtime_collection.wdl" as rtc
import "coverage_workflow.wdl" as cov
import "cnv_workflow.wdl" as cnv
import "snv_workflow.wdl" as snv
import "clonal_analysis_workflow.wdl" as clone


workflow MultiSampleSomaticWorkflow {
    input {
        # todo: add options for inputting cached output

        Patient patient = GetPatient.patient

        # This string is used to label the outputs of the workflow.
        String individual_id
        String? sex
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

        WorkflowArguments? args
        WorkflowResources? resources
        RuntimeCollection? runtime_collection
    }

    if (!defined(runtime_collection)) {
        call rtc.DefineRuntimeCollection as RuntimeParameters {
            input:
                num_bams = length(bams),
                bam_size = ceil(size(bams, "GB") + size(bais, "GB")),
                scatter_count = scatter_count,
        }
    }
    RuntimCollection this_runtime_collection = select_first([runtime_collection, RuntimeParameters.rtc])

    if (!defined(resources)) {
        call wfres.DefineWorkflowResources as Files
    }
    WorkflowResources this_resources = select_first([resources, Files.resources])

    if (!defined(args)) {
        call wfargs.DefineWorkflowArguments as Parameters {
            input:
                scatter_count = scatter_count,
                resources = this_resources,
                runtime_collection = this_runtime_collection,
        }
    }
    WorkflowArguments this_args = select_first([args, Parameters.arguments])

    call p_def.DefinePatient as GetPatient {
        input:
            individual_id = individual_id,
            sex = sex,

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

            runtime_collection = this_runtime_collection,
    }

    call cov.CoverageWorkflow {
        input:
            args = this_args,
            patient = patient,
            runtime_collection = this_runtime_collection,
    }

    call snv.SNVWorkflow {
        input:
            args = this_args,
            patient = CoverageWorkflow.updated_patient,
            runtime_collection = this_runtime_collection,
    }

    call cnv.CNVWorkflow {
        input:
            args = this_args,
            patient = SNVWorkflow.updated_patient,
            runtime_collection = this_runtime_collection,
    }

    call clone.ClonalAnalysisWorkflow {
        input:
            args = this_args,
            patient = CNVWorkflow.updated_patient,
            runtime_collection = this_runtime_collection,
    }

    output {
        Array[File?]? target_read_counts = CoverageWorkflow.target_read_counts
        Array[File?]? denoised_copy_ratios = CoverageWorkflow.denoised_copy_ratios
        Array[File?]? covered_regions_interval_list = CoverageWorkflow.covered_regions_interval_list

        File? unfiltered_vcf = SNVWorkflow.unfiltered_vcf
        File? unfiltered_vcf_idx = SNVWorkflow.unfiltered_vcf_idx
        File? mutect_stats = SNVWorkflow.mutect_stats
        File? somatic_calls_bam = SNVWorkflow.somatic_calls_bam
        File? somatic_calls_bai = SNVWorkflow.somatic_calls_bai
        File? orientation_bias = SNVWorkflow.orientation_bias
        File? filtered_vcf = SNVWorkflow.filtered_vcf
        File? filtered_vcf_idx = SNVWorkflow.filtered_vcf_idx
        File? somatic_vcf = SNVWorkflow.somatic_vcf
        File? somatic_vcf_idx = SNVWorkflow.somatic_vcf_idx
        File? filtering_stats = SNVWorkflow.filtering_stats
        Array[File]? annotated_variants = SNVWorkflow.annotated_variants
        Array[File?]? annotated_variants_idx = SNVWorkflow.annotated_variants_idx
        Array[File?]? somatic_allelic_counts = SNVWorkflow.somatic_allelic_counts

        File? germline_vcf = CNVWorkflow.snppanel_genotyped_vcf
        File? germline_vcf_idx = CNVWorkflow.snppanel_genotyped_vcf_idx
        File? germline_ref_counts = CNVWorkflow.snppanel_ref_counts
        File? germline_alt_counts = CNVWorkflow.snppanel_alt_counts
        File? germline_other_alt_counts = CNVWorkflow.snppanel_other_alt_counts
        File? germline_sample_correlation = CNVWorkflow.snppanel_sample_correlation
        Array[File]? germline_sample_genotype_likelihoods = CNVWorkflow.snppanel_sample_genotype_likelihoods
        Array[File]? contamination_tables = select_first([CNVWorkflow.contamination_tables, CoverageWorkflow.contamination_tables])
        Array[File]? segmentation_tables = select_first([CNVWorkflow.segmentation_tables, CoverageWorkflow.segmentation_tables])
        Array[File]? germline_pileups = select_first([CNVWorkflow.snppanel_pileups, CoverageWorkflow.snppanel_pileups])
        Array[File]? germline_allelic_counts = CNVWorkflow.snppanel_allelic_counts

        File? modeled_segments = CNVWorkflow.modeled_segments
        Array[File]? cr_segmentations = CNVWorkflow.cr_segmentations
        Array[File]? cr_plots = CNVWorkflow.cr_plots
        Array[File]? af_model_parameters = CNVWorkflow.af_model_parameters
        Array[File]? cr_model_parameters = CNVWorkflow.cr_model_parameters

        Array[File]? absolute_acr_plots = ClonalAnalysisWorkflow.absolute_acr_plots
        Array[File]? absolute_acr_rdata = ClonalAnalysisWorkflow.absolute_acr_rdata
#        Array[File]? absolute_tcr_plots = ClonalAnalysisWorkflow.absolute_tcr_plots
#        Array[File]? absolute_tcr_rdata = ClonalAnalysisWorkflow.absolute_tcr_rdata
    }
}
