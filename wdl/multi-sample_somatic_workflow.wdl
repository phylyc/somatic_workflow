version development

import "patient.wdl" as p
import "patient.define.wdl" as p_def
import "patient.out.wdl" as p_out
import "workflow_arguments.wdl" as wfargs
import "workflow_resources.wdl" as wfres
import "runtime_collection.wdl" as rtc
import "coverage_workflow.wdl" as cov
import "cnv_workflow.wdl" as cnv
import "snv_workflow.wdl" as snv
import "clonal_analysis_workflow.wdl" as clone


workflow MultiSampleSomaticWorkflow {
    input {
        # This string is used to label the outputs of the workflow.
        String patient_id
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
        # The target_intervals annotated with gc content, mappability, and segmental duplications.
        Array[File]? annotated_target_intervals
        # If a panel of normals is not available for the sequencing platform of a sample,
        # its corresponding path must point to an empty file (of size 0B). The
        # annotated_target_intervals will instead be used for denoising.
        Array[File]? cnv_panel_of_normals
        # Setting this avoids double counting evidence from paired-end reads. This is
        # particularly important for cell-free DNA samples, where the majority of
        # templates is shorter than twice the read length.
        Array[Boolean]? is_paired_end
        # Whether to use the sample for total copy ratio (tCR) and/or allelic copy ratio
        # (aCR) estimation. If not provided, all samples will be used.
        Array[Boolean]? use_sample_for_tCR  # Boolean inputs, ensure correct format!
        Array[Boolean]? use_sample_for_aCR  # Boolean inputs, ensure correct format!

        # A list of normal sample names. If not provided, all samples will be treated as
        # tumor samples.
        Array[String]? normal_sample_names

        Int scatter_count = 10

        Patient? patient
        WorkflowArguments? args
        WorkflowResources? resources
        RuntimeCollection? runtime_collection
    }

    if (!defined(runtime_collection)) {
        call rtc.DefineRuntimeCollection as RuntimeParameters {
            input:
                num_bams = length(bams),
                scatter_count = scatter_count,
        }
    }
    RuntimeCollection this_runtime_collection = select_first([runtime_collection, RuntimeParameters.rtc])

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

    if (!defined(patient)) {
        call p_def.DefinePatient as Cache {
            input:
                name = patient_id,
                sex = sex,
                sample_names = sample_names,
                bams = bams,
                bais = bais,
                target_intervals = target_intervals,
                annotated_target_intervals = annotated_target_intervals,
                cnv_panel_of_normals = cnv_panel_of_normals,
                is_paired_end = is_paired_end,
                use_for_tCR = use_sample_for_tCR,
                use_for_aCR = use_sample_for_aCR,
                normal_sample_names = normal_sample_names,
                scattered_intervals = this_args.files.scattered_intervals,
                runtime_collection = this_runtime_collection,
        }
    }
    Patient this_patient = select_first([patient, Cache.patient])

    # TODO: add parse_input task to check for validity, then add "after parse_input" to all calls

    call cov.CoverageWorkflow {
        input:
            args = this_args,
            patient = this_patient,
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

    Patient out_patient = ClonalAnalysisWorkflow.updated_patient

    call p_out.Output {
        input:
            patient = out_patient
    }

    output {
        # for each sequencing run:
        # CACHE (as returned by the workflow)
        Array[Array[File]]? callable_loci = Output.callable_loci
        Array[Array[File]]? total_read_counts = Output.total_read_counts
        Array[Array[File]]? denoised_total_copy_ratios = Output.denoised_total_copy_ratios
        Array[Array[File]]? snppanel_allelic_pileup_summaries = Output.snppanel_allelic_pileup_summaries

        # for each sample:
        # CACHE (as returned by the workflow)
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
        Boolean? samples_are_from_same_patient = out_patient.samples_are_from_same_patient
        File? modeled_segments = out_patient.modeled_segments
    }
}
