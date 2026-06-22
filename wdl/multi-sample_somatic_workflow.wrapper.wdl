version development

# Curated entry-point wrapper (hoisted form).
#
# Terra lists workflow inputs required-first, then alphabetically by the
# `Workflow.callAlias.var` key, and womtool-85 exposes only ONE level of call
# nesting. This wrapper exploits both: it owns the four setup calls directly (so
# their knobs are the only nested inputs Terra shows) and delegates the actual
# pipeline to MultiSampleSomaticWorkflow via `zzz_run` (whose internals are two
# levels down and therefore hidden). The `z1_`..`z4_` alias prefixes sort the
# setup groups after the lowercase top-level inputs and in a chosen order:
#
#   1. MultiSampleSomaticWorkflowWrapper.*   curated user inputs
#   2. z1_Files.*                            WorkflowResources (reference files)
#   3. z2_Cache.*                            DefinePatient cache-injection fields
#   4. z3_Parameters.*                       WorkflowArguments (run_*, scatter, ...)
#   5. z4_RuntimeParameters.*                RuntimeCollection (mem_*, time_*, ...)
#   (no group 6: the implementation's own inputs are bound here, so they are hidden)
#
# NOTE: while MultiSampleSomaticWorkflow is left unchanged, its identical setup
# block (Define*/reorder/DefinePatient) is duplicated below and runs here instead
# of inside the implementation (which receives the resulting structs and skips it).
# The implementation still runs its own ParseInput gate on the raw inputs. When the
# wrapper is promoted to the sole entry point, that setup should be removed from the
# implementation so it is defined in exactly one place.

import "multi-sample_somatic_workflow.wdl" as impl
import "patient.wdl" as p
import "patient.define.wdl" as p_def
import "workflow_arguments.wdl" as wfargs
import "workflow_resources.wdl" as wfres
import "runtime_collection.wdl" as rtc
import "tasks.wdl"


workflow MultiSampleSomaticWorkflowWrapper {
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
        # Ordinal timepoints of sample collection for phylogenetic inference.
        Array[Int]? timepoints
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

        Patient? input_patient
        WorkflowArguments? input_args
        WorkflowResources? input_resources
        RuntimeCollection? input_runtime_collection
    }

    if (!defined(input_resources)) {
        call wfres.DefineWorkflowResources as z1_Files
    }
    WorkflowResources resources = select_first([input_resources, z1_Files.resources])

    if (!defined(input_runtime_collection)) {
        call rtc.DefineRuntimeCollection as z4_RuntimeParameters {
            input:
                num_bams = length(bams),
        }
    }
    RuntimeCollection runtime_collection = select_first([input_runtime_collection, z4_RuntimeParameters.rtc])

    if (!defined(input_args)) {
        call wfargs.DefineWorkflowArguments as z3_Parameters {
            input:
                resources = resources,
                runtime_collection = runtime_collection,
        }
    }
    WorkflowArguments args = select_first([input_args, z3_Parameters.arguments])

    # Pre-flight gate: validate the raw inputs before any expensive task. ParseInput
    # exits non-zero (failing the workflow) on a mis-shaped input or an unmet implicit
    # requirement. Each downstream work block is guarded by `&& ParseInput.validated`,
    # which makes the block depend on ParseInput so nothing heavy starts until it passes.
    call tasks.ParseInput {
        input:
            script = args.script_validate_inputs,
            n_bams = length(bams),
            n_bais = length(bais),
            sample_names = sample_names,
            normal_sample_names = normal_sample_names,
            is_paired_end = is_paired_end,
            use_for_tCR = use_sample_for_tCR,
            use_for_aCR = use_sample_for_aCR,
            timepoints = timepoints,
            target_intervals = target_intervals,
            annotated_target_intervals = annotated_target_intervals,
            cnv_panel_of_normals = cnv_panel_of_normals,
            has_common_germline_alleles = defined(args.files.common_germline_alleles),
            has_common_germline_alleles_idx = defined(args.files.common_germline_alleles_idx),
            has_realignment_image = defined(args.files.realignment_bwa_mem_index_image),
            has_germline_resource = defined(args.files.germline_resource),
            has_snv_panel_of_normals = defined(args.files.snv_panel_of_normals),
            has_germline_resource_v4_1 = defined(args.files.germline_resource_v4_1),
            has_snv_panel_of_normals_v4_1 = defined(args.files.snv_panel_of_normals_v4_1),
            has_funcotator_sources = defined(args.files.funcotator_data_sources_tar_gz),
            has_sex = defined(sex),
            run_collect_total_read_counts = args.run_collect_total_read_counts,
            run_model_segments = args.run_model_segments,
            run_collect_allelic_read_counts = args.run_collect_allelic_read_counts,
            run_contamination_model = args.run_contamination_model,
            run_realignment_filter = args.run_realignment_filter,
            run_variant_calling_mutect1 = args.run_variant_calling_mutect1,
            run_variant_calling = args.run_variant_calling,
            run_variant_annotation = args.run_variant_annotation,
            run_clonal_decomposition = args.run_clonal_decomposition,
            deep = args.run_input_validation_deep,
            bams = if args.run_input_validation_deep then bams else [],
            ref_dict = args.files.ref_dict,
            runtime_params = runtime_collection.parse_input,
    }

    if (!defined(input_patient) && ParseInput.validated) {
        # Only necessary if bams contig order does not match reference contig
        if (args.run_reorder_bam_contigs) {
            scatter (pair in zip(bams, bais)) {
                call tasks.ReorderSam as ReorderSam {
                    input:
                        ref_fasta = args.files.ref_fasta,
                        ref_fasta_index = args.files.ref_fasta_index,
                        ref_dict = args.files.ref_dict,
                        bam = pair.left,
                        bai = pair.right,
                        runtime_params = runtime_collection.reorder_sam
                }
            }
            Array[File] reordered_bams = select_all(ReorderSam.reordered_bam)
            Array[File] reordered_bais = select_all(ReorderSam.reordered_bai)
        }
        call p_def.DefinePatient as z2_Cache {
            input:
                name = patient_id,
                sex = sex,
                sample_names = sample_names,
                timepoints = timepoints,
                bams = select_first([reordered_bams, bams]),
                bais = select_first([reordered_bais, bais]),
                target_intervals = target_intervals,
                annotated_target_intervals = annotated_target_intervals,
                cnv_panel_of_normals = cnv_panel_of_normals,
                is_paired_end = is_paired_end,
                use_for_tCR = use_sample_for_tCR,
                use_for_aCR = use_sample_for_aCR,
                normal_sample_names = normal_sample_names,
                scattered_intervals_for_variant_calling = args.files.scattered_intervals_for_variant_calling_m2,
                runtime_collection = runtime_collection,
        }
    }

    Patient patient = select_first([input_patient, z2_Cache.patient])

    # --- Delegate the pipeline; its internals are bound here, so Terra hides them.
    call impl.MultiSampleSomaticWorkflow as zzz_run {
        input:
            ### REMOVE
            patient_id = patient_id,
            sex = sex,
            sample_names = sample_names,
            timepoints = timepoints,
            bams = bams,
            bais = bais,
            target_intervals = target_intervals,
            annotated_target_intervals = annotated_target_intervals,
            cnv_panel_of_normals = cnv_panel_of_normals,
            is_paired_end = is_paired_end,
            use_sample_for_tCR = use_sample_for_tCR,
            use_sample_for_aCR = use_sample_for_aCR,
            normal_sample_names = normal_sample_names,
            ### END REMOVE
            input_patient = patient,
            input_args = args,
            input_resources = resources,
            input_runtime_collection = runtime_collection,
    }

    output {
        Array[Array[File]]? callable_loci = zzz_run.callable_loci
        Array[Array[File]]? total_read_counts = zzz_run.total_read_counts
        Array[Array[File]]? denoised_total_copy_ratios = zzz_run.denoised_total_copy_ratios
        Array[Array[File]]? snppanel_allelic_pileup_summaries = zzz_run.snppanel_allelic_pileup_summaries

        Array[String]? sample_names_ordered = zzz_run.sample_names_ordered
        Array[File]? harmonized_callable_loci = zzz_run.harmonized_callable_loci
        Array[File]? harmonized_denoised_total_copy_ratios = zzz_run.harmonized_denoised_total_copy_ratios
        Array[File]? harmonized_snppanel_allelic_pileup_summaries = zzz_run.harmonized_snppanel_allelic_pileup_summaries
        Array[File]? contamination_table = zzz_run.contamination_table
        Array[File]? af_segmentation_table = zzz_run.af_segmentation_table
        Array[File]? allelic_pileup_summaries = zzz_run.allelic_pileup_summaries
        Array[File]? aggregated_allelic_read_counts = zzz_run.aggregated_allelic_read_counts
        Array[Float]? genotype_error_probabilities = zzz_run.genotype_error_probabilities
        Array[File]? af_model_parameters = zzz_run.af_model_parameters
        Array[File]? cr_model_parameters = zzz_run.cr_model_parameters
        Array[File]? called_copy_ratio_segmentation = zzz_run.called_copy_ratio_segmentation
        Array[File]? cr_plot = zzz_run.cr_plot
        Array[File]? annotated_somatic_variants = zzz_run.annotated_somatic_variants
        Array[File?]? annotated_somatic_variants_idx = zzz_run.annotated_somatic_variants_idx
        Array[File]? acs_copy_ratio_segmentation = zzz_run.acs_copy_ratio_segmentation
        Array[Float]? acs_copy_ratio_skew = zzz_run.acs_copy_ratio_skew

        Array[File]? absolute_snv_maf = zzz_run.absolute_snv_maf
        Array[File]? absolute_indel_maf = zzz_run.absolute_indel_maf

        Array[File]? absolute_acr_rdata = zzz_run.absolute_acr_rdata
        Array[File]? absolute_acr_plot = zzz_run.absolute_acr_plot
        Array[Int]? absolute_acr_solution = zzz_run.absolute_acr_solution
        Array[File]? absolute_acr_maf = zzz_run.absolute_acr_maf
        Array[File]? absolute_acr_segtab = zzz_run.absolute_acr_segtab
        Array[File]? absolute_acr_segtab_igv = zzz_run.absolute_acr_segtab_igv
        Array[File]? absolute_acr_table = zzz_run.absolute_acr_table
        Array[Float]? absolute_acr_purity = zzz_run.absolute_acr_purity
        Array[Float]? absolute_acr_ploidy = zzz_run.absolute_acr_ploidy

        Array[File]? absolute_tcr_rdata = zzz_run.absolute_tcr_rdata
        Array[File]? absolute_tcr_plot = zzz_run.absolute_tcr_plot
        Array[Int]? absolute_tcr_solution = zzz_run.absolute_tcr_solution
        Array[File]? absolute_tcr_maf = zzz_run.absolute_tcr_maf
        Array[File]? absolute_tcr_segtab = zzz_run.absolute_tcr_segtab
        Array[File]? absolute_tcr_segtab_igv = zzz_run.absolute_tcr_segtab_igv
        Array[File]? absolute_tcr_table = zzz_run.absolute_tcr_table
        Array[Float]? absolute_tcr_purity = zzz_run.absolute_tcr_purity
        Array[Float]? absolute_tcr_ploidy = zzz_run.absolute_tcr_ploidy

        Array[Float]? user_purity = zzz_run.user_purity
        Array[Float]? user_ploidy = zzz_run.user_ploidy

        Array[File]? first_pass_cr_segmentations = zzz_run.first_pass_cr_segmentations
        Array[File]? first_pass_cr_plots = zzz_run.first_pass_cr_plots
        Array[File]? first_pass_af_model_parameters = zzz_run.first_pass_af_model_parameters
        Array[File]? first_pass_cr_model_parameters = zzz_run.first_pass_cr_model_parameters

        Array[File]? raw_calls_mutect2_vcf_scattered = zzz_run.raw_calls_mutect2_vcf_scattered
        Array[File]? raw_calls_mutect2_vcf_idx_scattered = zzz_run.raw_calls_mutect2_vcf_idx_scattered
        Array[File]? raw_mutect2_stats_scattered = zzz_run.raw_mutect2_stats_scattered
        Array[File]? raw_mutect2_bam_out_scattered = zzz_run.raw_mutect2_bam_out_scattered
        Array[File]? raw_mutect2_bai_out_scattered = zzz_run.raw_mutect2_bai_out_scattered
        Array[File]? raw_mutect2_artifact_priors_scattered = zzz_run.raw_mutect2_artifact_priors_scattered

        File? raw_snv_calls_vcf = zzz_run.raw_snv_calls_vcf
        File? raw_snv_calls_vcf_idx = zzz_run.raw_snv_calls_vcf_idx
        File? mutect2_stats = zzz_run.mutect2_stats
        File? orientation_bias = zzz_run.orientation_bias
        File? filtered_vcf = zzz_run.filtered_vcf
        File? filtered_vcf_idx = zzz_run.filtered_vcf_idx
        File? filtering_stats = zzz_run.filtering_stats
        File? somatic_vcf = zzz_run.somatic_vcf
        File? somatic_vcf_idx = zzz_run.somatic_vcf_idx
        Int? num_somatic_variants = zzz_run.num_somatic_variants
        File? germline_vcf = zzz_run.germline_vcf
        File? germline_vcf_idx = zzz_run.germline_vcf_idx
        Int? num_germline_variants = zzz_run.num_germline_variants
        File? rare_germline_alleles = zzz_run.rare_germline_alleles
        File? rare_germline_alleles_idx = zzz_run.rare_germline_alleles_idx
        File? somatic_calls_bam = zzz_run.somatic_calls_bam
        File? somatic_calls_bai = zzz_run.somatic_calls_bai
        File? gvcf = zzz_run.gvcf
        File? gvcf_idx = zzz_run.gvcf_idx
        File? snp_ref_counts = zzz_run.snp_ref_counts
        File? snp_alt_counts = zzz_run.snp_alt_counts
        File? snp_other_alt_counts = zzz_run.snp_other_alt_counts
        File? snp_sample_correlation = zzz_run.snp_sample_correlation
        Float? snp_sample_correlation_min = zzz_run.snp_sample_correlation_min
        File? modeled_segments = zzz_run.modeled_segments
        File? phylogic_sif_file = zzz_run.phylogic_sif_file
        File? phylogic_report = zzz_run.phylogic_report
        File? phylogic_ccfs_cnvs = zzz_run.phylogic_ccfs_cnvs
        File? phylogic_ccfs_snvs = zzz_run.phylogic_ccfs_snvs
        File? phylogic_constrained_ccf = zzz_run.phylogic_constrained_ccf
        File? phylogic_cluster_ccfs = zzz_run.phylogic_cluster_ccfs
        File? phylogic_build_tree_posteriors = zzz_run.phylogic_build_tree_posteriors
        File? phylogic_growth_rates = zzz_run.phylogic_growth_rates
        File? phylogic_growth_rate_plot = zzz_run.phylogic_growth_rate_plot
        File? phylogic_timing_report = zzz_run.phylogic_timing_report
        File? phylogic_timing_wgd_supporting_events = zzz_run.phylogic_timing_wgd_supporting_events
        File? phylogic_timing_graph = zzz_run.phylogic_timing_graph
        File? phylogic_timing_comparison = zzz_run.phylogic_timing_comparison
        File? phylogic_timing_table = zzz_run.phylogic_timing_table

        Patient output_patient = zzz_run.output_patient
        WorkflowArguments output_args = zzz_run.output_args
        WorkflowResources output_resources = zzz_run.output_resources
        RuntimeCollection output_runtime_collection = zzz_run.output_runtime_collection
    }
}
