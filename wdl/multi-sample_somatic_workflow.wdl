version development

# Public entry point. This workflow owns all preprocessing — resource/argument/
# runtime construction, input validation, and patient assembly — then delegates the
# compute pipeline to MultiSampleSomaticWorkflowRun. Splitting it this way keeps the
# user-facing input surface clean on Terra.
#
# Terra lists workflow inputs required-first, then alphabetically by the
# `Workflow.callAlias.var` key, and womtool exposes only ONE level of call nesting.
# Both are leveraged here: the four setup calls are made directly from this workflow,
# so their knobs are the only nested inputs Terra surfaces; everything inside the
# delegated run is two levels down and therefore hidden. The `z1_`..`z4_` alias
# prefixes sort after the lowercase top-level inputs, fixing the displayed order:
#
#   1. MultiSampleSomaticWorkflow.*   curated user inputs
#   2. z1_Files.*                     WorkflowResources (reference files)
#   3. z2_Cache.*                     DefinePatient cache-injection fields
#   4. z3_Parameters.*                WorkflowArguments (run_*, scatter, ...)
#   5. z4_RuntimeParameters.*         RuntimeCollection (mem_*, time_*, ...)
#
# The setup is defined here only; MultiSampleSomaticWorkflowRun receives the finished
# structs and patient and contains no preprocessing of its own.

import "multi-sample_somatic_workflow.run.wdl" as mssw_run
import "patient.wdl" as p
import "patient.out.wdl" as p_out
import "patient.define.wdl" as p_def
import "workflow_arguments.wdl" as wfargs
import "workflow_resources.wdl" as wfres
import "runtime_collection.wdl" as rtc
import "tasks.wdl"


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

        # Used to split the genome for SNV calling so each shard on average completes
        # in ~2-3h. The default of 500 reads is chosen for a NGS WES paired tumor-normal
        # setting where each sample contributes ~250x read depth. Increase if more
        # samples are supplied.
        Int total_mean_read_depth = 500

        Patient? input_patient
        WorkflowArguments? input_args
        WorkflowResources? input_resources
        RuntimeCollection? input_runtime_collection
    }


###############################################################################
#                                                                             #
#                              PREPROCESSING                                  #
#                                                                             #
###############################################################################


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
                total_mean_read_depth = total_mean_read_depth,
                resources = resources,
                runtime_collection = runtime_collection,
        }
    }
    WorkflowArguments args = select_first([input_args, z3_Parameters.arguments])

    # Pre-flight gate: validate the raw inputs before any expensive task. ParseInput
    # exits non-zero (failing the workflow) on a mis-shaped input or an unmet implicit
    # requirement. Patient assembly and the delegated run are both guarded by
    # `ParseInput.validated`, which makes them depend on this task so that nothing
    # heavy starts until validation passes.
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


###############################################################################
#                                                                             #
#                                  RUN                                        #
#                                                                             #
###############################################################################


    # --- Delegate the pipeline; its internals are bound here, so Terra hides them.
    if (ParseInput.validated) {
        call mssw_run.MultiSampleSomaticWorkflowRun as run {
            input:
                patient = patient,
                args = args,
                resources = resources,
                runtime_collection = runtime_collection,
        }
    }


###############################################################################
#                                                                             #
#                                OUTPUT                                       #
#                                                                             #
###############################################################################


    # `run` is gated on ParseInput.validated, so its outputs are optional; the gate
    # only fails by aborting the workflow, so when we reach here `run` has produced a
    # patient. Coerce it once to a non-optional Patient for member access below.
    Patient out_patient = select_first([run.output_patient])

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
        Array[File]? annotated_somatic_variants = Output.annotated_somatic_variants
        Array[File?]? annotated_somatic_variants_idx = Output.annotated_somatic_variants_idx
        Array[File]? acs_copy_ratio_segmentation = Output.acs_copy_ratio_segmentation
        Array[Float]? acs_copy_ratio_skew = Output.acs_copy_ratio_skew

        Array[File]? absolute_snv_maf = Output.absolute_snv_maf
        Array[File]? absolute_indel_maf = Output.absolute_indel_maf

        Array[File]? absolute_acr_rdata = Output.absolute_acr_rdata
        Array[File]? absolute_acr_plot = Output.absolute_acr_plot
        Array[Int]? absolute_acr_solution = Output.absolute_acr_solution
        Array[File]? absolute_acr_maf = Output.absolute_acr_maf
        Array[File]? absolute_acr_segtab = Output.absolute_acr_segtab
        Array[File]? absolute_acr_segtab_igv = Output.absolute_acr_segtab_igv
        Array[File]? absolute_acr_table = Output.absolute_acr_table
        Array[Float]? absolute_acr_purity = Output.absolute_acr_purity
        Array[Float]? absolute_acr_ploidy = Output.absolute_acr_ploidy

        Array[File]? absolute_tcr_rdata = Output.absolute_tcr_rdata
        Array[File]? absolute_tcr_plot = Output.absolute_tcr_plot
        Array[Int]? absolute_tcr_solution = Output.absolute_tcr_solution
        Array[File]? absolute_tcr_maf = Output.absolute_tcr_maf
        Array[File]? absolute_tcr_segtab = Output.absolute_tcr_segtab
        Array[File]? absolute_tcr_segtab_igv = Output.absolute_tcr_segtab_igv
        Array[File]? absolute_tcr_table = Output.absolute_tcr_table
        Array[Float]? absolute_tcr_purity = Output.absolute_tcr_purity
        Array[Float]? absolute_tcr_ploidy = Output.absolute_tcr_ploidy

        # user-supplied purity/ploidy override, echoed back for round-tripping:
        Array[Float]? user_purity = Output.user_purity
        Array[Float]? user_ploidy = Output.user_ploidy

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
        Int? num_somatic_variants = out_patient.num_somatic_variants
        File? germline_vcf = out_patient.germline_vcf
        File? germline_vcf_idx = out_patient.germline_vcf_idx
        Int? num_germline_variants = out_patient.num_germline_variants
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
        Float? snp_sample_correlation_min = out_patient.snp_sample_correlation_min
        String? ancestry_pred = out_patient.ancestry_pred
        Float? ancestry_prob = out_patient.ancestry_prob
        File? ancestry_background_pca_table = out_patient.ancestry_background_pca_table
        File? ancestry_pca_plot = out_patient.ancestry_pca_plot
        File? modeled_segments = out_patient.modeled_segments
        File? phylogic_sif_file = out_patient.phylogic_sif_file
        File? phylogic_report = out_patient.phylogic_report
        File? phylogic_ccfs_cnvs = out_patient.phylogic_ccfs_cnvs
        File? phylogic_ccfs_snvs = out_patient.phylogic_ccfs_snvs
        File? phylogic_constrained_ccf = out_patient.phylogic_constrained_ccf
        File? phylogic_cluster_ccfs = out_patient.phylogic_cluster_ccfs
        File? phylogic_build_tree_posteriors = out_patient.phylogic_build_tree_posteriors
        File? phylogic_growth_rates = out_patient.phylogic_growth_rates
        File? phylogic_growth_rate_plot = out_patient.phylogic_growth_rate_plot
        File? phylogic_timing_report = out_patient.phylogic_timing_report
        File? phylogic_timing_wgd_supporting_events = out_patient.phylogic_timing_wgd_supporting_events
        File? phylogic_timing_graph = out_patient.phylogic_timing_graph
        File? phylogic_timing_comparison = out_patient.phylogic_timing_comparison
        File? phylogic_timing_table = out_patient.phylogic_timing_table

        # composite cache
        Patient output_patient = out_patient
        WorkflowArguments output_args = args
        WorkflowResources output_resources = resources
        RuntimeCollection output_runtime_collection = runtime_collection
    }
}
