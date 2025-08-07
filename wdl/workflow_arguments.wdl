version development

import "runtime_collection.wdl" as rtc
import "tasks.wdl"
import "workflow_resources.wdl" as wfres
import "workflow_resources.update.wdl" as wfres_update


struct WorkflowArguments {
    WorkflowResources files

    String analyst_id

    Int scatter_count_for_variant_calling
    Int scatter_count_for_pileups
    Int variants_per_scatter

    Boolean run_reorder_bam_contigs
    Boolean run_collect_callable_loci
    Boolean run_collect_total_read_counts
    Boolean run_collect_allelic_read_counts
    Boolean run_contamination_model
    Boolean run_orientation_bias_mixture_model
    Boolean run_variant_calling
    Boolean run_variant_calling_mutect1
    Boolean run_variant_filter
    Boolean run_realignment_filter
    Boolean run_variant_annotation
    Boolean run_variant_annotation_scattered
    Boolean run_model_segments
    Boolean run_clonal_decomposition

    Boolean keep_germline
    Boolean compress_output
    Boolean make_bamout

    # arguments
    # CNV WORKFLOW
    Int preprocess_intervals_bin_length
    Int preprocess_intervals_padding
    Int collect_read_counts_max_soft_clipped_bases
    Float min_snppanel_pop_af
    Float max_snppanel_pop_af
    Int min_snppanel_read_depth
    Float genotype_variants_normal_to_tumor_weight
    Float genotype_variants_min_genotype_likelihood
    Float genotype_variants_outlier_prior
    Int genotype_variants_overdispersion
    Float genotype_variants_ref_bias
    Int harmonize_min_target_length
    Int het_to_interval_mapping_max_distance
    Int model_segments_max_number_of_segments_per_chromosome
    Array[Int] model_segments_window_sizes
    Int model_segments_kernel_approximation_dimension
    Boolean model_segments_use_multi_sample_cr_segmentation
    Float model_segments_smoothing_credible_interval_threshold
    Float call_copy_ratios_neutral_segment_copy_ratio_lower_bound
    Float call_copy_ratios_neutral_segment_copy_ratio_upper_bound
    Float call_copy_ratios_outlier_neutral_segment_copy_ratio_z_score_threshold
    Float call_copy_ratios_z_score_threshold
    Int filter_germline_cnvs_min_segment_length

    String script_acs_conversion
    String script_genotype_variants
    String script_harmonize_copy_ratios
    String script_map_to_absolute_copy_number
    String script_merge_pileups
    String script_pileup_to_allelic_counts
    String script_phylogicndt_create_sif

    Int absolute_min_hets
    Int absolute_min_probes
    Float absolute_maf90_threshold
    String absolute_genome_build

    Boolean phylogic_use_segtab

    # SNV WORKFLOW
    Int min_read_depth
    Float mutect1_initial_tumor_lod
    Float mutect1_tumor_lod_to_emit
    Float mutect2_initial_tumor_lod
    Float mutect2_tumor_lod_to_emit
    Float mutect2_high_mem_factor
    Boolean mutect2_native_pair_hmm_use_double_precision
    Boolean mutect2_dont_use_soft_clipped_bases
    Boolean mutect2_use_linked_de_bruijn_graph
    Boolean mutect2_recover_all_dangling_branches
    Boolean mutect2_pileup_detection
    Boolean mutect2_genotype_germline_sites
    Int mutect2_downsampling_stride
    Int mutect2_max_reads_per_alignment_start
    Int mutect2_pcr_snv_qual
    Int mutect2_pcr_indel_qual
    Int filter_mutect2_max_median_fragment_length_difference
    Int filter_alignment_artifacts_max_reasonable_fragment_length
    Int hard_filter_min_base_quality
    Int hard_filter_min_mapping_quality
    Int hard_filter_min_fragment_length
    Int hard_filter_min_total_depth
    Int hard_filter_min_total_alt_count
    Int hard_filter_min_position_from_end_of_read
    Int hard_filter_min_read_orientation_quality
    Float hard_filter_germline_min_population_af
    Array[String] hard_filter_expressions
    Array[String] hard_filter_names
    String somatic_filter_whitelist
    String germline_filter_whitelist
    String funcotator_reference_version
    String funcotator_output_format
    String funcotator_variant_type
    String funcotator_transcript_selection_mode
    Boolean funcotator_prefer_mane_transcripts
    Boolean funcotator_use_gnomad
    Array[String]? funcotator_data_sources_paths
    Array[String]? funcotator_annotation_defaults
    Array[String]? funcotator_annotation_overrides
    Array[String]? funcotator_exclude_fields

    # expose extra arguments for import of this workflow
    String? split_intervals_extra_args
    String? getpileupsummaries_extra_args
    String? mutect2_extra_args
    String? filter_mutect2_extra_args
    String? variant_filtration_extra_args
    String? left_align_and_trim_variants_extra_args
    String? select_variants_extra_args
    String? select_high_conficence_variants_jexl_arg
    String? realignment_extra_args
    String? funcotate_extra_args
}


workflow DefineWorkflowArguments {
    input {
        WorkflowResources resources

        String analyst_id = "PH"

        Int scatter_count_base_for_variant_calling = 25
        Int scatter_count_for_pileups = 1
        Int total_mean_read_depth = 500
        Int total_mean_read_depth_per_scatter = 500
        Int variants_per_scatter = 50

        # workflow options
        Boolean run_reorder_bam_contigs = false
        Boolean run_collect_callable_loci = false
        Boolean run_collect_total_read_counts = true
        Boolean run_collect_allelic_read_counts = true
        Boolean run_contamination_model = true
        Boolean run_orientation_bias_mixture_model = true
        Boolean run_variant_calling = true
        Boolean run_variant_calling_mutect1 = true
        Boolean run_variant_filter = true
        Boolean run_realignment_filter = true
        Boolean run_variant_annotation = true
        Boolean run_variant_annotation_scattered = false
        Boolean run_model_segments = true
        Boolean run_clonal_decomposition = true

        Boolean keep_germline = true
        Boolean compress_output = true
        Boolean make_bamout = true

        # arguments
        Int preprocess_intervals_bin_length = 0
        Int preprocess_intervals_padding = 0

        # CNV WORKFLOW
        Int collect_read_counts_max_soft_clipped_bases = 0
        Float min_snppanel_pop_af = 0.01
        Float max_snppanel_pop_af = 1.0  # default: 0.2
        Int min_snppanel_read_depth = 10
        Float genotype_variants_normal_to_tumor_weight = 10.0
        Float genotype_variants_min_genotype_likelihood = 0.995  # = LOD threshold 5.3
        Float genotype_variants_outlier_prior = 0.00001
        Int genotype_variants_overdispersion = 10
        Float genotype_variants_ref_bias = 1.05
        Int harmonize_min_target_length = 20
        Int het_to_interval_mapping_max_distance = 250
        Int model_segments_max_number_of_segments_per_chromosome = 10000
        Array[Int] model_segments_window_sizes = [4, 8, 16, 32, 64, 128, 256, 512, 1024]
        Int model_segments_kernel_approximation_dimension = 200
        Boolean model_segments_use_multi_sample_cr_segmentation = true
        Float model_segments_smoothing_credible_interval_threshold = 3.0
        Float call_copy_ratios_neutral_segment_copy_ratio_lower_bound = 0.9
        Float call_copy_ratios_neutral_segment_copy_ratio_upper_bound = 1.1
        Float call_copy_ratios_outlier_neutral_segment_copy_ratio_z_score_threshold = 2.0
        Float call_copy_ratios_z_score_threshold = 2.0
        Int filter_germline_cnvs_min_segment_length = 100

        String script_acs_conversion =              "https://github.com/phylyc/somatic_workflow/raw/master/python/acs_conversion.py"
        String script_genotype_variants =           "https://github.com/phylyc/somatic_workflow/raw/master/python/genotype.py"
        String script_harmonize_copy_ratios =       "https://github.com/phylyc/somatic_workflow/raw/master/python/harmonize_copy_ratios.py"
        String script_map_to_absolute_copy_number = "https://github.com/phylyc/somatic_workflow/raw/master/python/map_to_absolute_copy_number.py"
        String script_merge_pileups =               "https://github.com/phylyc/somatic_workflow/raw/master/python/merge_pileups.py"
        String script_pileup_to_allelic_counts =    "https://github.com/phylyc/somatic_workflow/raw/master/python/pileup_to_allelic_counts.py"
        String script_phylogicndt_create_sif =      "https://github.com/phylyc/somatic_workflow/raw/master/python/create_patient_sif.py"

        Int absolute_min_hets = 0
        Int absolute_min_probes = 2
        Float absolute_maf90_threshold = 0.485
        String absolute_genome_build = "hg19"

        Boolean phylogic_use_segtab = true

        # SNV WORKFLOW
        Int min_read_depth = 4
        Float mutect1_initial_tumor_lod = 4.0
        Float mutect1_tumor_lod_to_emit = 6.0
        Float mutect2_initial_tumor_lod = 2.0
        Float mutect2_tumor_lod_to_emit = 3.0
        Float mutect2_high_mem_factor = 1.5
        # This is essentially a custom implementation of the mitochondiral model:
        Boolean mutect2_native_pair_hmm_use_double_precision = true
        Boolean mutect2_dont_use_soft_clipped_bases = false
        Boolean mutect2_use_linked_de_bruijn_graph = true
        Boolean mutect2_recover_all_dangling_branches = true
        Boolean mutect2_pileup_detection = true
        Boolean mutect2_genotype_germline_sites = true
        # The stride is the window in which the AVERAGE depth is required to meet
        # the max_reads_per_alignment_start. Usually a good idea to have a value of 20-50.
        Int mutect2_downsampling_stride = 1  # default: 1
        # This guards against amplicons.
        Int mutect2_max_reads_per_alignment_start = 100  # default: 50
        # Increase for high quality (de-duplexed, high-depth) panel sequencing data
        Int mutect2_pcr_snv_qual = 40 # default: 40
        Int mutect2_pcr_indel_qual = 40  # default: 40
        Int filter_mutect2_max_median_fragment_length_difference = 10000  # default: 10000
        Int filter_alignment_artifacts_max_reasonable_fragment_length = 10000 # default: 100000
        Int hard_filter_min_base_quality = 20
        Int hard_filter_min_mapping_quality = 20
        Int hard_filter_min_fragment_length = 18
        Int hard_filter_min_total_depth = 10
        Int hard_filter_min_total_alt_count = 3
        Int hard_filter_min_position_from_end_of_read = 6
        Int hard_filter_min_read_orientation_quality = 10
        Float hard_filter_germline_min_population_af = 3
        Array[String] hard_filter_expressions = []
        Array[String] hard_filter_names = []
        String somatic_filter_whitelist = "PASS,normal_artifact"
        String germline_filter_whitelist = "normal_artifact,panel_of_normals"
        String funcotator_reference_version = "hg19"
        String funcotator_output_format = "MAF"
        String funcotator_variant_type = "somatic"  # alternative: germline
        String funcotator_transcript_selection_mode = "CANONICAL"  # GATK default: "CANONICAL"
        Boolean funcotator_prefer_mane_transcripts = true
        Boolean funcotator_use_gnomad = false
        Array[String]? funcotator_data_sources_paths
        Array[String]? funcotator_annotation_defaults
        Array[String]? funcotator_annotation_overrides
        Array[String]? funcotator_exclude_fields

        # expose extra arguments for import of this workflow
        String? split_intervals_extra_args
        String? getpileupsummaries_extra_args
        String? mutect2_extra_args
        String? filter_mutect2_extra_args
        String? variant_filtration_extra_args
        String? left_align_and_trim_variants_extra_args
        String? select_variants_extra_args
        String? select_high_conficence_variants_jexl_arg = "GERMQ > 30"
        String? realignment_extra_args
        String? funcotate_extra_args

        RuntimeCollection runtime_collection
    }

    if (!defined(resources.preprocessed_intervals)) {
        call tasks.PreprocessIntervals {
            input:
                interval_list = resources.interval_list,
                interval_blacklist = resources.interval_blacklist,
                interval_lists = resources.interval_lists,
                ref_fasta = resources.ref_fasta,
                ref_fasta_index = resources.ref_fasta_index,
                ref_dict = resources.ref_dict,
                bin_length = preprocess_intervals_bin_length,
                padding = preprocess_intervals_padding,
                runtime_params = runtime_collection.preprocess_intervals,
        }
    }

    # The shards for variant calling always need to be defined.
    if (!defined(resources.scattered_intervals_for_variant_calling)) {
        Int good_scatter_count = ceil(scatter_count_base_for_variant_calling * (total_mean_read_depth + 1) / (total_mean_read_depth_per_scatter + 1))
        call tasks.SplitIntervals as VariantCallingSplitIntervals {
            input:
                interval_list = select_first([resources.preprocessed_intervals, PreprocessIntervals.preprocessed_interval_list]),
                ref_fasta = resources.ref_fasta,
                ref_fasta_index = resources.ref_fasta_index,
                ref_dict = resources.ref_dict,
                scatter_count = good_scatter_count,
                split_intervals_extra_args = split_intervals_extra_args,
                runtime_params = runtime_collection.split_intervals,
        }
    }

    if (!defined(resources.scattered_intervals_for_pileups) && run_collect_allelic_read_counts && (scatter_count_for_pileups > 1)) {
        call tasks.SplitIntervals as CollectAllelicCountsSplitIntervals {
            input:
                interval_list = select_first([resources.preprocessed_intervals, PreprocessIntervals.preprocessed_interval_list]),
                ref_fasta = resources.ref_fasta,
                ref_fasta_index = resources.ref_fasta_index,
                ref_dict = resources.ref_dict,
                scatter_count = scatter_count_for_pileups,
                split_intervals_extra_args = split_intervals_extra_args,
                runtime_params = runtime_collection.split_intervals,
        }
    }

    call wfres_update.UpdateWorkflowResources {
        input:
            resources = resources,
            preprocessed_intervals = select_first([resources.preprocessed_intervals, PreprocessIntervals.preprocessed_interval_list]),
            scattered_intervals_for_variant_calling = select_first([resources.scattered_intervals_for_variant_calling, VariantCallingSplitIntervals.interval_files]),
            scattered_intervals_for_pileups = CollectAllelicCountsSplitIntervals.interval_files
    }

    WorkflowArguments args = object {
        files: UpdateWorkflowResources.updated_resources,

        analyst_id: analyst_id,

        scatter_count_for_variant_calling: length(select_first([resources.scattered_intervals_for_variant_calling, VariantCallingSplitIntervals.interval_files])),
        scatter_count_for_pileups: length(select_first([resources.scattered_intervals_for_pileups, CollectAllelicCountsSplitIntervals.interval_files, ["ONE"]])),
        variants_per_scatter: variants_per_scatter,

        run_reorder_bam_contigs: run_reorder_bam_contigs,
        run_collect_callable_loci: run_collect_callable_loci,
        run_collect_total_read_counts: run_collect_total_read_counts,
        run_collect_allelic_read_counts: run_collect_allelic_read_counts,
        run_contamination_model: run_contamination_model,
        run_orientation_bias_mixture_model: run_orientation_bias_mixture_model,
        run_variant_calling: run_variant_calling,
        run_variant_calling_mutect1: run_variant_calling_mutect1,
        run_variant_filter: run_variant_filter,
        run_realignment_filter: run_realignment_filter,
        run_variant_annotation: run_variant_annotation,
        run_variant_annotation_scattered: run_variant_annotation_scattered,
        run_model_segments: run_model_segments,
        run_clonal_decomposition: run_clonal_decomposition,

        keep_germline: keep_germline,
        compress_output: compress_output,
        make_bamout: make_bamout,

        preprocess_intervals_bin_length: preprocess_intervals_bin_length,
        preprocess_intervals_padding: preprocess_intervals_padding,
        collect_read_counts_max_soft_clipped_bases: collect_read_counts_max_soft_clipped_bases,
        min_snppanel_pop_af: min_snppanel_pop_af,
        max_snppanel_pop_af: max_snppanel_pop_af,
        min_snppanel_read_depth: min_snppanel_read_depth,
        genotype_variants_normal_to_tumor_weight: genotype_variants_normal_to_tumor_weight,
        genotype_variants_min_genotype_likelihood: genotype_variants_min_genotype_likelihood,
        genotype_variants_outlier_prior: genotype_variants_outlier_prior,
        genotype_variants_overdispersion: genotype_variants_overdispersion,
        genotype_variants_ref_bias: genotype_variants_ref_bias,
        harmonize_min_target_length: harmonize_min_target_length,
        het_to_interval_mapping_max_distance: het_to_interval_mapping_max_distance,
        model_segments_max_number_of_segments_per_chromosome: model_segments_max_number_of_segments_per_chromosome,
        model_segments_window_sizes: model_segments_window_sizes,
        model_segments_kernel_approximation_dimension: model_segments_kernel_approximation_dimension,
        model_segments_use_multi_sample_cr_segmentation: model_segments_use_multi_sample_cr_segmentation,
        model_segments_smoothing_credible_interval_threshold: model_segments_smoothing_credible_interval_threshold,
        call_copy_ratios_neutral_segment_copy_ratio_lower_bound: call_copy_ratios_neutral_segment_copy_ratio_lower_bound,
        call_copy_ratios_neutral_segment_copy_ratio_upper_bound: call_copy_ratios_neutral_segment_copy_ratio_upper_bound,
        call_copy_ratios_outlier_neutral_segment_copy_ratio_z_score_threshold: call_copy_ratios_outlier_neutral_segment_copy_ratio_z_score_threshold,
        call_copy_ratios_z_score_threshold: call_copy_ratios_z_score_threshold,
        filter_germline_cnvs_min_segment_length: filter_germline_cnvs_min_segment_length,

        script_acs_conversion: script_acs_conversion,
        script_genotype_variants: script_genotype_variants,
        script_harmonize_copy_ratios: script_harmonize_copy_ratios,
        script_map_to_absolute_copy_number: script_map_to_absolute_copy_number,
        script_merge_pileups: script_merge_pileups,
        script_pileup_to_allelic_counts: script_pileup_to_allelic_counts,
        script_phylogicndt_create_sif: script_phylogicndt_create_sif,

        absolute_min_hets: absolute_min_hets,
        absolute_min_probes: absolute_min_probes,
        absolute_maf90_threshold: absolute_maf90_threshold,
        absolute_genome_build: absolute_genome_build,

        phylogic_use_segtab: phylogic_use_segtab,

        min_read_depth: min_read_depth,
        mutect1_initial_tumor_lod: mutect1_initial_tumor_lod,
        mutect1_tumor_lod_to_emit: mutect1_tumor_lod_to_emit,
        mutect2_initial_tumor_lod: mutect2_initial_tumor_lod,
        mutect2_tumor_lod_to_emit: mutect2_tumor_lod_to_emit,
        mutect2_high_mem_factor: mutect2_high_mem_factor,
        mutect2_native_pair_hmm_use_double_precision: mutect2_native_pair_hmm_use_double_precision,
        mutect2_dont_use_soft_clipped_bases: mutect2_dont_use_soft_clipped_bases,
        mutect2_use_linked_de_bruijn_graph: mutect2_use_linked_de_bruijn_graph,
        mutect2_recover_all_dangling_branches: mutect2_recover_all_dangling_branches,
        mutect2_pileup_detection: mutect2_pileup_detection,
        mutect2_genotype_germline_sites: mutect2_genotype_germline_sites,
        mutect2_downsampling_stride: mutect2_downsampling_stride,
        mutect2_max_reads_per_alignment_start: mutect2_max_reads_per_alignment_start,
        mutect2_pcr_snv_qual: mutect2_pcr_snv_qual,
        mutect2_pcr_indel_qual: mutect2_pcr_indel_qual,
        filter_mutect2_max_median_fragment_length_difference: filter_mutect2_max_median_fragment_length_difference,
        filter_alignment_artifacts_max_reasonable_fragment_length: filter_alignment_artifacts_max_reasonable_fragment_length,
        hard_filter_min_base_quality: hard_filter_min_base_quality,
        hard_filter_min_mapping_quality: hard_filter_min_mapping_quality,
        hard_filter_min_fragment_length: hard_filter_min_fragment_length,
        hard_filter_min_total_depth: hard_filter_min_total_depth,
        hard_filter_min_total_alt_count: hard_filter_min_total_alt_count,
        hard_filter_min_position_from_end_of_read: hard_filter_min_position_from_end_of_read,
        hard_filter_min_read_orientation_quality: hard_filter_min_read_orientation_quality,
        hard_filter_germline_min_population_af: hard_filter_germline_min_population_af,
        hard_filter_expressions: hard_filter_expressions,
        hard_filter_names: hard_filter_names,
        somatic_filter_whitelist: somatic_filter_whitelist,
        germline_filter_whitelist: germline_filter_whitelist,
        funcotator_reference_version: funcotator_reference_version,
        funcotator_output_format: funcotator_output_format,
        funcotator_variant_type: funcotator_variant_type,
        funcotator_transcript_selection_mode: funcotator_transcript_selection_mode,
        funcotator_prefer_mane_transcripts: funcotator_prefer_mane_transcripts,
        funcotator_use_gnomad: funcotator_use_gnomad,
        funcotator_data_sources_paths: funcotator_data_sources_paths,
        funcotator_annotation_defaults: funcotator_annotation_defaults,
        funcotator_annotation_overrides: funcotator_annotation_overrides,
        funcotator_exclude_fields: funcotator_exclude_fields,

        split_intervals_extra_args: split_intervals_extra_args,
        getpileupsummaries_extra_args: getpileupsummaries_extra_args,
        mutect2_extra_args: mutect2_extra_args,
        filter_mutect2_extra_args: filter_mutect2_extra_args,
        variant_filtration_extra_args: variant_filtration_extra_args,
        left_align_and_trim_variants_extra_args: left_align_and_trim_variants_extra_args,
        select_variants_extra_args: select_variants_extra_args,
        select_high_conficence_variants_jexl_arg: select_high_conficence_variants_jexl_arg,
        realignment_extra_args: realignment_extra_args,
        funcotate_extra_args: funcotate_extra_args
    }

    output {
        WorkflowArguments arguments = args
    }
}