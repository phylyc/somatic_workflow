version development

import "runtime_collection.wdl" as rtc
import "tasks.wdl"
import "workflow_resources.wdl"


struct WorkflowArguments {
    WorkflowResources files

    Int scatter_count
    Int total_mean_read_depth
    Int variants_per_scatter

    File preprocessed_interval_list
    Array[File] scattered_interval_list

    Boolean run_collect_covered_regions
    Boolean run_collect_target_coverage
    Boolean run_collect_allelic_coverage
    Boolean run_contamination_model
    Boolean run_orientation_bias_mixture_model
    Boolean run_variant_calling
    Boolean run_variant_filter
    Boolean run_variant_hard_filter
    Boolean run_realignment_filter
    Boolean run_realignment_filter_only_on_high_confidence_variants
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
    Float genotype_variants_min_genotype_likelihood
    Int genotype_variants_overdispersion
    Float genotype_variants_ref_bias
    Int harmonize_min_target_length
    Array[Int] model_segments_window_sizes
    Float call_copy_ratios_neutral_segment_copy_ratio_lower_bound
    Float call_copy_ratios_neutral_segment_copy_ratio_upper_bound
    Float call_copy_ratios_outlier_neutral_segment_copy_ratio_z_score_threshold
    Float call_copy_ratios_z_score_threshold
    Int filter_germline_cnvs_min_segment_length

    String genotype_variants_script
    String harmonize_copy_ratios_script
    String merge_pileups_script
    String acs_conversion_script

    Int absolute_min_hets
    Int absolute_min_probes
    Float absolute_maf90_threshold

    # SNV WORKFLOW
    Int min_read_depth
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
    Int filter_mutect2_min_alt_median_base_quality
    Int filter_mutect2_min_alt_median_mapping_quality
    Int filter_mutect2_min_median_read_position
    Int filter_alignment_artifacts_max_reasonable_fragment_length
    Array[String] hard_filter_expressions
    Array[String] hard_filter_names
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
    String? select_low_conficence_variants_jexl_arg
    String? realignment_extra_args
    String? funcotate_extra_args
}


workflow DefineWorkflowArguments {
    input {
        WorkflowResources resources

        Int scatter_count = 10
        Int total_mean_read_depth = 500
        Int total_mean_read_depth_per_scatter = 500
        Int variants_per_scatter = 50

        # workflow options
        Boolean run_collect_covered_regions = false
        Boolean run_collect_target_coverage = true
        Boolean run_collect_allelic_coverage = true
        Boolean run_contamination_model = true
        Boolean run_model_segments = true
        Boolean run_orientation_bias_mixture_model = true
        Boolean run_variant_calling = true
        Boolean run_variant_filter = true
        Boolean run_variant_hard_filter = true
        Boolean run_realignment_filter = true
        Boolean run_realignment_filter_only_on_high_confidence_variants = true
        Boolean run_variant_annotation = true
        Boolean run_variant_annotation_scattered = false
        Boolean run_clonal_decomposition = true

        Boolean keep_germline = true
        Boolean compress_output = true
        Boolean make_bamout = false

        # arguments
        Int preprocess_intervals_bin_length = 0
        Int preprocess_intervals_padding = 0

        # CNV WORKFLOW
        Int collect_read_counts_max_soft_clipped_bases = 0
        Float min_snppanel_pop_af = 0.01
        Float max_snppanel_pop_af = 1.0  # default: 0.2
        Int min_snppanel_read_depth = 10
        Float genotype_variants_min_genotype_likelihood = 0.9
        Int genotype_variants_overdispersion = 50
        Float genotype_variants_ref_bias = 1.05
        Int harmonize_min_target_length = 100
        Array[Int] model_segments_window_sizes = [8, 16, 32, 64, 128, 256, 512]
        Float call_copy_ratios_neutral_segment_copy_ratio_lower_bound = 0.9
        Float call_copy_ratios_neutral_segment_copy_ratio_upper_bound = 1.1
        Float call_copy_ratios_outlier_neutral_segment_copy_ratio_z_score_threshold = 2.0
        Float call_copy_ratios_z_score_threshold = 2.0
        Int filter_germline_cnvs_min_segment_length = 100

        String genotype_variants_script =       "https://github.com/phylyc/somatic_workflow/raw/master/python/genotype.py"
        String harmonize_copy_ratios_script =   "https://github.com/phylyc/somatic_workflow/raw/master/python/harmonize_copy_ratios.py"
        String merge_pileups_script =           "https://github.com/phylyc/somatic_workflow/raw/master/python/merge_pileups.py"
        String acs_conversion_script =          "https://github.com/phylyc/somatic_workflow/raw/master/python/acs_conversion.py"

        Int absolute_min_hets = 10
        Int absolute_min_probes = 4
        Float absolute_maf90_threshold = 0.485

        # SNV WORKFLOW
        Int min_read_depth = 4
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
        Int filter_mutect2_min_alt_median_base_quality = 20  # default: 20
        Int filter_mutect2_min_alt_median_mapping_quality = 20  # default: -1
        Int filter_mutect2_min_median_read_position = 5  # default: 1
        Int filter_alignment_artifacts_max_reasonable_fragment_length = 10000 # default: 100000
        Array[String] hard_filter_expressions = [
            "DP < 10",
            "MBQ.0 < 20", "MBQ.1 < 20",
            "MMQ.0 < 20", "MMQ.1 < 20",
            "MFRL.0 < 18", "MFRL.1 < 18",
            "MPOS.0 < 6",
            "ROQ < 10"
        ]
        Array[String] hard_filter_names = [
            "lowDP",
            "lowMBQ.0", "lowMBQ.1",
            "lowMMQ.0", "lowMMQ.1",
            "lowMFRL.0", "lowMFRL.1",
            "lowMPOS",
            "lowROQ"
        ]
        String germline_filter_whitelist = "normal_artifact,panel_of_normals"
        String funcotator_reference_version = "hg19"
        String funcotator_output_format = "MAF"
        String funcotator_variant_type = "somatic"  # alternative: germline
        String funcotator_transcript_selection_mode = "CANONICAL"  # GATK default: "CANONICAL"
        Boolean funcotator_prefer_mane_transcripts = true
        Boolean funcotator_use_gnomad = true
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
        String? select_low_conficence_variants_jexl_arg = "GERMQ < 30"
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

    if (!defined(resources.scattered_intervals)) {
        Int good_scatter_count = ceil(scatter_count * (total_mean_read_depth + 1) / (total_mean_read_depth_per_scatter + 1))
        call tasks.SplitIntervals {
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

    WorkflowArguments args = object {
        files: resources,

        scatter_count: length(select_first([resources.scattered_intervals, SplitIntervals.interval_files])),
        total_mean_read_depth: total_mean_read_depth,
        variants_per_scatter: variants_per_scatter,

        preprocessed_interval_list: select_first([resources.preprocessed_intervals, PreprocessIntervals.preprocessed_interval_list]),
        scattered_interval_list: select_first([resources.scattered_intervals, SplitIntervals.interval_files]),

        run_collect_covered_regions: run_collect_covered_regions,
        run_collect_target_coverage: run_collect_target_coverage,
        run_collect_allelic_coverage: run_collect_allelic_coverage,
        run_contamination_model: run_contamination_model,
        run_model_segments: run_model_segments,
        run_orientation_bias_mixture_model: run_orientation_bias_mixture_model,
        run_variant_calling: run_variant_calling,
        run_variant_filter: run_variant_filter,
        run_variant_hard_filter: run_variant_hard_filter,
        run_realignment_filter: run_realignment_filter,
        run_realignment_filter_only_on_high_confidence_variants: run_realignment_filter_only_on_high_confidence_variants,
        run_variant_annotation: run_variant_annotation,
        run_variant_annotation_scattered: run_variant_annotation_scattered,
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
        genotype_variants_min_genotype_likelihood: genotype_variants_min_genotype_likelihood,
        genotype_variants_overdispersion: genotype_variants_overdispersion,
        genotype_variants_ref_bias: genotype_variants_ref_bias,
        harmonize_min_target_length: harmonize_min_target_length,
        model_segments_window_sizes: model_segments_window_sizes,
        call_copy_ratios_neutral_segment_copy_ratio_lower_bound: call_copy_ratios_neutral_segment_copy_ratio_lower_bound,
        call_copy_ratios_neutral_segment_copy_ratio_upper_bound: call_copy_ratios_neutral_segment_copy_ratio_upper_bound,
        call_copy_ratios_outlier_neutral_segment_copy_ratio_z_score_threshold: call_copy_ratios_outlier_neutral_segment_copy_ratio_z_score_threshold,
        call_copy_ratios_z_score_threshold: call_copy_ratios_z_score_threshold,
        filter_germline_cnvs_min_segment_length: filter_germline_cnvs_min_segment_length,

        genotype_variants_script: genotype_variants_script,
        harmonize_copy_ratios_script: harmonize_copy_ratios_script,
        merge_pileups_script: merge_pileups_script,
        acs_conversion_script: acs_conversion_script,

        absolute_min_hets: absolute_min_hets,
        absolute_min_probes: absolute_min_probes,
        absolute_maf90_threshold: absolute_maf90_threshold,

        min_read_depth: min_read_depth,
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
        filter_mutect2_min_alt_median_base_quality: filter_mutect2_min_alt_median_base_quality,
        filter_mutect2_min_alt_median_mapping_quality: filter_mutect2_min_alt_median_mapping_quality,
        filter_mutect2_min_median_read_position: filter_mutect2_min_median_read_position,
        filter_alignment_artifacts_max_reasonable_fragment_length: filter_alignment_artifacts_max_reasonable_fragment_length,
        hard_filter_expressions: hard_filter_expressions,
        hard_filter_names: hard_filter_names,
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
        select_low_conficence_variants_jexl_arg: select_low_conficence_variants_jexl_arg,
        realignment_extra_args: realignment_extra_args,
        funcotate_extra_args: funcotate_extra_args
    }

    output {
        WorkflowArguments arguments = args
    }
}