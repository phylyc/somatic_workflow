version development

import "runtime_collection.wdl" as rtc
import "tasks.wdl"
import "workflow_resources.wdl"


struct WorkflowArguments {
    WorkflowResources files

    Int scatter_count

    File preprocessed_interval_list
    Array[File] scattered_interval_list

    Boolean run_collect_covered_regions
    Boolean run_collect_target_coverage
    Boolean run_collect_allelic_coverage
    Boolean run_contamination_model
    Boolean run_model_segments
    Boolean run_orientation_bias_mixture_model
    Boolean run_variant_calling
    Boolean run_variant_filter
    Boolean run_realignment_filter
    Boolean run_realignment_filter_only_on_high_confidence_variants
    Boolean run_collect_called_variants_allelic_coverage
    Boolean run_variant_annotation
    Boolean run_variant_annotation_scattered
    Boolean run_clonal_decomposition

    Boolean keep_germline
    Boolean compress_output
    Boolean make_bamout

    # arguments
    # CNV WORKFLOW
    Int preprocess_intervals_bin_length
    Int preprocess_intervals_padding
    Float min_snp_array_pop_af
    Float max_snp_array_pop_af
    Int min_snp_array_read_depth
    Array[Int] model_segments_window_sizes
    Float call_copy_ratios_neutral_segment_copy_ratio_lower_bound
    Float call_copy_ratios_neutral_segment_copy_ratio_upper_bound
    Float call_copy_ratios_outlier_neutral_segment_copy_ratio_z_score_threshold
    Float call_copy_ratios_z_score_threshold

    String genotype_variants_script
    String harmonize_copy_ratios_script
    String merge_pileups_script
    String acs_conversion_script

    # SNV WORKFLOW
    Boolean mutect2_native_pair_hmm_use_double_precision
    Boolean mutect2_use_linked_de_bruijn_graph
    Boolean mutect2_recover_all_dangling_branches
    Boolean mutect2_pileup_detection
    Boolean mutect2_genotype_germline_sites
    Int mutect2_downsampling_stride
    Int mutect2_max_reads_per_alignment_start
    Int filter_mutect2_max_median_fragment_length_difference
    Int filter_mutect2_min_alt_median_base_quality
    Int filter_mutect2_min_alt_median_mapping_quality
    Int filter_mutect2_min_median_read_position
    Int filter_alignment_artifacts_max_reasonable_fragment_length
    String funcotator_reference_version
    String funcotator_output_format
    String funcotator_variant_type
    String funcotator_transcript_selection_mode
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
    String? select_variants_extra_args
    String? select_low_conficence_variants_jexl_arg
    String? realignment_extra_args
    String? funcotate_extra_args
}


workflow DefineWorkflowArguments {
    input {
        WorkflowResources resources

        Int scatter_count = 10

        # workflow options
        Boolean run_collect_covered_regions = false
        Boolean run_collect_target_coverage = true
        Boolean run_collect_allelic_coverage = true
        Boolean run_contamination_model = true
        Boolean run_model_segments = true
        Boolean run_orientation_bias_mixture_model = true
        Boolean run_variant_calling = true
        Boolean run_variant_filter = true
        Boolean run_realignment_filter = true
        Boolean run_realignment_filter_only_on_high_confidence_variants = true
        Boolean run_collect_called_variants_allelic_coverage = true
        Boolean run_variant_annotation = true
        Boolean run_variant_annotation_scattered = false
        Boolean run_clonal_decomposition = true

        Boolean keep_germline = false
        Boolean compress_output = true
        Boolean make_bamout = false

        # arguments
        Int preprocess_intervals_bin_length = 0
        Int preprocess_intervals_padding = 0

        # CNV WORKFLOW
        Float min_snp_array_pop_af = 0.01
        Float max_snp_array_pop_af = 1.0  # default: 0.2
        Int min_snp_array_read_depth = 10
        Array[Int] model_segments_window_sizes = [2, 4, 8, 16, 32, 64, 128, 256, 512]
        Float call_copy_ratios_neutral_segment_copy_ratio_lower_bound = 0.9
        Float call_copy_ratios_neutral_segment_copy_ratio_upper_bound = 1.1
        Float call_copy_ratios_outlier_neutral_segment_copy_ratio_z_score_threshold = 2.0
        Float call_copy_ratios_z_score_threshold = 2.0

        String genotype_variants_script = "https://github.com/phylyc/genomics_workflows/raw/master/python/genotype.py"
        String harmonize_copy_ratios_script = "https://github.com/phylyc/genomics_workflows/raw/master/python/harmonize_copy_ratios.py"
        String merge_pileups_script = "https://github.com/phylyc/genomics_workflows/raw/master/python/merge_pileups.py"
        String acs_conversion_script = "https://github.com/phylyc/genomics_workflows/raw/master/python/acs_conversion.py"

        # SNV WORKFLOW
        Int min_read_depth = 4
        Boolean mutect2_native_pair_hmm_use_double_precision = true
        Boolean mutect2_use_linked_de_bruijn_graph = true
        Boolean mutect2_recover_all_dangling_branches = true
        Boolean mutect2_pileup_detection = true
        Boolean mutect2_genotype_germline_sites = false  # use with care! (see above)
        Int mutect2_downsampling_stride = 1  # default: 1
        Int mutect2_max_reads_per_alignment_start = 100  # default: 50
        Int filter_mutect2_max_median_fragment_length_difference = 10000  # default: 10000
        Int filter_mutect2_min_alt_median_base_quality = 20  # default: 20
        Int filter_mutect2_min_alt_median_mapping_quality = 20  # default: -1
        Int filter_mutect2_min_median_read_position = 5  # default: 1
        Int filter_alignment_artifacts_max_reasonable_fragment_length = 10000 # default: 100000
        String funcotator_reference_version = "hg19"
        String funcotator_output_format = "MAF"
        String funcotator_variant_type = "somatic"  # alternative: germline
        String funcotator_transcript_selection_mode = "CANONICAL"  # GATK default: "CANONICAL"
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
        String? select_variants_extra_args = "-select '(vc.getAttribute(\"DP\") > 3) && (vc.getAttribute(\"MBQ\").0 > 19) && (vc.getAttribute(\"MBQ\").1 > 19) && (vc.getAttribute(\"MMQ\").0 > 19) && (vc.getAttribute(\"MMQ\").1 > 19) && (vc.getAttribute(\"MFRL\").0 > 17) && (vc.getAttribute(\"MFRL\").1 > 17)'"
        String? select_low_conficence_variants_jexl_arg = "'(vc.getAttribute(\"GERMQ\") < 30)'"
        String? realignment_extra_args
        String? funcotate_extra_args

        RuntimeCollection runtime_collection
    }

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

    call tasks.SplitIntervals {
    	input:
            interval_list = PreprocessIntervals.preprocessed_interval_list,
            ref_fasta = resources.ref_fasta,
            ref_fasta_index = resources.ref_fasta_index,
            ref_dict = resources.ref_dict,
            scatter_count = scatter_count,
            split_intervals_extra_args = split_intervals_extra_args,
            runtime_params = runtime_collection.split_intervals,
    }

    WorkflowArguments args = object {
        files: resources,

        scatter_count: scatter_count,

        preprocessed_interval_list: PreprocessIntervals.preprocessed_interval_list,
        scattered_interval_list: SplitIntervals.interval_files,

        run_collect_covered_regions: run_collect_covered_regions,
        run_collect_target_coverage: run_collect_target_coverage,
        run_collect_allelic_coverage: run_collect_allelic_coverage,
        run_contamination_model: run_contamination_model,
        run_model_segments: run_model_segments,
        run_orientation_bias_mixture_model: run_orientation_bias_mixture_model,
        run_variant_calling: run_variant_calling,
        run_variant_filter: run_variant_filter,
        run_realignment_filter: run_realignment_filter,
        run_realignment_filter_only_on_high_confidence_variants: run_realignment_filter_only_on_high_confidence_variants,
        run_collect_called_variants_allelic_coverage: run_collect_called_variants_allelic_coverage,
        run_variant_annotation: run_variant_annotation,
        run_variant_annotation_scattered: run_variant_annotation_scattered,
        run_clonal_decomposition: run_clonal_decomposition,

        keep_germline: keep_germline,
        compress_output: compress_output,
        make_bamout: make_bamout,

        preprocess_intervals_bin_length: preprocess_intervals_bin_length,
        preprocess_intervals_padding: preprocess_intervals_padding,
        min_snp_array_pop_af: min_snp_array_pop_af,
        max_snp_array_pop_af: max_snp_array_pop_af,
        min_snp_array_read_depth: min_snp_array_read_depth,
        model_segments_window_sizes: model_segments_window_sizes,
        call_copy_ratios_neutral_segment_copy_ratio_lower_bound: call_copy_ratios_neutral_segment_copy_ratio_lower_bound,
        call_copy_ratios_neutral_segment_copy_ratio_upper_bound: call_copy_ratios_neutral_segment_copy_ratio_upper_bound,
        call_copy_ratios_outlier_neutral_segment_copy_ratio_z_score_threshold: call_copy_ratios_outlier_neutral_segment_copy_ratio_z_score_threshold,
        call_copy_ratios_z_score_threshold: call_copy_ratios_z_score_threshold,

        genotype_variants_script: genotype_variants_script,
        harmonize_copy_ratios_script: harmonize_copy_ratios_script,
        merge_pileups_script: merge_pileups_script,
        acs_conversion_script: acs_conversion_script,

        min_read_depth: min_read_depth,
        mutect2_native_pair_hmm_use_double_precision: mutect2_native_pair_hmm_use_double_precision,
        mutect2_use_linked_de_bruijn_graph: mutect2_use_linked_de_bruijn_graph,
        mutect2_recover_all_dangling_branches: mutect2_recover_all_dangling_branches,
        mutect2_pileup_detection: mutect2_pileup_detection,
        mutect2_genotype_germline_sites: mutect2_genotype_germline_sites,
        mutect2_downsampling_stride: mutect2_downsampling_stride,
        mutect2_max_reads_per_alignment_start: mutect2_max_reads_per_alignment_start,
        filter_mutect2_max_median_fragment_length_difference: filter_mutect2_max_median_fragment_length_difference,
        filter_mutect2_min_alt_median_base_quality: filter_mutect2_min_alt_median_base_quality,
        filter_mutect2_min_alt_median_mapping_quality: filter_mutect2_min_alt_median_mapping_quality,
        filter_mutect2_min_median_read_position: filter_mutect2_min_median_read_position,
        filter_alignment_artifacts_max_reasonable_fragment_length: filter_alignment_artifacts_max_reasonable_fragment_length,
        funcotator_reference_version: funcotator_reference_version,
        funcotator_output_format: funcotator_output_format,
        funcotator_variant_type: funcotator_variant_type,
        funcotator_transcript_selection_mode: funcotator_transcript_selection_mode,
        funcotator_use_gnomad: funcotator_use_gnomad,
        funcotator_data_sources_paths: funcotator_data_sources_paths,
        funcotator_annotation_defaults: funcotator_annotation_defaults,
        funcotator_annotation_overrides: funcotator_annotation_overrides,
        funcotator_exclude_fields: funcotator_exclude_fields,

        split_intervals_extra_args: split_intervals_extra_args,
        getpileupsummaries_extra_args: getpileupsummaries_extra_args,
        mutect2_extra_args: mutect2_extra_args,
        filter_mutect2_extra_args: filter_mutect2_extra_args,
        select_variants_extra_args: select_variants_extra_args,
        select_low_conficence_variants_jexl_arg: select_low_conficence_variants_jexl_arg,
        realignment_extra_args: realignment_extra_args,
        funcotate_extra_args: funcotate_extra_args
    }

    output {
        WorkflowArguments arguments = args
    }
}