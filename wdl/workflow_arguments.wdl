version development

import "runtimes.wdl"
import "tasks.wdl"


struct WorkflowArguments {
    Int scatter_count
    File? interval_list
    File? interval_blacklist
    Array[File]? interval_lists

    File preprocessed_interval_list
    Array[File] scattered_interval_list

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    File? force_call_alleles
    File? force_call_alleles_idx
    File? snv_panel_of_normals
    File? snv_panel_of_normals_idx
    File? germline_resource
    File? germline_resource_tbi
    File? common_germline_alleles
    File? common_germline_alleles_idx
    File? realignment_bwa_mem_index_image
    File? funcotator_transcript_list
    File? funcotator_data_sources_tar_gz

    Boolean run_collect_covered_regions
    Boolean run_collect_target_coverage
    Boolean run_collect_allelic_coverage
    Boolean run_contamination_model
    Boolean run_orientation_bias_mixture_model
    Boolean run_variant_calling
    Boolean run_variant_filter
    Boolean run_realignment_filter
    Boolean run_realignment_filter_only_on_high_confidence_variants
    Boolean run_collect_called_variants_allelic_coverage
    Boolean run_variant_annotation
    Boolean run_variant_annotation_scattered

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
    String genotype_variants_script
    Boolean genotype_variants_save_sample_genotype_likelihoods

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
        Int scatter_count

        File? interval_list
        File? interval_blacklist
        Array[File]? interval_lists

        File ref_fasta
        File ref_fasta_index
        File ref_dict

        # resources
        File? force_call_alleles
        File? force_call_alleles_idx
        File? snv_panel_of_normals
        File? snv_panel_of_normals_idx
        File? germline_resource
        File? germline_resource_tbi
        File? common_germline_alleles
        File? common_germline_alleles_idx
        File? realignment_bwa_mem_index_image
        File? funcotator_transcript_list
        File? funcotator_data_sources_tar_gz

        # workflow options
        Boolean run_collect_covered_regions = false
        Boolean run_collect_target_coverage = true
        Boolean run_collect_allelic_coverage = true
        Boolean run_contamination_model = true
        Boolean run_orientation_bias_mixture_model = true
        Boolean run_variant_calling = true
        Boolean run_variant_filter = true
        Boolean run_realignment_filter = true
        Boolean run_realignment_filter_only_on_high_confidence_variants = false
        Boolean run_collect_called_variants_allelic_coverage = true
        Boolean run_variant_annotation = true
        Boolean run_variant_annotation_scattered = false

        Boolean keep_germline = false
        Boolean compress_output = true
        Boolean make_bamout = false

        # arguments
        # CNV WORKFLOW
        Int preprocess_intervals_bin_length = 0
        Int preprocess_intervals_padding = 0
        Float min_snp_array_pop_af = 0.01
        Float max_snp_array_pop_af = 1.0  # default: 0.2
        Int min_snp_array_read_depth = 10
        String genotype_variants_script = "https://github.com/phylyc/genomics_workflows/raw/master/python/genotype.py"
        Boolean genotype_variants_save_sample_genotype_likelihoods = false

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
        String? select_variants_extra_args
        String? select_low_conficence_variants_jexl_arg = "'(vc.getAttribute(\"GERMQ\") < 30) || (vc.getAttribute(\"DP\") < 4) || (vc.getAttribute(\"MBQ\").0 == 0) || (vc.getAttribute(\"MFRL\").0 == 0)'"
        String? realignment_extra_args
        String? funcotate_extra_args

        RuntimeCollection runtime_collection
    }

    call tasks.PreprocessIntervals {
        input:
            interval_list = interval_list,
            interval_blacklist = interval_blacklist,
            interval_lists = interval_lists,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            bin_length = preprocess_intervals_bin_length,
            padding = preprocess_intervals_padding,
            runtime_params = runtime_collection.preprocess_intervals,
    }

    call tasks.SplitIntervals {
    	input:
            interval_list = PreprocessIntervals.preprocessed_interval_list,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            scatter_count = scatter_count,
            split_intervals_extra_args = split_intervals_extra_args,
            runtime_params = runtime_collection.split_intervals,
    }

    WorkflowArguments args = object {
        scatter_count: scatter_count,

        interval_list: interval_list,
        interval_blacklist: interval_blacklist,
        interval_lists: interval_lists,

        preprocessed_interval_list: PreprocessIntervals.preprocessed_interval_list,
        scattered_interval_list: SplitIntervals.interval_files,

        ref_fasta: ref_fasta,
        ref_fasta_index: ref_fasta_index,
        ref_dict: ref_dict,

        force_call_alleles: force_call_alleles,
        force_call_alleles_idx: force_call_alleles_idx,
        snv_panel_of_normals: snv_panel_of_normals,
        snv_panel_of_normals_idx: snv_panel_of_normals_idx,
        germline_resource: germline_resource,
        germline_resource_tbi: germline_resource_tbi,
        common_germline_alleles: common_germline_alleles,
        common_germline_alleles_idx: common_germline_alleles_idx,
        realignment_bwa_mem_index_image: realignment_bwa_mem_index_image,
        funcotator_transcript_list: funcotator_transcript_list,
        funcotator_data_sources_tar_gz: funcotator_data_sources_tar_gz,

        run_collect_covered_regions: run_collect_covered_regions,
        run_collect_target_coverage: run_collect_target_coverage,
        run_collect_allelic_coverage: run_collect_allelic_coverage,
        run_contamination_model: run_contamination_model,
        run_orientation_bias_mixture_model: run_orientation_bias_mixture_model,
        run_variant_calling: run_variant_calling,
        run_variant_filter: run_variant_filter,
        run_realignment_filter: run_realignment_filter,
        run_realignment_filter_only_on_high_confidence_variants: run_realignment_filter_only_on_high_confidence_variants,
        run_collect_called_variants_allelic_coverage: run_collect_called_variants_allelic_coverage,
        run_variant_annotation: run_variant_annotation,
        run_variant_annotation_scattered: run_variant_annotation_scattered,

        keep_germline: keep_germline,
        compress_output: compress_output,
        make_bamout: make_bamout,

        preprocess_intervals_bin_length: preprocess_intervals_bin_length,
        preprocess_intervals_padding: preprocess_intervals_padding,
        min_snp_array_pop_af: min_snp_array_pop_af,
        max_snp_array_pop_af: max_snp_array_pop_af,
        min_snp_array_read_depth: min_snp_array_read_depth,
        genotype_variants_script: genotype_variants_script,
        genotype_variants_save_sample_genotype_likelihoods: genotype_variants_save_sample_genotype_likelihoods,

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