version development

import "patient.wdl" as p
import "patient.define.wdl" as p_def
import "workflow_arguments.wdl" as wfargs
import "runtime_collection.wdl" as rtc
import "tasks.wdl"
import "cnv_workflow.wdl" as cnv
import "snv_workflow.wdl" as snv
import "absolute.wdl" as abs


workflow MultiSampleSomaticWorkflow {
    input {
        # todo: add options for inputting cached output

        Patient patient = GetPatient.patient

        String individual_id
        Array[String]? tumor_sample_names
        Array[File]+ tumor_bams
        Array[File]+ tumor_bais
        Array[File]+ tumor_target_intervals
        Array[File]? tumor_annotated_target_intervals
        Array[File]? tumor_cnv_panel_of_normals
        Array[String]? normal_sample_names
        Array[File]? normal_bams
        Array[File]? normal_bais
        Array[File]? normal_target_intervals
        Array[File]? normal_annotated_target_intervals
        Array[File]? normal_cnv_panel_of_normals

        WorkflowArguments args = GetWorkflowArguments.arguments

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
        # Mutect2
        Boolean mutect2_native_pair_hmm_use_double_precision = true
        Boolean mutect2_use_linked_de_bruijn_graph = true
        Boolean mutect2_recover_all_dangling_branches = true
        Boolean mutect2_pileup_detection = true
        Boolean mutect2_genotype_germline_sites = false  # use with care! (see above)
        Int mutect2_downsampling_stride = 1  # default: 1
        Int mutect2_max_reads_per_alignment_start = 100  # default: 50
        # FilterMutectCalls
        Int filter_mutect2_max_median_fragment_length_difference = 10000  # default: 10000
        Int filter_mutect2_min_alt_median_base_quality = 20  # default: 20
        Int filter_mutect2_min_alt_median_mapping_quality = 20  # default: -1
        Int filter_mutect2_min_median_read_position = 5  # default: 1
        # FilterAlignmentArtifacts
        Int filter_alignment_artifacts_max_reasonable_fragment_length = 10000 # default: 100000
        # Funcotator
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
        String? cnn_score_variants_extra_args
        String? funcotate_extra_args

        RuntimeCollection runtime_collection = GetRTC.rtc

        Int scatter_count = 10
        # Needs docker image with bedtools, samtools, and gatk
        # todo: find smaller image. This one takes ~13 mins to spin up.
        String jupyter_docker = "us.gcr.io/broad-dsp-gcr-public/terra-jupyter-gatk"  # 27.5GB
        String gatk_docker = "broadinstitute/gatk:4.3.0.0"
        String bcftools_docker = "staphb/bcftools"
        String genotype_docker = "civisanalytics/datascience-python:latest"
        String ubuntu_docker = "ubuntu"
        File? gatk_override
        Int preemptible = 1
        Int max_retries = 1
        Int disk_sizeGB = 1

        # memory assignments in MB
        Int mem_machine_overhead = 512
        Int mem_additional_per_sample = 256  # this depends on bam size (WES vs WGS)
        Int mem_get_sample_name = 256
        Int mem_preprocess_intervals = 3072
        Int mem_split_intervals = 1024
        Int mem_collect_covered_regions = 8192
        Int mem_collect_read_counts = 2048
        Int mem_denoise_read_counts = 2048
        Int mem_mutect2_base = 3072
        Int mem_learn_read_orientation_model_base = 4096
        Int mem_get_pileup_summaries = 4096  # needs at least 2G
        Int mem_gather_pileup_summaries = 512  # 64
        Int mem_select_pileup_summaries = 512  # 64
        Int mem_harmonize_copy_ratios = 1024
        Int mem_merge_allelic_counts = 4096
        Int mem_calculate_contamination = 4096  # depends on the common_germline_alleles resource
        Int mem_genotype_variants = 4096
        Int mem_filter_mutect_calls = 4096
        Int mem_select_variants = 3072
        Int mem_filter_alignment_artifacts_base = 3072 # needs to be increased in some cases
        Int mem_merge_vcfs = 2048
        Int mem_merge_mafs = 512
        Int mem_merge_mutect_stats = 512 # 64
        Int mem_merge_bams = 8192  # wants at least 6G
        Int mem_funcotate = 6144

        # runtime assignments in minutes (for HPC cluster)
        # Derived from reasonable maximum values amongst >1000 patients.
        Int time_startup = 10
        Int time_get_sample_name = 1
        Int time_preprocess_intervals = 60
        Int time_split_intervals = 1
        Int time_collect_covered_regions = 300
        Int time_collect_read_counts = 300
        Int time_denoise_read_counts = 120
        Int time_mutect2_total = 10000  # 6 d / scatter_count
        Int time_learn_read_orientation_model = 180  # 3 h
        Int time_get_pileup_summaries = 4500  # 3 d / scatter_count
        Int time_gather_pileup_summaries = 5
        Int time_select_pileup_summaries = 5
        Int time_harmonize_copy_ratios = 1440  # 24 h
        Int time_merge_allelic_counts = 10
        Int time_calculate_contamination = 10
        Int time_genotype_variants = 30
        Int time_filter_mutect_calls = 800  # 13 h
        Int time_select_variants = 5
        Int time_filter_alignment_artifacts_total = 10000  # 12 d / scatter_count
        Int time_merge_vcfs = 10
        Int time_merge_mafs = 5
        Int time_merge_mutect_stats = 1
        Int time_merge_bams = 60
        Int time_funcotate = 1440  # 24 h

        Int cpu = 1
        Int cpu_mutect2 = 1  # good for PairHMM: 2
        Int cpu_filter_alignment_artifacts = 1  # good for PairHMM: 4
        Int cpu_harmonize_copy_ratios = 1  # 4 is worthwhile
    }

    call rtc.DefineRuntimeCollection as GetRTC {
        input:
            num_bams = length(tumor_bams) + length(select_first([normal_bams, []])),
            scatter_count = scatter_count,
            jupyter_docker = jupyter_docker,
            gatk_docker = gatk_docker,
            genotype_docker = genotype_docker,
            bcftools_docker = bcftools_docker,
            ubuntu_docker = ubuntu_docker,
            gatk_override = gatk_override,
            preemptible = preemptible,
            max_retries = max_retries,
            cpu = cpu,
            disk_sizeGB = disk_sizeGB,

            run_variant_anntation_scattered = run_variant_annotation_scattered,

            mem_machine_overhead = mem_machine_overhead,
            mem_additional_per_sample = mem_additional_per_sample,
            mem_get_sample_name = mem_get_sample_name,
            mem_preprocess_intervals = mem_preprocess_intervals,
            mem_split_intervals = mem_split_intervals,
            mem_collect_covered_regions = mem_collect_covered_regions,
            mem_collect_read_counts = mem_collect_read_counts,
            mem_denoise_read_counts = mem_denoise_read_counts,
            mem_mutect2_base = mem_mutect2_base,
            mem_learn_read_orientation_model_base = mem_learn_read_orientation_model_base,
            mem_get_pileup_summaries = mem_get_pileup_summaries,
            mem_gather_pileup_summaries = mem_gather_pileup_summaries,
            mem_select_pileup_summaries = mem_select_pileup_summaries,
            mem_harmonize_copy_ratios = mem_harmonize_copy_ratios,
            mem_merge_allelic_counts = mem_merge_allelic_counts,
            mem_calculate_contamination = mem_calculate_contamination,
            mem_genotype_variants = mem_genotype_variants,
            mem_filter_mutect_calls = mem_filter_mutect_calls,
            mem_select_variants = mem_select_variants,
            mem_filter_alignment_artifacts_base = mem_filter_alignment_artifacts_base,
            mem_merge_vcfs = mem_merge_vcfs,
            mem_merge_mafs = mem_merge_mafs,
            mem_merge_mutect_stats = mem_merge_mutect_stats,
            mem_merge_bams = mem_merge_bams,
            mem_funcotate = mem_funcotate,

            time_startup = time_startup,
            time_get_sample_name = time_get_sample_name,
            time_preprocess_intervals = time_preprocess_intervals,
            time_split_intervals = time_split_intervals,
            time_collect_covered_regions = time_collect_covered_regions,
            time_collect_read_counts = time_collect_read_counts,
            time_denoise_read_counts = time_denoise_read_counts,
            time_mutect2_total = time_mutect2_total,
            time_learn_read_orientation_model = time_learn_read_orientation_model,
            time_get_pileup_summaries = time_get_pileup_summaries,
            time_gather_pileup_summaries = time_gather_pileup_summaries,
            time_select_pileup_summaries = time_select_pileup_summaries,
            time_harmonize_copy_ratios = time_harmonize_copy_ratios,
            time_merge_allelic_counts = time_merge_allelic_counts,
            time_calculate_contamination = time_calculate_contamination,
            time_genotype_variants = time_genotype_variants,
            time_filter_mutect_calls = time_filter_mutect_calls,
            time_select_variants = time_select_variants,
            time_filter_alignment_artifacts_total = time_filter_alignment_artifacts_total,
            time_merge_vcfs = time_merge_vcfs,
            time_merge_mafs = time_merge_mafs,
            time_merge_mutect_stats = time_merge_mutect_stats,
            time_merge_bams = time_merge_bams,
            time_funcotate = time_funcotate,

            cpu_mutect2 = cpu_mutect2,
            cpu_filter_alignment_artifacts = cpu_filter_alignment_artifacts,
            cpu_harmonize_copy_ratios = cpu_harmonize_copy_ratios,
    }

    call wfargs.DefineWorkflowArguments as GetWorkflowArguments {
        input:
            scatter_count = scatter_count,
            interval_list = interval_list,
            interval_blacklist = interval_blacklist,
            interval_lists = interval_lists,

            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,

            force_call_alleles = force_call_alleles,
            force_call_alleles_idx = force_call_alleles_idx,
            snv_panel_of_normals = snv_panel_of_normals,
            snv_panel_of_normals_idx = snv_panel_of_normals_idx,
            germline_resource = germline_resource,
            germline_resource_tbi = germline_resource_tbi,
            common_germline_alleles = common_germline_alleles,
            common_germline_alleles_idx = common_germline_alleles_idx,
            realignment_bwa_mem_index_image = realignment_bwa_mem_index_image,
            funcotator_transcript_list = funcotator_transcript_list,
            funcotator_data_sources_tar_gz = funcotator_data_sources_tar_gz,

            run_collect_target_coverage = run_collect_target_coverage,
            run_collect_allelic_coverage = run_collect_allelic_coverage,
            run_contamination_model = run_contamination_model,
            run_orientation_bias_mixture_model = run_orientation_bias_mixture_model,
            run_variant_calling = run_variant_calling,
            run_variant_filter = run_variant_filter,
            run_realignment_filter = run_realignment_filter,
            run_realignment_filter_only_on_high_confidence_variants = run_realignment_filter_only_on_high_confidence_variants,
            run_collect_called_variants_allelic_coverage = run_collect_called_variants_allelic_coverage,
            run_variant_annotation = run_variant_annotation,
            run_variant_annotation_scattered = run_variant_annotation_scattered,

            keep_germline = keep_germline,
            compress_output = compress_output,
            make_bamout = make_bamout,

            preprocess_intervals_bin_length = preprocess_intervals_bin_length,
            preprocess_intervals_padding = preprocess_intervals_padding,
            min_snp_array_pop_af = min_snp_array_pop_af,
            max_snp_array_pop_af = max_snp_array_pop_af,
            min_snp_array_read_depth = min_snp_array_read_depth,
            genotype_variants_script = genotype_variants_script,
            genotype_variants_save_sample_genotype_likelihoods = genotype_variants_save_sample_genotype_likelihoods,

            min_read_depth = min_read_depth,
            mutect2_native_pair_hmm_use_double_precision = mutect2_native_pair_hmm_use_double_precision,
            mutect2_use_linked_de_bruijn_graph = mutect2_use_linked_de_bruijn_graph,
            mutect2_recover_all_dangling_branches = mutect2_recover_all_dangling_branches,
            mutect2_pileup_detection = mutect2_pileup_detection,
            mutect2_genotype_germline_sites = mutect2_genotype_germline_sites,
            mutect2_downsampling_stride = mutect2_downsampling_stride,
            mutect2_max_reads_per_alignment_start = mutect2_max_reads_per_alignment_start,
            filter_mutect2_max_median_fragment_length_difference = filter_mutect2_max_median_fragment_length_difference,
            filter_mutect2_min_alt_median_base_quality = filter_mutect2_min_alt_median_base_quality,
            filter_mutect2_min_alt_median_mapping_quality = filter_mutect2_min_alt_median_mapping_quality,
            filter_mutect2_min_median_read_position = filter_mutect2_min_median_read_position,
            filter_alignment_artifacts_max_reasonable_fragment_length = filter_alignment_artifacts_max_reasonable_fragment_length,
            funcotator_reference_version = funcotator_reference_version,
            funcotator_output_format = funcotator_output_format,
            funcotator_variant_type = funcotator_variant_type,
            funcotator_transcript_selection_mode = funcotator_transcript_selection_mode,
            funcotator_use_gnomad = funcotator_use_gnomad,
            funcotator_data_sources_paths = funcotator_data_sources_paths,
            funcotator_annotation_defaults = funcotator_annotation_defaults,
            funcotator_annotation_overrides = funcotator_annotation_overrides,
            funcotator_exclude_fields = funcotator_exclude_fields,

            split_intervals_extra_args = split_intervals_extra_args,
            getpileupsummaries_extra_args = getpileupsummaries_extra_args,
            mutect2_extra_args = mutect2_extra_args,
            filter_mutect2_extra_args = filter_mutect2_extra_args,
            select_variants_extra_args = select_variants_extra_args,
            select_low_conficence_variants_jexl_arg = select_low_conficence_variants_jexl_arg,
            realignment_extra_args = realignment_extra_args,
            funcotate_extra_args = funcotate_extra_args,

            runtime_collection = runtime_collection,
    }

    call p_def.DefinePatient as GetPatient {
        input:
            individual_id = individual_id,

            tumor_sample_names = tumor_sample_names,
            tumor_bams = tumor_bams,
            tumor_bais = tumor_bais,
            tumor_target_intervals = tumor_target_intervals,
            tumor_annotated_target_intervals = tumor_annotated_target_intervals,
            tumor_cnv_panel_of_normals = tumor_cnv_panel_of_normals,

            normal_sample_names = normal_sample_names,
            normal_bams = normal_bams,
            normal_bais = normal_bais,
            normal_target_intervals = normal_target_intervals,
            normal_annotated_target_intervals = normal_annotated_target_intervals,
            normal_cnv_panel_of_normals = normal_cnv_panel_of_normals,

            runtime_collection = runtime_collection,
    }

#    if (run_collect_covered_regions) {
#        scatter (sample in patient.samples) {
#            scatter (sequencing_run in sample.sequencing_runs) {
#                call ccr.CollectCoveredRegions {
#                    input:
#                        ref_fasta = args.ref_fasta,
#                        ref_fasta_index = args.ref_fasta_index,
#                        ref_dict = args.ref_dict,
#                        sample_name = sequencing_run.name,
#                        bam = sequencing_run.bam,
#                        bai = sequencing_run.bai,
#                        interval_list = sequencing_run.target_intervals,
#                        paired_end = sequencing_run.paired_end,
#                        min_read_depth_threshold = args.min_read_depth_threshold,
#                        output_format = "bam",
#                        runtime_collection = runtime_collection,
#                }
#            }
#        }
#    }

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

    call abs.Absolute {

    }

    output {
        # todo: flatten
#        Array[File]? covered_regions_bed = CollectCoveredRegions.regions_bed
#        Array[File?]? covered_regions_bam = CollectCoveredRegions.regions_bam
#        Array[File?]? covered_regions_bai = CollectCoveredRegions.regions_bai
#        Array[File?]? covered_regions_interval_list = CollectCoveredRegions.regions_interval_list

        File? unfiltered_vcf = SNVWorkflow.unfiltered_vcf
        File? unfiltered_vcf_idx = SNVWorkflow.unfiltered_vcf_idx
        File? mutect_stats = SNVWorkflow.mutect_stats
        File? bam = SNVWorkflow.bam
        File? bai = SNVWorkflow.bai
        File? orientation_bias = SNVWorkflow.orientation_bias
        File? filtered_vcf = SNVWorkflow.filtered_vcf
        File? filtered_vcf_idx = SNVWorkflow.filtered_vcf_idx
        File? somatic_vcf = SNVWorkflow.somatic_vcf
        File? somatic_vcf_idx = SNVWorkflow.somatic_vcf_idx
        File? germline_vcf = SNVWorkflow.germline_vcf
        File? germline_vcf_idx = SNVWorkflow.germline_vcf_idx
        File? filtering_stats = SNVWorkflow.filtering_stats
        Array[File?]? called_germline_allelic_counts = SNVWorkflow.called_germline_allelic_counts
        Array[File?]? called_somatic_allelic_counts = SNVWorkflow.called_somatic_allelic_counts
        Array[File]? annotated_variants = SNVWorkflow.annotated_variants
        Array[File?]? annotated_variants_idx = SNVWorkflow.annotated_variants_idx

        File? genotyped_snparray_vcf = CNVWorkflow.genotyped_snparray_vcf
        File? genotyped_snparray_vcf_idx = CNVWorkflow.genotyped_snparray_vcf_idx
        File? snparray_ref_counts = CNVWorkflow.snparray_ref_counts
        File? snparray_alt_counts = CNVWorkflow.snparray_alt_counts
        File? snparray_other_alt_counts = CNVWorkflow.snparray_other_alt_counts
        File? sample_snp_correlation = CNVWorkflow.sample_snp_correlation
        Array[File]? sample_snparray_genotype_likelihoods = CNVWorkflow.sample_snparray_genotype_likelihoods
        Array[File]? snparray_allelic_counts = CNVWorkflow.snparray_allelic_counts
        Array[File]? contamination_table = CNVWorkflow.contamination_tables
        Array[File]? segmentation_table = CNVWorkflow.segmentation_tables
        Array[File?]? target_read_counts = CNVWorkflow.target_read_counts
        Array[File?]? denoised_copy_ratios = CNVWorkflow.denoised_copy_ratios
        Array[File?]? standardized_copy_ratios = CNVWorkflow.standardized_copy_ratios
    }
}
