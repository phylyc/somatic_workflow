version development

import "runtimes.wdl"
import "tasks.wdl"
import "patient.wdl"
import "calculate_contamination.wdl" as cc
import "collect_allelic_counts.wdl" as cac
import "collect_covered_regions.wdl" as ccr
import "collect_read_counts.wdl" as crc
#import "genotype_variants.wdl" as gv
#import "phase_variants.wdl" as ph
import "call_variants.wdl" as cv
import "filter_variants.wdl" as fv
import "annotate_variants.wdl" as av


workflow MultiSampleSomaticWorkflow {
    input {
        File? interval_list
        File? interval_blacklist  # For CalculateContamination
        Array[File]? interval_lists
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        String individual_id
        Array[File]+ tumor_bams
        Array[File]+ tumor_bais
        Array[File]+ tumor_target_intervals
        Array[String]+ tumor_sample_names
        Array[File]? normal_bams
        Array[File]? normal_bais
        Array[File]? normal_target_intervals
        Array[String]? normal_sample_names

        # resources
        File? force_call_alleles
        File? force_call_alleles_idx
        File? panel_of_normals
        File? panel_of_normals_idx
        File? read_count_panel_of_normals
        File? germline_resource
        File? germline_resource_tbi
        File? variants_for_contamination
        File? variants_for_contamination_idx
        File? bwa_mem_index_image
        File? funcotator_transcript_list
        File? funcotator_data_sources_tar_gz

        # workflow options
        Boolean run_collect_coverage = false
        Boolean run_contamination_model = true
        Boolean run_phase_hets = false
        Boolean run_orientation_bias_mixture_model = true
        Boolean run_variant_filter = true
        Boolean run_realignment_filter = true
        Boolean run_realignment_filter_only_on_high_confidence_variants = false
        Boolean run_final_pileup_summaries = true
        Boolean run_annotate_variants = true

        Boolean keep_germline = true
        Boolean compress_output = true
        Boolean make_bamout = false

        # arguments
        Int preprocess_intervals_bin_length = 0
        Int preprocess_intervals_padding = 0
        Int min_read_depth_threshold = 1
        Boolean mutect2_native_pair_hmm_use_double_precision = true
        Boolean mutect2_use_linked_de_bruijn_graph = true
        Boolean mutect2_recover_all_dangling_branches = true
        Boolean mutect2_pileup_detection = true
        Boolean mutect2_genotype_germline_sites = false  # use with care! (see above)
        Int mutect2_downsampling_stride = 1
        Int mutect2_max_reads_per_alignment_start = 50
        Int filter_mutect2_max_median_fragment_length_difference = 10000  # default: 10000
        Int filter_mutect2_min_alt_median_base_quality = 20  # default: 20
        Int filter_mutect2_min_alt_median_mapping_quality = 20  # default: -1
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
        String? cnn_score_variants_extra_args
        String? funcotate_extra_args

        # runtime
        Int scatter_count = 10
        # Needs docker image with bedtools, samtools, and gatk
        # todo: find smaller image. This one takes ~13 mins to spin up.
        String jupyter_docker = "us.gcr.io/broad-dsp-gcr-public/terra-jupyter-gatk"  # 27.5GB
        String gatk_docker = "broadinstitute/gatk"
        String bcf_tools_docker = "dceoy/bcftools"
        String ubuntu_docker = "ubuntu"
        File? gatk_override
        Int preemptible = 1
        Int max_retries = 1
        Int disk_sizeGB = 1

        # memory assignments in MB
        Int mem_additional_per_sample = 256  # this depends on bam size (WES vs WGS)
        Int mem_get_sample_name = 512
        Int mem_preprocess_intervals = 2048
        Int mem_split_intervals = 512
        Int mem_collect_covered_regions = 8192
        Int mem_collect_read_counts = 2048
        Int mem_denoise_read_counts = 2048
        Int mem_variant_call_base = 4096
        Int mem_learn_read_orientation_model_base = 6144
        Int mem_get_pileup_summaries = 4096  # needs at least 2G
        Int mem_gather_pileup_summaries = 512  # 64
        Int mem_select_pileup_summaries = 256  # 64
        Int mem_calculate_contamination = 8192  # depends on the variants_for_contamination resource
        Int mem_filter_mutect_calls = 4096
        Int mem_select_variants = 2048
        Int mem_filter_alignment_artifacts_base = 2048  # needs to be increased in some cases
        Int mem_merge_vcfs = 512
        Int mem_merge_mutect_stats = 512 # 64
        Int mem_merge_bams = 8192  # wants at least 6G
        Int mem_cnn_scoring = 4096
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
        Int time_variant_call_total = 10000  # 6 d / scatter_count
        Int time_learn_read_orientation_model = 180  # 3 h
        Int time_get_pileup_summaries = 4500  # 3 d / scatter_count
        Int time_gather_pileup_summaries = 5
        Int time_select_pileup_summaries = 5
        Int time_calculate_contamination = 10
        Int time_filter_mutect_calls = 800  # 13 h
        Int time_select_variants = 5
        Int time_filter_alignment_artifacts_total = 10000  # 12 d / scatter_count
        Int time_merge_vcfs = 10
        Int time_merge_mutect_stats = 1
        Int time_merge_bams = 60
        Int time_cnn_scoring = 10
        Int time_funcotate = 500  # 8 h

        # Increasing cpus likely increases costs by the same factor.
        Int cpu = 1
        Int cpu_variant_call = 1  # good for PairHMM: 2
        Int cpu_filter_alignment_artifacts = 1  # good for PairHMM: 4
        Int cpu_cnn_scoring = 1
    }

    Array[File] non_optional_normal_bams = select_first([normal_bams, []])
    Array[File] non_optional_normal_bais = select_first([normal_bais, []])
    Array[File] non_optional_normal_target_intervals = select_first([normal_target_intervals, []])

    Boolean normal_is_present = defined(normal_bams) && (length(non_optional_normal_bams) > 0)

    call runtimes.DefineRuntimes as Runtimes {
        input:
            num_bams = length(tumor_bams) + length(non_optional_normal_bams),
            scatter_count = scatter_count,
            jupyter_docker = jupyter_docker,
            gatk_docker = gatk_docker,
            bcf_tools_docker = bcf_tools_docker,
            ubuntu_docker = ubuntu_docker,
            gatk_override = gatk_override,
            preemptible = preemptible,
            max_retries = max_retries,
            cpu = cpu,
            disk_sizeGB = disk_sizeGB,

            mem_additional_per_sample = mem_additional_per_sample,
            mem_get_sample_name = mem_get_sample_name,
            mem_preprocess_intervals = mem_preprocess_intervals,
            mem_split_intervals = mem_split_intervals,
            mem_collect_covered_regions = mem_collect_covered_regions,
            mem_collect_read_counts = mem_collect_read_counts,
            mem_denoise_read_counts = mem_denoise_read_counts,
            mem_variant_call_base = mem_variant_call_base,
            mem_learn_read_orientation_model_base = mem_learn_read_orientation_model_base,
            mem_get_pileup_summaries = mem_get_pileup_summaries,
            mem_gather_pileup_summaries = mem_gather_pileup_summaries,
            mem_select_pileup_summaries = mem_select_pileup_summaries,
            mem_calculate_contamination = mem_calculate_contamination,
            mem_filter_mutect_calls = mem_filter_mutect_calls,
            mem_select_variants = mem_select_variants,
            mem_filter_alignment_artifacts_base = mem_filter_alignment_artifacts_base,
            mem_merge_vcfs = mem_merge_vcfs,
            mem_merge_mutect_stats = mem_merge_mutect_stats,
            mem_merge_bams = mem_merge_bams,
            mem_cnn_scoring = mem_cnn_scoring,
            mem_funcotate = mem_funcotate,

            time_startup = time_startup,
            time_get_sample_name = time_get_sample_name,
            time_preprocess_intervals = time_preprocess_intervals,
            time_split_intervals = time_split_intervals,
            time_collect_covered_regions = time_collect_covered_regions,
            time_collect_read_counts = time_collect_read_counts,
            time_denoise_read_counts = time_denoise_read_counts,
            time_variant_call_total = time_variant_call_total,
            time_learn_read_orientation_model = time_learn_read_orientation_model,
            time_get_pileup_summaries = time_get_pileup_summaries,
            time_gather_pileup_summaries = time_gather_pileup_summaries,
            time_select_pileup_summaries = time_select_pileup_summaries,
            time_calculate_contamination = time_calculate_contamination,
            time_filter_mutect_calls = time_filter_mutect_calls,
            time_select_variants = time_select_variants,
            time_filter_alignment_artifacts_total = time_filter_alignment_artifacts_total,
            time_merge_vcfs = time_merge_vcfs,
            time_merge_mutect_stats = time_merge_mutect_stats,
            time_merge_bams = time_merge_bams,
            time_cnn_scoring = time_cnn_scoring,
            time_funcotate = time_funcotate,

            cpu_variant_call = cpu_variant_call,
            cpu_filter_alignment_artifacts = cpu_filter_alignment_artifacts,
            cpu_cnn_scoring = cpu_cnn_scoring
    }

    scatter (tumor_bam in tumor_bams) {
        call tasks.GetSampleName as GetTumorSampleName {
            input:
                bam = tumor_bam,
                runtime_params = Runtimes.get_sample_name_runtime,
        }
    }

    if (normal_is_present) {
        scatter (normal_bam in non_optional_normal_bams) {
            call tasks.GetSampleName as GetNormalSampleName {
                input:
                    bam = normal_bam,
                    runtime_params = Runtimes.get_sample_name_runtime,
            }
        }
    }

    call patient.Patient {
        input:
            individual_id = individual_id,
            tumor_bams = tumor_bams,
            tumor_bais = tumor_bais,
            tumor_target_intervals = tumor_target_intervals,
            tumor_bam_names = GetTumorSampleName.sample_name,
            tumor_sample_names = tumor_sample_names,
            normal_bams = normal_bams,
            normal_bais = normal_bais,
            normal_target_intervals = normal_target_intervals,
            normal_bam_names = GetNormalSampleName.sample_name,
            normal_sample_names = normal_sample_names,
    }

    call tasks.PreprocessIntervals {
        input:
            interval_list = interval_list,
            interval_lists = interval_lists,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            bin_length = preprocess_intervals_bin_length,
            padding = preprocess_intervals_padding,
            runtime_params = Runtimes.preprocess_intervals_runtime,
    }

    call tasks.SplitIntervals {
    	input:
            interval_list = PreprocessIntervals.preprocessed_interval_list,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            scatter_count = scatter_count,
            split_intervals_extra_args = split_intervals_extra_args,
            runtime_params = Runtimes.split_intervals_runtime,
    }

    if (run_collect_coverage) {
        scatter (sample in Patient.samples) {
            call ccr.CollectCoveredRegions {
                input:
                    interval_list = sample.target_intervals,
                    ref_fasta = ref_fasta,
                    ref_fasta_index = ref_fasta_index,
                    ref_dict = ref_dict,
                    bam = sample.bam,
                    bai = sample.bai,
                    sample_name = sample.assigned_sample_name,
                    min_read_depth_threshold = min_read_depth_threshold,
                    paired_end = sample.paired_end,
                    output_format = "bam",
                    collect_covered_regions_runtime = Runtimes.collect_covered_regions_runtime,
            }

            call crc.CollectReadCounts {
                input:
                    interval_list = sample.target_intervals,
                    ref_fasta = ref_fasta,
                    ref_fasta_index = ref_fasta_index,
                    ref_dict = ref_dict,
                    bam = sample.bam,
                    bai = sample.bai,
                    sample_name = sample.assigned_sample_name,
                    # annotated_interval_list = annotated_interval_list,
                    read_count_panel_of_normals = read_count_panel_of_normals,

                    collect_read_counts_runtime = Runtimes.collect_read_counts_runtime,
                    denoise_read_counts_runtime = Runtimes.denoise_read_counts_runtime,
            }
        }
    }

    if (run_contamination_model) {
        scatter (normal in Patient.normal_samples) {
            call cc.CalculateContamination as CalculateNormalContamination {
                input:
                    interval_list = normal.target_intervals,
                    interval_blacklist = interval_blacklist,
                    scattered_interval_list = SplitIntervals.interval_files,
                    ref_dict = ref_dict,
                    tumor_bam = normal.bam,
                    tumor_bai = normal.bai,
                    tumor_sample_name = normal.assigned_sample_name,
                    variants = variants_for_contamination,
                    variants_idx = variants_for_contamination_idx,
                    getpileupsummaries_extra_args = getpileupsummaries_extra_args,
                    get_pileup_summaries_runtime = Runtimes.get_pileup_summaries_runtime,
                    gather_pileup_summaries_runtime = Runtimes.gather_pileup_summaries_runtime,
                    select_pileup_summaries_runtime = Runtimes.select_pileup_summaries_runtime,
                    calculate_contamination_runtime = Runtimes.calculate_contamination_runtime
            }
        }

        # If multiple normals are present, it is not that important which of them
        # to choose for the contamination workflow.
        # todo: Choose the normal with the greatest sequencing depth.
        if (Patient.has_normal) {
            File? normal_pileups = select_first(CalculateNormalContamination.tumor_pileup_summaries)
        }

        scatter (tumor in Patient.tumor_samples) {
            call cc.CalculateContamination as CalculateTumorContamination {
                input:
                    interval_list = tumor.target_intervals,
                    interval_blacklist = interval_blacklist,
                    scattered_interval_list = SplitIntervals.interval_files,
                    ref_dict = ref_dict,
                    tumor_bam = tumor.bam,
                    tumor_bai = tumor.bai,
                    tumor_sample_name = tumor.assigned_sample_name,
                    normal_pileups = normal_pileups,
                    variants = variants_for_contamination,
                    variants_idx = variants_for_contamination_idx,
                    getpileupsummaries_extra_args = getpileupsummaries_extra_args,
                    get_pileup_summaries_runtime = Runtimes.get_pileup_summaries_runtime,
                    gather_pileup_summaries_runtime = Runtimes.gather_pileup_summaries_runtime,
                    select_pileup_summaries_runtime = Runtimes.select_pileup_summaries_runtime,
                    calculate_contamination_runtime = Runtimes.calculate_contamination_runtime
            }
        }

        Array[File] panel_pileup_summaries = flatten([
            CalculateTumorContamination.tumor_pileup_summaries,
            select_first([CalculateNormalContamination.tumor_pileup_summaries, []])
        ])

        Array[File] contamination_tables = flatten([
            CalculateTumorContamination.contamination_table,
            select_first([CalculateNormalContamination.contamination_table, []])
        ])
        Array[File] segmentation_tables = flatten([
            CalculateTumorContamination.segmentation,
            select_first([CalculateNormalContamination.segmentation, []])
        ])

#        call gv.GenotypeVariants {
#            input:
#                pileups = panel_pileup_summaries,
#                contamination_tables = contamination_tables,
#                segmentation_tables = segmentation_tables,
#                runtime_params = standard_runtime,
#        }
#
#        call SelectVariants as SelectHETs {
#            input:
#                input_vcf = GenotypeVariants.vcf,
#                input_vcf_idx = GenotypeVariants.vcf_idx,
#                select_variants_extra_args = select_variants_extra_args,
#                runtime_params = standard_runtime,
#                mem_select_variants = mem_select_variants,
#                time_select_variants = time_startup + time_select_variants
#        }
#
#        if (run_phase_hets) {
#            call pv.PhaseVariants as PhaseHETs {
#                input:
#                    bams = Patient.bams,
#                    bais = Patient.bais,
#                    vcf = SelectHETs.vcf,
#                    vcf_idx = SelectHETs.vcf_idx,
#            }
#        }
    }

    call cv.CallVariants {
        input:
            scattered_interval_list = SplitIntervals.interval_files,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            individual_id = individual_id,
            tumor_bams = tumor_bams,
            tumor_bais = tumor_bais,
            normal_bams = normal_bams,
            normal_bais = normal_bais,
            normal_bam_names = GetNormalSampleName.sample_name,

            force_call_alleles = force_call_alleles,
            force_call_alleles_idx = force_call_alleles_idx,
            panel_of_normals = panel_of_normals,
            panel_of_normals_idx = panel_of_normals_idx,
            germline_resource = germline_resource,
            germline_resource_tbi = germline_resource_tbi,
            make_bamout = make_bamout,
            run_orientation_bias_mixture_model = run_orientation_bias_mixture_model,
            compress_output = compress_output,

            genotype_germline_sites = mutect2_genotype_germline_sites,
            native_pair_hmm_use_double_precision = mutect2_native_pair_hmm_use_double_precision,
            use_linked_de_bruijn_graph = mutect2_use_linked_de_bruijn_graph,
            recover_all_dangling_branches = mutect2_recover_all_dangling_branches,
            pileup_detection = mutect2_pileup_detection,
            downsampling_stride = mutect2_downsampling_stride,
            max_reads_per_alignment_start = mutect2_max_reads_per_alignment_start,
            mutect2_extra_args = mutect2_extra_args,

            variant_call_runtime = Runtimes.variant_call_runtime,
            merge_vcfs_runtime = Runtimes.merge_vcfs_runtime,
            merge_mutect_stats_runtime = Runtimes.merge_mutect_stats_runtime,
            merge_bams_runtime = Runtimes.merge_bams_runtime,
            learn_read_orientation_model_runtime = Runtimes.learn_read_orientation_model_runtime,
    }

    if (run_variant_filter) {
        call fv.FilterVariants {
            input:
                scattered_interval_list = SplitIntervals.interval_files,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                ref_dict = ref_dict,
                tumor_bams = tumor_bams,
                tumor_bais = tumor_bais,
                vcf = CallVariants.vcf,
                vcf_idx = CallVariants.vcf_idx,
                mutect_stats = CallVariants.mutect_stats,
                orientation_bias = CallVariants.orientation_bias,
                contamination_tables = contamination_tables,
                segmentation_tables = segmentation_tables,
                bwa_mem_index_image = bwa_mem_index_image,
                run_realignment_filter = run_realignment_filter,
                run_realignment_filter_only_on_high_confidence_variants = run_realignment_filter_only_on_high_confidence_variants,
                keep_germline = keep_germline,
                compress_output = compress_output,
                max_median_fragment_length_difference = filter_mutect2_max_median_fragment_length_difference,
                min_alt_median_base_quality = filter_mutect2_min_alt_median_base_quality,
                min_alt_median_mapping_quality = filter_mutect2_min_alt_median_mapping_quality,
                max_reasonable_fragment_length = filter_alignment_artifacts_max_reasonable_fragment_length,
                filter_mutect2_extra_args = filter_mutect2_extra_args,
                select_variants_extra_args = select_variants_extra_args,
                select_low_conficence_variants_jexl_arg = select_low_conficence_variants_jexl_arg,
                realignment_extra_args = realignment_extra_args,

                filter_mutect_calls_runtime = Runtimes.filter_mutect_calls_runtime,
                filter_alignment_artifacts_runtime = Runtimes.filter_alignment_artifacts_runtime,
                select_variants_runtime = Runtimes.select_variants_runtime,
                merge_vcfs_runtime = Runtimes.merge_vcfs_runtime,
        }

        if (run_final_pileup_summaries) {
            scatter (sample in Patient.samples) {
                call cac.CollectAllelicCounts as GermlineAllelicCounts {
                    input:
                        scattered_interval_list = SplitIntervals.interval_files,
                        bam = sample.bam,
                        bai = sample.bai,
                        ref_dict = ref_dict,
                        vcf = FilterVariants.germline_vcf,
                        vcf_idx = FilterVariants.germline_vcf_idx,
                        getpileupsummaries_extra_args = getpileupsummaries_extra_args,
                        output_base_name = sample.assigned_sample_name + ".germline",
                        vcf_to_pileup_variants_runtime = Runtimes.vcf_to_pileup_variants_runtime,
                        get_pileup_summaries_runtime = Runtimes.get_pileup_summaries_runtime,
                        gather_pileup_summaries_runtime = Runtimes.gather_pileup_summaries_runtime,
                }

                call cac.CollectAllelicCounts as SomaticAllelicCounts {
                    input:
                        scattered_interval_list = SplitIntervals.interval_files,
                        bam = sample.bam,
                        bai = sample.bai,
                        ref_dict = ref_dict,
                        vcf = FilterVariants.somatic_vcf,
                        vcf_idx = FilterVariants.somatic_vcf_idx,
                        getpileupsummaries_extra_args = getpileupsummaries_extra_args,
                        output_base_name = sample.assigned_sample_name + ".somatic",
                        vcf_to_pileup_variants_runtime = Runtimes.vcf_to_pileup_variants_runtime,
                        get_pileup_summaries_runtime = Runtimes.get_pileup_summaries_runtime,
                        gather_pileup_summaries_runtime = Runtimes.gather_pileup_summaries_runtime,
                }
            }
        }
    }

    if (run_annotate_variants) {
        call av.AnnotateVariants {
            input:
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                ref_dict = ref_dict,
                interval_list = PreprocessIntervals.preprocessed_interval_list,
                vcf = select_first([FilterVariants.filtered_vcf, CallVariants.vcf]),
                vcf_idx = select_first([FilterVariants.filtered_vcf_idx, CallVariants.vcf_idx]),
                individual_id = individual_id,
                samples = Patient.samples,

                reference_version = funcotator_reference_version,
                output_format = funcotator_output_format,
                variant_type = funcotator_variant_type,
                transcript_selection_mode = funcotator_transcript_selection_mode,
                transcript_list = funcotator_transcript_list,
                data_sources_tar_gz = funcotator_data_sources_tar_gz,
                use_gnomad = funcotator_use_gnomad,
                compress_output = compress_output,
                data_sources_paths = funcotator_data_sources_paths,
                annotation_defaults = funcotator_annotation_defaults,
                annotation_overrides = funcotator_annotation_overrides,
                exclude_fields = funcotator_exclude_fields,
                select_variants_extra_args = select_variants_extra_args,
                funcotate_extra_args = funcotate_extra_args,

                select_variants_runtime = Runtimes.select_variants_runtime,
                funcotate_runtime = Runtimes.funcotate_runtime,
        }
    }

    output {
        File unfiltered_vcf = CallVariants.vcf
        File unfiltered_vcf_idx = CallVariants.vcf_idx
        File mutect_stats = CallVariants.mutect_stats
        File? bam = CallVariants.bam
        File? bai = CallVariants.bai
        File? orientation_bias = CallVariants.orientation_bias
        File? filtered_vcf = FilterVariants.filtered_vcf
        File? filtered_vcf_idx = FilterVariants.filtered_vcf_idx
        File? somatic_vcf = FilterVariants.somatic_vcf
        File? somatic_vcf_idx = FilterVariants.somatic_vcf_idx
        File? germline_vcf = FilterVariants.germline_vcf
        File? germline_vcf_idx = FilterVariants.germline_vcf_idx
        File? filtering_stats = FilterVariants.filtering_stats
        Array[File]? annotated_variants = AnnotateVariants.annotated_variants
        Array[File?]? annotated_variants_idx = AnnotateVariants.annotated_variants_idx
        Array[File]? somatic_allelic_counts = SomaticAllelicCounts.pileup_summaries
        Array[File]? germline_allelic_counts = GermlineAllelicCounts.pileup_summaries
        Array[File]? panel_allelic_counts = panel_pileup_summaries
        Array[File]? contamination_table = contamination_tables
        Array[File]? segmentation_table = segmentation_tables
        Array[File]? read_counts = CollectReadCounts.read_counts
        Array[File?]? denoised_copy_ratios = CollectReadCounts.denoised_copy_ratios
        Array[File?]? standardized_copy_ratios = CollectReadCounts.standardized_copy_ratios
        Array[File]? covered_regions_bed = CollectCoveredRegions.regions_bed
        Array[File?]? covered_regions_bam = CollectCoveredRegions.regions_bam
        Array[File?]? covered_regions_bai = CollectCoveredRegions.regions_bai
        Array[File?]? covered_regions_interval_list = CollectCoveredRegions.regions_interval_list
    }
}
