version development

import "runtimes.wdl" as rt


struct RuntimeCollection {
    Runtime get_tumor_sample_names
    Runtime get_sample_name
    Runtime annotate_intervals
    Runtime preprocess_intervals
    Runtime split_intervals
    Runtime collect_covered_regions
    Runtime collect_read_counts
    Runtime denoise_read_counts
    Runtime vcf_to_pileup_variants
    Runtime get_pileup_summaries
    Runtime gather_pileup_summaries
    Runtime select_pileup_summaries
    Runtime pileup_to_allelic_counts
    Runtime harmonize_copy_ratios
    Runtime merge_allelic_counts
    Runtime calculate_contamination
    Runtime genotype_variants
    Runtime model_segments
    Runtime call_copy_ratio_segments
    Runtime merge_calls_with_modeled_segments
    Runtime plot_modeled_segments
    Runtime filter_germline_cnvs
    Runtime recount_markers
    Runtime model_segments_to_acs_conversion
    Runtime process_maf_for_absolute
    Runtime absolute
    Runtime absolute_extract
    Runtime mutect2
    Runtime learn_read_orientation_model
    Runtime merge_vcfs
    Runtime merge_mafs
    Runtime merge_mutect_stats
    Runtime print_reads
    Runtime filter_mutect_calls
    Runtime variant_filtration
    Runtime left_align_and_trim_variants
    Runtime filter_alignment_artifacts
    Runtime select_variants
    Runtime funcotate
    Runtime create_empty_annotation
    Runtime create_cnv_panel
    Runtime create_mutect2_panel
    Runtime select_af_only_from_vcf
}


workflow DefineRuntimeCollection {
    input {
        Int num_bams = 1
        Int bam_size = 0

        Int scatter_count = 10
        String gatk_docker = "broadinstitute/gatk:4.6.1.0"
        # Needs docker image with bedtools, samtools, and gatk
        String jupyter_docker = "us.gcr.io/broad-dsp-gcr-public/terra-jupyter-gatk"  # 27.5GB todo: find smaller image. This one takes ~13 mins to spin up.
        String absolute_docker = "phylyc/absolute:1.6"
        String ubuntu_docker = "ubuntu"
        String bcftools_docker = "staphb/bcftools:1.21"  # @sha256:176f4c7c10e57c8c3e2d26f0f105bd680e9ddff65c9e20dd4d3ebff228f17188
        String python_docker = "civisanalytics/datascience-python:8.0.1"  # @sha256:3482b19792546214a6952b369472c9d4d50d60b3a38300127ce346b7bab5fd51
        File? gatk_override
        Int preemptible = 1
        Int max_retries = 1
        Int cpu = 1
        Int disk_sizeGB = 1
        Int boot_disk_size = 12  # needs to be > 10

        # memory assignments in MB
        # runtime assignments in minutes (for HPC cluster)
        # Derived from reasonable maximum values amongst >1000 patients.
        # Increasing cpus likely increases costs by the same factor.

        Int mem_machine_overhead = 512

        Int time_startup = 10

        #######################################################################
        ### Preprocessing
        #######################################################################

        # GetTumorSampleNames
        Int mem_get_tumor_sample_names = 256
        Int time_get_tumor_sample_names = 1

        # gatk: GetSampleName
        Int mem_get_sample_name = 256
        Int time_get_sample_name = 1

        # gatk: AnnotateIntervals
        Int mem_annotate_intervals = 2048
        Int time_annotate_intervals = 60

        # gatk: PreprocessIntervals
        Int mem_preprocess_intervals = 3072
        Int time_preprocess_intervals = 60

        # gatk: SplitIntervals
        Int mem_split_intervals = 1024
        Int time_split_intervals = 1

        #######################################################################
        # CNV workflow
        #######################################################################

        # CollectCoveredRegions
        Int mem_collect_covered_regions = 8192
        Int time_collect_covered_regions = 300

        # gatk: CollectReadCounts
        Int mem_collect_read_counts = 2048
        Int time_collect_read_counts = 300

        # gatk: DenoiseReadCounts
        Int mem_denoise_read_counts = 2048
        Int time_denoise_read_counts = 120

        # VcfToPileupVariants
        Int mem_vcf_to_pileup_variants = 512
        Int time_vcf_to_pileup_variants = 5

        # gatk: GetPileupSummaries
        Int mem_get_pileup_summaries = 2560  # needs at least 2G
        Int time_get_pileup_summaries = 4500  # 3 d / scatter_count

        # gatk: GatherPileupSummaries
        Int mem_gather_pileup_summaries = 512  # 64
        Int time_gather_pileup_summaries = 5

        # SelectPileupSummaries
        Int mem_select_pileup_summaries = 512
        Int time_select_pileup_summaries = 5

        # PileupToAllelicCounts
        Int mem_pileup_to_allelic_counts = 2048
        Int time_pileup_to_allelic_counts = 5

        # HarmonizeCopyRatios
        Int cpu_harmonize_copy_ratios = 4
        Int mem_harmonize_copy_ratios_base = 8192
        Int mem_harmonize_copy_ratios_additional_per_sample = 1280
        Int time_harmonize_copy_ratios = 60

        # MergeAllelicCounts
        Int mem_merge_allelic_counts = 4096
        Int time_merge_allelic_counts = 10

        # gatk: CalculateContamination
        # Memory depends on size of SNP array; gnomad v2.1.1 in WES target regions gives ~50k variants, which uses ~120 MB
        Int mem_calculate_contamination = 3072
        Int time_calculate_contamination = 10

        # custom genotyping script based on CalculateContamination model
        Int mem_genotype_variants = 8192
        Int time_genotype_variants = 30

        # gatk: ModelSegments
        Int mem_model_segments_base = 1024
        Int mem_model_segments_additional_per_sample = 128
        Int time_model_segments = 60

        # gatk: CallCopyRatioSegments
        Int mem_call_copy_ratio_segments = 1024
        Int time_call_copy_ratio_segments = 10

        # custom task to merge calls with modeled segments
        Int mem_merge_calls_with_modeled_segments = 512
        Int time_merge_calls_with_modeled_segments = 1

        # gatk: PlotModeledSegments
        Int mem_plot_modeled_segments = 1024
        Int time_plot_modeled_segments = 10

        # custom filtering script
        Int mem_filter_germline_cnvs = 2048
        Int time_filter_germline_cnvs = 10

        # custom recount marker script
        Int mem_recount_markers = 2048
        Int time_recount_markers = 10

        #######################################################################
        ### SNV workflow
        #######################################################################

        # gatk: Mutect2
        # The GATK only parallelizes a few parts of the computation, so any extra cores would be idle for a large fraction of time.
        Int cpu_mutect2 = 1
        Int mem_mutect2_base = 3072
        Int mem_mutect2_additional_per_sample = 512
        Int mem_mutect2_overhead = 1024  # needs to be at least 1GB to run decently
        Int time_mutect2_total = 10000  # 6 d / scatter_count
        Int preemptible_mutect2 = 1
        Int max_retries_mutect2 = 2
        Int disk_mutect2_total = bam_size

        # gatk: MergeVCFs
        Int mem_merge_vcfs = 2048
        Int time_merge_vcfs = 10

        # gatk: MergeMAFs
        Int mem_merge_mafs = 512
        Int time_merge_mafs = 5

        # gatk: MergeMutectStats
        Int mem_merge_mutect_stats = 512 # 64
        Int time_merge_mutect_stats = 1

        # gatk: PrintReads
        Int mem_print_reads = 32768
        Int time_print_reads = 60
        Int disk_print_reads = bam_size

        # gatk: LearnReadOrientationModel
        Int mem_learn_read_orientation_model_base = 8192
        Int mem_learn_read_orientation_model_additional_per_sample = 1024
        Int time_learn_read_orientation_model = 180  # 3 h

        # gatk: FilterMutectCalls
        Int mem_filter_mutect_calls = 1024
        Int time_filter_mutect_calls = 800  # 13 h

        # gatk: VariantFiltration
        Int mem_variant_filtration = 512
        Int time_variant_filtration = 5

        # gatk: LeftAlignAndTrimVariants
        Int mem_left_align_and_trim_variants = 1024
        Int time_left_align_and_trim_variants = 60

        # gatk: FilterAlignmentArtifacts
        Int cpu_filter_alignment_artifacts = 1
        Int mem_filter_alignment_artifacts_base = 1024
        Int mem_filter_alignment_artifacts_additional_per_sample = 192
        Int time_filter_alignment_artifacts_total = 10000  # 12 d / scatter_count

        # gatk: SelectVariants
        Int mem_select_variants = 3072
        Int time_select_variants = 5

        # gatk: Funcotator
        Int mem_funcotate = 6144
        Int time_funcotate = 1440  # 24 h

        # CreateEmptyAnnotations
        Int mem_create_empty_annotations = 128
        Int time_create_empty_annotations = 1


        #######################################################################
        ### Clonal Analysis workflow
        #######################################################################

        # ModelSegmentsToACSConversion
        Int mem_model_segments_to_acs_conversion = 1024
        Int time_model_segments_to_acs_conversion = 10

        # ProcessMafForAbsolute
        Int mem_process_maf_for_absolute = 2048
        Int time_process_maf_for_absolute = 10

        # Absolute
        Int mem_absolute = 6144
        Int time_absolute = 60

        #######################################################################
        ### Assorted
        #######################################################################

        # gatk: CreateReadCountPanelOfNormals
        Int mem_create_cnv_panel = 16384
        Int time_create_cnv_panel = 1200  # 20 h
        Int disk_create_cnv_panel = 10

        # gatk: GenomicsDBImport / CreateSomaticPanelOfNormals
        Int mem_create_mutect2_panel = 16384
        Int time_create_mutect2_panel = 1200  # 20 h
        Int disk_create_mutect2_panel = 10

        Int mem_select_af_only_from_vcf = 2024
        Int time_select_af_only_from_vcf = 1440
    }

    Int gatk_override_size = ceil(size(gatk_override, "GB"))
    Int disk = 2 + gatk_override_size + disk_sizeGB

    Runtime get_tumor_sample_names = {
        "docker": ubuntu_docker,
        "jar_override": gatk_override,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu,
        "machine_mem": mem_get_tumor_sample_names + mem_machine_overhead,
        "command_mem": mem_get_tumor_sample_names,
        "runtime_minutes": time_startup + time_get_tumor_sample_names,
        "disk": disk,
        "boot_disk_size": boot_disk_size
    }

    Runtime get_sample_name = {
        "docker": gatk_docker,
        "jar_override": gatk_override,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu,
        "machine_mem": mem_get_sample_name + mem_machine_overhead,
        "command_mem": mem_get_sample_name,
        "runtime_minutes": time_startup + time_get_sample_name,
        "disk": disk,
        "boot_disk_size": boot_disk_size
    }

    Runtime annotate_intervals = {
        "docker": gatk_docker,
        "jar_override": gatk_override,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu,
        "machine_mem": mem_annotate_intervals + mem_machine_overhead,
        "command_mem": mem_annotate_intervals,
        "runtime_minutes": time_startup + time_annotate_intervals,
        "disk": disk,
        "boot_disk_size": boot_disk_size
    }

    Runtime preprocess_intervals = {
        "docker": gatk_docker,
        "jar_override": gatk_override,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu,
        "machine_mem": mem_preprocess_intervals + mem_machine_overhead,
        "command_mem": mem_preprocess_intervals,
        "runtime_minutes": time_startup + time_preprocess_intervals,
        "disk": disk,
        "boot_disk_size": boot_disk_size
    }

    Runtime split_intervals = {
        "docker": gatk_docker,
        "jar_override": gatk_override,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu,
        "machine_mem": mem_split_intervals + mem_machine_overhead,
        "command_mem": mem_split_intervals,
        "runtime_minutes": time_startup + time_split_intervals,
        "disk": disk,
        "boot_disk_size": boot_disk_size
    }

    Runtime collect_covered_regions = {
        "docker": jupyter_docker,
        "jar_override": gatk_override,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu,
        "machine_mem": mem_collect_covered_regions + mem_machine_overhead,
        "command_mem": mem_collect_covered_regions,
        "runtime_minutes": time_startup + time_collect_covered_regions,
        "disk": disk,
        "boot_disk_size": boot_disk_size
    }

    Runtime collect_read_counts = {
        "docker": gatk_docker,
        "jar_override": gatk_override,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu,
        "machine_mem": mem_collect_read_counts + mem_machine_overhead,
        "command_mem": mem_collect_read_counts,
        "runtime_minutes": time_startup + time_collect_read_counts,
        "disk": disk,
        "boot_disk_size": boot_disk_size
    }

    Runtime denoise_read_counts = {
        "docker": gatk_docker,
        "jar_override": gatk_override,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu,
        "machine_mem": mem_denoise_read_counts + mem_machine_overhead,
        "command_mem": mem_denoise_read_counts,
        "runtime_minutes": time_startup + time_denoise_read_counts,
        "disk": disk,
        "boot_disk_size": boot_disk_size
    }

    Runtime vcf_to_pileup_variants = {
        "docker": bcftools_docker,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu,
        "machine_mem": mem_vcf_to_pileup_variants + mem_machine_overhead,
        "command_mem": mem_vcf_to_pileup_variants,
        "runtime_minutes": time_startup + time_vcf_to_pileup_variants,
        "disk": disk,
        "boot_disk_size": boot_disk_size
    }

    Runtime get_pileup_summaries = {
        "docker": gatk_docker,
        "jar_override": gatk_override,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu,
        "machine_mem": mem_get_pileup_summaries + mem_machine_overhead,
        "command_mem": mem_get_pileup_summaries,
        "runtime_minutes": time_startup + ceil(time_get_pileup_summaries / scatter_count),
        "disk": disk,
        "boot_disk_size": boot_disk_size
    }

    Runtime gather_pileup_summaries = {
        "docker": gatk_docker,
        "jar_override": gatk_override,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu,
        "machine_mem": mem_gather_pileup_summaries + mem_machine_overhead,
        "command_mem": mem_gather_pileup_summaries,
        "runtime_minutes": time_startup + time_gather_pileup_summaries,
        "disk": disk,
        "boot_disk_size": boot_disk_size
    }

    Runtime select_pileup_summaries = {
        "docker": gatk_docker,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu,
        "machine_mem": mem_select_pileup_summaries + mem_machine_overhead,
        "command_mem": mem_select_pileup_summaries,
        "runtime_minutes": time_startup + time_select_pileup_summaries,
        "disk": disk,
        "boot_disk_size": boot_disk_size
    }

    Runtime pileup_to_allelic_counts = {
        "docker": python_docker,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu,
        "machine_mem": mem_pileup_to_allelic_counts + mem_machine_overhead,
        "command_mem": mem_pileup_to_allelic_counts,
        "runtime_minutes": time_startup + time_pileup_to_allelic_counts,
        "disk": disk,
        "boot_disk_size": boot_disk_size
    }

    Int mem_harmonize_copy_ratios = mem_harmonize_copy_ratios_base + num_bams * mem_harmonize_copy_ratios_additional_per_sample
    Runtime harmonize_copy_ratios = {
        "docker": python_docker,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu_harmonize_copy_ratios,
        "machine_mem": mem_harmonize_copy_ratios + mem_machine_overhead,
        "command_mem": mem_harmonize_copy_ratios,
        "runtime_minutes": time_startup + time_harmonize_copy_ratios,
        "disk": disk,
        "boot_disk_size": boot_disk_size
    }

    Runtime merge_allelic_counts = {
        "docker": python_docker,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu,
        "machine_mem": mem_merge_allelic_counts + mem_machine_overhead,
        "command_mem": mem_merge_allelic_counts,
        "runtime_minutes": time_startup + time_merge_allelic_counts,
        "disk": disk,
        "boot_disk_size": boot_disk_size
    }

    Runtime calculate_contamination = {
        "docker": gatk_docker,
        "jar_override": gatk_override,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu,
        "machine_mem": mem_calculate_contamination + mem_machine_overhead,
        "command_mem": mem_calculate_contamination,
        "runtime_minutes": time_startup + time_calculate_contamination,
        "disk": disk,
        "boot_disk_size": boot_disk_size
    }

    Runtime genotype_variants = {
        "docker": python_docker,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu,
        "machine_mem": mem_genotype_variants + mem_machine_overhead,
        "command_mem": mem_genotype_variants,
        "runtime_minutes": time_startup + time_genotype_variants,
        "disk": disk,
        "boot_disk_size": boot_disk_size
    }

    Int mem_model_segments = mem_model_segments_base + num_bams * mem_model_segments_additional_per_sample
    Runtime model_segments = {
        "docker": gatk_docker,
        "jar_override": gatk_override,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu,
        "machine_mem": mem_model_segments + mem_machine_overhead,
        "command_mem": mem_model_segments,
        "runtime_minutes": time_startup + time_model_segments,
        "disk": disk,
        "boot_disk_size": boot_disk_size
    }

    Runtime call_copy_ratio_segments = {
        "docker": gatk_docker,
        "jar_override": gatk_override,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu,
        "machine_mem": mem_call_copy_ratio_segments + mem_machine_overhead,
        "command_mem": mem_call_copy_ratio_segments,
        "runtime_minutes": time_startup + time_call_copy_ratio_segments,
        "disk": disk,
        "boot_disk_size": boot_disk_size
    }

    Runtime merge_calls_with_modeled_segments = {
        "docker": ubuntu_docker,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu,
        "machine_mem": mem_merge_calls_with_modeled_segments + mem_machine_overhead,
        "command_mem": mem_merge_calls_with_modeled_segments,
        "runtime_minutes": time_startup + time_merge_calls_with_modeled_segments,
        "disk": disk,
        "boot_disk_size": boot_disk_size
    }

    Runtime plot_modeled_segments = {
        "docker": gatk_docker,
        "jar_override": gatk_override,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu,
        "machine_mem": mem_plot_modeled_segments + mem_machine_overhead,
        "command_mem": mem_plot_modeled_segments,
        "runtime_minutes": time_startup + time_plot_modeled_segments,
        "disk": disk,
        "boot_disk_size": boot_disk_size
    }

    Runtime filter_germline_cnvs = {
        "docker": python_docker,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu,
        "machine_mem": mem_filter_germline_cnvs + mem_machine_overhead,
        "command_mem": mem_filter_germline_cnvs,
        "runtime_minutes": time_startup + time_filter_germline_cnvs,
        "disk": disk,
        "boot_disk_size": boot_disk_size
    }

    Runtime recount_markers = {
        "docker": python_docker,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu,
        "machine_mem": mem_recount_markers + mem_machine_overhead,
        "command_mem": mem_recount_markers,
        "runtime_minutes": time_startup + time_recount_markers,
        "disk": disk,
        "boot_disk_size": boot_disk_size
    }

    Runtime model_segments_to_acs_conversion = {
        "docker": python_docker,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu,
        "machine_mem": mem_model_segments_to_acs_conversion + mem_machine_overhead,
        "command_mem": mem_model_segments_to_acs_conversion,
        "runtime_minutes": time_startup + time_model_segments_to_acs_conversion,
        "disk": disk,
        "boot_disk_size": boot_disk_size
    }

    Runtime process_maf_for_absolute = {
        "docker": python_docker,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu,
        "machine_mem": mem_process_maf_for_absolute + mem_machine_overhead,
        "command_mem": mem_process_maf_for_absolute,
        "runtime_minutes": time_startup + time_process_maf_for_absolute,
        "disk": disk,
        "boot_disk_size": boot_disk_size
    }

    Runtime absolute = {
        "docker": absolute_docker,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu,
        "machine_mem": mem_absolute + mem_machine_overhead,
        "command_mem": mem_absolute,
        "runtime_minutes": time_startup + time_absolute,
        "disk": disk,
        "boot_disk_size": boot_disk_size
    }

    Runtime absolute_extract = {
        "docker": absolute_docker,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu,
        "machine_mem": mem_absolute + mem_machine_overhead,
        "command_mem": mem_absolute,
        "runtime_minutes": time_startup + time_absolute,
        "disk": disk,
        "boot_disk_size": boot_disk_size
    }

    Int mem_mutect2 = mem_mutect2_base + num_bams * mem_mutect2_additional_per_sample
    Runtime mutect2 = {
        "docker": gatk_docker,
        "jar_override": gatk_override,
        "preemptible": preemptible_mutect2,
        "max_retries": max_retries_mutect2,
        "cpu": cpu_mutect2,
        "machine_mem": mem_mutect2 + mem_mutect2_overhead,
        "command_mem": mem_mutect2,
        "runtime_minutes": time_startup + ceil(time_mutect2_total / scatter_count),
        "disk": disk + ceil(disk_mutect2_total / scatter_count),
        "boot_disk_size": boot_disk_size
    }

    Int mem_learn_read_orientation_model = mem_learn_read_orientation_model_base + num_bams * mem_learn_read_orientation_model_additional_per_sample
    Runtime learn_read_orientation_model = {
        "docker": gatk_docker,
        "jar_override": gatk_override,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu,
        "machine_mem": mem_learn_read_orientation_model + mem_machine_overhead,
        "command_mem": mem_learn_read_orientation_model,
        "runtime_minutes": time_startup + time_learn_read_orientation_model,
        "disk": disk,
        "boot_disk_size": boot_disk_size
    }

    Runtime merge_vcfs = {
        "docker": gatk_docker,
        "jar_override": gatk_override,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu,
        "machine_mem": mem_merge_vcfs + mem_machine_overhead,
        "command_mem": mem_merge_vcfs,
        "runtime_minutes": time_startup + time_merge_vcfs,
        "disk": 1 + disk,
        "boot_disk_size": boot_disk_size
    }

    Runtime merge_mafs = {
        "docker": ubuntu_docker,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu,
        "machine_mem": mem_merge_mafs + mem_machine_overhead,
        "command_mem": mem_merge_mafs,
        "runtime_minutes": time_startup + time_merge_mafs,
        "disk": disk,
        "boot_disk_size": boot_disk_size
    }

    Runtime merge_mutect_stats = {
        "docker": gatk_docker,
        "jar_override": gatk_override,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu,
        "machine_mem": mem_merge_mutect_stats + mem_machine_overhead,
        "command_mem": mem_merge_mutect_stats,
        "runtime_minutes": time_startup + time_merge_mutect_stats,
        "disk": disk,
        "boot_disk_size": boot_disk_size
    }

    Runtime print_reads = {
        "docker": gatk_docker,
        "jar_override": gatk_override,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu,
        "machine_mem": mem_print_reads + mem_machine_overhead,
        "command_mem": mem_print_reads,
        "runtime_minutes": time_startup + time_print_reads,
        "disk": disk + disk_print_reads,
        "boot_disk_size": boot_disk_size
    }

    Runtime filter_mutect_calls = {
        "docker": gatk_docker,
        "jar_override": gatk_override,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu,
        "machine_mem": mem_filter_mutect_calls + mem_machine_overhead,
        "command_mem": mem_filter_mutect_calls,
        "runtime_minutes": time_startup + time_filter_mutect_calls,
        "disk": disk,
        "boot_disk_size": boot_disk_size
    }

    Runtime variant_filtration = {
        "docker": gatk_docker,
        "jar_override": gatk_override,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu,
        "machine_mem": mem_variant_filtration + mem_machine_overhead,
        "command_mem": mem_variant_filtration,
        "runtime_minutes": time_startup + time_variant_filtration,
        "disk": disk,
        "boot_disk_size": boot_disk_size
    }

    Runtime left_align_and_trim_variants = {
        "docker": gatk_docker,
        "jar_override": gatk_override,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu,
        "machine_mem": mem_left_align_and_trim_variants + mem_machine_overhead,
        "command_mem": mem_left_align_and_trim_variants,
        "runtime_minutes": time_startup + time_left_align_and_trim_variants,
        "disk": disk,
        "boot_disk_size": boot_disk_size
    }

    Int mem_filter_alignment_artifacts = mem_filter_alignment_artifacts_base + num_bams * mem_filter_alignment_artifacts_additional_per_sample
    Runtime filter_alignment_artifacts = {
        "docker": gatk_docker,
        "jar_override": gatk_override,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu_filter_alignment_artifacts,
        "machine_mem": mem_filter_alignment_artifacts + mem_machine_overhead,
        "command_mem": mem_filter_alignment_artifacts,
        "runtime_minutes": time_startup + ceil(time_filter_alignment_artifacts_total / scatter_count),
        "disk": disk,
        "boot_disk_size": boot_disk_size
    }

    Runtime select_variants = {
        "docker": gatk_docker,
        "jar_override": gatk_override,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu,
        "machine_mem": mem_select_variants + mem_machine_overhead,
        "command_mem": mem_select_variants,
        "runtime_minutes": time_startup + time_select_variants,
        "disk": disk,
        "boot_disk_size": boot_disk_size
    }

    Runtime funcotate = {
        "docker": gatk_docker,
        "jar_override": gatk_override,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu,
        "machine_mem": mem_funcotate + mem_machine_overhead,
        "command_mem": mem_funcotate,
#        "runtime_minutes": time_startup + if run_variant_anntation_scattered then ceil(time_funcotate / scatter_count) else time_funcotate,
        "runtime_minutes": time_startup + time_funcotate,
        "disk": disk,
        "boot_disk_size": boot_disk_size
    }

    Runtime create_empty_annotation = {
        "docker": ubuntu_docker,
        "preemptible": preemptible,
        "max_retries": max_retries,
        'cpu': cpu,
        "machine_mem": mem_create_empty_annotations + mem_machine_overhead,
        "command_mem": mem_create_empty_annotations,
        "runtime_minutes": time_startup + time_create_empty_annotations,
        "disk": disk_sizeGB,
        "boot_disk_size": boot_disk_size
    }

    Runtime create_cnv_panel = {
        "docker": gatk_docker,
        "jar_override": gatk_override,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu,
        "machine_mem": mem_create_cnv_panel + mem_machine_overhead,
        "command_mem": mem_create_cnv_panel,
        "runtime_minutes": time_startup + time_create_cnv_panel,
        "disk": disk + disk_create_cnv_panel,
        "boot_disk_size": boot_disk_size
    }

    Runtime create_mutect2_panel = {
        "docker": gatk_docker,
        "jar_override": gatk_override,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu,
        "machine_mem": mem_create_mutect2_panel + mem_machine_overhead,
        "command_mem": mem_create_mutect2_panel,
        "runtime_minutes": time_startup + time_create_mutect2_panel,
        "disk": disk + disk_create_mutect2_panel,
        "boot_disk_size": boot_disk_size
    }

    Runtime select_af_only_from_vcf = {
        "docker": bcftools_docker,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu,
        "machine_mem": mem_select_af_only_from_vcf + mem_machine_overhead,
        "command_mem": mem_select_af_only_from_vcf,
        "runtime_minutes": time_startup + time_select_af_only_from_vcf,
        "disk": disk,
        "boot_disk_size": boot_disk_size
    }

    RuntimeCollection runtime_collection = {
        "get_tumor_sample_names": get_tumor_sample_names,
        "get_sample_name": get_sample_name,
        "annotate_intervals": annotate_intervals,
        "preprocess_intervals": preprocess_intervals,
        "split_intervals": split_intervals,

        "collect_read_counts": collect_read_counts,
        "denoise_read_counts": denoise_read_counts,
        "vcf_to_pileup_variants": vcf_to_pileup_variants,
        "get_pileup_summaries": get_pileup_summaries,
        "gather_pileup_summaries": gather_pileup_summaries,
        "select_pileup_summaries": select_pileup_summaries,
        "pileup_to_allelic_counts": pileup_to_allelic_counts,
        "harmonize_copy_ratios": harmonize_copy_ratios,
        "merge_allelic_counts": merge_allelic_counts,
        "calculate_contamination": calculate_contamination,
        "genotype_variants": genotype_variants,
        "model_segments": model_segments,
        "call_copy_ratio_segments": call_copy_ratio_segments,
        "merge_calls_with_modeled_segments": merge_calls_with_modeled_segments,
        "plot_modeled_segments": plot_modeled_segments,
        "filter_germline_cnvs": filter_germline_cnvs,
        "recount_markers": recount_markers,

        "model_segments_to_acs_conversion": model_segments_to_acs_conversion,
        "process_maf_for_absolute": process_maf_for_absolute,
        "absolute": absolute,
        "absolute_extract": absolute_extract,

        "mutect2": mutect2,
        "learn_read_orientation_model": learn_read_orientation_model,
        "merge_vcfs": merge_vcfs,
        "merge_mafs": merge_mafs,
        "merge_mutect_stats": merge_mutect_stats,
        "print_reads": print_reads,
        "filter_mutect_calls": filter_mutect_calls,
        "variant_filtration": variant_filtration,
        "left_align_and_trim_variants": left_align_and_trim_variants,
        "filter_alignment_artifacts": filter_alignment_artifacts,
        "select_variants": select_variants,
        "funcotate": funcotate,
        "create_empty_annotation": create_empty_annotation,
        "create_cnv_panel": create_cnv_panel,
        "create_mutect2_panel": create_mutect2_panel,

        "collect_covered_regions": collect_covered_regions,
        "select_af_only_from_vcf": select_af_only_from_vcf,
    }

    output {
        RuntimeCollection rtc = runtime_collection
    }
}