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
    Runtime model_segments_to_acs_conversion
    Runtime process_maf_for_absolute
    Runtime absolute
    Runtime mutect2
    Runtime learn_read_orientation_model
    Runtime merge_vcfs
    Runtime merge_mafs
    Runtime merge_mutect_stats
    Runtime merge_bams
    Runtime filter_mutect_calls
    Runtime variant_filtration
    Runtime filter_alignment_artifacts
    Runtime select_variants
    Runtime funcotate
    Runtime create_cnv_panel
    Runtime create_mutect2_panel
    Runtime whatshap
    Runtime shapeit4
    Runtime shapeit5
}


workflow DefineRuntimeCollection {
    input {
        Int num_bams = 1

        Int scatter_count = 10
        String gatk_docker = "broadinstitute/gatk:4.5.0.0"
        # Needs docker image with bedtools, samtools, and gatk
        String jupyter_docker = "us.gcr.io/broad-dsp-gcr-public/terra-jupyter-gatk"  # 27.5GB todo: find smaller image. This one takes ~13 mins to spin up.
        String tag_cga_pipline_docker = "us.gcr.io/tag-team-160914/neovax-tag-cga-pipeline:v1"
        String ubuntu_docker = "ubuntu"
        String bcftools_docker = "staphb/bcftools"
        String python_docker = "civisanalytics/datascience-python:latest"
        String whatshap_docker = "hangsuunc/whatshap:v1"
        String shapeit4_docker = "yussab/shapeit4:4.2.2"
        String shapeit5_docker = "lindonkambule/shapeit5_2023-05-05_d6ce1e2"
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
        Int mem_additional_per_sample = 384  # this depends on bam size (WES vs WGS)

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
        Int mem_get_pileup_summaries = 4096  # needs at least 2G
        Int time_get_pileup_summaries = 4500  # 3 d / scatter_count

        # gatk: GatherPileupSummaries
        Int mem_gather_pileup_summaries = 512  # 64
        Int time_gather_pileup_summaries = 5

        # SelectPileupSummaries
        Int mem_select_pileup_summaries = 512
        Int time_select_pileup_summaries = 5

        # PileupToAllelicCounts
        Int mem_pileup_to_allelic_counts = 512
        Int time_pileup_to_allelic_counts = 5

        # HarmonizeCopyRatios
        Int cpu_harmonize_copy_ratios = 4
        Int mem_harmonize_copy_ratios_base = 8192
        Int time_harmonize_copy_ratios = 60

        # MergeAllelicCounts
        Int mem_merge_allelic_counts = 4096
        Int time_merge_allelic_counts = 10

        # gatk: CalculateContamination
        Int mem_calculate_contamination = 3072  # depends on the variants_for_contamination resource
        Int time_calculate_contamination = 10

        # custom genotyping script based on CalculateContamination model
        Int mem_genotype_variants = 8192
        Int time_genotype_variants = 30

        # gatk: ModelSegments
        Int mem_model_segments = 2048
        Int time_model_segments = 60

        # gatk: CallCopyRatioSegments
        Int mem_call_copy_ratio_segments = 2048
        Int time_call_copy_ratio_segments = 10

        # custom task to merge calls with modeled segments
        Int mem_merge_calls_with_modeled_segments = 512
        Int time_merge_calls_with_modeled_segments = 1

        # gatk: PlotModeledSegments
        Int mem_plot_modeled_segments = 4096
        Int time_plot_modeled_segments = 10

        # custom filtering script
        Int mem_filter_germline_cnvs = 2048
        Int time_filter_germline_cnvs = 10

        #######################################################################
        ### SNV workflow
        #######################################################################

        # gatk: Mutect2
        Int cpu_mutect2 = 4
        Int mem_mutect2_base = 3072
        Int time_mutect2_total = 10000  # 6 d / scatter_count
        Int disk_mutect2 = 0

        # gatk: MergeVCFs
        Int mem_merge_vcfs = 2048
        Int time_merge_vcfs = 10

        # gatk: MergeMAFs
        Int mem_merge_mafs = 512
        Int time_merge_mafs = 5

        # gatk: MergeMutectStats
        Int mem_merge_mutect_stats = 512 # 64
        Int time_merge_mutect_stats = 1

        # gatk: MergeBams
        Int mem_merge_bams = 8192  # wants at least 6G
        Int time_merge_bams = 60

        # gatk: LearnReadOrientationModel
        Int mem_learn_read_orientation_model_base = 4096
        Int time_learn_read_orientation_model = 180  # 3 h

        # gatk: FilterMutectCalls
        Int mem_filter_mutect_calls = 4096
        Int time_filter_mutect_calls = 800  # 13 h

        # gatk: VariantFiltration
        Int mem_variant_filtration = 2048
        Int time_variant_filtration = 5

        # gatk: FilterAlignmentArtifacts
        Int cpu_filter_alignment_artifacts = 4
        Int mem_filter_alignment_artifacts_base = 3072  # needs to be increased in some cases
        Int time_filter_alignment_artifacts_total = 10000  # 12 d / scatter_count

        # gatk: SelectVariants
        Int mem_select_variants = 3072
        Int time_select_variants = 5

        # gatk: Funcotator
        Int mem_funcotate = 6144
        Int time_funcotate = 1440  # 24 h

        #######################################################################
        ### Clonal Analysis worklflow
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

        # Whatshap
        Int cpu_whatshap = 1
        Int mem_whatshap = 4096
        Int time_whatshap = 60

        # Shapeit
        Int cpu_shapeit = 4
        Int mem_shapeit = 4096
        Int time_shapeit = 60
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
        "docker": ubuntu_docker,
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
        "docker": ubuntu_docker,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu,
        "machine_mem": mem_pileup_to_allelic_counts + mem_machine_overhead,
        "command_mem": mem_pileup_to_allelic_counts,
        "runtime_minutes": time_startup + time_pileup_to_allelic_counts,
        "disk": disk,
        "boot_disk_size": boot_disk_size
    }

    Int mem_harmonize_copy_ratios = mem_harmonize_copy_ratios_base + num_bams * mem_additional_per_sample
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
        "docker": tag_cga_pipline_docker,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu,
        "machine_mem": mem_absolute + mem_machine_overhead,
        "command_mem": mem_absolute,
        "runtime_minutes": time_startup + time_absolute,
        "disk": disk,
        "boot_disk_size": boot_disk_size
    }

    Int mem_mutect2 = mem_mutect2_base + num_bams * mem_additional_per_sample
    Runtime mutect2 = {
        "docker": gatk_docker,
        "jar_override": gatk_override,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu_mutect2,
        "machine_mem": mem_mutect2 + mem_machine_overhead,
        "command_mem": mem_mutect2,
        "runtime_minutes": time_startup + ceil(time_mutect2_total / scatter_count),
        "disk": disk + disk_mutect2,
        "boot_disk_size": boot_disk_size
    }

    Int mem_learn_read_orientation_model = mem_learn_read_orientation_model_base + num_bams * mem_additional_per_sample
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

    Runtime merge_bams = {
        "docker": gatk_docker,
        "jar_override": gatk_override,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu,
        "machine_mem": mem_merge_bams + mem_machine_overhead,
        "command_mem": mem_merge_bams,
        "runtime_minutes": time_startup + time_merge_bams,
        "disk": disk,
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

    Int mem_filter_alignment_artifacts = mem_filter_alignment_artifacts_base + num_bams * mem_additional_per_sample
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

    Runtime whatshap = {
        "docker": whatshap_docker,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu_whatshap,
        "machine_mem": mem_whatshap + mem_machine_overhead,
        "command_mem": mem_whatshap,
        "runtime_minutes": time_startup + time_whatshap,
        "disk": disk,
        "boot_disk_size": boot_disk_size
    }

    Runtime shapeit4 = {
        "docker": shapeit4_docker,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu_shapeit,
        "machine_mem": mem_shapeit + mem_machine_overhead,
        "command_mem": mem_shapeit,
        "runtime_minutes": time_startup + time_shapeit,
        "disk": disk,
        "boot_disk_size": boot_disk_size
    }

    Runtime shapeit5 = {
        "docker": shapeit5_docker,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu_shapeit,
        "machine_mem": mem_shapeit + mem_machine_overhead,
        "command_mem": mem_shapeit,
        "runtime_minutes": time_startup + time_shapeit,
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

        "model_segments_to_acs_conversion": model_segments_to_acs_conversion,
        "process_maf_for_absolute": process_maf_for_absolute,
        "absolute": absolute,

        "mutect2": mutect2,
        "learn_read_orientation_model": learn_read_orientation_model,
        "merge_vcfs": merge_vcfs,
        "merge_mafs": merge_mafs,
        "merge_mutect_stats": merge_mutect_stats,
        "merge_bams": merge_bams,
        "filter_mutect_calls": filter_mutect_calls,
        "variant_filtration": variant_filtration,
        "filter_alignment_artifacts": filter_alignment_artifacts,
        "select_variants": select_variants,
        "funcotate": funcotate,
        "create_cnv_panel": create_cnv_panel,
        "create_mutect2_panel": create_mutect2_panel,

        "collect_covered_regions": collect_covered_regions,

        "whatshap": whatshap,
        "shapeit4": shapeit4,
        "shapeit5": shapeit5
    }

    output {
        RuntimeCollection rtc = runtime_collection
    }
}