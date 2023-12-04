version development

# Asking for a machine on gcp with specified cpu and memory will make gcp provide
# a machine with AT LEAST that much cpu and memory. E.g. a request for a machine
# with 4 cpus and 64GB mem can result in a delivered machine with 10 cpus. Thus,
# tasks should not necessarily request a specific number of threads, but use
# $(nproc) instead. That said, tasks with more threads might need higher memory.

    
struct Runtime {
    String docker
    File? jar_override
    Int preemptible
    Int max_retries
    Int cpu
    Int machine_mem
    Int command_mem
    Int runtime_minutes
    Int disk
    Int boot_disk_size
}


workflow UpdateRuntimeParameters {
    input {
        Runtime runtime_params
        String? docker
        File? jar_override
        Int? preemptible
        Int? max_retries
        Int? cpu
        Int? machine_mem
        Int? command_mem
        Int? runtime_minutes
        Int? disk
        Int? boot_disk_size
    }

    Runtime updated_runtime = {
        "docker": select_first([docker, runtime_params.docker]),
        "jar_override": select_first([jar_override, runtime_params.jar_override]),
        "preemptible": select_first([preemptible, runtime_params.preemptible]),
        "max_retries": select_first([max_retries, runtime_params.max_retries]),
        "cpu": select_first([cpu, runtime_params.cpu]),
        "machine_mem": select_first([machine_mem, runtime_params.machine_mem]),
        "command_mem": select_first([command_mem, runtime_params.command_mem]),
        "runtime_minutes": select_first([runtime_minutes, runtime_params.runtime_minutes]),
        "disk": select_first([disk, runtime_params.disk]),
        "boot_disk_size": select_first([boot_disk_size, runtime_params.boot_disk_size])
    }

    output {
        Runtime params = updated_runtime
    }
}


workflow DefineRuntimes {
    input {
        Int num_bams = 1

        Int scatter_count = 10
        String gatk_docker = "broadinstitute/gatk"
        # Needs docker image with bedtools, samtools, and gatk
        String jupyter_docker = "us.gcr.io/broad-dsp-gcr-public/terra-jupyter-gatk"  # 27.5GB todo: find smaller image. This one takes ~13 mins to spin up.
        String ubuntu_docker = "ubuntu"
        String bcftools_docker = "staphb/bcftools"
        String genotype_docker = "civisanalytics/datascience-python:latest"
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
        Int mem_additional_per_sample = 256  # this depends on bam size (WES vs WGS)

        Int time_startup = 10

        # gatk: GetSampleName
        Int mem_get_sample_name = 256
        Int time_get_sample_name = 1

        # gatk: PreprocessIntervals
        Int mem_preprocess_intervals = 3072
        Int time_preprocess_intervals = 60

        # gatk: SplitIntervals
        Int mem_split_intervals = 1024
        Int time_split_intervals = 1

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

        # gatk: CalculateContamination
        Int mem_calculate_contamination = 3072  # depends on the variants_for_contamination resource
        Int time_calculate_contamination = 10

        # custom genotyping script based on CalculateContamination model
        Int mem_genotype_variants = 8192
        Int time_genotype_variants = 30

        # gatk: Mutect2
        Int cpu_variant_call = 1  # good for PairHMM: 2
        Int mem_variant_call_base = 3072
        Int time_variant_call_total = 10000  # 6 d / scatter_count
        Int disk_variant_call = 0

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

        # gatk: FilterAlignmentArtifacts
        Int cpu_filter_alignment_artifacts = 1  # good for PairHMM: 4
        Int mem_filter_alignment_artifacts_base = 3072  # needs to be increased in some cases
        Int time_filter_alignment_artifacts_total = 10000  # 12 d / scatter_count

        # gatk: SelectVariants
        Int mem_select_variants = 3072
        Int time_select_variants = 5

        # gatk: CNNScoreVariants
        Int cpu_cnn_scoring = 1
        Int mem_cnn_scoring = 4096
        Int time_cnn_scoring = 10

        # gatk: Funcotator
        Int mem_funcotate = 6144
        Int time_funcotate = 1440  # 24 h

        # gatk: GenomicsDBImport / CreateSomaticPanelOfNormals
        Int mem_create_panel = 16384
        Int time_create_panel = 1200  # 20 h
        Int disk_create_panel = 0

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
        "docker": genotype_docker,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu,
        "machine_mem": mem_genotype_variants + mem_machine_overhead,
        "command_mem": mem_genotype_variants,
        "runtime_minutes": time_startup + time_genotype_variants,
        "disk": disk,
        "boot_disk_size": boot_disk_size
    }

    Int mem_variant_call = mem_variant_call_base + num_bams * mem_additional_per_sample
    Runtime variant_call = {
        "docker": gatk_docker,
        "jar_override": gatk_override,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu_variant_call,
        "machine_mem": mem_variant_call + mem_machine_overhead,
        "command_mem": mem_variant_call,
        "runtime_minutes": time_startup + ceil(time_variant_call_total / scatter_count),
        "disk": disk + disk_variant_call,
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

    Runtime cnn_scoring = {
        "docker": gatk_docker,
        "jar_override": gatk_override,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu_cnn_scoring,
        "machine_mem": mem_cnn_scoring + mem_machine_overhead,
        "command_mem": mem_cnn_scoring,
        "runtime_minutes": time_startup + time_cnn_scoring,
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
        "runtime_minutes": time_startup + ceil(time_funcotate / scatter_count),
        "disk": disk,
        "boot_disk_size": boot_disk_size
    }

    Runtime create_panel = {
        "docker": gatk_docker,
        "jar_override": gatk_override,
        "preemptible": preemptible,
        "max_retries": max_retries,
        "cpu": cpu,
        "machine_mem": mem_create_panel + mem_machine_overhead,
        "command_mem": mem_create_panel,
        "runtime_minutes": time_startup + time_create_panel,
        "disk": disk + disk_create_panel,
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

    output {
        Runtime get_sample_name_runtime = get_sample_name
        Runtime preprocess_intervals_runtime = preprocess_intervals
        Runtime split_intervals_runtime = split_intervals
        Runtime collect_covered_regions_runtime = collect_covered_regions
        Runtime collect_read_counts_runtime = collect_read_counts
        Runtime denoise_read_counts_runtime = denoise_read_counts
        Runtime vcf_to_pileup_variants_runtime = vcf_to_pileup_variants
        Runtime get_pileup_summaries_runtime = get_pileup_summaries
        Runtime gather_pileup_summaries_runtime = gather_pileup_summaries
        Runtime select_pileup_summaries_runtime = select_pileup_summaries
        Runtime calculate_contamination_runtime = calculate_contamination
        Runtime genotype_variants_runtime = genotype_variants
        Runtime variant_call_runtime = variant_call
        Runtime learn_read_orientation_model_runtime = learn_read_orientation_model
        Runtime merge_vcfs_runtime = merge_vcfs
        Runtime merge_mafs_runtime = merge_mafs
        Runtime merge_mutect_stats_runtime = merge_mutect_stats
        Runtime merge_bams_runtime = merge_bams
        Runtime filter_mutect_calls_runtime = filter_mutect_calls
        Runtime filter_alignment_artifacts_runtime = filter_alignment_artifacts
        Runtime select_variants_runtime = select_variants
        Runtime cnn_scoring_runtime = cnn_scoring
        Runtime funcotate_runtime = funcotate
        Runtime create_panel_runtime = create_panel
        Runtime whatshap_runtime = whatshap
        Runtime shapeit4_runtime = shapeit4
        Runtime shapeit5_runtime = shapeit5
    }
}