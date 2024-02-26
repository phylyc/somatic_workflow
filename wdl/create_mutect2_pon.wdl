version development

import "runtimes.wdl"
import "tasks.wdl"
import "multi-sample_somatic_workflow.wdl" as mssw


workflow CreateMutect2PanelOfNormals {
    input {
        File target_interval_list
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        Array[File]? normal_bams
        Array[File]? normal_bais
        File? normal_bams_file
        File? normal_bais_file

        File germline_resource
        File germline_resource_idx

        Boolean compress_output = true
        String mutect2_extra_args = ""
        String pon_name

        Int min_contig_size = 1000000
        Int num_contigs = 24

        RuntimeCollection runtime_collection = GetRTC.rtc

        # runtime
        Int scatter_count = 10
        String gatk_docker = "broadinstitute/gatk"
        File? gatk_override
        Int preemptible = 1
        Int max_retries = 1
        Int disk_sizeGB = 1
        Int cpu = 1

        Int mem_machine_overhead = 2048
        Int mem_create_mutect2_panel = 16384

        Int disk_create_mutect2_panel = 10

        Int time_startup = 10
        Int time_create_mutect2_panel = 1200  # 20 h
    }

    call runtimes.DefineRuntimeCollection as GetRTC {
        input:
            scatter_count = scatter_count,
            gatk_docker = gatk_docker,
            gatk_override = gatk_override,
            max_retries = max_retries,
            preemptible = preemptible,
            cpu = cpu,
            disk_sizeGB = disk_sizeGB,
            time_startup = time_startup,
            mem_machine_overhead = mem_machine_overhead,
            mem_create_mutect2_panel = mem_create_mutect2_panel,
            time_create_mutect2_panel = time_create_mutect2_panel,
            disk_create_mutect2_panel = disk_create_mutect2_panel
    }

    # todo: assert either normal_bams or normal_bams_file is defined

    Array[File] non_optional_normal_bams = if defined(normal_bams) then select_first([normal_bams, []]) else read_lines(select_first([normal_bams_file, ""]))
    Array[File] non_optional_normal_bais = if defined(normal_bais) then select_first([normal_bais, []]) else read_lines(select_first([normal_bais_file, ""]))

    scatter (normal in zip(non_optional_normal_bams, non_optional_normal_bais)) {
        call tasks.GetSampleName {
            input:
                bam = normal.left,
                runtime_params = runtime_collection.get_sample_name,
        }

        call mssw.MultiSampleSomaticWorkflow {
            input:
                interval_list = target_interval_list,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                ref_dict = ref_dict,

                individual_id = GetSampleName.sample_name,
                tumor_bams = [normal.left],
                tumor_bais = [normal.right],
                tumor_target_intervals = [target_interval_list],

                run_collect_covered_regions = false,
                run_collect_target_coverage = false,
                run_collect_allelic_coverage = false,
                run_contamination_model = false,
                run_orientation_bias_mixture_model = false,
                run_variant_calling = true,
                run_variant_filter = false,
                run_realignment_filter = false,
                run_realignment_filter_only_on_high_confidence_variants = false,
                run_collect_called_variants_allelic_coverage = false,
                run_variant_annotation = false,

                keep_germline = false,
                compress_output = true,
                make_bamout = false,

                compress_output = compress_output,
                scatter_count = scatter_count,
                mutect2_extra_args = mutect2_extra_args + " --max-mnp-distance 0",  # GenpmicsDBImport can't handle MNPs

                runtime_collection = runtime_collection,
        }
    }

    String split_intervals_extra_args = (
        "--subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION"
        + " --min-contig-size " + min_contig_size
    )
    call tasks.SplitIntervals {
        input:
            interval_list = target_interval_list,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            scatter_count = num_contigs,
            split_intervals_extra_args = split_intervals_extra_args,
            runtime_params = runtime_collection.split_intervals
    }

    call runtimes.UpdateRuntimeParameters as CreateMutect2PanelRuntime {
        input:
            runtime_params = runtime_collection.create_mutect2_panel,
            disk = runtime_collection.create_mutect2_panel.disk + 3 * ceil(size(MultiSampleSomaticWorkflow.unfiltered_vcf, "GB")) + ceil(length(MultiSampleSomaticWorkflow.unfiltered_vcf) / 10)
    }

    scatter (scattered_intervals in SplitIntervals.interval_files) {
        call CreateMutect2Panel {
            input:
                input_vcfs = select_all(MultiSampleSomaticWorkflow.unfiltered_vcf),
                input_vcf_indices = select_all(MultiSampleSomaticWorkflow.unfiltered_vcf_idx),
                interval_list = scattered_intervals,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                ref_dict = ref_dict,
                compress_output = compress_output,
                gnomad = germline_resource,
                gnomad_idx = germline_resource_idx,
                output_vcf_name = pon_name,
                runtime_params = CreateMutect2PanelRuntime.params
        }
    }

    call tasks.MergeVCFs {
        input:
            vcfs = CreateMutect2Panel.output_vcf,
            vcfs_idx = CreateMutect2Panel.output_vcf_index,
            output_name = pon_name,
            compress_output = compress_output,
            runtime_params = runtime_collection.merge_vcfs
    }

    output {
        File pon = MergeVCFs.merged_vcf
        File pon_idx = MergeVCFs.merged_vcf_idx
        Array[File] normal_calls = select_all(MultiSampleSomaticWorkflow.unfiltered_vcf)
        Array[File] normal_calls_idx = select_all(MultiSampleSomaticWorkflow.unfiltered_vcf_idx)
    }
}

task CreateMutect2Panel {
    input {
        File interval_list
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        Array[File]+ input_vcfs
        Array[File]+ input_vcf_indices
        String output_vcf_name

        Boolean compress_output = true
        File gnomad
        File gnomad_idx

        # UMCCR found that 5 is optimizing the F2 score, but not by much compared to 2:
        # https://umccr.org/blog/panel-of-normals/
        Int min_sample_count = 2
        String? create_pon_extra_args

        # runtime
        Runtime runtime_params
    }

    String pon_file = "pon.vcf.gz"
    String output_file = output_vcf_name + if compress_output then ".vcf.gz" else ".vcf"
    String output_file_idx = output_file + if compress_output then ".tbi" else ".idx"

    command {
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.jar_override}

        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            GenomicsDBImport \
            --genomicsdb-workspace-path pon_db \
            -R '~{ref_fasta}' \
            ~{sep="' " prefix("-V '", input_vcfs)}' \
            -L '~{interval_list}'

        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            CreateSomaticPanelOfNormals \
            -R '~{ref_fasta}' \
            --germline-resource '~{gnomad}' \
            -V gendb://pon_db \
            -O '~{pon_file}' \
            --min-sample-count ~{min_sample_count} \
            ~{create_pon_extra_args}

        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            SortVcf \
            -I '~{pon_file}' \
            -O '~{output_file}'
    }

    output {
        File output_vcf = output_file
        File output_vcf_index = output_file_idx
    }

    runtime {
        docker: runtime_params.docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: runtime_params.machine_mem + " MB"
        runtime_minutes: runtime_params.runtime_minutes
        disks: "local-disk " + runtime_params.disk + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

    # GenomicsDB requires that the reference be a local file.
    parameter_meta{
        input_vcfs: {localization_optional: true}
        input_vcf_indices: {localization_optional: true}
        gnomad: {localization_optional: true}
        gnomad_idx: {localization_optional: true}
    }
}