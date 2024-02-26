version development

import "runtimes.wdl"
import "tasks.wdl"
import "collect_read_counts.wdl" as crc


workflow CreateReadCountPanelOfNormals {
    input {
        File interval_list
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        Array[File]? normal_bams
        Array[File]? normal_bais
        File? normal_bams_file
        File? normal_bais_file

        String pon_name

        File? annotated_interval_list
        # If annotated_interval_list is specified, those args are ignored:
        File? mappability_track
        File? mappability_track_idx
        File? segmental_duplication_track
        File? segmental_duplication_track_idx

        RuntimeCollection runtime_collection = GetRTC.rtc

        String gatk_docker = "broadinstitute/gatk"
        File? gatk_override
        Int preemptible = 1
        Int max_retries = 1
        Int disk_sizeGB = 1

        # memory assignments in MB
        Int mem_annotate_intervals = 2048
        Int mem_collect_read_counts = 2048
        Int mem_create_cnv_panel = 16384

        Int disk_create_cnv_panel = 10

        Int time_annotate_intervals = 60
        Int time_collect_read_counts = 300
        Int time_create_cnv_panel = 1200
    }

    # todo: assert either normal_bams or normal_bams_file is defined

    Array[File] non_optional_normal_bams = if defined(normal_bams) then select_first([normal_bams, []]) else read_lines(select_first([normal_bams_file, ""]))
    Array[File] non_optional_normal_bais = if defined(normal_bais) then select_first([normal_bais, []]) else read_lines(select_first([normal_bais_file, ""]))

    call runtimes.DefineRuntimeCollection as GetRTC {
        input:
            gatk_docker = gatk_docker,
            gatk_override = gatk_override,
            preemptible = preemptible,
            max_retries = max_retries,
            disk_sizeGB = disk_sizeGB,
            mem_annotate_intervals = mem_annotate_intervals,
            mem_collect_read_counts = mem_collect_read_counts,
            mem_create_cnv_panel = mem_create_cnv_panel,
            disk_create_cnv_panel = disk_create_cnv_panel,
            time_annotate_intervals = time_annotate_intervals,
            time_collect_read_counts = time_collect_read_counts,
            time_create_cnv_panel = time_create_cnv_panel,
    }

    scatter (normal_bam in zip(non_optional_normal_bams, non_optional_normal_bais)) {
        call crc.CollectReadCounts {
            input:
                interval_list = interval_list,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                ref_dict = ref_dict,
                bam = normal_bam.left,
                bai = normal_bam.right,
                format = "HDF5",
                gatk_override = gatk_override,
                gatk_docker = gatk_docker,
                preemptible = preemptible,
                max_retries = max_retries,
                runtime_collection = runtime_collection,
        }
    }

    if (!defined(annotated_interval_list)) {
        call tasks.AnnotateIntervals {
            input:
                interval_list = interval_list,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                ref_dict = ref_dict,
                mappability_track = mappability_track,
                mappability_track_idx = mappability_track_idx,
                segmental_duplication_track = segmental_duplication_track,
                segmental_duplication_track_idx = segmental_duplication_track_idx,
                runtime_params = runtime_collection.annotate_intervals,
        }
    }
    File this_annotated_interval_list = select_first([annotated_interval_list, AnnotateIntervals.annotated_interval_list])

    call runtimes.UpdateRuntimeParameters as CreateCNVPanelRuntime {
        input:
            runtime_params = runtime_collection.create_cnv_panel,
            disk = runtime_collection.create_cnv_panel.disk + ceil(size(CollectReadCounts.read_counts, "GB"))
    }

    call CreateReadCountPanelOfNormals {
        input:
            input_counts = CollectReadCounts.read_counts,
            output_name = pon_name,
            annotated_interval_list = this_annotated_interval_list,
            runtime_params = CreateCNVPanelRuntime.params,
	}

    output {
        File pon = CreateReadCountPanelOfNormals.cnv_pon
        File? new_annotated_interval_list = AnnotateIntervals.annotated_interval_list
    }
}

task CreateReadCountPanelOfNormals {
    input {
        Array[File]+ input_counts
        String output_name

        File? annotated_interval_list
        Int number_of_eigensamples = 20

        Runtime runtime_params
    }

    String output_pon = output_name + ".hdf5"

	command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.jar_override}
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            CreateReadCountPanelOfNormals \
            ~{sep="' " prefix("-I '", input_counts)}' \
            -O '~{output_pon}' \
            ~{"--annotated-intervals '" + annotated_interval_list + "'"} \
            --number-of-eigensamples ~{number_of_eigensamples}
	>>>

	output {
        File cnv_pon = output_pon
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

    parameter_meta {
#        input_counts: {localization_optional: true}
#        annotated_interval_list: {localization_optional: true}
    }
}