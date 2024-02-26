version development

import "runtimes.wdl"
import "tasks.wdl"


workflow CollectReadCounts {
    input {
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        String format = "TSV"

        # SequencingRun sequencing_run = GetSeqRun.sequencing_run
        String? sample_name
        File bam
        File bai
        File interval_list
        File? annotated_interval_list
        File? read_count_panel_of_normals
        Boolean paired_end = false

        RuntimeCollection runtime_collection = GetRTC.rtc
        String gatk_docker = "broadinstitute/gatk"
        File? gatk_override
        Int preemptible = 1
        Int max_retries = 1
        # memory assignments in MB
        Int mem_get_sample_name = 512  # 256
        Int mem_collect_read_counts = 2048
        Int mem_denoise_read_counts = 2048
        Int time_get_sample_name = 1
        Int time_collect_read_counts = 300
        Int time_denoise_read_counts = 120
    }

    if (!defined(sample_name)) {
        call tasks.GetSampleName {
            input:
                bam = bam,
                runtime_params = runtime_collection.get_sample_name
        }
    }

#    call sequencing_run.DefineSequencingRun as GetSeqRun {
#        input:
#            name = select_first([sample_name, GetSampleName.sample_name]),
#            bam = bam,
#            bai = bai,
#            interval_list = interval_list,
#            annotated_interval_list = annotated_interval_list,
#            read_count_panel_of_normals = read_count_panel_of_normals,
#            paired_end = paired_end
#    }

    call runtimes.DefineRuntimeCollection as GetRTC {
        input:
            gatk_docker = gatk_docker,
            gatk_override = gatk_override,
            max_retries = max_retries,
            preemptible = preemptible,
            mem_get_sample_name = mem_get_sample_name,
            mem_collect_read_counts = mem_collect_read_counts,
            mem_denoise_read_counts = mem_denoise_read_counts,
            time_get_sample_name = time_get_sample_name,
            time_collect_read_counts = time_collect_read_counts,
            time_denoise_read_counts = time_denoise_read_counts
    }

    String this_sample_name = select_first([sample_name, GetSampleName.sample_name])

	call CollectReadCounts {
		input:
			interval_list = interval_list,
            ref_fasta = ref_fasta,
			ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            bam = bam,
            bai = bai,
            format = format,
            sample_name = this_sample_name,
            paired_end = paired_end,
            runtime_params = runtime_collection.collect_read_counts
	}

    if (defined(annotated_interval_list) || defined(read_count_panel_of_normals)) {
        call DenoiseReadCounts {
            input:
                read_counts = CollectReadCounts.read_counts,
                sample_name = this_sample_name,
                annotated_interval_list = annotated_interval_list,
                count_panel_of_normals = read_count_panel_of_normals,
                runtime_params = runtime_collection.denoise_read_counts
        }
    }

    output {
        String sample_name = this_sample_name
        File read_counts = CollectReadCounts.read_counts
        File? denoised_copy_ratios = DenoiseReadCounts.denoised_copy_ratios
        File? standardized_copy_ratios = DenoiseReadCounts.standardized_copy_ratios
    }
}

task CollectReadCounts {
    input {
        File interval_list
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        File bam
        File bai
        String sample_name
        String interval_merging_rule = "OVERLAPPING_ONLY"
        String format = "TSV"
        Boolean paired_end = false

        Runtime runtime_params
    }

    parameter_meta {
        interval_list: {localization_optional: true}
        ref_fasta: {localization_optional: true}
        ref_fasta_index: {localization_optional: true}
        ref_dict: {localization_optional: true}
        bam: {localization_optional: true}
        bai: {localization_optional: true}
    }

    String output_name = sample_name + ".read_counts.tsv"

	command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.jar_override}
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            CollectReadCounts \
            -I '~{bam}' \
            -L '~{interval_list}' \
            -R '~{ref_fasta}' \
            -O '~{output_name}' \
            --interval-merging-rule ~{interval_merging_rule} \
            --format ~{format} \
            ~{true="--read-filter MateUnmappedAndUnmappedReadFilter " false="" paired_end} \
            ~{true="--read-filter PairedReadFilter " false="" paired_end}
	>>>

	output {
		File read_counts = output_name
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
}

task DenoiseReadCounts {
    input {
        File read_counts
        String sample_name
        File? annotated_interval_list
        File? count_panel_of_normals
        Int? number_of_eigensamples

        Runtime runtime_params
    }

    parameter_meta {
        # read_counts: {localization_optional: true}
        # annotated_interval_list: {localization_optional: true}
        # count_panel_of_normals: {localization_optional: true}
    }

    String output_denoised_copy_ratios = sample_name + ".denoised_CR.tsv"
    String output_standardized_copy_ratios = sample_name + ".standardized_CR.tsv"

	command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.jar_override}
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            DenoiseReadCounts \
            -I '~{read_counts}' \
            --denoised-copy-ratios '~{output_denoised_copy_ratios}' \
            --standardized-copy-ratios '~{output_standardized_copy_ratios}' \
            ~{"--number-of-eigensamples " + number_of_eigensamples} \
            ~{"--annotated-intervals '" + annotated_interval_list + "'"} \
            ~{"--count-panel-of-normals '" + count_panel_of_normals + "'"}
	>>>

	output {
        File denoised_copy_ratios = output_denoised_copy_ratios
        File standardized_copy_ratios = output_standardized_copy_ratios
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
}