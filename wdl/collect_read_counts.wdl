version development

import "runtime_collection.wdl" as rtc
import "runtimes.wdl" as rt
import "sequencing_run.define.wdl" as seqrun_def


workflow CollectReadCounts {
    input {
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        String format = "TSV"

        String sample_name
        File bam
        File bai
        File target_intervals
        File? annotated_target_intervals
        File? read_count_panel_of_normals
        Boolean is_paired_end = false

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

    call rtc.DefineRuntimeCollection as GetRTC {
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

	call CollectReadCountsTask {
		input:
			interval_list = target_intervals,
            ref_fasta = ref_fasta,
			ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            bam = bam,
            bai = bai,
            format = format,
            sample_name = sample_name,
            is_paired_end = is_paired_end,
            runtime_params = runtime_collection.collect_read_counts
	}

    if (defined(annotated_target_intervals) || defined(read_count_panel_of_normals)) {
        call DenoiseReadCounts {
            input:
                read_counts = CollectReadCountsTask.read_counts,
                sample_name = sample_name,
                annotated_interval_list = annotated_target_intervals,
                count_panel_of_normals = read_count_panel_of_normals,
                runtime_params = runtime_collection.denoise_read_counts
        }
    }

    output {
        String sample_name = sample_name
        File read_counts = CollectReadCountsTask.read_counts
        File? denoised_copy_ratios = DenoiseReadCounts.denoised_copy_ratios
        File? standardized_copy_ratios = DenoiseReadCounts.standardized_copy_ratios
    }
}

task CollectReadCountsTask {
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
        Boolean is_paired_end = false

        Runtime runtime_params
    }

    String output_name = sample_name + ".read_counts.tsv"

	command <<<
        set -e
        export GATK_LOCAL_JAR=~{select_first([runtime_params.jar_override, "/root/gatk.jar"])}
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            CollectReadCounts \
            -I '~{bam}' \
            -L '~{interval_list}' \
            -R '~{ref_fasta}' \
            -O '~{output_name}' \
            --interval-merging-rule ~{interval_merging_rule} \
            --format ~{format} \
            ~{if is_paired_end then "--read-filter FirstOfPairReadFilter " else ""} \
            ~{if is_paired_end then "--read-filter PairedReadFilter " else ""} \
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

    parameter_meta {
        interval_list: {localization_optional: true}
        ref_fasta: {localization_optional: true}
        ref_fasta_index: {localization_optional: true}
        ref_dict: {localization_optional: true}
        bam: {localization_optional: true}
        bai: {localization_optional: true}
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

    String output_denoised_copy_ratios = sample_name + ".denoised_CR.tsv"
    String output_standardized_copy_ratios = sample_name + ".standardized_CR.tsv"

	command <<<
        set -e
        export GATK_LOCAL_JAR=~{select_first([runtime_params.jar_override, "/root/gatk.jar"])}
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

    parameter_meta {
        # read_counts: {localization_optional: true}
        # annotated_interval_list: {localization_optional: true}
        # count_panel_of_normals: {localization_optional: true}
    }
}