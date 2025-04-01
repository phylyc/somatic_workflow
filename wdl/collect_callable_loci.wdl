version development

import "runtime_collection.wdl" as rtc
import "runtimes.wdl" as rt
import "sequencing_run.define.wdl" as seqrun_def


workflow CollectCallableLoci {
    input {
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        String sample_name = basename(bam, ".bam")
        File bam
        File bai
        Boolean is_paired_end = false

        Boolean compress_output = false

        RuntimeCollection runtime_collection = RuntimeParameters.rtc
        String gatk_docker = "broadinstitute/gatk"
        File? gatk_override
        Int preemptible = 1
        Int max_retries = 1
        # memory assignments in MB
        Int mem_callable_loci = 2048
        Int time_callable_loci = 300
    }

    call rtc.DefineRuntimeCollection as RuntimeParameters {
        input:
            gatk_docker = gatk_docker,
            gatk_override = gatk_override,
            max_retries = max_retries,
            preemptible = preemptible,
            mem_callable_loci = mem_callable_loci,
            time_callable_loci = time_callable_loci
    }

	call CollectCallableLociTask {
		input:
            ref_fasta = ref_fasta,
			ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            bam = bam,
            bai = bai,
            sample_name = sample_name,
            is_paired_end = is_paired_end,
            compress_output = compress_output,
            runtime_params = runtime_collection.collect_callable_loci
	}

    output {
        File bed = CollectCallableLociTask.bed
        File summary = CollectCallableLociTask.summary
    }
}

task CollectCallableLociTask {
    input {
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        File bam
        File bai
        String sample_name
        Boolean is_paired_end = false

        Boolean compress_output = false

        Runtime runtime_params
    }

    String bed_out = sample_name + ".callable.bed"
    String bed_out_gz = sample_name + ".callable.bed" + if compress_output then ".gz" else ""
    String summary_out = sample_name + ".callable.summary.txt"

	command <<<
        set -e
        export GATK_LOCAL_JAR=~{select_first([runtime_params.jar_override, "/root/gatk.jar"])}
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            CallableLoci \
            -I '~{bam}' \
            --read-index '~{bai}' \
            -R '~{ref_fasta}' \
            -O '~{bed_out}' \
            --summary '~{summary_out}' \
            ~{if is_paired_end then "--read-filter FirstOfPairReadFilter " else ""} \
            ~{if is_paired_end then "--read-filter PairedReadFilter " else ""} \
            --seconds-between-progress-updates 60

        if [ "~{compress_output}" == "true" ] ; then
            bgzip -c '~{bed_out}' > '~{bed_out_gz}'
        fi
	>>>

	output {
		File bed = bed_out_gz
        File summary = summary_out
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
        ref_fasta: {localization_optional: true}
        ref_fasta_index: {localization_optional: true}
        ref_dict: {localization_optional: true}
        bam: {localization_optional: true}
        bai: {localization_optional: true}
    }
}