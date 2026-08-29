version development

import "runtime_collection.wdl" as rtc
import "tasks.wdl"


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
        File? gatk_override
        Int max_retries = 1
    }

    call rtc.DefineRuntimeCollection as RuntimeParameters {
        input:
            gatk_override = gatk_override,
            max_retries = max_retries,
    }

	call tasks.CollectCallableLoci as CollectCallableLociTask {
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
