version development

import "runtimes.wdl"
import "tasks.wdl"


workflow SelectVariants {
    input {
        File? interval_list
        File? ref_fasta
        File? ref_fasta_index
        File? ref_dict
        File vcf
        File vcf_idx
        Boolean select_passing = false
        Boolean keep_germline = false
        Boolean compress_output = false
        String? tumor_sample_name
        String? normal_sample_name
        String? select_variants_extra_args

        RuntimeCollection runtime_collection = GetRTC.rtc

        String gatk_docker = "broadinstitute/gatk"
        Int preemptible = 1
        Int max_retries = 1
        Int cpu = 1
        Int disk_sizeGB = 1
        Int boot_disk_size = 12  # needs to be > 10
        Int mem_machine_overhead = 512
        Int time_startup = 10

        Int mem_select_variants = 2048
        Int time_select_variants = 5

    }

    call runtimes.DefineRuntimeCollection as GetRTC {
        input:
            gatk_docker = gatk_docker,
            preemptible = preemptible,
            max_retries = max_retries,
            cpu = cpu,
            disk_sizeGB = disk_sizeGB,
            boot_disk_size = boot_disk_size,
            mem_machine_overhead = mem_machine_overhead,
            time_startup = time_startup,
            mem_select_variants = mem_select_variants,
            time_select_variants = time_select_variants
    }

    call tasks.SelectVariants {
        input:
            interval_list = interval_list,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            vcf = vcf,
            vcf_idx = vcf_idx,
            select_passing = select_passing,
            keep_germline = keep_germline,
            compress_output = compress_output,
            tumor_sample_name = tumor_sample_name,
            normal_sample_name = normal_sample_name,
            select_variants_extra_args = select_variants_extra_args,
            runtime_params = runtime_collection.select_variants
    }

    output {
        File selected_vcf = SelectVariants.selected_vcf
        File selected_vcf_idx = SelectVariants.selected_vcf_idx
        Int num_selected_variants = SelectVariants.num_selected_variants
    }
}