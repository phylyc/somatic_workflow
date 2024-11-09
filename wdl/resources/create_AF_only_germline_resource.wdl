version development

import "runtime_collection.wdl" as rtc
import "runtimes.wdl" as rt
import "tasks.wdl"


workflow CreateAFonlyVCF {
    input {
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        File vcf
        File vcf_idx

        Boolean compress_output = true
        Boolean create_biallelic = false
        Boolean create_multiallelic = false

        RuntimeCollection runtime_collection = RuntimeParameters.rtc
    }

    call rtc.DefineRuntimeCollection as RuntimeParameters

    call SelectAFonly {
        input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            compress_output = true,  # need to compress as index can only be created for compressed vcf
            runtime_params = runtime_collection.select_af_only_from_vcf
    }

    call SelectVariants {
        input:
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            vcf = SelectAFonly.af_only_vcf,
            vcf_idx = SelectAFonly.af_only_vcf_idx,
            compress_output = compress_output,
            create_biallelic = create_biallelic,
            create_multiallelic = create_multiallelic,
            runtime_params = runtime_collection.select_variants
    }

    output {
        File af_only_vcf = SelectVariants.filtered_vcf
        File af_only_vcf_idx = SelectVariants.filtered_vcf_idx
  		File? biallelic_af_only_vcf = SelectVariants.biallelic_af_only_vcf
  		File? biallelic_af_only_vcf_idx = SelectVariants.biallelic_af_only_vcf_idx
        File? multiallelic_af_only_vcf = SelectVariants.multiallelic_af_only_vcf
        File? multiallelic_af_only_vcf_idx = SelectVariants.multiallelic_af_only_vcf_idx
    }
}

task SelectAFonly {
    input {
        File vcf
        File vcf_idx
        Boolean compress_output = true
        String? select_variants_extra_args

        Runtime runtime_params
    }

    parameter_meta {
#        vcf: {localization_optional: true}
#        vcf_idx: {localization_optional: true}
    }

    Int diskGB = runtime_params.disk + ceil(1.5 * size(vcf, "GB"))

    String output_base_name = basename(basename(basename(vcf, ".gz"), ".bgz"), ".vcf")
    String af_only_vcf_ = output_base_name + ".af_only" + if compress_output then ".vcf.gz" else ".vcf"
    String af_only_vcf_idx_ = af_only_vcf_ + if compress_output then ".tbi" else ".idx"

	command <<<
        set -e
        echo "Subsetting to AF only data ... "
        bcftools annotate -x FORMAT,^INFO/AC,INFO/AF -o ~{af_only_vcf_} ~{vcf}
        echo "Creating index for vcf file ... "
        bcftools index -t -o ~{af_only_vcf_idx_} ~{af_only_vcf_}
	>>>

	output {
        File af_only_vcf = af_only_vcf_
        File af_only_vcf_idx = af_only_vcf_idx_
	}

    runtime {
        docker: runtime_params.docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: runtime_params.machine_mem + " MB"
        runtime_minutes: runtime_params.runtime_minutes
        disks: "local-disk " + diskGB + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }
}

task SelectVariants {
    input {
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        File vcf
        File vcf_idx
        Boolean compress_output = true
        Boolean create_biallelic = false
        Boolean create_multiallelic = false
        String? select_variants_extra_args

        Runtime runtime_params
    }

    parameter_meta {
        ref_fasta: {localization_optional: true}
        ref_fasta_index: {localization_optional: true}
        ref_dict: {localization_optional: true}
        vcf: {localization_optional: true}
        vcf_idx: {localization_optional: true}
    }

    Int diskGB = 3 * ceil(size(vcf, "GB"))

    String output_base_name = basename(basename(basename(vcf, ".gz"), ".bgz"), ".vcf")
    String filtered_vcf_ = output_base_name + ".filtered" + if compress_output then ".vcf.gz" else ".vcf"
    String filtered_vcf_idx_ = filtered_vcf_ + if compress_output then ".tbi" else ".idx"
    String biallelic_af_only_vcf_ = output_base_name + ".filtered.biallelic" + if compress_output then ".vcf.gz" else ".vcf"
    String biallelic_af_only_vcf_idx_ = biallelic_af_only_vcf_ + if compress_output then ".tbi" else ".idx"
    String multiallelic_af_only_vcf_ = output_base_name + ".filtered.multiallelic" + if compress_output then ".vcf.gz" else ".vcf"
    String multiallelic_af_only_vcf_idx_ = multiallelic_af_only_vcf_ + if compress_output then ".tbi" else ".idx"

	command <<<
        set -e
        export GATK_LOCAL_JAR=~{select_first([runtime_params.jar_override, "/root/gatk.jar"])}
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            SelectVariants \
            -R ~{ref_fasta} \
            -V ~{vcf} \
            --output ~{filtered_vcf_} \
            --exclude-filtered true \
            ~{select_variants_extra_args}

        if [ "~{create_biallelic}" == "true" ] ; then
            gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
                SelectVariants \
                -R ~{ref_fasta} \
                -V ~{filtered_vcf_} \
                --output ~{biallelic_af_only_vcf_} \
                --restrict-alleles-to BIALLELIC \
                ~{select_variants_extra_args}
        fi

        if [ "~{create_multiallelic}" == "true" ] ; then
            gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
                SelectVariants \
                -R ~{ref_fasta} \
                -V ~{filtered_vcf_} \
                --output ~{multiallelic_af_only_vcf_} \
                --restrict-alleles-to MULTIALLELIC \
                ~{select_variants_extra_args}
        fi
	>>>

	output {
        File filtered_vcf = filtered_vcf_
        File filtered_vcf_idx = filtered_vcf_idx_
  		File? biallelic_af_only_vcf = biallelic_af_only_vcf_
  		File? biallelic_af_only_vcf_idx = biallelic_af_only_vcf_idx_
        File? multiallelic_af_only_vcf = multiallelic_af_only_vcf_
        File? multiallelic_af_only_vcf_idx = multiallelic_af_only_vcf_idx_
	}

    runtime {
        docker: runtime_params.docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: runtime_params.machine_mem + " MB"
        runtime_minutes: runtime_params.runtime_minutes
        disks: "local-disk " + diskGB + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }
}