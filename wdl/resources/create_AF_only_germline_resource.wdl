version development

import "../runtime_collection.wdl" as rtc
import "../runtimes.wdl" as rt
import "../tasks.wdl"


workflow CreateAFonlyVCF {
    input {
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        File vcf
        File vcf_idx

        Float minAF = 0.05

        Boolean compress_output = true

        RuntimeCollection runtime_collection = RuntimeParameters.rtc
    }

    call rtc.DefineRuntimeCollection as RuntimeParameters

    call SelectAFonly {
        input:
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            vcf = vcf,
            vcf_idx = vcf_idx,
            minAF = minAF,
            compress_output = compress_output,
            runtime_params = runtime_collection.select_af_only_from_vcf
    }

    output {
        File af_only_vcf = SelectAFonly.af_only_vcf
        File af_only_vcf_idx = SelectAFonly.af_only_vcf_idx
        File common_af_only_vcf = SelectAFonly.common_af_only_vcf
        File common_af_only_vcf_idx = SelectAFonly.common_af_only_vcf_idx
    }
}

task SelectAFonly {
    input {
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        File vcf
        File vcf_idx
        Float minAF = 0.05

        Boolean compress_output = true

        Runtime runtime_params
    }

    Int diskGB = runtime_params.disk + ceil(1.5 * size(vcf, "GB"))

    String output_base_name = basename(basename(basename(vcf, ".gz"), ".bgz"), ".vcf")
    String vcf_af_only = output_base_name + ".af_only.vcf" + if compress_output then ".gz" else ""
    String vcf_af_only_idx = vcf_af_only + if compress_output then ".tbi" else ".idx"
    String vcf_af_only_filtered = output_base_name + ".af_only.filtered.vcf" + if compress_output then ".gz" else ""
    String vcf_af_only_filtered_idx = vcf_af_only_filtered + if compress_output then ".tbi" else ".idx"
    String vcf_af_only_filtered_AFgtX = output_base_name + ".af_only.AFgt" + minAF + ".filtered.vcf" + if compress_output then ".gz" else ""
    String vcf_af_only_filtered_AFgtX_idx = vcf_af_only_filtered_AFgtX + if compress_output then ".tbi" else ".idx"
    String vcf_af_only_filtered_AFgtX_onlySNP = output_base_name + ".af_only.AFgt" + minAF + ".filtered.onlySNP.vcf" + if compress_output then ".gz" else ""
    String vcf_af_only_filtered_AFgtX_onlySNP_idx = vcf_af_only_filtered_AFgtX_onlySNP + if compress_output then ".tbi" else ".idx"
    String vcf_af_only_filtered_AFgtX_onlySNP_biallelic = output_base_name + ".af_only.AFgt" + minAF + ".filtered.onlySNP.biallelic.vcf" + if compress_output then ".gz" else ""
    String vcf_af_only_filtered_AFgtX_onlySNP_biallelic_idx = vcf_af_only_filtered_AFgtX_onlySNP_biallelic + if compress_output then ".tbi" else ".idx"

	command <<<
        set -e
        export GATK_LOCAL_JAR=~{select_first([runtime_params.jar_override, "/root/gatk.jar"])}

        echo "Subsetting to AF only data ... "
        bcftools annotate \
            -x FORMAT,^INFO/AC,INFO/AF \
            -o '~{vcf_af_only}' \
            '~{vcf}'

        echo "Creating index for vcf_af_only file ... "
        bcftools index -t \
            -o '~{vcf_af_only_idx}' \
            '~{vcf_af_only}'

        # Also outputs idx file
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            SelectVariants \
            -R '~{ref_fasta}' \
            -V '~{vcf_af_only}' \
            --output '~{vcf_af_only_filtered}' \
            --exclude-filtered true

        # make space
        rm '~{vcf_af_only}'
        rm '~{vcf_af_only_idx}'

        # --- Transformations for common germline resource ---
        echo "Subsetting to AF > ~{minAF} ... "
        bcftools view \
            -i "INFO/AF > ~{minAF}" \
            -o '~{vcf_af_only_filtered_AFgtX}' \
            '~{vcf_af_only_filtered}'

        echo "Creating index for vcf file ... "
        bcftools index -t \
            -o '~{vcf_af_only_filtered_AFgtX_idx}' \
            '~{vcf_af_only_filtered_AFgtX}'
        
        # Subset to only SNPs 
        echo "Subsetting to only SNPs ... "
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            SelectVariants \
            -R '~{ref_fasta}' \
            -V '~{vcf_af_only_filtered_AFgtX}' \
            --select-type-to-include SNP \
            -O '~{vcf_af_only_filtered_AFgtX_onlySNP}'
        
        # For multiallelic sites (>1 rows), keep first record
        echo "Making biallelic ... "
        bcftools sort '~{vcf_af_only_filtered_AFgtX_onlySNP}' | \
            bcftools norm --rm-dup both -Oz > '~{vcf_af_only_filtered_AFgtX_onlySNP_biallelic}'

        # Index the biallelic only SNP file
        echo "Creating index for biallelic only SNP vcf file ... "
        bcftools index -t \
            -o '~{vcf_af_only_filtered_AFgtX_onlySNP_biallelic_idx}' \
            '~{vcf_af_only_filtered_AFgtX_onlySNP_biallelic}'
	>>>

	output {
        File af_only_vcf = vcf_af_only_filtered
        File af_only_vcf_idx = vcf_af_only_filtered_idx
        File common_af_only_vcf = vcf_af_only_filtered_AFgtX_onlySNP_biallelic
        File common_af_only_vcf_idx = vcf_af_only_filtered_AFgtX_onlySNP_biallelic_idx
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

    parameter_meta {
        ref_fasta: {localization_optional: true}
        ref_fasta_index: {localization_optional: true}
        ref_dict: {localization_optional: true}
    }
}