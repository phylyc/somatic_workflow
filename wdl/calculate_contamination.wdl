version development

## Calculate percent contamination of a sample by another sample from another patient.
##
## Known issues:
## The contamination estimates are based on allelic pileups at biallelic SNP panel
## sites. It is important to subset the SNP panel by supplying an interval list
## (scattered or not) that contains only the sites that are callable in the sample.
## This can be approximated by the target intervals.
## Additionally, GetPileupSummaries allows for the selection of a read depth minimum
## for the returned pileups.
## Too many low coverage sites will lead to inaccurate estimates of contamination.

import "runtime_collection.wdl" as rtc
import "runtimes.wdl" as rt
import "collect_allelic_counts.wdl" as cac


workflow CalculateContamination {
    input {
        File? interval_list
        File? interval_blacklist
        Array[File]? scattered_interval_list

        File? ref_dict

        File? tumor_bam
        File? tumor_bai
        String? tumor_sample_name
        File? tumor_pileups

        File? normal_bam
        File? normal_bai
        String? normal_sample_name
        File? normal_pileups

        File? common_germline_alleles  # AF-only biallelic VCF file with germline alleles
        File? common_germline_alleles_idx
        File? vcf  # VCF file with germline alleles
        File? vcf_idx
        String? getpileupsummaries_extra_args

        # parameters for GetPileupSummaries to select loci of common_germline_alleles for contamination
        Float minimum_population_allele_frequency = 0.01
        Float maximum_population_allele_frequency = 0.2
        Int minimum_read_depth = 10

        RuntimeCollection runtime_collection = GetRTC.rtc

        String bcftools_docker = "stephb/bcftools"
        String gatk_docker = "broadinstitute/gatk"
        String ubuntu_docker = "ubuntu"
        File? gatk_override
        Int max_retries = 1
        Int preemptible = 1
        Int cpu = 1
        Int disk_sizeGB = 1

        Int mem_vcf_to_pileup_variants = 512  # 64
        Int mem_get_pileup_summaries = 4096  # needs at least 2G
        Int mem_gather_pileup_summaries = 512  # 64
        Int mem_select_pileup_summaries = 256  # 64
        Int mem_calculate_contamination = 8192  # depends on the common_germline_alleles resource
        Int time_startup = 10
        Int time_vcf_to_pileup_variants = 5
        Int time_get_pileup_summaries = 90  # 1.5 h
        Int time_gather_pileup_summaries = 5
        Int time_select_pileup_summaries = 5
        Int time_calculate_contamination = 300
    }

    # todo: assert tumor_bam and tumor_bai or tumor_pileups is defined.

    Int scatter_count = if defined(scattered_interval_list) then length(select_first([scattered_interval_list])) else 1

    call rtc.DefineRuntimeCollection as GetRTC {
        input:
            scatter_count = scatter_count,
            bcftools_docker = bcftools_docker,
            gatk_docker = gatk_docker,
            ubuntu_docker = ubuntu_docker,
            gatk_override = gatk_override,
            max_retries = max_retries,
            preemptible = preemptible,
            cpu = cpu,
            disk_sizeGB = disk_sizeGB,
            mem_vcf_to_pileup_variants = mem_vcf_to_pileup_variants,
            mem_get_pileup_summaries = mem_get_pileup_summaries,
            mem_gather_pileup_summaries = mem_gather_pileup_summaries,
            mem_select_pileup_summaries = mem_select_pileup_summaries,
            mem_calculate_contamination = mem_calculate_contamination,
            time_startup = time_startup,
            time_vcf_to_pileup_variants = time_vcf_to_pileup_variants,
            time_get_pileup_summaries = time_get_pileup_summaries,
            time_gather_pileup_summaries = time_gather_pileup_summaries,
            time_select_pileup_summaries = time_select_pileup_summaries,
            time_calculate_contamination = time_calculate_contamination,
    }

    if (!defined(tumor_pileups)) {
        # todo: assert ref_dict and sample_name is defined
        # todo: assert common_germline_alleles or vcf is defined
        call cac.CollectAllelicCounts as TumorPileupSummaries {
            input:
                interval_list = interval_list,
                interval_blacklist = interval_blacklist,
                scattered_interval_list = scattered_interval_list,
                bam = select_first(select_all([tumor_bam])),
                bai = select_first(select_all([tumor_bai])),
                ref_dict = select_first(select_all([ref_dict])),
                sample_name = select_first([tumor_sample_name, "tumor"]),
                common_germline_alleles = common_germline_alleles,
                common_germline_alleles_idx = common_germline_alleles_idx,
                vcf = vcf,
                vcf_idx = vcf_idx,
                minimum_population_allele_frequency = minimum_population_allele_frequency,
                maximum_population_allele_frequency = maximum_population_allele_frequency,
                minimum_read_depth = minimum_read_depth,
                getpileupsummaries_extra_args = getpileupsummaries_extra_args,
                runtime_collection = runtime_collection,
        }
    }

    if (!defined(normal_pileups) && defined(normal_bam)) {
        call cac.CollectAllelicCounts as NormalPileupSummaries {
            input:
                interval_list = interval_list,
                interval_blacklist = interval_blacklist,
                scattered_interval_list = scattered_interval_list,
                bam = select_first(select_all([normal_bam])),
                bai = select_first(select_all([normal_bai])),
                ref_dict = select_first(select_all([ref_dict])),
                sample_name = select_first([normal_sample_name, "normal"]),
                common_germline_alleles = common_germline_alleles,
                common_germline_alleles_idx = common_germline_alleles_idx,
                vcf = vcf,
                vcf_idx = vcf_idx,
                minimum_population_allele_frequency = minimum_population_allele_frequency,
                maximum_population_allele_frequency = maximum_population_allele_frequency,
                minimum_read_depth = minimum_read_depth,
                getpileupsummaries_extra_args = getpileupsummaries_extra_args,
                runtime_collection = runtime_collection,
        }
    }

    File non_optional_tumor_pileups = select_first([tumor_pileups, TumorPileupSummaries.pileup_summaries])
    if (defined(normal_pileups) || defined(NormalPileupSummaries.pileup_summaries)) {
        File optional_normal_pileup_summaries = select_first([normal_pileups, NormalPileupSummaries.pileup_summaries])
    }

    call CalculateContaminationTask {
        input:
            tumor_pileups = non_optional_tumor_pileups,
            normal_pileups = optional_normal_pileup_summaries,
            runtime_params = runtime_collection.calculate_contamination,
    }

    output {
        File tumor_pileup_summaries = non_optional_tumor_pileups
        File? normal_pileup_summaries = optional_normal_pileup_summaries
        File contamination_table = CalculateContaminationTask.contamination_table
        File segmentation = CalculateContaminationTask.segmentation
    }
}

task CalculateContaminationTask {
    input {
        File tumor_pileups
        File? normal_pileups

        Runtime runtime_params
    }

    String tumor_sample_id = basename(tumor_pileups, ".pileup")
    String output_contamination = tumor_sample_id + ".contamination"
    String output_segments = tumor_sample_id + ".segments"

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{select_first([runtime_params.jar_override, "/root/gatk.jar"])}
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            CalculateContamination \
            --input '~{tumor_pileups}' \
            ~{"--matched-normal '" + normal_pileups + "'"} \
            --output '~{output_contamination}' \
            --tumor-segmentation '~{output_segments}'
    >>>

    output {
        File contamination_table = output_contamination
        File segmentation = output_segments
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

    # Optional localization leads to cromwell error.
    # parameter_meta {
    #     tumor_pileups: {localization_optional: true}
    #     normal_pileups: {localization_optional: true}
    # }
}