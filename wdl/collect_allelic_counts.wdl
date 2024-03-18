version development

## Get pileup summaries of allelic read counts for a bam file at variant sites.
##
## Known issues:
## - GetPileupSummaries ignores certain MNVs and INDELs which do not appear in
##      the pileup output. This tool should eventually be replaced by CollectByBaseCounts.
##      https://github.com/broadinstitute/gatk/pull/6545
## - @Robert Klein, 2020: "Paired-end reads that overlap at some sites of interest
##      lead to double counting. This process in general is more of an issue in
##      cell-free DNA where the vast majority of templates are ~166bp long, which
##      is shorter than twice the read length. The result is that many bases on
##      a given template are reported twice, once from each paired-end read."
##      -> Solve by using FirstOfPairReadFilter?

import "runtime_collection.wdl" as rtc
import "runtimes.wdl" as rt


workflow CollectAllelicCounts {
	input {
        File ref_dict

        # SequencingRun sequencing_run = GetSeqRun.sequencing_run
        String sample_name
        File bam
        File bai
        Boolean? is_paired_end
        File? interval_list
        File? interval_blacklist
        Array[File]? scattered_interval_list

        File? common_germline_alleles
        File? common_germline_alleles_idx
        File? vcf
        File? vcf_idx
        String? getpileupsummaries_extra_args

        Float minimum_population_allele_frequency = 0.01
        Float maximum_population_allele_frequency = 0.2
        Int minimum_read_depth = 0

        RuntimeCollection runtime_collection = RuntimeParameters.rtc

        String bcftools_docker = "stephb/bcftools"
        String gatk_docker = "broadinstitute/gatk"
        String ubuntu_docker = "ubuntu"
        File? gatk_override
        Int preemptible = 1
        Int max_retries = 1
        Int disk_sizeGB = 1

        Int mem_vcf_to_pileup_variants = 512  # 64
        Int mem_get_pileup_summaries = 4096  # needs at least 2G
        Int mem_gather_pileup_summaries = 512  # 64
        Int mem_select_pileup_summaries = 256  # 64
        Int time_startup = 10
        Int time_vcf_to_pileup_variants = 5
        Int time_get_pileup_summaries = 90  # 1.5 h
        Int time_gather_pileup_summaries = 5
        Int time_select_pileup_summaries = 5
	}

    call rtc.DefineRuntimeCollection as RuntimeParameters {
        input:
            bcftools_docker = bcftools_docker,
            gatk_docker = gatk_docker,
            ubuntu_docker = ubuntu_docker,
            gatk_override = gatk_override,
            preemptible = preemptible,
            max_retries = max_retries,
            disk_sizeGB = disk_sizeGB,
            mem_vcf_to_pileup_variants = mem_vcf_to_pileup_variants,
            mem_get_pileup_summaries = mem_get_pileup_summaries,
            mem_gather_pileup_summaries = mem_gather_pileup_summaries,
            mem_select_pileup_summaries = mem_select_pileup_summaries,
            time_startup = time_startup,
            time_vcf_to_pileup_variants = time_vcf_to_pileup_variants,
            time_get_pileup_summaries = time_get_pileup_summaries,
            time_gather_pileup_summaries = time_gather_pileup_summaries,
            time_select_pileup_summaries = time_select_pileup_summaries,
    }

    if (defined(vcf) && !defined(common_germline_alleles)) {
        call VcfToPileupVariants {
            input:
                vcf = select_first([vcf]),
                vcf_idx = select_first([vcf_idx]),
                runtime_params = runtime_collection.vcf_to_pileup_variants,
        }
    }

    if (defined(scattered_interval_list)) {
        scatter (scattered_intervals in select_all(select_first([scattered_interval_list]))) {
            call GetPileupSummaries as ScatteredGetPileupSummaries {
                input:
                    input_bam = bam,
                    input_bai = bai,
                    is_paired_end = is_paired_end,
                    interval_list = interval_list,
                    interval_blacklist = interval_blacklist,
                    scattered_intervals = scattered_intervals,
                    common_germline_alleles = select_first([common_germline_alleles, VcfToPileupVariants.common_germline_alleles]),
                    common_germline_alleles_idx = select_first([common_germline_alleles_idx, VcfToPileupVariants.common_germline_alleles_idx]),
                    getpileupsummaries_extra_args = getpileupsummaries_extra_args,
                    minimum_population_allele_frequency = minimum_population_allele_frequency,
                    maximum_population_allele_frequency = maximum_population_allele_frequency,
                    runtime_params = runtime_collection.get_pileup_summaries,
            }
        }

        call GatherPileupSummaries {
            input:
                input_tables = ScatteredGetPileupSummaries.pileup_summaries,
                ref_dict = ref_dict,
                sample_name = sample_name,
                runtime_params = runtime_collection.gather_pileup_summaries,
        }
    }
    # else
    if (!defined(scattered_interval_list)) {
        call GetPileupSummaries {
            input:
                input_bam = bam,
                input_bai = bai,
                is_paired_end = is_paired_end,
                interval_list = interval_list,
                interval_blacklist = interval_blacklist,
                sample_name = sample_name,
                common_germline_alleles = select_first([common_germline_alleles, VcfToPileupVariants.common_germline_alleles]),
                common_germline_alleles_idx = select_first([common_germline_alleles_idx, VcfToPileupVariants.common_germline_alleles_idx]),
                getpileupsummaries_extra_args = getpileupsummaries_extra_args,
                minimum_population_allele_frequency = minimum_population_allele_frequency,
                maximum_population_allele_frequency = maximum_population_allele_frequency,
                runtime_params = runtime_collection.get_pileup_summaries,
        }
    }

    if (minimum_read_depth > 0) {
        call SelectPileups {
            input:
                pileup_summaries = select_first([GatherPileupSummaries.merged_pileup_summaries, GetPileupSummaries.pileup_summaries]),
                sample_name = sample_name,
                minimum_read_depth = minimum_read_depth,
                runtime_params = runtime_collection.select_pileup_summaries,
        }
    }

    output {
        File pileup_summaries = select_first([SelectPileups.selected_pileup_summaries, GatherPileupSummaries.merged_pileup_summaries, GetPileupSummaries.pileup_summaries])
    }
}

task VcfToPileupVariants {
    # Input: a (multi-sample) VCF, e.g. from Mutect2
    # create a VCF with AF for all common_germline_alleles in the input VCF, dropping all
    # samples, resulting in a gnomad-style VCF with only AF in the INFO field.

    input {
        File vcf
        File vcf_idx
        Float AF = 0.1  # AF that GetPileupSummaries will consider is by default between (0.01, 0.2)

        Runtime runtime_params
    }

    Int diskGB = ceil(2 * size(vcf, "GB")) + runtime_params.disk

    String sample_name = basename(basename(basename(vcf, ".gz"), ".bgz"), ".vcf")
    String tmp_vcf = sample_name + ".tmp.vcf"
    String uncompressed_vcf = sample_name + ".af_only.vcf"
    String af_only_vcf = sample_name + ".af_only.vcf.gz"
    String af_only_vcf_idx = af_only_vcf + ".tbi"

    command <<<
        set -euxo pipefail

        # Filter the VCF file to retain only rows with genotypes
        # Remove FORMAT field and retain only INFO/AF field
        bcftools view -G '~{vcf}' \
            | bcftools annotate -x FORMAT,^INFO/AF \
            > '~{tmp_vcf}'

        set +e  # grep returns 1 if no lines are found
        # Separate header lines (lines starting with '#') into a separate file
        grep "^#" '~{tmp_vcf}' > '~{uncompressed_vcf}'

        # Filter and modify non-header lines to include AF information
        # Use AWK to set the 'AF' INFO field to the specified value
        grep -v "^#" '~{tmp_vcf}' \
            | awk 'BEGIN {OFS="\t"} {$8 = "AF=~{AF}"; print}' \
            >> '~{uncompressed_vcf}'
        set -e

        # Compress the modified VCF file (bgzip is not available)
        bcftools convert -O z -o '~{af_only_vcf}' '~{uncompressed_vcf}'

        # Index the compressed VCF file
        bcftools index -t -o '~{af_only_vcf_idx}' '~{af_only_vcf}'

        # Clean up temporary files
        rm -f '~{tmp_vcf}' '~{uncompressed_vcf}'
    >>>

    output {
        File common_germline_alleles = af_only_vcf
        File common_germline_alleles_idx = af_only_vcf_idx
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

task GetPileupSummaries {
    # If the common_germline_alleles for contamination and the intervals for this scatter don't
    # intersect, GetPileupSummaries throws an error. However, there is nothing wrong
    # with an empty intersection for our purposes; it simply doesn't contribute to the
    # merged pileup summaries that we create downstream. We implement this by creating
    # an empty pileup file that is overwritten by GetPileupSummaries and force an
    # exit code of 0.

	input {
        File? interval_list
        File? interval_blacklist
        File? scattered_intervals
        File input_bam
        File input_bai
        File common_germline_alleles
        File common_germline_alleles_idx
        Boolean is_paired_end = false
        String? sample_name
        String? getpileupsummaries_extra_args

        Float minimum_population_allele_frequency = 0.01
        Float maximum_population_allele_frequency = 0.2
        Int min_mapping_quality = 50

        Runtime runtime_params
	}

    String sample_id = if defined(sample_name) then sample_name else basename(input_bam, ".bam")
    String output_file = sample_id + ".pileup"

    command <<<
        set +e
        export GATK_LOCAL_JAR=~{select_first([runtime_params.jar_override, "/root/gatk.jar"])}

        # Create an empty pileup file if there are no common_germline_alleles in the intersection
        # between the common_germline_alleles and the intervals. Will be overwritten by GetPileupSummaries
        printf "#<METADATA>SAMPLE=~{sample_id}\n" > '~{output_file}'
        printf "contig\tposition\tref_count\talt_count\tother_alt_count\tallele_frequency\n" >> '~{output_file}'

        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            GetPileupSummaries \
            --input '~{input_bam}' \
            ~{"--intervals '" +  interval_list + "'"} \
            ~{"--intervals '" +  scattered_intervals + "'"} \
            ~{"--exclude-intervals '" +  interval_blacklist + "'"} \
            --intervals '~{common_germline_alleles}' \
            --interval-set-rule INTERSECTION \
            --variant '~{common_germline_alleles}' \
            -min-af '~{minimum_population_allele_frequency}' \
            -max-af '~{maximum_population_allele_frequency}' \
            --min-mapping-quality ~{min_mapping_quality} \
            --output '~{output_file}' \
            ~{if is_paired_end then "--read-filter FirstOfPairReadFilter " else ""} \
            ~{if is_paired_end then "--read-filter PairedReadFilter " else ""} \
            --seconds-between-progress-updates 60 \
            ~{getpileupsummaries_extra_args}

        # It only fails due to empty intersection between common_germline_alleles and intervals, which is ok.
        exit 0
    >>>

    output {
        File pileup_summaries = output_file
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
        input_bam: {localization_optional: true}
        input_bai: {localization_optional: true}
        common_germline_alleles: {localization_optional: true}
        common_germline_alleles_idx: {localization_optional: true}
    }
}

task GatherPileupSummaries {
    input {
        Array[File] input_tables
        File ref_dict
        String sample_name

        Runtime runtime_params
    }

    String output_file = sample_name + ".pileup"

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{select_first([runtime_params.jar_override, "/root/gatk.jar"])}
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            GatherPileupSummaries \
            --sequence-dictionary '~{ref_dict}' \
            ~{sep="' " prefix("-I '", input_tables)}' \
            -O '~{output_file}'
    >>>

    output {
        File merged_pileup_summaries = output_file
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
    #     input_tables: {localization_optional: true}
    #     ref_dict: {localization_optional: true}
    # }
}

task SelectPileups {
    input {
        File pileup_summaries
        String sample_name
        Int minimum_read_depth

        Runtime runtime_params
    }

    String output_file = sample_name + ".pileup"
    String tmp_output_file = "tmp." + output_file

    command <<<
        set -uxo pipefail
        # Extract leading comment lines
        grep '^#' '~{pileup_summaries}' > '~{tmp_output_file}'

        # Extract column headers
        grep -v '^#' '~{pileup_summaries}' | head -n 1 >> '~{tmp_output_file}'

        # Count the number of lines that are not comments (headers)
        num_variants_plus_one=$(grep -vc '^#' '~{pileup_summaries}')

        if [ "$num_variants_plus_one" -gt 1 ]; then
            # Extract table and select lines with read depth >= min_read_depth
            grep -v '^#' '~{pileup_summaries}' | tail -n +2 \
                | awk -F"\t" '$3 + $4 + $5 >= ~{minimum_read_depth}' \
                >> '~{tmp_output_file}'
        fi

        mv '~{tmp_output_file}' '~{output_file}'
    >>>

    output {
        File selected_pileup_summaries = output_file
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