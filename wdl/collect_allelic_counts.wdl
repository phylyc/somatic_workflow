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
        String sample_name = basename(bam, ".bam")
        File bam
        File bai
        Boolean? is_paired_end
        File? interval_list
        File? interval_list_idx
        File? interval_blacklist
        File? interval_blacklist_idx
        Array[File]? scattered_interval_list

        File? variants
        File? variants_idx
        File? vcf
        File? vcf_idx
        String? getpileupsummaries_extra_args

        Float minimum_population_allele_frequency = 0.0
        Float maximum_population_allele_frequency = 1.0
        Int minimum_read_depth = 0

        Boolean compress_output = false

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

    if (defined(vcf) && !defined(variants)) {
        call VcfToPileupVariants {
            input:
                vcf = select_first([vcf]),
                vcf_idx = select_first([vcf_idx]),
                interval_list = interval_list,
                interval_list_idx = interval_list_idx,
                interval_blacklist = interval_blacklist,
                interval_blacklist_idx = interval_blacklist_idx,
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
                    interval_list_idx = interval_list_idx,
                    interval_blacklist = interval_blacklist,
                    interval_blacklist_idx = interval_blacklist_idx,
                    scattered_intervals = scattered_intervals,
                    variants = select_first([variants, VcfToPileupVariants.variants]),
                    variants_idx = select_first([variants_idx, VcfToPileupVariants.variants_idx]),
                    getpileupsummaries_extra_args = getpileupsummaries_extra_args,
                    minimum_population_allele_frequency = minimum_population_allele_frequency,
                    maximum_population_allele_frequency = maximum_population_allele_frequency,
                    compress_output = false,
                    runtime_params = runtime_collection.get_pileup_summaries,
            }
        }

        call GatherPileupSummaries {
            input:
                input_tables = ScatteredGetPileupSummaries.pileup_summaries,
                ref_dict = ref_dict,
                sample_name = sample_name,
                compress_output = compress_output,
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
                interval_list_idx = interval_list_idx,
                interval_blacklist = interval_blacklist,
                interval_blacklist_idx = interval_blacklist_idx,
                sample_name = sample_name,
                variants = select_first([variants, VcfToPileupVariants.variants]),
                variants_idx = select_first([variants_idx, VcfToPileupVariants.variants_idx]),
                getpileupsummaries_extra_args = getpileupsummaries_extra_args,
                minimum_population_allele_frequency = minimum_population_allele_frequency,
                maximum_population_allele_frequency = maximum_population_allele_frequency,
                compress_output = compress_output,
                runtime_params = runtime_collection.get_pileup_summaries,
        }
    }

    if (minimum_read_depth > 0) {
        call SelectPileups {
            input:
                pileup_summaries = select_first([GatherPileupSummaries.merged_pileup_summaries, GetPileupSummaries.pileup_summaries]),
                sample_name = sample_name,
                minimum_read_depth = minimum_read_depth,
                compress_output = compress_output,
                runtime_params = runtime_collection.select_pileup_summaries,
        }
    }

    output {
        File? variants_from_vcf = VcfToPileupVariants.variants
        File? variants_from_vcf_idx = VcfToPileupVariants.variants_idx
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
        File? interval_list
        File? interval_list_idx
        File? interval_blacklist
        File? interval_blacklist_idx
        Array[String]? sample_names
        Float AF = 0.00000007
        Boolean compress_output = false

        Runtime runtime_params
    }

    Int diskGB = ceil(2 * size(vcf, "GB")) + runtime_params.disk

    String sample_name = basename(basename(basename(vcf, ".gz"), ".bgz"), ".vcf")
    String tmp_vcf = sample_name + ".tmp.vcf"
    String uncompressed_vcf = sample_name + ".af_only.vcf"
    String af_only_vcf = sample_name + ".af_only.vcf.gz"
    String af_only_vcf_idx = af_only_vcf + ".tbi"

    String dollar = "$"

    command <<<
        set -euxo pipefail

        # Prepare the AF INFO field for the header
        echo '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">' > "header_file"

        # Convert interval_list to BED if needed
        if [[ -n "~{interval_list}" ]]; then
            if [[ "~{interval_list}" == *.interval_list ]]; then
                grep -v "^@" "~{interval_list}" | awk '{print $1, $2-1, $3}' OFS='\t' > intervals.bed
                interval_bed="intervals.bed"
            else
                interval_bed="~{interval_list}"
            fi
        fi

        # Convert interval_blacklist to BED if needed
        if [[ -n "~{interval_blacklist}" ]]; then
            if [[ "~{interval_blacklist}" == *.interval_list ]]; then
                grep -v "^@" "~{interval_blacklist}" | awk '{print $1, $2-1, $3}' OFS='\t' > blacklist.bed
                blacklist_bed="blacklist.bed"
            else
                blacklist_bed="~{interval_blacklist}"
            fi
        fi

        # Build bcftools view filter arguments
        filter_args=()
        if [[ -n "~{dollar}{interval_bed:-}" ]]; then
            filter_args+=("-R" "$interval_bed")
        fi
        if [[ -n "~{dollar}{blacklist_bed:-}" ]]; then
            filter_args+=("-T" "^$blacklist_bed")
        fi

        # Filter the VCF file to retain only rows with genotypes
        # Remove FORMAT field and retain only INFO/POPAF field
        # Add header for AF field
        bcftools view -G "~{dollar}{filter_args[@]}" '~{vcf}' \
            | bcftools annotate -x FORMAT,^INFO/POPAF \
            | bcftools annotate -h "header_file" \
            > '~{tmp_vcf}'

        set +e  # grep returns 1 if no lines are found
        grep "^#" '~{tmp_vcf}' > '~{uncompressed_vcf}'

        # Set AF INFO field based on POPAF or default task variable AF
        grep -v "^#" '~{tmp_vcf}' \
            | awk -v default_AF='~{AF}' '
                BEGIN {OFS="\t"}
                {
                    # Check if POPAF is present in the INFO field
                    if ($8 ~ /^POPAF=/) {
                        # Extract POPAF value and calculate AF as 10^(-POPAF)
                        popaf_value = substr($8, 7)  # Extract value after 'POPAF='
                        AF_value = sprintf("%.6f", exp(log(10) * -popaf_value))
                    } else {
                        # Use default AF if POPAF is not found
                        AF_value = default_AF
                    }
                    # Replace INFO field with AF=calculated_value
                    $8 = "AF=" AF_value
                    print
                }' >> '~{uncompressed_vcf}'
        set -e

        # Compress the modified VCF file (bgzip is not available in docker)
        bcftools convert -O z -o '~{af_only_vcf}' '~{uncompressed_vcf}'

        # Index the compressed VCF file
        bcftools index -t -o '~{af_only_vcf_idx}' '~{af_only_vcf}'

        # Generate pileup tables for each sample
        if [ "~{defined(sample_names)}" == "true" ]; then
            for sample in ~{sep=" " sample_names}; do
                printf "#<METADATA>SAMPLE=$sample\n" > "$sample.pileup"
                bcftools query -s "$sample" -f '%CHROM\t%POS\t%INFO/POPAF\t[%DP\t%AD]\n' "~{vcf}" \
                    | awk -v default_AF='~{AF}' '
                        BEGIN {OFS="\t"; print "contig", "position", "ref_count", "alt_count", "other_alt_count", "allele_frequency"}
                        {
                            split($5, ad, ",");
                            ref_count = ad[1];
                            alt_count = ad[2];
                            total_depth = $4;
                            other_alt_count = total_depth - ref_count - alt_count;
                            allele_frequency = ($3 != "" ? exp(log(10) * -$3) : default_AF);  # Use default AF if missing
                            print $1, $2, ref_count, alt_count, other_alt_count, allele_frequency;
                        }' >> "$sample.pileup"
                if [ "~{compress_output}" == "true" ]; then
                    gzip -c "$sample.pileup" > "$sample.pileup.gz"
                    rm -f "$sample.pileup"
                fi
            done
        fi

        # Clean up temporary files
        rm -f '~{tmp_vcf}' '~{uncompressed_vcf}' 'header_file'
    >>>

    output {
        File variants = af_only_vcf
        File variants_idx = af_only_vcf_idx
        Array[File] pileups = glob("*.pileup" + (if compress_output then ".gz" else ""))
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
	input {
        File? interval_list
        File? interval_list_idx
        File? interval_blacklist
        File? interval_blacklist_idx
        File? scattered_intervals
        File input_bam
        File input_bai
        File variants
        File variants_idx
        Boolean is_paired_end = false
        String? sample_name
        String? getpileupsummaries_extra_args

        Float minimum_population_allele_frequency = 0.0
        Float maximum_population_allele_frequency = 1.0

        Boolean compress_output

        Runtime runtime_params
	}

    String sample_id = if defined(sample_name) then sample_name else basename(input_bam, ".bam")
    String pileup_file = sample_id + ".pileup"
    String output_file = pileup_file + if compress_output then ".gz" else ""

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{select_first([runtime_params.jar_override, "/root/gatk.jar"])}

        select_variants() {
            gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
                SelectVariants \
                "$@" \
                -O tmp.selected_loci.vcf
            mv tmp.selected_loci.vcf selected_loci.vcf
            mv tmp.selected_loci.vcf.idx selected_loci.vcf.idx

            set +e   # grep returns 1 if no lines are found
            num_loci=$(grep -v "^#" selected_loci.vcf | wc -l || echo 0)
            set -e
            echo ">> Selected $num_loci loci."
        }

        if [ "~{defined(scattered_intervals)}" == "true" ]; then
            select_variants \
                -V '~{variants}' \
                -L '~{scattered_intervals}'
        fi
        if [ "~{defined(interval_list)}" == "true" ]; then
            select_variants \
                -V '~{if defined(scattered_intervals) then "selected_loci.vcf" else variants}' \
                -L '~{interval_list}'
        fi
        if [ "~{defined(interval_blacklist)}" == "true" ]; then
            select_variants \
                -V '~{if defined(scattered_intervals) || defined(interval_list) then "selected_loci.vcf" else variants}' \
                -XL '~{interval_blacklist}'
        fi

        if [ -f selected_loci.vcf ] ; then
            set +e  # grep returns 1 if no lines are found
            num_loci=$(grep -v "^#" selected_loci.vcf | wc -l || echo 0)
            set -e
        else
            num_loci=1
        fi

        if [ "$num_loci" -eq 0 ] ; then
            # Create an empty pileup file if there are no variants in the intersection
            # between the variants and the intervals.
            printf "#<METADATA>SAMPLE=~{sample_id}\n" > '~{pileup_file}'
            printf "contig\tposition\tref_count\talt_count\tother_alt_count\tallele_frequency\n" >> '~{pileup_file}'
        else
            gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
                GetPileupSummaries \
                --input '~{input_bam}' \
                --read-index '~{input_bai}' \
                --intervals '~{if defined(scattered_intervals) || defined(interval_list) || defined(interval_blacklist) then "selected_loci.vcf" else variants}' \
                --variant '~{if defined(scattered_intervals) || defined(interval_list) || defined(interval_blacklist) then "selected_loci.vcf" else variants}' \
                -min-af '~{minimum_population_allele_frequency}' \
                -max-af '~{maximum_population_allele_frequency}' \
                --output '~{pileup_file}' \
                ~{if is_paired_end then "--read-filter FirstOfPairReadFilter " else ""} \
                ~{if is_paired_end then "--read-filter PairedReadFilter " else ""} \
                --seconds-between-progress-updates 60 \
                ~{getpileupsummaries_extra_args}
        fi

        if [ "~{compress_output}" == "true" ] ; then
            bgzip -c '~{pileup_file}' > '~{output_file}'
        fi
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
        interval_list_idx: {localization_optional: true}
        interval_blacklist: {localization_optional: true}
        interval_blacklist_idx: {localization_optional: true}
        scattered_intervals: {localization_optional: true}
        input_bam: {localization_optional: true}
        input_bai: {localization_optional: true}
        variants: {localization_optional: true}
        variants_idx: {localization_optional: true}
    }
}

task GatherPileupSummaries {
    input {
        Array[File] input_tables
        File ref_dict
        String sample_name

        Boolean compress_output

        Runtime runtime_params
    }

    # TODO: handle compressed input tables

    String pileup_file = sample_name + ".pileup"
    String output_file = pileup_file + if compress_output then ".gz" else ""

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{select_first([runtime_params.jar_override, "/root/gatk.jar"])}
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            GatherPileupSummaries \
            --sequence-dictionary '~{ref_dict}' \
            ~{sep="' " prefix("-I '", input_tables)}' \
            -O '~{pileup_file}'

        if [ "~{compress_output}" == "true" ] ; then
            bgzip -c '~{pileup_file}' > '~{output_file}'
        fi
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

        Boolean compress_output

        Runtime runtime_params
    }

    String uncompressed_pileup_summaries = basename(pileup_summaries, ".gz")
    Boolean is_compressed = uncompressed_pileup_summaries != basename(pileup_summaries)

    String pileup_file = sample_name + ".pileup"
    String output_file = pileup_file + if compress_output then ".gz" else ""
    String tmp_pileup_file = "tmp." + pileup_file

    command <<<
        set -uxo pipefail

        if [ "~{is_compressed}" == "true" ] ; then
            bgzip -cd '~{pileup_summaries}' > '~{uncompressed_pileup_summaries}'
        else
            mv '~{pileup_summaries}' '~{uncompressed_pileup_summaries}'
        fi

        # Extract leading comment lines
        grep '^#' '~{uncompressed_pileup_summaries}' > '~{tmp_pileup_file}'

        # Extract column headers
        grep -v '^#' '~{uncompressed_pileup_summaries}' | head -n 1 >> '~{tmp_pileup_file}'

        # Count the number of lines that are not comments (headers) or column headers
        num_variants=$(( $(grep -vc '^#' '~{uncompressed_pileup_summaries}') - 1 ))

        if [ "$num_variants" -gt 0 ]; then
            # Extract table and select lines with read depth >= min_read_depth
            grep -v '^#' '~{uncompressed_pileup_summaries}' | tail -n +2 \
                | awk -F"\t" '$3 + $4 + $5 >= ~{minimum_read_depth}' \
                >> '~{tmp_pileup_file}'
        fi

        mv '~{tmp_pileup_file}' '~{pileup_file}'

        # Count the number of lines that are not comments (headers) or column headers
        num_selected_variants=$(( $(grep -vc '^#' '~{pileup_file}') - 1 ))

        echo ">> Selected $num_selected_variants loci out of $num_variants."

        if [ "~{compress_output}" == "true" ] ; then
            bgzip -c '~{pileup_file}' > '~{output_file}'
        fi
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