version development

## Collection of Tasks

import "runtimes.wdl" as rt


task GetSampleName {
    input {
        File bam
        File bai

        Runtime runtime_params
    }

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{select_first([runtime_params.jar_override, "/root/gatk.jar"])}
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            GetSampleName \
            -I '~{bam}' \
            --read-index '~{bai}' \
            -O bam_name.txt
    >>>

    output {
        String sample_name = read_string("bam_name.txt")
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
        bam: {localization_optional: true}
        bai: {localization_optional: true}
    }
}

task AnnotateIntervals {
    input {
        File interval_list

        File ref_fasta
        File ref_fasta_index
        File ref_dict

        File? mappability_track
        File? mappability_track_idx
        File? segmental_duplication_track
        File? segmental_duplication_track_idx

        Runtime runtime_params
    }

    String output_file = basename(interval_list, ".interval_list") + ".annotated.interval_list"

	command <<<
        set -e
        export GATK_LOCAL_JAR=~{select_first([runtime_params.jar_override, "/root/gatk.jar"])}
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            AnnotateIntervals \
            -R '~{ref_fasta}' \
            -L '~{interval_list}' \
            -O '~{output_file}' \
            --interval-merging-rule OVERLAPPING_ONLY \
            ~{"--mappability-track " + mappability_track} \
            ~{"--segmental-duplication-track " + segmental_duplication_track}
	>>>

	output {
		File annotated_interval_list = output_file
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
        mappability_track: {localization_optional: true}
        mappability_track_idx: {localization_optional: true}
        segmental_duplication_track: {localization_optional: true}
        segmental_duplication_track_idx: {localization_optional: true}
    }
}

task PreprocessIntervals {
    input {
        File? interval_list
        File? interval_blacklist
        Array[File]? interval_lists
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        Int bin_length = 0
        Int padding = 0
        String? preprocess_intervals_extra_args

        Runtime runtime_params
    }

    String preprocessed_intervals = "preprocessed.interval_list"

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{select_first([runtime_params.jar_override, "/root/gatk.jar"])}
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            PreprocessIntervals \
            -R '~{ref_fasta}' \
            ~{"-L '" + interval_list + "'"} \
            ~{"-XL '" + interval_blacklist + "'"} \
            ~{true="-L '" false="" defined(interval_lists)}~{default="" sep="' -L '" interval_lists}~{true="'" false="" defined(interval_lists)} \
            --bin-length ~{bin_length} \
            --padding ~{padding} \
            --interval-merging-rule OVERLAPPING_ONLY \
            -O '~{preprocessed_intervals}' \
            ~{preprocess_intervals_extra_args}
    >>>

    output {
        File preprocessed_interval_list = preprocessed_intervals
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
        interval_lists: {localization_optional: true}
        ref_fasta: {localization_optional: true}
        ref_fasta_index: {localization_optional: true}
    }
}

task SplitIntervals {
    input {
        File? interval_list
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        Int scatter_count
        String? split_intervals_extra_args

        Runtime runtime_params
    }

    String extra_args = (
        select_first([split_intervals_extra_args, ""])
        # to avoid splitting intervals:
        # + " --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW"
        # Applied after inital scatter, so leads to more scattered intervals.
        # + " --dont-mix-contigs"
    )

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{select_first([runtime_params.jar_override, "/root/gatk.jar"])}
        mkdir interval-files
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            SplitIntervals \
            -R '~{ref_fasta}' \
            ~{"-L '" + interval_list + "'"} \
            -scatter ~{scatter_count} \
            -O interval-files \
            ~{extra_args}
    >>>

    output {
        Array[File] interval_files = glob("interval-files/*.interval_list")
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
    }
}

#task VariantFiltration {
#    input {
#        File? interval_list
#        File? ref_fasta
#        File? ref_fasta_index
#        File? ref_dict
#        File vcf
#        File vcf_idx
#
#        Boolean compress_output = false
#
#        Array[String] filter_expressions
#        Array[String] filter_names
#        String? variant_filtration_extra_args
#
#        Runtime runtime_params
#    }
#
#    String output_vcf = basename(basename(vcf, ".gz"), ".vcf") + ".hard_filtered.vcf" + if compress_output then ".gz" else ""
#    String output_vcf_idx = output_vcf + if compress_output then ".tbi" else ".idx"
#
#    command <<<
#        set -e
#        export GATK_LOCAL_JAR=~{select_first([runtime_params.jar_override, "/root/gatk.jar"])}
#
#        echo ""
#        # Some variants don't have certain INFO fields, so we suppress the warning messages.
#        echo "Suppressing the following warning message: 'WARN  JexlEngine - '"
#        echo ""
#
#        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
#            VariantFiltration \
#            ~{"-R '" + ref_fasta + "'"} \
#            ~{"-L '" + interval_list + "'"} \
#            -V '~{vcf}' \
#            ~{if (length(filter_names) > 0) then " --filter-name '" else ""}~{default="" sep="' --filter-name '" filter_names}~{if (length(filter_names) > 0) then "'" else ""} \
#            ~{if (length(filter_expressions) > 0) then " --filter-expression '" else ""}~{default="" sep="' --filter-expression '" filter_expressions}~{if (length(filter_expressions) > 0) then "'" else ""} \
#            --output '~{output_vcf}' \
#            ~{variant_filtration_extra_args} \
#            2> >(grep -v "WARN  JexlEngine - " >&2)
#    >>>
#
#    output {
#        File filtered_vcf = output_vcf
#        File filtered_vcf_idx = output_vcf_idx
#    }
#
#    runtime {
#        docker: runtime_params.docker
#        bootDiskSizeGb: runtime_params.boot_disk_size
#        memory: runtime_params.machine_mem + " MB"
#        runtime_minutes: runtime_params.runtime_minutes
#        disks: "local-disk " + runtime_params.disk + " HDD"
#        preemptible: runtime_params.preemptible
#        maxRetries: runtime_params.max_retries
#        cpu: runtime_params.cpu
#    }
#
#    parameter_meta {
#        interval_list: {localization_optional: true}
#        ref_fasta: {localization_optional: true}
#        ref_fasta_index: {localization_optional: true}
#        ref_dict: {localization_optional: true}
#        vcf: {localization_optional: true}
#        vcf_idx: {localization_optional: true}
#    }
#}
#
#task LeftAlignAndTrimVariants {
#    input {
#        File? ref_fasta
#        File? ref_fasta_index
#        File? ref_dict
#        File vcf
#        File vcf_idx
#        Int max_indel_length = 200
#        Boolean dont_trim_alleles = false
#        Boolean split_multi_allelics = false
#
#        Boolean compress_output = false
#        String? left_align_and_trim_variants_extra_args
#
#        Runtime runtime_params
#    }
#
#    String output_vcf_ = basename(basename(vcf, ".gz"), ".vcf") + ".split.trimmed.vcf" + if compress_output then ".gz" else ""
#    String output_vcf_idx_ = output_vcf_ + if compress_output then ".tbi" else ".idx"
#
#    command <<<
#        set -e
#        export GATK_LOCAL_JAR=~{select_first([runtime_params.jar_override, "/root/gatk.jar"])}
#        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
#            LeftAlignAndTrimVariants \
#            -R '~{ref_fasta}' \
#            -V '~{vcf}' \
#            --output '~{output_vcf_}' \
#            --max-indel-length ~{max_indel_length} \
#            ~{if (dont_trim_alleles) then " --dont-trim-alleles " else ""} \
#            ~{if (split_multi_allelics) then " --split-multi-allelics " else ""} \
#            ~{left_align_and_trim_variants_extra_args}
#    >>>
#
#    output {
#        File output_vcf = output_vcf_
#        File output_vcf_idx = output_vcf_idx_
#    }
#
#    runtime {
#        docker: runtime_params.docker
#        bootDiskSizeGb: runtime_params.boot_disk_size
#        memory: runtime_params.machine_mem + " MB"
#        runtime_minutes: runtime_params.runtime_minutes
#        disks: "local-disk " + runtime_params.disk + " HDD"
#        preemptible: runtime_params.preemptible
#        maxRetries: runtime_params.max_retries
#        cpu: runtime_params.cpu
#    }
#
#    parameter_meta {
#        ref_fasta: {localization_optional: true}
#        ref_fasta_index: {localization_optional: true}
#        ref_dict: {localization_optional: true}
#        vcf: {localization_optional: true}
#        vcf_idx: {localization_optional: true}
#    }
#}

task SelectVariants {
    input {
        File? ref_fasta
        File? ref_fasta_index
        File? ref_dict
        File? interval_list
        File? interval_blacklist
        File? interval_blacklist_idx
        File vcf
        File vcf_idx
        Boolean compress_output = false
        Boolean select_somatic = false
        Boolean select_germline = false
        String somatic_filter_whitelist = "PASS,normal_artifact"
        String germline_filter_whitelist = "normal_artifact,panel_of_normals"
        String suffix = ""
        String? tumor_sample_name
        String? normal_sample_name
        String? select_variants_extra_args

        Runtime runtime_params
    }

    String uncompressed_input_vcf = basename(vcf, ".gz")
    Boolean is_compressed = uncompressed_input_vcf != basename(vcf)
    String base_name = if defined(tumor_sample_name) then sub(select_first([tumor_sample_name, ""]), " ", "+") else basename(uncompressed_input_vcf, ".vcf")
    String output_base_name = base_name + ".selected" + suffix
    
    String select_variants_output_vcf = output_base_name + ".tmp.vcf"
    String select_variants_output_vcf_idx = select_variants_output_vcf + ".idx"
    String uncompressed_selected_vcf = output_base_name + ".vcf"
    String uncompressed_selected_vcf_idx = uncompressed_selected_vcf + ".idx"
    String output_vcf = uncompressed_selected_vcf + if compress_output then ".gz" else ""
    String output_vcf_idx = output_vcf + if compress_output then ".tbi" else ".idx"

    String output_not_selected_vcf = output_base_name + ".not_selected.vcf" + if compress_output then ".gz" else ""
    String output_not_selected_vcf_idx = output_not_selected_vcf + if compress_output then ".tbi" else ".idx"

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{select_first([runtime_params.jar_override, "/root/gatk.jar"])}
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            SelectVariants \
            ~{"-R '" + ref_fasta + "'"} \
            ~{"-L '" + interval_list + "'"} \
            ~{"-XL '" + interval_blacklist + "'"} \
            -V '~{vcf}' \
            --output '~{select_variants_output_vcf}' \
            --exclude-filtered false \
            ~{"--sample-name '" + tumor_sample_name + "'"} \
            ~{"--sample-name '" + normal_sample_name + "'"} \
            ~{select_variants_extra_args}

        set -uo pipefail
        set +e  # grep returns 1 if no lines are found
        num_vars=$(grep -v "^#" '~{select_variants_output_vcf}' | wc -l)
        echo ">> Selected $num_vars variants."

        # =======================================
        # We do the selection step using grep|awk to also select germline variants.
        # ASSUMPTION: multi-allelic variants are split into one variant per row.
        # Otherwise passing variants that are accompanied by an artifactual other-allelic
        # variant will not be selected.

        if [ "$num_vars" -eq 0 ] || [ "~{select_somatic}" == "false" ] && [ "~{select_germline}" == "false" ] ; then
            cp '~{select_variants_output_vcf}' '~{uncompressed_selected_vcf}'
        else
            grep "^#" '~{select_variants_output_vcf}' > '~{uncompressed_selected_vcf}'
            if [ "~{select_somatic}" == "true" ] ; then
                echo ">> Selecting PASSing/whitelisted variants ... "
                # FilterMutectCalls assumes a normal sample with no tumor cell
                # contamination. If there is contamination from tumor cells,
                # somatic variants will be annotated as "normal_artifact", thus
                # it is desirable to whitelist them.
                mkdir -p tmp
                grep -v "^#" '~{select_variants_output_vcf}' \
                    | awk -F'\t' -v whitelist="~{somatic_filter_whitelist}" '
                        BEGIN {
                            n = split(whitelist, allowed_list, ",")
                            for (i = 1; i <= n; i++) {
                                allowed_tags[allowed_list[i]]
                            }
                        }
                        {
                            split($7, tags, ";")
                            is_somatic = 1
                            for (i in tags) {
                                if (tags[i] in allowed_tags) {
                                    continue
                                } else {
                                    is_somatic = 0
                                    break
                                }
                            }
                            if (is_somatic) {
                                print $0
                            }
                        }' \
                    > tmp/somatic.vcf
                cat tmp/somatic.vcf >> '~{uncompressed_selected_vcf}'
                num_selected_vars=$(cat tmp/somatic.vcf | wc -l)
                rm -rf tmp
                echo ">> Selected $num_selected_vars PASSing out of $num_vars variants."
            fi
            if [ "~{select_germline}" == "true" ] ; then
                echo ">> Selecting germline variants ... "
                # FilterMutectCalls does not distinguish well between germline
                # and artifacts; it's calibrated towards somatic vs non-somatic.
                # Thus, many good germline calls may also have an artifact flag.
                # Additional filtering of those may be necessary.
                mkdir -p tmp
                grep -v "^#" '~{select_variants_output_vcf}' \
                    | awk -F'\t' -v whitelist="~{germline_filter_whitelist}" '
                        BEGIN {
                            n = split(whitelist, allowed_list, ",")
                            for (i = 1; i <= n; i++) {
                                allowed_tags[allowed_list[i]]
                            }
                        }
                        {
                            split($7, tags, ";")
                            contains_germline = 0
                            valid = 1
                            for (i in tags) {
                                if (tags[i] == "germline") {
                                    contains_germline = 1
                                } else if (tags[i] in allowed_tags) {
                                    continue
                                } else {
                                    valid = 0
                                    break
                                }
                            }
                            if (contains_germline && valid) {
                                print $0
                            }
                        }' \
                    > tmp/germline.vcf
                cat tmp/germline.vcf >> '~{uncompressed_selected_vcf}'
                num_selected_vars=$(cat tmp/germline.vcf | wc -l)
                rm -rf tmp
                echo ">> Selected $num_selected_vars germline out of $num_vars variants."
            fi
        fi

        set -e

        rm -f '~{select_variants_output_vcf}' '~{select_variants_output_vcf_idx}'

        # =======================================
        # Hack to correct a SelectVariants output bug. When selecting for samples, this
        # task only retains the first sample annotation in the header. Those annotations
        # are important for Funcotator to fill the t_alt_count and t_ref_count coverage
        # columns. This hack assumes that only one tumor sample and/or only one normal
        # sample have been selected.

        if [ "~{defined(tumor_sample_name)}" == "true" ] ; then
            echo ">> Fixing tumor sample name in vcf header ... "
            input_header=$(grep "##tumor_sample=" '~{uncompressed_selected_vcf}')
            corrected_header="##tumor_sample=~{tumor_sample_name}"
            sed -i "s/$input_header/$corrected_header/g" '~{uncompressed_selected_vcf}'
        fi
        if [ "~{defined(normal_sample_name)}" == "true" ] ; then
            echo ">> Fixing normal sample name in vcf header ... "
            input_header=$(grep "##normal_sample=" '~{uncompressed_selected_vcf}')
            corrected_header="##normal_sample=~{normal_sample_name}"
            sed -i "s/$input_header/$corrected_header/g" '~{uncompressed_selected_vcf}'
        fi

        # Selecting both PASSing and germline variants can lead to unsorted vcf.
        mv '~{uncompressed_selected_vcf}' 'unsorted.~{uncompressed_selected_vcf}'
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            SortVcf \
            -I 'unsorted.~{uncompressed_selected_vcf}' \
            -O '~{uncompressed_selected_vcf}' \
            ~{"-SD '" +  ref_dict + "'"}
        rm -f 'unsorted.~{uncompressed_selected_vcf}'

        set +e  # grep returns 1 if no lines are found
        grep -v "^#" '~{uncompressed_selected_vcf}' | wc -l > num_selected_vars.txt
        set -e

        if [ "~{compress_output}" == "true" ] ; then
            echo ">> Compressing selected vcf."
            bgzip -c '~{uncompressed_selected_vcf}' > '~{output_vcf}'
            gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
                IndexFeatureFile \
                --input '~{output_vcf}' \
                --output '~{output_vcf_idx}'
            rm -f '~{uncompressed_selected_vcf}' '~{uncompressed_selected_vcf_idx}'
        fi

        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            SelectVariants \
            ~{"-R '" + ref_fasta + "'"} \
            -V '~{vcf}' \
            --discordance '~{output_vcf}' \
            --output '~{output_not_selected_vcf}'

        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            CountVariants \
            -V '~{output_not_selected_vcf}' \
            -O 'num_not_selected_vars.txt'
    >>>

    output {
        File selected_vcf = output_vcf
        File selected_vcf_idx = output_vcf_idx
        Int num_selected_variants = read_int("num_selected_vars.txt")
        File not_selected_vcf = output_not_selected_vcf
        File not_selected_vcf_idx = output_not_selected_vcf_idx
        Int num_not_selected_variants = read_int("num_not_selected_vars.txt")
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
        interval_blacklist: {localization_optional: true}
        interval_blacklist_idx: {localization_optional: true}
        ref_fasta: {localization_optional: true}
        ref_fasta_index: {localization_optional: true}
        # ref_dict: {localization_optional: true}  # needs to be localized for SortVcf
        vcf: {localization_optional: true}
        vcf_idx: {localization_optional: true}
    }
}

task GatherVCFs {
	input {
        File? ref_fasta
        File? ref_fasta_index
        File? ref_dict
        Array[File] vcfs
        Array[File] vcfs_idx
        String output_name
        Boolean compress_output = false
        Boolean drop_duplicate_sites = false

        Runtime runtime_params
    }

    Int diskGB = runtime_params.disk + ceil(2 * size(vcfs, "GB"))

    String output_vcf = output_name + ".vcf" + if compress_output then ".gz" else ""
    String output_vcf_idx = output_vcf + if compress_output then ".tbi" else ".idx"

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{select_first([runtime_params.jar_override, "/root/gatk.jar"])}
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            MergeVcfs \
            ~{sep="' " prefix("-I '", vcfs)}' \
            ~{"-R '" + ref_fasta + "'"} \
            ~{"-D '" + ref_dict + "'"} \
            -O 'tmp.~{output_vcf}'

        if [ "~{drop_duplicate_sites}" == "true" ]; then
            bcftools norm \
                -d exact \
                -o '~{output_vcf}' \
                'tmp.~{output_vcf}'

            bcftools index -t \
                -o '~{output_vcf_idx}' \
                '~{output_vcf}'
        else
            mv 'tmp.~{output_vcf}' '~{output_vcf}'
            mv 'tmp.~{output_vcf_idx}' '~{output_vcf_idx}'
        fi
    >>>

    output {
    	File merged_vcf = output_vcf
        File merged_vcf_idx = output_vcf_idx
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

    # Optional localization leads to cromwell error.
    # parameter_meta {
    #     vcfs: {localization_optional: true}
    #     vcfs_idx: {localization_optional: true}
    # }
}

task MergeMAFs {
    # This tasks weakly assumes that all mafs have the same header
    # and stronly assumes the same column order.

	input {
        Array[File] mafs  # assumes uncompressed
        String output_name
        Boolean compress_output = false

        Runtime runtime_params
    }

    String uncompressed_output_maf = output_name + ".maf"
    String output_maf = output_name + ".maf" + if compress_output then ".gz" else ""
    String dollar = "$"

    command <<<
        set -euxo pipefail

        # Convert WDL array to a temporary file
        printf "~{sep='\n' mafs}" > temp_mafs.txt

        # Read temporary file into a shell array
        mapfile -t maf_files < temp_mafs.txt

        # Extract leading comment lines from first file
        grep "^#" "~{dollar}{maf_files[0]}" > '~{uncompressed_output_maf}'

        # Extract column headers from first file
        # (|| true is necessary since either grep or head return non-zero exit code; don't understand why.)
        grep -v "^#" "~{dollar}{maf_files[0]}" | head -n 1 >> '~{uncompressed_output_maf}' || true

        # Extract variants
        for maf in "~{dollar}{maf_files[@]}" ; do
            grep -v "^#" "$maf" | tail -n +2 >> '~{uncompressed_output_maf}' || true
        done

        if [ "~{compress_output}" == "true" ] ; then
            echo ">> Compressing merged MAF."
            gzip -c '~{uncompressed_output_maf}' > '~{output_maf}'
            rm -f '~{uncompressed_output_maf}'
        fi
        # else: uncompressed_output_maf == output_maf by design

        rm -f temp_mafs.txt
    >>>

    output {
    	File merged_maf = output_maf
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

task PrintReads {
    input {
        File? ref_fasta
        File? ref_fasta_index
        File? ref_dict

        String prefix
        Array[File] bams
        Array[File] bais
        File? interval_list
        File? vcf
        File? vcf_idx

        Runtime runtime_params
    }

    String output_file = prefix + ".bam"
    String output_index = prefix + ".bai"

    # Even if we subset the bams to the interval_list, this task is very short,
    # so we won't spend a lot of money on it.
    Int diskGB = runtime_params.disk + ceil(size(bams, "GB"))

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{select_first([runtime_params.jar_override, "/root/gatk.jar"])}
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            PrintReads \
            ~{sep="' " prefix("-I '", bams)}' \
            ~{sep="' " prefix("--read-index '", bais)}' \
            -O '~{output_file}' \
            ~{"-R '" + ref_fasta + "'"} \
            ~{"-L '" + interval_list + "'"} \
            ~{"-L '" + vcf + "'"}
    >>>

    output {
        File output_bam = output_file
        File output_bai = output_index
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
        interval_list: {localization_optional: true}
        bams: {localization_optional: true}
        bais: {localization_optional: true}
        vcf: {localization_optional: true}
        vcf_idx: {localization_optional: true}
    }
}

# Ensure ordering of bam contigs are in the same order as ref_dict
task ReorderSam {
    input {
        File? ref_fasta
        File? ref_fasta_index
        File ref_dict

        File bam
        File bai
        Runtime runtime_params
    }
    
    String prefix = basename(bam, ".bam")
    String output_bam = prefix + ".reordered.bam"
    String output_bai = prefix + ".reordered.bai"
    Int diskGB = runtime_params.disk + ceil(size(bam, "GB") * 2) + ceil(size(ref_dict, "GB"))

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{select_first([runtime_params.jar_override, "/root/gatk.jar"])}
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            ReorderSam \
            -I '~{bam}' \
            -O '~{output_bam}' \
            -SD '~{ref_dict}' \
            --CREATE_INDEX true
    >>>

    output {
        File reordered_bam = output_bam
        File reordered_bai = output_bai
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
        # Picard ReorderSam requires ref_dict to be localized
        # ref_dict: {localization_optional: true} 
        # bam: {localization_optional: true}
        # bai: {localization_optional: true}
    }
}