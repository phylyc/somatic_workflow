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

task IndexFeatureFile {
    input {
        File vcf

        Runtime runtime_params
    }

    String uncompressed_vcf = basename(vcf, ".gz")
    Boolean is_compressed = uncompressed_vcf != basename(vcf)
    String output_vcf_idx = basename(vcf) + if is_compressed then ".tbi" else ".idx"
    Int diskGB = runtime_params.disk + ceil(size(vcf, "GB")) + 1

    command <<<
        set -euxo pipefail
        export GATK_LOCAL_JAR=~{select_first([runtime_params.jar_override, "/root/gatk.jar"])}
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            IndexFeatureFile \
            --input '~{vcf}' \
            --output '~{output_vcf_idx}'
    >>>

    output {
        File vcf_idx = output_vcf_idx
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

task ParseInput {
    # Pre-flight validation gate. Fails the workflow before any expensive task if
    # the inputs are mis-shaped or violate an implicit downstream requirement. The
    # heavy reference resources are tested for presence via the has_* booleans, so
    # this task localizes nothing in the default (non-deep) mode; only the small
    # per-run interval / panel files are localized for content checks. The deep
    # mode additionally localizes the bams for samtools integrity + contig-order
    # checks (see run_input_validation_deep).
    input {
        String script = "https://github.com/phylyc/somatic_workflow/raw/master/python/validate_inputs.py"

        Int n_bams
        Int n_bais
        Array[String]? sample_names
        Array[String]? normal_sample_names
        Array[Boolean]? is_paired_end
        Array[Boolean]? use_for_tCR
        Array[Boolean]? use_for_aCR
        Array[Int]? timepoints
        Array[File]? target_intervals
        Array[File]? annotated_target_intervals
        Array[File]? cnv_panel_of_normals

        Boolean has_common_germline_alleles
        Boolean has_common_germline_alleles_idx
        Boolean has_realignment_image
        Boolean has_germline_resource
        Boolean has_snv_panel_of_normals
        Boolean has_germline_resource_v4_1
        Boolean has_snv_panel_of_normals_v4_1
        Boolean has_funcotator_sources
        Boolean has_sex

        Boolean run_collect_total_read_counts
        Boolean run_model_segments
        Boolean run_collect_allelic_read_counts
        Boolean run_contamination_model
        Boolean run_realignment_filter
        Boolean run_variant_calling_mutect1
        Boolean run_variant_calling
        Boolean run_variant_annotation
        Boolean run_clonal_decomposition

        Boolean deep = false
        Array[File] bams = []
        File? ref_dict

        Runtime runtime_params
    }

    # In deep mode the bams are localized for samtools, so the request must cover
    # them; bams is empty in the default mode, leaving just the base disk.
    Int disk = runtime_params.disk + if deep then ceil(size(bams, "GB")) + 1 else 0

    # Terra may pass an empty optional array as defined([]). Normalize optional
    # arrays and gate command-line expansion on non-empty values so repeated
    # argparse flags are omitted instead of receiving a single empty string.
    Array[String] sample_names_args = select_first([sample_names, []])
    Array[String] normal_sample_names_args = select_first([normal_sample_names, []])
    Array[Boolean] is_paired_end_args = select_first([is_paired_end, []])
    Array[Boolean] use_for_tCR_args = select_first([use_for_tCR, []])
    Array[Boolean] use_for_aCR_args = select_first([use_for_aCR, []])
    Array[Int] timepoints_args = select_first([timepoints, []])
    Array[File] target_intervals_args = select_first([target_intervals, []])
    Array[File] annotated_target_intervals_args = select_first([annotated_target_intervals, []])
    Array[File] cnv_panel_of_normals_args = select_first([cnv_panel_of_normals, []])

    Boolean has_sample_names_args = length(sample_names_args) > 0
    Boolean has_normal_sample_names_args = length(normal_sample_names_args) > 0
    Boolean has_is_paired_end_args = length(is_paired_end_args) > 0
    Boolean has_use_for_tCR_args = length(use_for_tCR_args) > 0
    Boolean has_use_for_aCR_args = length(use_for_aCR_args) > 0
    Boolean has_timepoints_args = length(timepoints_args) > 0
    Boolean has_target_intervals_args = length(target_intervals_args) > 0
    Boolean has_annotated_target_intervals_args = length(annotated_target_intervals_args) > 0
    Boolean has_cnv_panel_of_normals_args = length(cnv_panel_of_normals_args) > 0
    Boolean has_bams_args = length(bams) > 0

    command <<<
        set -euxo pipefail
        wget -O validate_inputs.py ~{script}
        python validate_inputs.py \
            --n_bams ~{n_bams} \
            --n_bais ~{n_bais} \
            ~{true="--sample_names '" false="" has_sample_names_args}~{default="" sep="' --sample_names '" sample_names_args}~{true="'" false="" has_sample_names_args} \
            ~{true="--normal_sample_names '" false="" has_normal_sample_names_args}~{default="" sep="' --normal_sample_names '" normal_sample_names_args}~{true="'" false="" has_normal_sample_names_args} \
            ~{true="--is_paired_end '" false="" has_is_paired_end_args}~{default="" sep="' --is_paired_end '" is_paired_end_args}~{true="'" false="" has_is_paired_end_args} \
            ~{true="--use_for_tCR '" false="" has_use_for_tCR_args}~{default="" sep="' --use_for_tCR '" use_for_tCR_args}~{true="'" false="" has_use_for_tCR_args} \
            ~{true="--use_for_aCR '" false="" has_use_for_aCR_args}~{default="" sep="' --use_for_aCR '" use_for_aCR_args}~{true="'" false="" has_use_for_aCR_args} \
            ~{true="--timepoints '" false="" has_timepoints_args}~{default="" sep="' --timepoints '" timepoints_args}~{true="'" false="" has_timepoints_args} \
            ~{true="--target_intervals '" false="" has_target_intervals_args}~{default="" sep="' --target_intervals '" target_intervals_args}~{true="'" false="" has_target_intervals_args} \
            ~{true="--annotated_target_intervals '" false="" has_annotated_target_intervals_args}~{default="" sep="' --annotated_target_intervals '" annotated_target_intervals_args}~{true="'" false="" has_annotated_target_intervals_args} \
            ~{true="--cnv_panel_of_normals '" false="" has_cnv_panel_of_normals_args}~{default="" sep="' --cnv_panel_of_normals '" cnv_panel_of_normals_args}~{true="'" false="" has_cnv_panel_of_normals_args} \
            --has_common_germline_alleles ~{has_common_germline_alleles} \
            --has_common_germline_alleles_idx ~{has_common_germline_alleles_idx} \
            --has_realignment_image ~{has_realignment_image} \
            --has_germline_resource ~{has_germline_resource} \
            --has_snv_panel_of_normals ~{has_snv_panel_of_normals} \
            --has_germline_resource_v4_1 ~{has_germline_resource_v4_1} \
            --has_snv_panel_of_normals_v4_1 ~{has_snv_panel_of_normals_v4_1} \
            --has_funcotator_sources ~{has_funcotator_sources} \
            --has_sex ~{has_sex} \
            --run_collect_total_read_counts ~{run_collect_total_read_counts} \
            --run_model_segments ~{run_model_segments} \
            --run_collect_allelic_read_counts ~{run_collect_allelic_read_counts} \
            --run_contamination_model ~{run_contamination_model} \
            --run_realignment_filter ~{run_realignment_filter} \
            --run_variant_calling_mutect1 ~{run_variant_calling_mutect1} \
            --run_variant_calling ~{run_variant_calling} \
            --run_variant_annotation ~{run_variant_annotation} \
            --run_clonal_decomposition ~{run_clonal_decomposition} \
            ~{if deep then "--deep" else ""} \
            ~{true="--bams '" false="" deep && has_bams_args}~{default="" sep="' --bams '" bams}~{true="'" false="" deep && has_bams_args} \
            ~{if (deep && defined(ref_dict)) then "--ref_dict '" + ref_dict + "'" else ""}
    >>>

    output {
        Boolean validated = true
    }

    runtime {
        docker: runtime_params.docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: runtime_params.machine_mem + " MB"
        runtime_minutes: runtime_params.runtime_minutes
        disks: "local-disk " + disk + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
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
    Array[File] interval_lists_args = select_first([interval_lists, []])
    Boolean has_interval_lists_args = length(interval_lists_args) > 0

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{select_first([runtime_params.jar_override, "/root/gatk.jar"])}
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            PreprocessIntervals \
            -R '~{ref_fasta}' \
            ~{"-L '" + interval_list + "'"} \
            ~{"-XL '" + interval_blacklist + "'"} \
            ~{true="-L '" false="" has_interval_lists_args}~{default="" sep="' -L '" interval_lists_args}~{true="'" false="" has_interval_lists_args} \
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
        String somatic_filter_whitelist = "PASS,normal_artifact,RESCUED"
        String germline_filter_whitelist = "normal_artifact,panel_of_normals,lowGERMQ"
        String suffix = ""
        String? tumor_sample_name
        String? normal_sample_name
        String? select_variants_extra_args

        Runtime runtime_params
    }

    Int diskGB = runtime_params.disk + ceil(10 * size(vcf, "GB"))

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
        disks: "local-disk " + diskGB + " HDD"
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
    Int diskGB = runtime_params.disk + ceil(size(bam, "GB") * 5) + ceil(size(ref_dict, "GB"))

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

task RunPeddy {
    input {
        String patient_id
        File gvcf
        File gvcf_idx
        String genome_build

        Runtime runtime_params
    }

    command <<<
        set -e

        # Filter the VCF to only retain autosomes
        if [[ "~{genome_build}" == "hg38" ]]; then
            regions=$(echo chr{1..22} | tr ' ' ',')
        else
            regions=$(echo {1..22} | tr ' ' ',')
        fi
        bcftools view -r "$regions" "~{gvcf}" -Oz -o "~{patient_id}.filtered.vcf.gz"
        bcftools index --tbi "~{patient_id}.filtered.vcf.gz"

        # Create dummy ped file
        printf "~{patient_id}\t~{patient_id}\t0\t0\t0\t-9\n" > "~{patient_id}.ped"

        # Run peddy
        if [ "~{genome_build}" == "hg38" ]; then
            peddy --plot -p 4 --sites ~{genome_build} --prefix ~{patient_id} "~{patient_id}.filtered.vcf.gz" "~{patient_id}.ped"
        else
            peddy --plot -p 4 --prefix ~{patient_id} "~{patient_id}.filtered.vcf.gz" "~{patient_id}.ped"
        fi

        # Getting the prediction and associated probability
        awk -F',' 'NR==2 {print $12}' "~{patient_id}.het_check.csv" > pred.txt
        awk -F',' 'NR==2 {print $13}' "~{patient_id}.het_check.csv" > prob.txt
    >>>

    output {
        String ancestry_pred = read_string("pred.txt")
        Float ancestry_prob = read_float("prob.txt")
        File ancestry_background_pca_table = "~{patient_id}.background_pca.json"
        File ancestry_pca_plot = "~{patient_id}.pca_check.png"
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
