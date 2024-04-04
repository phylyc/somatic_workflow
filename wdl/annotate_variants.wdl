version development

import "workflow_arguments.wdl" as wfargs
import "runtime_collection.wdl" as rtc
import "runtimes.wdl" as rt
import "tasks.wdl"


workflow AnnotateVariants {
    input {
        File vcf
        File vcf_idx

        String? individual_id
        String tumor_bam_name
        String tumor_sample_name
        String? normal_bam_name
        String? normal_sample_name

        WorkflowArguments args
        RuntimeCollection runtime_collection
    }

    Array[File] scattered_interval_list = if args.run_variant_annotation_scattered then args.scattered_interval_list else [args.preprocessed_interval_list]

    scatter (intervals in scattered_interval_list) {
        call tasks.SelectVariants as SelectSampleVariants {
            input:
                ref_fasta = args.files.ref_fasta,
                ref_fasta_index = args.files.ref_fasta_index,
                ref_dict = args.files.ref_dict,
                interval_list = intervals,
                vcf = vcf,
                vcf_idx = vcf_idx,
                tumor_sample_name = tumor_bam_name,
                normal_sample_name = normal_bam_name,
                compress_output = args.compress_output,
                select_variants_extra_args = args.select_variants_extra_args,
                runtime_params = runtime_collection.select_variants
        }

        call Funcotate {
            input:
                ref_fasta = args.files.ref_fasta,
                ref_fasta_index = args.files.ref_fasta_index,
                ref_dict = args.files.ref_dict,
                interval_list = intervals,
                vcf = SelectSampleVariants.selected_vcf,
                vcf_idx = SelectSampleVariants.selected_vcf_idx,
                individual_id = individual_id,
                tumor_sample_name = tumor_sample_name,
                normal_sample_name = normal_sample_name,
                output_base_name = tumor_sample_name + ".annotated",
                reference_version = args.funcotator_reference_version,
                output_format = args.funcotator_output_format,
                variant_type = args.funcotator_variant_type,
                transcript_selection_mode = args.funcotator_transcript_selection_mode,
                transcript_list = args.files.funcotator_transcript_list,
                data_sources_tar_gz = args.files.funcotator_data_sources_tar_gz,
                use_gnomad = args.funcotator_use_gnomad,
                compress_output = args.compress_output,
                data_sources_paths = args.funcotator_data_sources_paths,
                annotation_defaults = args.funcotator_annotation_defaults,
                annotation_overrides = args.funcotator_annotation_overrides,
                exclude_fields = args.funcotator_exclude_fields,
                funcotate_extra_args = args.funcotate_extra_args,
                runtime_params = runtime_collection.funcotate
        }
    }

    if (args.funcotator_output_format == "VCF") {
        call tasks.MergeVCFs {
            input:
                vcfs = Funcotate.annotations,
                vcfs_idx = select_all(Funcotate.annotations_idx),
                output_name = tumor_sample_name + ".annotated",
                compress_output = args.compress_output,
                runtime_params = runtime_collection.merge_vcfs
        }
    }
    if (args.funcotator_output_format == "MAF") {
        call tasks.MergeMAFs {
            input:
                mafs = Funcotate.annotations,
                output_name = tumor_sample_name + ".annotated",
                compress_output = args.compress_output,
                runtime_params = runtime_collection.merge_mafs
        }
    }

    output {
        File annotated_variants = select_first(select_all([MergeVCFs.merged_vcf, MergeMAFs.merged_maf]))
        File? annotated_variants_idx = MergeVCFs.merged_vcf_idx
    }
}

task Funcotate {
    input {
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        File? interval_list

        File vcf
        File vcf_idx
        String? individual_id = "Unknown"
        String? tumor_sample_name = "Unknown"
        String? normal_sample_name = "Unknown"

        String reference_version = "hg19"
        String output_base_name
        String output_format = "MAF"
        String variant_type = "somatic"  # alternative: germline
        String transcript_selection_mode = "CANONICAL"  # GATK default: "CANONICAL"
        File? transcript_list
        File? data_sources_tar_gz  # most recent version is downloaded if not chosen
        Boolean use_gnomad = false
        Boolean compress_output = false  # ignored if output_format == "MAF"
        Array[String]? data_sources_paths
        Array[String]? annotation_defaults
        Array[String]? annotation_overrides
        Array[String]? exclude_fields
        String? funcotate_extra_args

        Runtime runtime_params
    }

    # TODO: Annotate phase

    # ==============
    # Process input args:
    String output_maf = output_base_name + ".maf"
    String output_maf_index = output_maf + ".idx"
    String output_vcf = output_base_name + if compress_output then ".vcf.gz" else ".vcf"
    String output_vcf_idx = output_vcf +  if compress_output then ".tbi" else ".idx"
    String output_file = if output_format == "MAF" then output_maf else output_vcf
    String output_file_index = if output_format == "MAF" then output_maf_index else output_vcf_idx

    # Calculate disk size:
    Int funco_tar_sizeGB = if defined(data_sources_paths) then 0 else (if defined(data_sources_tar_gz) then 4 * ceil(size(data_sources_tar_gz, "GB")) else 100)
    Int disk = runtime_params.disk + funco_tar_sizeGB

    String annotation_default_arg = if defined(annotation_defaults) then sep(" ", prefix("--annotation-default ", select_first([annotation_defaults, ["none"]]))) else ""
    String annotation_override_arg = if defined(annotation_overrides) then sep(" ", prefix("--annotation-override ", select_first([annotation_overrides, ["none"]]))) else ""
    String exclude_fields_arg = if defined(exclude_fields) then sep(" ", prefix("--exclude-field ", select_first([exclude_fields, ["none"]]))) else ""

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{select_first([runtime_params.jar_override, "/root/gatk.jar"])}

        if [ "~{output_format}" != "MAF" ] && [ "~{output_format}" != "VCF" ] ; then
            echo "ERROR: Output format must be MAF or VCF."
            false
        fi

        if [ "~{defined(data_sources_paths)}" == "false" ] ; then
            mkdir datasources_dir
            DATA_SOURCES_FOLDER="$PWD/datasources_dir"
            echo "Obtaining Funcotator data sources..."
            if [[ ! -z "~{data_sources_tar_gz}" ]]; then
                data_sources_tar_gz=~{data_sources_tar_gz}
            else
                data_sources_tar_gz="funcotator_datasources.tar.gz"
                gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
                    FuncotatorDataSourceDownloader \
                    --~{variant_type} \
                    --validate-integrity \
                    --output $data_sources_tar_gz
            fi
            echo "Unzipping Funcotator data sources..."
            tar zxvf $data_sources_tar_gz -C datasources_dir --strip-components 1

            if [ "~{use_gnomad}" == "true" ] ; then
                echo "Enabling gnomAD..."
                for potential_gnomad_gz in gnomAD_exome.tar.gz gnomAD_genome.tar.gz ; do
                    if [[ -f $DATA_SOURCES_FOLDER/$potential_gnomad_gz ]] ; then
                        cd $DATA_SOURCES_FOLDER
                        tar -zvxf $potential_gnomad_gz
                        cd -
                    else
                        echo "ERROR: Cannot find gnomAD folder: $potential_gnomad_gz" 1>&2
                        false
                    fi
                done
            fi
        fi

        echo ""
        # using a custom list of transcripts gives this INFO message:
        printf "Suppressing the following INFO messages: 'Adding transcript ID to transcript set:'\n" >&2
        echo ""

        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            Funcotator \
            ~{sep=" " prefix("--data-sources-path ", select_first([data_sources_paths, ["$DATA_SOURCES_FOLDER"]]))} \
            --ref-version ~{reference_version} \
            --output-file-format ~{output_format} \
            -R '~{ref_fasta}' \
            -V '~{vcf}' \
            -O '~{output_file}' \
            ~{"-L '" + interval_list + "'"} \
            ~{"--transcript-selection-mode " + transcript_selection_mode} \
            ~{"--transcript-list '" + transcript_list + "'"} \
            --annotation-default 'individual_id:~{individual_id}' \
            --annotation-default 'tumor_barcode:~{tumor_sample_name}' \
            --annotation-default 'normal_barcode:~{normal_sample_name}' \
            --seconds-between-progress-updates 60 \
            ~{annotation_default_arg} \
            ~{annotation_override_arg} \
            ~{exclude_fields_arg} \
            ~{funcotate_extra_args} \
            2> >( \
                grep -v 'Adding transcript ID to transcript set' \
                >&2 \
            )

        rm -rf $DATA_SOURCES_FOLDER
    >>>

    output {
        File annotations = output_file
        File? annotations_idx = output_file_index
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

    parameter_meta {
        ref_fasta: {localization_optional: true}
        ref_fasta_index: {localization_optional: true}
        ref_dict: {localization_optional: true}
        interval_list: {localization_optional: true}
        vcf: {localization_optional: true}
        vcf_idx: {localization_optional: true}
        # transcript_selection_file: {localization_optional: true}
    }
}