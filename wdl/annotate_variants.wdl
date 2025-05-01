version development

import "workflow_arguments.wdl" as wfargs
import "runtime_collection.wdl" as rtc
import "runtimes.wdl" as rt
import "tasks.wdl"
import "workflow_arguments.wdl" as wfargs
import "workflow_resources.wdl" as wfres


workflow AnnotateVariants {
    input {
        File? ref_fasta
        File? ref_fasta_index
        File? ref_dict

        File vcf
        File vcf_idx
        Int? num_variants = 40

        String? individual_id
        String tumor_bam_name
        String tumor_sample_name
        String? normal_bam_name
        String? normal_sample_name

        WorkflowArguments? args
        RuntimeCollection runtime_collection = RuntimeParameters.rtc
    }

    # Define for standalone workflow
    call rtc.DefineRuntimeCollection as RuntimeParameters

    if (!defined(args)) {
        call wfres.DefineWorkflowResources as Files {
            input:
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                ref_dict = ref_dict,
        }

        call wfargs.DefineWorkflowArguments as Parameters {
            input:
                resources = Files.resources,
                runtime_collection = runtime_collection,
        }
    }
    WorkflowArguments this_args = select_first([args, Parameters.arguments])

    if (this_args.run_variant_annotation_scattered) {
        # Due to its long runtime, we scatter the annotation task over intervals.
        Int scatter_count = ceil((num_variants + 1) / this_args.variants_per_scatter)
        call tasks.SplitIntervals {
            input:
                interval_list = this_args.files.preprocessed_intervals,
                ref_fasta = this_args.files.ref_fasta,
                ref_fasta_index = this_args.files.ref_fasta_index,
                ref_dict = this_args.files.ref_dict,
                scatter_count = scatter_count,
                split_intervals_extra_args = this_args.split_intervals_extra_args,
                runtime_params = runtime_collection.split_intervals,
        }
    }

    scatter (intervals in select_first([SplitIntervals.interval_files, [this_args.files.preprocessed_intervals]])) {
        call tasks.SelectVariants as SelectSampleVariants {
            input:
                ref_fasta = this_args.files.ref_fasta,
                ref_fasta_index = this_args.files.ref_fasta_index,
                ref_dict = this_args.files.ref_dict,
                interval_list = intervals,
                vcf = vcf,
                vcf_idx = vcf_idx,
                tumor_sample_name = tumor_bam_name,
                normal_sample_name = normal_bam_name,
                compress_output = this_args.compress_output,
                select_variants_extra_args = this_args.select_variants_extra_args,
                runtime_params = runtime_collection.select_variants
        }

        if (SelectSampleVariants.num_selected_variants > 0) {
            call Funcotate {
                input:
                    ref_fasta = this_args.files.ref_fasta,
                    ref_fasta_index = this_args.files.ref_fasta_index,
                    ref_dict = this_args.files.ref_dict,
                    interval_list = intervals,
                    vcf = SelectSampleVariants.selected_vcf,
                    vcf_idx = SelectSampleVariants.selected_vcf_idx,
                    individual_id = individual_id,
                    tumor_sample_name = tumor_sample_name,
                    normal_sample_name = normal_sample_name,
                    output_base_name = tumor_sample_name + ".annotated",
                    reference_version = this_args.funcotator_reference_version,
                    output_format = this_args.funcotator_output_format,
                    variant_type = this_args.funcotator_variant_type,
                    transcript_selection_mode = this_args.funcotator_transcript_selection_mode,
                    transcript_list = this_args.files.funcotator_transcript_list,
                    data_sources_tar_gz = this_args.files.funcotator_data_sources_tar_gz,
                    use_gnomad = this_args.funcotator_use_gnomad,
                    compress_output = this_args.compress_output,
                    data_sources_paths = this_args.funcotator_data_sources_paths,
                    annotation_defaults = this_args.funcotator_annotation_defaults,
                    annotation_overrides = this_args.funcotator_annotation_overrides,
                    exclude_fields = this_args.funcotator_exclude_fields,
                    funcotate_extra_args = this_args.funcotate_extra_args,
                    runtime_params = runtime_collection.funcotate
            }
        }
        if (SelectSampleVariants.num_selected_variants <= 0) {
            call CreateEmptyAnnotation {
                input:
                    selected_vcf = SelectSampleVariants.selected_vcf,
                    output_base_name = tumor_sample_name + ".annotated",
                    output_format = this_args.funcotator_output_format,
                    compress_output = this_args.compress_output,
                    runtime_params = runtime_collection.create_empty_annotation
            }
        }
        # TODO: Add CADD annotation

        File this_annotation = select_first([Funcotate.annotations, CreateEmptyAnnotation.empty_annotation])
        File this_annotation_idx = select_first([Funcotate.annotations_idx, CreateEmptyAnnotation.empty_annotation_idx])
    }

    if (this_args.funcotator_output_format == "VCF") {
        call tasks.GatherVCFs {
            input:
                vcfs = this_annotation,
                vcfs_idx = this_annotation_idx,
                output_name = tumor_sample_name + ".annotated",
                compress_output = this_args.compress_output,
                runtime_params = runtime_collection.gather_vcfs
        }
    }
    if (this_args.funcotator_output_format == "MAF") {
        call tasks.MergeMAFs {
            input:
                mafs = this_annotation,
                output_name = tumor_sample_name + ".annotated",
                compress_output = this_args.compress_output,
                runtime_params = runtime_collection.merge_mafs
        }
    }

    output {
        File annotated_variants = select_first(select_all([GatherVCFs.merged_vcf, MergeMAFs.merged_maf]))
        File? annotated_variants_idx = GatherVCFs.merged_vcf_idx
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
        Boolean prefer_mane_transcripts = true
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

    # TODO: Annotate haplotype phase

    # ==============
    # Process input args:
    String output_maf = output_base_name + ".maf"
    String output_vcf = output_base_name + if compress_output then ".vcf.gz" else ".vcf"
    String output_file = if output_format == "MAF" then output_maf else output_vcf
    String output_file_index = output_file +  if compress_output then ".tbi" else ".idx"

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

        echo "Creating dummy index file."
        touch '~{output_file_index}'

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
            ~{if prefer_mane_transcripts then "--prefer-mane-transcripts true" else ""} \
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
        File annotations_idx = output_file_index
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

task CreateEmptyAnnotation {
    input {
        File selected_vcf
        String output_base_name
        String output_format  # "VCF" or "MAF"
        Boolean compress_output

        Runtime runtime_params
    }

    # Determine if vcf is gzipped
    String uncompressed_vcf = basename(selected_vcf, ".gz")
    Boolean is_compressed = (uncompressed_vcf != basename(selected_vcf))

    String output_file = if output_format == "VCF" then output_base_name + ".vcf" else output_base_name + ".maf"
    String output_file_idx = output_file +  if compress_output then ".tbi" else ".idx"

    command <<<
        # Uncompress the selected_vcf if it is gzipped
        if [[ "~{is_compressed}" == "true" ]]; then
            gunzip -c '~{selected_vcf}' > '~{uncompressed_vcf}'
        else
            mv '~{selected_vcf}' '~{uncompressed_vcf}'
        fi

        if [ "~{output_format}" == "VCF" ]; then
            # Copy all headers from selected_vcf to create an empty VCF
            grep "^#" '~{uncompressed_vcf}' > '~{output_base_name}.vcf'

        elif [ "~{output_format}" == "MAF" ]; then
            # Copy header lines except VCF schema to create an empty MAF
            grep "^##" '~{uncompressed_vcf}' > '~{output_base_name}.maf'

            # Add MAF schema to the empty MAF
            echo -e "Hugo_Symbol\tEntrez_Gene_Id\tCenter\tNCBI_Build\tChromosome\tStart_Position\tEnd_Position\tStrand\t"\
            "Variant_Classification\tVariant_Type\tReference_Allele\tTumor_Seq_Allele1\tTumor_Seq_Allele2\tdbSNP_RS\t"\
            "dbSNP_Val_Status\tTumor_Sample_Barcode\tMatched_Norm_Sample_Barcode\tMatch_Norm_Seq_Allele1\t"\
            "Match_Norm_Seq_Allele2\tTumor_Validation_Allele1\tTumor_Validation_Allele2\t"\
            "Match_Norm_Validation_Allele1\tMatch_Norm_Validation_Allele2\tVerification_Status\tValidation_Status\t"\
            "Mutation_Status\tSequencing_Phase\tSequence_Source\tValidation_Method\tScore\tBAM_File\tSequencer\t"\
            "Tumor_Sample_UUID\tMatched_Norm_Sample_UUID\tHGVSc\tHGVSp\tHGVSp_Short\tTranscript_ID\tExon_Number\t"\
            "t_depth\tt_ref_count\tt_alt_count\tn_depth\tn_ref_count\tn_alt_count\tProtein_Change\tUniProt_AApos" >> '~{output_base_name}.maf'
        fi

        # Create an empty index file for VCF based on compress_output
        touch '~{output_file_idx}'
    >>>

    output {
        File empty_annotation = output_file
        File empty_annotation_idx = output_file_idx
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