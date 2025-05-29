version development

import "../runtime_collection.wdl" as rtc
import "../runtimes.wdl" as rt
import "../workflow_arguments.wdl" as wfargs
import "../workflow_resources.wdl" as wfres
import "../tasks.wdl"
import "../multi-sample_somatic_workflow.wdl" as mssw


workflow CreateSNVPanelOfNormals {
    input {
        String pon_name

        File ref_fasta
        File ref_fasta_index
        File ref_dict

        File germline_resource
        File germline_resource_idx

        File interval_list
        Int scatter_count = 10

        Array[File]? normal_bams
        Array[File]? normal_bais
        File? normal_bams_file
        File? normal_bais_file

        String mutect2_extra_args = ""

        Int min_contig_size = 1000000
        Int num_contigs = 24

        WorkflowArguments args = Parameters.arguments
        WorkflowResources resources = Files.resources
        RuntimeCollection runtime_collection = RuntimeParameters.rtc
    }

    call rtc.DefineRuntimeCollection as RuntimeParameters {
        input:
            scatter_count_for_variant_calling = scatter_count,
    }

    call wfres.DefineWorkflowResources as Files {
        input:
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            germline_resource = germline_resource,
            germline_resource_idx = germline_resource_idx,
            interval_list = interval_list
    }

    call wfargs.DefineWorkflowArguments as Parameters {
        input:
            scatter_count_base_for_variant_calling = scatter_count,
            resources = resources,

            run_collect_callable_loci = false,
            run_collect_total_read_counts = false,
            run_collect_allelic_read_counts = false,
            run_contamination_model = false,
            run_model_segments = false,

            run_orientation_bias_mixture_model = false,
            run_variant_calling = true,
            run_variant_filter = false,
            run_realignment_filter = false,
            run_variant_annotation = false,
            run_clonal_decomposition = false,

            keep_germline = false,
            make_bamout = false,

            mutect2_extra_args = mutect2_extra_args + " --max-mnp-distance 0",  # GenpmicsDBImport can't handle MNPs

            runtime_collection = runtime_collection,
    }

    # todo: assert either normal_bams or normal_bams_file is defined

    Array[File] non_optional_normal_bams = if defined(normal_bams) then select_first([normal_bams, []]) else read_lines(select_first([normal_bams_file, ""]))
    Array[File] non_optional_normal_bais = if defined(normal_bais) then select_first([normal_bais, []]) else read_lines(select_first([normal_bais_file, ""]))

    scatter (normal in zip(non_optional_normal_bams, non_optional_normal_bais)) {
        call tasks.GetSampleName {
            input:
                bam = normal.left,
                bai = normal.right,
                runtime_params = runtime_collection.get_sample_name,
        }

        call mssw.MultiSampleSomaticWorkflow {
            input:
                patient_id = GetSampleName.sample_name,
                bams = [normal.left],
                bais = [normal.right],
                target_intervals = select_all([args.files.preprocessed_intervals]),
                input_args = args,
                input_resources = resources,
                input_runtime_collection = runtime_collection,
        }
    }

    String split_intervals_extra_args = (
        "--subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION"
        + " --min-contig-size " + min_contig_size
    )
    call tasks.SplitIntervals {
        input:
            interval_list = args.files.preprocessed_intervals,
            ref_fasta = args.files.ref_fasta,
            ref_fasta_index = args.files.ref_fasta_index,
            ref_dict = args.files.ref_dict,
            scatter_count = num_contigs,
            split_intervals_extra_args = split_intervals_extra_args,
            runtime_params = runtime_collection.split_intervals
    }

    call rt.UpdateRuntimeParameters as CreateSNVPanelRuntime {
        input:
            runtime_params = runtime_collection.create_mutect2_panel,
            disk = runtime_collection.create_mutect2_panel.disk + 10 * ceil(size(MultiSampleSomaticWorkflow.raw_snv_calls_vcf, "GB")) + ceil(length(MultiSampleSomaticWorkflow.raw_snv_calls_vcf_idx) / 10)
    }

    scatter (scattered_intervals in SplitIntervals.interval_files) {
        call CreateSNVPanel {
            input:
                input_vcfs = select_all(MultiSampleSomaticWorkflow.raw_snv_calls_vcf),
                input_vcf_indices = select_all(MultiSampleSomaticWorkflow.raw_snv_calls_vcf_idx),
                interval_list = scattered_intervals,
                ref_fasta = args.files.ref_fasta,
                ref_fasta_index = args.files.ref_fasta_index,
                ref_dict = args.files.ref_dict,
                compress_output = args.compress_output,
                gnomad = args.files.germline_resource,
                gnomad_idx = args.files.germline_resource_idx,
                output_vcf_name = pon_name,
                runtime_params = CreateSNVPanelRuntime.params
        }
    }

    call tasks.GatherVCFs {
        input:
            ref_fasta = args.files.ref_fasta,
            ref_fasta_index = args.files.ref_fasta_index,
            ref_dict = args.files.ref_dict,
            vcfs = CreateSNVPanel.output_vcf,
            vcfs_idx = CreateSNVPanel.output_vcf_index,
            output_name = pon_name,
            compress_output = args.compress_output,
            runtime_params = runtime_collection.gather_vcfs
    }

    output {
        File pon = GatherVCFs.merged_vcf
        File pon_idx = GatherVCFs.merged_vcf_idx
        Array[File] normal_calls = select_all(MultiSampleSomaticWorkflow.raw_snv_calls_vcf)
        Array[File] normal_calls_idx = select_all(MultiSampleSomaticWorkflow.raw_snv_calls_vcf_idx)
    }
}

task CreateSNVPanel {
    input {
        File interval_list
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        Array[File]+ input_vcfs
        Array[File]+ input_vcf_indices
        String output_vcf_name

        Boolean compress_output = true
        File? gnomad
        File? gnomad_idx

        # UMCCR found that 5 is optimizing the F2 score, but not by much compared to 2:
        # https://umccr.org/blog/panel-of-normals/
        Int min_sample_count = 2
        String? create_pon_extra_args

        # runtime
        Runtime runtime_params
    }

    String pon_file = "pon.vcf.gz"
    String output_file = output_vcf_name + if compress_output then ".vcf.gz" else ".vcf"
    String output_file_idx = output_file + if compress_output then ".tbi" else ".idx"

    command {
        set -e
        export GATK_LOCAL_JAR=~{select_first([runtime_params.jar_override, "/root/gatk.jar"])}

        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            GenomicsDBImport \
            --genomicsdb-workspace-path pon_db \
            -R '~{ref_fasta}' \
            ~{sep="' " prefix("-V '", input_vcfs)}' \
            -L '~{interval_list}'

        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            CreateSomaticPanelOfNormals \
            -R '~{ref_fasta}' \
            ~{"--germline-resource '" + gnomad + "'"} \
            -V gendb://pon_db \
            -O '~{pon_file}' \
            --min-sample-count ~{min_sample_count} \
            ~{create_pon_extra_args}

        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            SortVcf \
            -I '~{pon_file}' \
            -O '~{output_file}'
    }

    output {
        File output_vcf = output_file
        File output_vcf_index = output_file_idx
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

    # GenomicsDB requires that the reference be a local file.
    parameter_meta{
        input_vcfs: {localization_optional: true}
        input_vcf_indices: {localization_optional: true}
        gnomad: {localization_optional: true}
        gnomad_idx: {localization_optional: true}
    }
}