version development

import "patient.wdl" as p
import "workflow_arguments.wdl" as wfargs
import "runtime_collection.wdl" as rtc
import "runtimes.wdl" as rt
import "tasks.wdl"


workflow FilterVariants {
    input {
        File vcf
        File vcf_idx
        File? mutect_stats
        File? orientation_bias

        Patient patient
        WorkflowArguments args
        RuntimeCollection runtime_collection
    }

    String vcf_name = basename(basename(vcf, ".gz"), ".vcf")

    scatter (sample in patient.samples) {
        File? c = sample.contamination
        File? s = sample.af_segmentation
    }
    if (length(select_all(c)) > 0) {
        Array[File] contamination_tables = select_all(c)
    }
    if (length(select_all(s)) > 0) {
        Array[File] segmentation_tables = select_all(s)
    }

    # From the documentation: "FilterMutectCalls goes over an unfiltered vcf in
    # three passes, two to learn any unknown parameters of the filters' models and
    # to set the threshold P(error), and one to apply the learned filters. [...]
    # [As such] it is critical to merge the unfiltered output of Mutect2 before
    # filtering."
    call FilterMutectCalls {
        input:
            ref_fasta = args.files.ref_fasta,
            ref_fasta_index = args.files.ref_fasta_index,
            ref_dict = args.files.ref_dict,
            vcf = vcf,
            vcf_idx = vcf_idx,
            orientation_bias = orientation_bias,
            contamination_tables = contamination_tables,
            tumor_segmentation = segmentation_tables,
            mutect_stats = mutect_stats,
            max_median_fragment_length_difference = args.filter_mutect2_max_median_fragment_length_difference,
            min_alt_median_base_quality = args.filter_mutect2_min_alt_median_base_quality,
            min_alt_median_mapping_quality = args.filter_mutect2_min_alt_median_mapping_quality,
            min_median_read_position = args.filter_mutect2_min_median_read_position,
            compress_output = args.compress_output,
            m2_filter_extra_args = args.filter_mutect2_extra_args,
            runtime_params = runtime_collection.filter_mutect_calls
    }

    call tasks.LeftAlignAndTrimVariants {
        input:
            ref_fasta = args.files.ref_fasta,
            ref_fasta_index = args.files.ref_fasta_index,
            ref_dict = args.files.ref_dict,
            vcf = FilterMutectCalls.filtered_vcf,
            vcf_idx = FilterMutectCalls.filtered_vcf_idx,
            split_multi_allelics = true,  # necessary for current SelectVariants implementation
            compress_output = args.compress_output,
            left_align_and_trim_variants_extra_args = args.left_align_and_trim_variants_extra_args,
            runtime_params = runtime_collection.left_align_and_trim_variants
    }

    if (args.run_variant_hard_filter) {
        call tasks.VariantFiltration as HardFilterVariants {
            input:
                ref_fasta = args.files.ref_fasta,
                ref_fasta_index = args.files.ref_fasta_index,
                ref_dict = args.files.ref_dict,
                vcf = LeftAlignAndTrimVariants.output_vcf,
                vcf_idx = LeftAlignAndTrimVariants.output_vcf_idx,
                compress_output = args.compress_output,
                filter_expressions = args.hard_filter_expressions,
                filter_names = args.hard_filter_names,
                variant_filtration_extra_args = args.variant_filtration_extra_args,
                runtime_params = runtime_collection.variant_filtration
        }
    }

    # TODO: add DeTiN

    call tasks.SelectVariants as SelectPassingVariants {
        input:
            ref_fasta = args.files.ref_fasta,
            ref_fasta_index = args.files.ref_fasta_index,
            ref_dict = args.files.ref_dict,
            vcf = select_first([HardFilterVariants.filtered_vcf, LeftAlignAndTrimVariants.output_vcf]),
            vcf_idx = select_first([HardFilterVariants.filtered_vcf_idx, LeftAlignAndTrimVariants.output_vcf_idx]),
            select_passing = true,
            keep_germline = args.keep_germline,
            germline_filter_whitelist = args.germline_filter_whitelist,
            compress_output = args.compress_output,
            select_variants_extra_args = args.select_variants_extra_args,
            runtime_params = runtime_collection.select_variants
    }

    if (args.run_realignment_filter) {
        # Realigning reads for variants can be expensive if many variants have been
        # called. Especially for tumor-only calling, plenty of variant calls are
        # sequencing artifacts or rare germline variants that have been missed by
        # FilterMutectCalls. In order to make the filter affordable, we divide
        # the called variants into groups of low and high somatic confidence.
        if (args.run_realignment_filter_only_on_high_confidence_variants) {
            call tasks.SelectVariants as SelectLowConfidenceVariants {
                input:
                    ref_fasta = args.files.ref_fasta,
                    ref_fasta_index = args.files.ref_fasta_index,
                    ref_dict = args.files.ref_dict,
                    vcf = SelectPassingVariants.selected_vcf,
                    vcf_idx = SelectPassingVariants.selected_vcf_idx,
                    compress_output = args.compress_output,
                    select_variants_extra_args = "-select '" + args.select_low_conficence_variants_jexl_arg + "'",
                    runtime_params = runtime_collection.select_variants
            }

            call tasks.SelectVariants as SelectHighConfidenceVariants {
                input:
                    ref_fasta = args.files.ref_fasta,
                    ref_fasta_index = args.files.ref_fasta_index,
                    ref_dict = args.files.ref_dict,
                    vcf = SelectPassingVariants.selected_vcf,
                    vcf_idx = SelectPassingVariants.selected_vcf_idx,
                    compress_output = args.compress_output,
                    select_variants_extra_args = "-select '" + args.select_low_conficence_variants_jexl_arg + "' -invertSelect true",
                    runtime_params = runtime_collection.select_variants
            }
        }

        File variants_to_realign = select_first([
            SelectHighConfidenceVariants.selected_vcf,
            SelectPassingVariants.selected_vcf
        ])
        File variants_to_realign_idx = select_first([
            SelectHighConfidenceVariants.selected_vcf_idx,
            SelectPassingVariants.selected_vcf_idx
        ])
        Int num_variants_to_realign = select_first([
            SelectHighConfidenceVariants.num_selected_variants,
            SelectPassingVariants.num_selected_variants
        ])
        scatter (tumor_sample in patient.tumor_samples) {
            scatter (seq_run in tumor_sample.sequencing_runs) {
                File seq_tumor_bams = seq_run.bam
                File seq_tumor_bais = seq_run.bai
            }
        }
        Array[File] tumor_bams = flatten(seq_tumor_bams)
        Array[File] tumor_bais = flatten(seq_tumor_bais)

        # Due to its long runtime, we scatter the realignment task over intervals.
        Int scatter_count = ceil((num_variants_to_realign + 1) / args.variants_per_scatter)
        call tasks.SplitIntervals {
            input:
                interval_list = args.preprocessed_interval_list,
                ref_fasta = args.files.ref_fasta,
                ref_fasta_index = args.files.ref_fasta_index,
                ref_dict = args.files.ref_dict,
                scatter_count = scatter_count,
                split_intervals_extra_args = args.split_intervals_extra_args,
                runtime_params = runtime_collection.split_intervals,
        }
        scatter (interval_list in SplitIntervals.interval_files) {
            call tasks.SelectVariants as SelectPreRealignmentVariants {
                input:
                    interval_list = interval_list,
                    ref_fasta = args.files.ref_fasta,
                    ref_fasta_index = args.files.ref_fasta_index,
                    ref_dict = args.files.ref_dict,
                    vcf = variants_to_realign,
                    vcf_idx = variants_to_realign_idx,
                    compress_output = args.compress_output,
                    runtime_params = runtime_collection.select_variants
            }

            if (SelectPreRealignmentVariants.num_selected_variants > 0) {
                call FilterAlignmentArtifacts {
                    input:
                        ref_fasta = args.files.ref_fasta,
                        ref_fasta_index = args.files.ref_fasta_index,
                        ref_dict = args.files.ref_dict,
                        tumor_bams = tumor_bams,
                        tumor_bais = tumor_bais,
                        vcf = SelectPreRealignmentVariants.selected_vcf,
                        vcf_idx = SelectPreRealignmentVariants.selected_vcf_idx,
                        bwa_mem_index_image = args.files.realignment_bwa_mem_index_image,
                        compress_output = args.compress_output,
                        max_reasonable_fragment_length = args.filter_alignment_artifacts_max_reasonable_fragment_length,
                        realignment_extra_args = args.realignment_extra_args,
                        runtime_params = runtime_collection.filter_alignment_artifacts
                }
            }
        }

        if (length(select_all(FilterAlignmentArtifacts.filtered_vcf)) > 0) {
            call tasks.MergeVCFs as MergeRealignmentFilteredVCFs {
                input:
                    vcfs = select_all(FilterAlignmentArtifacts.filtered_vcf),
                    vcfs_idx = select_all(FilterAlignmentArtifacts.filtered_vcf_idx),
                    output_name = vcf_name + ".filtered.selected.realignmentfiltered",
                    compress_output = args.compress_output,
                    runtime_params = runtime_collection.merge_vcfs
            }
        }

        call tasks.SelectVariants as SelectPostRealignmentVariants {
            input:
                ref_fasta = args.files.ref_fasta,
                ref_fasta_index = args.files.ref_fasta_index,
                ref_dict = args.files.ref_dict,
                vcf = select_first([MergeRealignmentFilteredVCFs.merged_vcf, variants_to_realign]),
                vcf_idx = select_first([MergeRealignmentFilteredVCFs.merged_vcf_idx, variants_to_realign_idx]),
                select_passing = true,
                keep_germline = args.keep_germline,
                germline_filter_whitelist = args.germline_filter_whitelist,
                compress_output = args.compress_output,
                select_variants_extra_args = args.select_variants_extra_args,
                runtime_params = runtime_collection.select_variants
        }

        if (args.run_realignment_filter_only_on_high_confidence_variants) {
            call tasks.MergeVCFs as MergeLowConfidenceAndRealignmentFilteredVCFs {
                input:
                    vcfs = select_all([
                        SelectLowConfidenceVariants.selected_vcf,
                        SelectPostRealignmentVariants.selected_vcf
                    ]),
                    vcfs_idx = select_all([
                        SelectLowConfidenceVariants.selected_vcf_idx,
                        SelectPostRealignmentVariants.selected_vcf_idx
                    ]),
                    output_name = vcf_name + ".filtered.selected.realignmentfiltered",
                    compress_output = args.compress_output,
                    runtime_params = runtime_collection.merge_vcfs
            }
        }
    }

    File selected_vcf = select_first([
        MergeLowConfidenceAndRealignmentFilteredVCFs.merged_vcf,
        SelectPostRealignmentVariants.selected_vcf,
        SelectPassingVariants.selected_vcf,
    ])
    File selected_vcf_idx = select_first([
        MergeLowConfidenceAndRealignmentFilteredVCFs.merged_vcf_idx,
        SelectPostRealignmentVariants.selected_vcf_idx,
        SelectPassingVariants.selected_vcf_idx,
    ])

    if (args.keep_germline) {
        call tasks.SelectVariants as SelectGermlineVariants {
            input:
                ref_fasta = args.files.ref_fasta,
                ref_fasta_index = args.files.ref_fasta_index,
                ref_dict = args.files.ref_dict,
                vcf = selected_vcf,
                vcf_idx = selected_vcf_idx,
                select_passing = false,
                keep_germline = true,
                suffix = ".germline",
                germline_filter_whitelist = args.germline_filter_whitelist,
                compress_output = args.compress_output,
                select_variants_extra_args = args.select_variants_extra_args,
                runtime_params = runtime_collection.select_variants
        }
    }

    call tasks.SelectVariants as SelectSomaticVariants {
        input:
            ref_fasta = args.files.ref_fasta,
            ref_fasta_index = args.files.ref_fasta_index,
            ref_dict = args.files.ref_dict,
            vcf = selected_vcf,
            vcf_idx = selected_vcf_idx,
            select_passing = true,
            keep_germline = false,
            suffix = ".somatic",
            compress_output = args.compress_output,
            select_variants_extra_args = args.select_variants_extra_args,
            runtime_params = runtime_collection.select_variants
    }

    output {
        File filtering_stats = FilterMutectCalls.filtering_stats
        File filtered_vcf = LeftAlignAndTrimVariants.output_vcf
        File filtered_vcf_idx = LeftAlignAndTrimVariants.output_vcf_idx
        File somatic_vcf = SelectSomaticVariants.selected_vcf
        File somatic_vcf_idx = SelectSomaticVariants.selected_vcf_idx
        Int num_somatic_variants = SelectSomaticVariants.num_selected_variants
        File? germline_vcf = SelectGermlineVariants.selected_vcf
        File? germline_vcf_idx = SelectGermlineVariants.selected_vcf_idx
        Int? num_germline_variants = SelectGermlineVariants.num_selected_variants
    }
}

task FilterMutectCalls {
    input {
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        File vcf
        File vcf_idx
        File? orientation_bias
        Array[File]? contamination_tables
        Array[File]? tumor_segmentation
        File? mutect_stats

        Int max_median_fragment_length_difference = 10000  # default: 10000
        Int min_alt_median_base_quality = 20  # default: 20
        Int min_alt_median_mapping_quality = 20  # default: -1
        Int min_median_read_position = 5  # default: 1

        Boolean compress_output = false
        String? m2_filter_extra_args

        Runtime runtime_params
    }

    Int disk = (
        ceil(size(vcf, "GB"))
        + ceil(size(orientation_bias, "GB"))
        + if defined(contamination_tables) then ceil(size(select_first([contamination_tables]), "GB")) else 0
        + if defined(tumor_segmentation) then ceil(size(select_first([tumor_segmentation]), "GB")) else 0
        + ceil(size(mutect_stats, "GB"))
        + runtime_params.disk
    )

    String output_base_name = basename(basename(vcf, ".gz"), ".vcf") + ".filtered"
    String output_vcf = output_base_name + if compress_output then ".vcf.gz" else ".vcf"
    String output_vcf_idx = output_vcf + if compress_output then ".tbi" else ".idx"

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{select_first([runtime_params.jar_override, "/root/gatk.jar"])}
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            FilterMutectCalls \
            --reference '~{ref_fasta}' \
            --variant '~{vcf}' \
            --output '~{output_vcf}' \
            ~{"--orientation-bias-artifact-priors '" + orientation_bias + "'"} \
            ~{true="--contamination-table '" false="" defined(contamination_tables)}~{default="" sep="' --contamination-table '" contamination_tables}~{true="'" false="" defined(contamination_tables)} \
            ~{true="--tumor-segmentation '" false="" defined(tumor_segmentation)}~{default="" sep="' --tumor-segmentation '" tumor_segmentation}~{true="'" false="" defined(tumor_segmentation)} \
            ~{"--stats '" + mutect_stats + "'"} \
            ~{"--max-median-fragment-length-difference " + max_median_fragment_length_difference} \
            ~{"--min-median-base-quality " + min_alt_median_base_quality} \
            ~{"--min-median-mapping-quality " + min_alt_median_mapping_quality} \
            ~{"--min-median-read-position " + min_median_read_position} \
            --filtering-stats '~{output_base_name}.stats' \
            --seconds-between-progress-updates 60 \
            ~{m2_filter_extra_args}
    >>>

    output {
        File filtered_vcf = output_vcf
        File filtered_vcf_idx = output_vcf_idx
        File filtering_stats = output_base_name + ".stats"
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

    # Optional localization leads to cromwell error.
    parameter_meta{
        ref_fasta: {localization_optional: true}
        ref_fasta_index: {localization_optional: true}
        ref_dict: {localization_optional: true}
        vcf: {localization_optional: true}
        vcf_idx: {localization_optional: true}
        # orientation_bias: {localization_optional: true}
        # contamination_tables: {localization_optional: true}
        # tumor_segmentation: {localization_optional: true}
        # mutect_stats: {localization_optional: true}
    }
}

task FilterAlignmentArtifacts {
    # Lifts the reads supporting the variant alleles from the tumor bams and aligns
    # these reads to the index image. If the reads supporting the variant allele have
    # ambiguous alignments, the tool will mark these variants as filtered in the VCF.

    input {
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        File vcf
        File vcf_idx
        Array[File] tumor_bams
        Array[File] tumor_bais

        File? bwa_mem_index_image

        Int max_reasonable_fragment_length = 10000 # default: 100000
        Boolean compress_output = false
        String? realignment_extra_args

        Runtime runtime_params
    }

    Int disk = ceil(size(bwa_mem_index_image, "GB")) + runtime_params.disk

    String output_base_name = basename(basename(vcf, ".gz"), ".vcf") + ".realignmentfiltered"
    String output_vcf = output_base_name + if compress_output then ".vcf.gz" else ".vcf"
    String output_vcf_idx = output_vcf + if compress_output then ".tbi" else ".idx"

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{select_first([runtime_params.jar_override, "/root/gatk.jar"])}
        # Not skipping filtered variants is important for keeping germline variants if requested.
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            FilterAlignmentArtifacts \
            ~{sep="' " prefix("-I '", tumor_bams)}' \
            ~{sep="' " prefix("--read-index '", tumor_bais)}' \
            --variant '~{vcf}' \
            --reference '~{ref_fasta}' \
            --bwa-mem-index-image '~{bwa_mem_index_image}' \
            --output '~{output_vcf}' \
            --max-reasonable-fragment-length ~{max_reasonable_fragment_length} \
            --dont-skip-filtered-variants true \
            --seconds-between-progress-updates 60 \
            ~{realignment_extra_args}
    >>>

    output {
        File filtered_vcf = output_vcf
        File filtered_vcf_idx = output_vcf_idx
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
        vcf: {localization_optional: true}
        vcf_idx: {localization_optional: true}
        tumor_bams: {localization_optional: true}
        tumor_bais: {localization_optional: true}
        ref_fasta: {localization_optional: true}
        ref_fasta_index: {localization_optional: true}
        ref_dict: {localization_optional: true}
        # bwa_mem_index_image: {localization_optional: true}  # needs to be localized
    }
}
