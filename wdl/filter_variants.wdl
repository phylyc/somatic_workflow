version development

import "runtimes.wdl"
import "tasks.wdl"


workflow FilterVariants {
    input {
        Array[File] scattered_interval_list
        File ref_fasta
        File ref_fasta_index

        Array[File]? tumor_bams
        Array[File]? tumor_bais
        File vcf
        File vcf_idx
        File? mutect_stats
        File? orientation_bias
        Array[File]? contamination_tables
        Array[File]? segmentation_tables

        # resources
        File? bwa_mem_index_image

        # workflow options
        Boolean run_realignment_filter = true
        Boolean run_realignment_filter_only_on_high_confidence_variants = false

        Boolean compress_output = true

        # arguments
        Int max_median_fragment_length_difference = 10000  # default: 10000
        Int min_alt_median_base_quality = 20  # default: 20
        Int min_alt_median_mapping_quality = 20  # default: -1
        Int max_reasonable_fragment_length = 10000 # default: 100000

        # expose extra arguments for import of this workflow
        String? filter_mutect2_extra_args
        String? select_variants_extra_args
        String? select_low_conficence_variants_jexl_arg = "'(vc.getAttribute(\"GERMQ\") < 30) || (vc.getAttribute(\"DP\") < 4) || (vc.getAttribute(\"MBQ\").0 == 0) || (vc.getAttribute(\"MFRL\").0 == 0)'"
        String? realignment_extra_args

        Runtime filter_mutect_calls_runtime = Runtimes.filter_mutect_calls_runtime
        Runtime filter_alignment_artifacts_runtime = Runtimes.filter_alignment_artifacts_runtime
        Runtime select_variants_runtime = Runtimes.select_variants_runtime
        Runtime merge_vcfs_runtime = Runtimes.merge_vcfs_runtime

        String gatk_docker = "broadinstitute/gatk"
        File? gatk_override
        Int preemptible = 1
        Int max_retries = 1
        Int cpu = 1
        Int disk_sizeGB = 1

        Int mem_filter_mutect_calls = 4096
        Int mem_filter_alignment_artifacts_base = 2048  # needs to be increased in some cases
        Int mem_select_variants = 1024
        Int mem_merge_vcfs = 512
        Int time_startup = 10
        Int time_filter_mutect_calls = 800  # 13 h
        Int time_filter_alignment_artifacts_total = 10000  # 12 d / scatter_count
        Int time_select_variants = 5
        Int time_merge_vcfs = 10
    }

    Int scatter_count = if defined(scattered_interval_list) then length(select_first([scattered_interval_list])) else 1

    call runtimes.DefineRuntimes as Runtimes {
        input:
            scatter_count = scatter_count,
            gatk_docker = gatk_docker,
            gatk_override = gatk_override,
            max_retries = max_retries,
            preemptible = preemptible,
            cpu = cpu,
            disk_sizeGB = disk_sizeGB,
            mem_filter_mutect_calls = mem_filter_mutect_calls,
            mem_select_variants = mem_select_variants,
            mem_filter_alignment_artifacts_base = mem_filter_alignment_artifacts_base,
            time_startup = time_startup,
            time_filter_mutect_calls = time_filter_mutect_calls,
            time_select_variants = time_select_variants,
            time_filter_alignment_artifacts_total = time_filter_alignment_artifacts_total,
    }

    String vcf_name = basename(basename(vcf, ".gz"), ".vcf")

    # From the documentation: "FilterMutectCalls goes over an unfiltered vcf in
    # three passes, two to learn any unknown parameters of the filters' models and
    # to set the threshold P(error), and one to apply the learned filters. [...]
    # [As such] it is critical to merge the unfiltered output of Mutect2 before
    # filtering."
    call FilterMutectCalls {
        input:
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            vcf = vcf,
            vcf_idx = vcf_idx,
            orientation_bias = orientation_bias,
            contamination_tables = contamination_tables,
            tumor_segmentation = segmentation_tables,
            mutect_stats = mutect_stats,
            max_median_fragment_length_difference = max_median_fragment_length_difference,
            min_alt_median_base_quality = min_alt_median_base_quality,
            min_alt_median_mapping_quality = min_alt_median_mapping_quality,
            compress_output = compress_output,
            m2_filter_extra_args = filter_mutect2_extra_args,
            runtime_params = filter_mutect_calls_runtime
    }

    call tasks.SelectVariants as SelectPassingVariants {
        input:
            vcf = FilterMutectCalls.filtered_vcf,
            vcf_idx = FilterMutectCalls.filtered_vcf_idx,
            select_passing = true,
            keep_germline = true,
            compress_output = compress_output,
            select_variants_extra_args = select_variants_extra_args,
            runtime_params = select_variants_runtime
    }

    if (run_realignment_filter && defined(tumor_bams) && defined(tumor_bais)) {
        # Realigning reads for variants can be expensive if many variants have been
        # called. Especially for tumor-only calling, plenty of variant calls are
        # still sequencing artifacts of sometimes obviously low quality that have
        # been missed by FilterMutectCalls. In order to make the filter affordable,
        # we divide the called variants into low and high confidence groups based
        # on read depth and VAF. Variants that come from reads that only support the
        # alternate allele are suspect. For those variants, the MBQ and MFRL are set
        # to zero.
        if (run_realignment_filter_only_on_high_confidence_variants) {
            call tasks.SelectVariants as SelectLowConfidenceVariants {
                input:
                    vcf = SelectPassingVariants.selected_vcf,
                    vcf_idx = SelectPassingVariants.selected_vcf_idx,
                    compress_output = compress_output,
                    select_variants_extra_args = "-select " + select_low_conficence_variants_jexl_arg,
                    runtime_params = select_variants_runtime
            }

            call tasks.SelectVariants as SelectHighConfidenceVariants {
                input:
                    vcf = SelectPassingVariants.selected_vcf,
                    vcf_idx = SelectPassingVariants.selected_vcf_idx,
                    compress_output = compress_output,
                    select_variants_extra_args = "-select " + select_low_conficence_variants_jexl_arg + " -invertSelect true",
                    runtime_params = select_variants_runtime
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

        # Due to its long runtime, we scatter the realignment task over intervals.
        scatter (interval_list in scattered_interval_list) {
            call tasks.SelectVariants as SelectPreRealignmentVariants {
                input:
                    interval_list = interval_list,
                    ref_fasta = ref_fasta,
                    ref_fasta_index = ref_fasta_index,
                    vcf = variants_to_realign,
                    vcf_idx = variants_to_realign_idx,
                    compress_output = compress_output,
                    runtime_params = select_variants_runtime
            }

            call FilterAlignmentArtifacts {
                input:
                    ref_fasta = ref_fasta,
                    ref_fasta_index = ref_fasta_index,
                    tumor_bams = select_first([tumor_bams]),
                    tumor_bais = select_first([tumor_bais]),
                    vcf = SelectPreRealignmentVariants.selected_vcf,
                    vcf_idx = SelectPreRealignmentVariants.selected_vcf_idx,
                    bwa_mem_index_image = bwa_mem_index_image,
                    compress_output = compress_output,
                    max_reasonable_fragment_length = max_reasonable_fragment_length,
                    realignment_extra_args = realignment_extra_args,
                    runtime_params = filter_alignment_artifacts_runtime
            }
        }

        call tasks.MergeVCFs as MergeRealignmentFilteredVCFs {
            input:
                vcfs = FilterAlignmentArtifacts.filtered_vcf,
                vcfs_idx = FilterAlignmentArtifacts.filtered_vcf_idx,
                output_name = vcf_name + ".filtered.selected.realignmentfiltered",
                compress_output = compress_output,
                runtime_params = select_variants_runtime
        }

        call tasks.SelectVariants as SelectPostRealignmentVariants {
            input:
                vcf = MergeRealignmentFilteredVCFs.merged_vcf,
                vcf_idx = MergeRealignmentFilteredVCFs.merged_vcf_idx,
                select_passing = true,
                keep_germline = true,
                compress_output = compress_output,
                select_variants_extra_args = select_variants_extra_args,
                runtime_params = select_variants_runtime
        }

        if (run_realignment_filter_only_on_high_confidence_variants) {
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
                    output_name = vcf_name + ".filtered.selected.realignmentfiltered.selected.merged",
                    compress_output = compress_output,
                    runtime_params = merge_vcfs_runtime
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

    call tasks.SelectVariants as SelectGermlineVariants {
        input:
            vcf = selected_vcf,
            vcf_idx = selected_vcf_idx,
            select_passing = false,
            keep_germline = true,
            compress_output = compress_output,
            select_variants_extra_args = select_variants_extra_args,
            runtime_params = select_variants_runtime
    }

    call tasks.SelectVariants as SelectSomaticVariants {
        input:
            vcf = selected_vcf,
            vcf_idx = selected_vcf_idx,
            select_passing = true,
            keep_germline = false,
            compress_output = compress_output,
            select_variants_extra_args = select_variants_extra_args,
            runtime_params = select_variants_runtime
    }

    output {
        File filtered_vcf = FilterMutectCalls.filtered_vcf
        File filtered_vcf_idx = FilterMutectCalls.filtered_vcf_idx
        File somatic_vcf = SelectSomaticVariants.selected_vcf
        File somatic_vcf_idx = SelectSomaticVariants.selected_vcf_idx
        File germline_vcf = SelectGermlineVariants.selected_vcf
        File germline_vcf_idx = SelectGermlineVariants.selected_vcf_idx
        File filtering_stats = FilterMutectCalls.filtering_stats

    }
}

task FilterMutectCalls {
    input {
        File ref_fasta
        File ref_fasta_index

        File vcf
        File vcf_idx
        File? orientation_bias
        Array[File]? contamination_tables
        Array[File]? tumor_segmentation
        File? mutect_stats

        Int max_median_fragment_length_difference = 10000  # default: 10000
        Int min_alt_median_base_quality = 20  # default: 20
        Int min_alt_median_mapping_quality = 20  # default: -1

        Boolean compress_output = false
        String? m2_filter_extra_args

        Runtime runtime_params
    }

    # Optional localization leads to cromwell error.
    parameter_meta{
        ref_fasta: {localization_optional: true}
        ref_fasta_index: {localization_optional: true}
        vcf: {localization_optional: true}
        vcf_idx: {localization_optional: true}
        # orientation_bias: {localization_optional: true}
        # contamination_tables: {localization_optional: true}
        # tumor_segmentation: {localization_optional: true}
        # mutect_stats: {localization_optional: true}
    }

    Int disk = (
        ceil(size(vcf, "GB"))
        + ceil(size(orientation_bias, "GB"))
        + ceil(size(contamination_tables, "GB"))
        + ceil(size(tumor_segmentation, "GB"))
        + ceil(size(mutect_stats, "GB"))
        + runtime_params.disk
    )

    String output_base_name = basename(basename(vcf, ".gz"), ".vcf") + ".filtered"
    String output_vcf = output_base_name + if compress_output then ".vcf.gz" else ".vcf"
    String output_vcf_idx = output_vcf + if compress_output then ".tbi" else ".idx"

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.jar_override}
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
            --filtering-stats '~{output_base_name}.stats' \
            --seconds-between-progress-updates 300 \
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
}

task FilterAlignmentArtifacts {
    # Lifts the reads supporting the variant alleles from the tumor bams and aligns
    # these reads to Hg38. If the reads supporting the variant allele have ambiguous
    # alignments, the tool will mark these variants as filtered in the VCF.

    input {
        File ref_fasta
        File ref_fasta_index
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

    parameter_meta {
        vcf: {localization_optional: true}
        vcf_idx: {localization_optional: true}
        tumor_bams: {localization_optional: true}
        tumor_bais: {localization_optional: true}
        ref_fasta: {localization_optional: true}
        ref_fasta_index: {localization_optional: true}
        # bwa_mem_index_image: {localization_optional: true}  # needs to be localized
    }

    Int disk = ceil(size(bwa_mem_index_image, "GB")) + runtime_params.disk

    String output_base_name = basename(basename(vcf, ".gz"), ".vcf") + ".realignmentfiltered"
    String output_vcf = output_base_name + if compress_output then ".vcf.gz" else ".vcf"
    String output_vcf_idx = output_vcf + if compress_output then ".tbi" else ".idx"

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.jar_override}

        if ~{!defined(bwa_mem_index_image)} ; then
            echo "ERROR: bwa_mem_index_image must be supplied."
            false
        fi

        # Not skipping filtered variants is important for keeping germline variants if requested.
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            FilterAlignmentArtifacts \
            ~{sep="' " prefix("-I '", tumor_bams)}' \
            --variant '~{vcf}' \
            --reference '~{ref_fasta}' \
            --bwa-mem-index-image '~{bwa_mem_index_image}' \
            --output '~{output_vcf}' \
            --max-reasonable-fragment-length ~{max_reasonable_fragment_length} \
            --dont-skip-filtered-variants true \
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
}
