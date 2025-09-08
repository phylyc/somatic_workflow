version development

import "patient.wdl" as p
import "workflow_arguments.wdl" as wfargs
import "runtime_collection.wdl" as rtc
import "runtimes.wdl" as rt
import "tasks.wdl"


workflow FilterVariants {
    input {
        Patient patient
        WorkflowArguments args
        RuntimeCollection runtime_collection
    }

    if (defined(patient.raw_snv_calls_vcf) && !defined(patient.filtered_vcf)) {
        scatter (sample in patient.samples) {
            File? c = sample.contamination_table
            File? s = sample.af_segmentation_table
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
        call FilterVariantCalls {
            input:
                ref_fasta = args.files.ref_fasta,
                ref_fasta_index = args.files.ref_fasta_index,
                ref_dict = args.files.ref_dict,
                vcf = select_first([patient.raw_snv_calls_vcf]),
                vcf_idx = select_first([patient.raw_snv_calls_vcf_idx]),
                orientation_bias = patient.orientation_bias,
                contamination_tables = contamination_tables,
                tumor_segmentation = segmentation_tables,
                mutect_stats = patient.mutect2_stats,
                max_median_fragment_length_difference = args.filter_mutect2_max_median_fragment_length_difference,
                min_base_quality = args.hard_filter_min_base_quality,
                min_mapping_quality = args.hard_filter_min_mapping_quality,
                min_fragment_length = args.hard_filter_min_fragment_length,
                min_total_depth = args.hard_filter_min_total_depth,
                min_total_alt_count = args.hard_filter_min_total_alt_count,
                min_position_from_end_of_read = args.hard_filter_min_position_from_end_of_read,
                min_read_orientation_quality = args.hard_filter_min_read_orientation_quality,
                min_not_germline_quality = args.hard_filter_min_not_germline_quality,
                germline_min_population_af = args.hard_filter_germline_min_population_af,
                split_multi_allelics = true,  # necessary for current SelectVariants implementation
                filter_expressions = args.hard_filter_expressions,
                filter_names = args.hard_filter_names,
                compress_output = args.compress_output,
                m2_filter_extra_args = args.filter_mutect2_extra_args,
                left_align_and_trim_variants_extra_args = args.left_align_and_trim_variants_extra_args,
                variant_filtration_extra_args = args.variant_filtration_extra_args,
                runtime_params = runtime_collection.filter_variant_calls
        }

        # TODO: add DeTiN
        # (not needed because we also collect variants flagged as normal_artifact)

        if (args.run_realignment_filter) {
            # Realigning reads for variants can be expensive if many variants have been
            # called. Especially for tumor-only calling, plenty of variant calls are
            # sequencing artifacts or rare germline variants that have been missed by
            # FilterMutectCalls. In order to make the filter affordable, we divide
            # the called variants into groups of low and high somatic confidence.
            call tasks.SelectVariants as SelectVariantsToRealign {
                input:
                    ref_fasta = args.files.ref_fasta,
                    ref_fasta_index = args.files.ref_fasta_index,
                    ref_dict = args.files.ref_dict,
                    vcf = FilterVariantCalls.filtered_vcf,
                    vcf_idx = FilterVariantCalls.filtered_vcf_idx,
                    compress_output = args.compress_output,
                    select_somatic = true,
                    somatic_filter_whitelist = args.somatic_filter_whitelist,
                    germline_filter_whitelist = args.germline_filter_whitelist,
                    select_variants_extra_args = "-select '" + args.select_high_conficence_variants_jexl_arg + "'",
                    runtime_params = runtime_collection.select_variants
            }

            if (SelectVariantsToRealign.num_selected_variants > 0) {
                scatter (tumor_sample in patient.tumor_samples) {
                    scatter (seq_run in tumor_sample.sequencing_runs) {
                        File seq_tumor_bams = seq_run.bam
                        File seq_tumor_bais = seq_run.bai
                    }
                }
                Array[File] tumor_bams = flatten(seq_tumor_bams)
                Array[File] tumor_bais = flatten(seq_tumor_bais)

                # Due to its long runtime, we scatter the realignment task over intervals.
                Int scatter_count = ceil(SelectVariantsToRealign.num_selected_variants / args.variants_per_scatter) + 1
                call tasks.SplitIntervals {
                    input:
                        interval_list = args.files.preprocessed_intervals,
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
                            vcf = SelectVariantsToRealign.selected_vcf,
                            vcf_idx = SelectVariantsToRealign.selected_vcf_idx,
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

                call tasks.GatherVCFs as MergeLowConfidenceAndRealignmentFilteredVCFs {
                    input:
                        vcfs = flatten([
                            [SelectVariantsToRealign.not_selected_vcf],
                            select_all(FilterAlignmentArtifacts.filtered_vcf)
                        ]),
                        vcfs_idx = flatten([
                            [SelectVariantsToRealign.not_selected_vcf_idx],
                            select_all(FilterAlignmentArtifacts.filtered_vcf_idx)
                        ]),
                        output_name = patient.name + ".filtered.realignmentfiltered",
                        compress_output = args.compress_output,
                        runtime_params = runtime_collection.gather_vcfs
                }
            }
        }
    }

    File filtered_vcf = select_first([
        MergeLowConfidenceAndRealignmentFilteredVCFs.merged_vcf,
        SelectVariantsToRealign.not_selected_vcf,
        FilterVariantCalls.filtered_vcf,
        patient.filtered_vcf
    ])
    File filtered_vcf_idx = select_first([
        MergeLowConfidenceAndRealignmentFilteredVCFs.merged_vcf_idx,
        SelectVariantsToRealign.not_selected_vcf_idx,
        FilterVariantCalls.filtered_vcf_idx,
        patient.filtered_vcf_idx
    ])

    if (!defined(patient.somatic_vcf)) {
        call tasks.SelectVariants as SelectSomaticVariants {
            input:
                ref_fasta = args.files.ref_fasta,
                ref_fasta_index = args.files.ref_fasta_index,
                ref_dict = args.files.ref_dict,
                vcf = filtered_vcf,
                vcf_idx = filtered_vcf,
                select_somatic = true,
                suffix = ".somatic",
                somatic_filter_whitelist = args.somatic_filter_whitelist,
                compress_output = args.compress_output,
                select_variants_extra_args = args.select_variants_extra_args,
                runtime_params = runtime_collection.select_variants
        }
    }

    if (args.keep_germline && (!defined(patient.germline_vcf))) {
        call tasks.SelectVariants as SelectGermlineVariants {
            input:
                ref_fasta = args.files.ref_fasta,
                ref_fasta_index = args.files.ref_fasta_index,
                ref_dict = args.files.ref_dict,
                vcf = filtered_vcf,
                vcf_idx = filtered_vcf,
                select_germline = true,
                suffix = ".germline",
                germline_filter_whitelist = args.germline_filter_whitelist,
                compress_output = args.compress_output,
                select_variants_extra_args = args.select_variants_extra_args,
                runtime_params = runtime_collection.select_variants
        }
    }

    scatter (shard in patient.shards) {
        File? raw_mutect2_bam_out_scattered = shard.raw_mutect2_bam_out
        File? raw_mutect2_bai_out_scattered = shard.raw_mutect2_bai_out
    }

    # subset CallVariants.bam to reads covering the somatic_vcf only
    if ((length(select_all(raw_mutect2_bam_out_scattered)) > 0) && (!defined(patient.somatic_calls_bam))) {
        call tasks.PrintReads as SomaticBam {
            input:
                ref_fasta = args.files.ref_fasta,
                ref_fasta_index = args.files.ref_fasta_index,
                ref_dict = args.files.ref_dict,
                prefix = patient.name + ".somatic",
                bams = select_all(raw_mutect2_bam_out_scattered),
                bais = select_all(raw_mutect2_bai_out_scattered),
                vcf = select_first([SelectSomaticVariants.selected_vcf, patient.somatic_vcf]),
                vcf_idx = select_first([SelectSomaticVariants.selected_vcf_idx, patient.somatic_vcf_idx]),
                runtime_params = runtime_collection.print_reads
        }
    }

    call p.UpdatePatient {
        input:
            patient = patient,
            filtered_vcf = filtered_vcf,
            filtered_vcf_idx = filtered_vcf_idx,
            filtering_stats = FilterVariantCalls.filtering_stats,
            somatic_vcf = SelectSomaticVariants.selected_vcf,
            somatic_vcf_idx = SelectSomaticVariants.selected_vcf_idx,
            num_somatic_variants = SelectSomaticVariants.num_selected_variants,
            germline_vcf = SelectGermlineVariants.selected_vcf,
            germline_vcf_idx = SelectGermlineVariants.selected_vcf_idx,
            num_germline_variants = SelectGermlineVariants.num_selected_variants,
            somatic_calls_bam = SomaticBam.output_bam,
            somatic_calls_bai = SomaticBam.output_bai,
    }

    output {
        Patient updated_patient = UpdatePatient.updated_patient
    }
}

task FilterVariantCalls {
    input {
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        File? interval_list

        File vcf
        File vcf_idx
        File? orientation_bias
        Array[File]? contamination_tables
        Array[File]? tumor_segmentation
        File? mutect_stats

        Int min_base_quality = 20
        Int min_mapping_quality = 20
        Int min_fragment_length = 18
        Int min_total_depth = 10
        Int min_total_alt_count = 4
        Int min_position_from_end_of_read = 6
        Int min_read_orientation_quality = 20
        Int min_not_germline_quality = 30
        Float germline_min_population_af = 3
        Int max_median_fragment_length_difference = 10000  # default: 10000
        Int max_indel_length = 1000  # default: 200
        Boolean dont_trim_alleles = false
        Boolean split_multi_allelics = false

        Array[String] filter_expressions = []
        Array[String] filter_names = []

        Boolean compress_output = false
        String? m2_filter_extra_args
        String? left_align_and_trim_variants_extra_args
        String? variant_filtration_extra_args

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
    String output_stats = output_base_name + ".stats"

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{select_first([runtime_params.jar_override, "/root/gatk.jar"])}
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            FilterMutectCalls \
            --reference '~{ref_fasta}' \
            --variant '~{vcf}' \
            --output 'filtered.vcf.gz' \
            ~{"--orientation-bias-artifact-priors '" + orientation_bias + "'"} \
            ~{true="--contamination-table '" false="" defined(contamination_tables)}~{default="" sep="' --contamination-table '" contamination_tables}~{true="'" false="" defined(contamination_tables)} \
            ~{true="--tumor-segmentation '" false="" defined(tumor_segmentation)}~{default="" sep="' --tumor-segmentation '" tumor_segmentation}~{true="'" false="" defined(tumor_segmentation)} \
            ~{"--stats '" + mutect_stats + "'"} \
            ~{"--max-median-fragment-length-difference " + max_median_fragment_length_difference} \
            ~{"--min-median-base-quality " + min_base_quality} \
            ~{"--min-median-mapping-quality " + min_mapping_quality} \
            ~{"--min-median-read-position " + min_position_from_end_of_read} \
            --filtering-stats '~{output_stats}' \
            --seconds-between-progress-updates 60 \
            ~{m2_filter_extra_args}

        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            LeftAlignAndTrimVariants \
            -R '~{ref_fasta}' \
            -V 'filtered.vcf.gz' \
            --output 'filtered.left_aligned_and_trimmed.vcf.gz' \
            --max-indel-length ~{max_indel_length} \
            ~{if (dont_trim_alleles) then " --dont-trim-alleles " else ""} \
            ~{if (split_multi_allelics) then " --split-multi-allelics " else ""} \
            ~{left_align_and_trim_variants_extra_args}

        rm -f 'filtered.vcf.gz' 'filtered.vcf.gz.tbi'

        #####
        # Add TAC field to INFO field that contains the sum of ALT allele counts
        # across samples to filter on.
        ####

        echo -e '##fileformat=VCFv4.2' > total_alt_count.vcf
        echo -e '##INFO=<ID=TAC,Number=1,Type=Integer,Description="Sum of ALT allele counts across samples">' >> total_alt_count.vcf
        echo -e '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO' >> total_alt_count.vcf

        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' \
            'filtered.left_aligned_and_trimmed.vcf.gz' \
            | awk 'BEGIN{OFS="\t"} {
                sum = 0
                for (i = 5; i <= NF; i++) {
                    n = split($i, ad, ",")
                    for (j = 2; j <= n; j++) {
                        sum += ad[j]  # sum only ALT counts
                    }
                }
                print $1, $2, ".", $3, $4, ".", ".", "TAC=" sum
            }' >> total_alt_count.vcf

        bgzip -f total_alt_count.vcf
        bcftools index total_alt_count.vcf.gz

        bcftools annotate \
            -o 'filtered.left_aligned_and_trimmed.TAC.vcf.gz' \
            -O z \
            -c CHROM,POS,INFO/TAC \
            -a total_alt_count.vcf.gz \
            'filtered.left_aligned_and_trimmed.vcf.gz'
        gatk IndexFeatureFile \
            -I 'filtered.left_aligned_and_trimmed.TAC.vcf.gz'

        rm -f total_alt_count.vcf.gz total_alt_count.vcf.gz.tbi total_alt_count.vcf.csi
        rm -f 'filtered.left_aligned_and_trimmed.vcf.gz' 'filtered.left_aligned_and_trimmed.vcf.gz.tbi'

        echo ""
        # Some variants don't have certain INFO fields, so we suppress the warning messages.
        echo "Suppressing the following warning message: 'WARN  JexlEngine - '"
        echo ""

        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            VariantFiltration \
            ~{"-R '" + ref_fasta + "'"} \
            ~{"-L '" + interval_list + "'"} \
            -V 'filtered.left_aligned_and_trimmed.TAC.vcf.gz' \
            --filter-name "lowDP" \
            --filter-expression 'DP < ~{min_total_depth}' \
            --filter-name "lowTAC" \
            --filter-expression 'TAC < ~{min_total_alt_count}' \
            --filter-name "lowMBQ.0" \
            --filter-expression 'MBQ.0 < ~{min_base_quality}' \
            --filter-name "lowMBQ.1" \
            --filter-expression 'MBQ.1 < ~{min_base_quality}' \
            --filter-name "lowMMQ.0" \
            --filter-expression 'MMQ.0 < ~{min_mapping_quality}' \
            --filter-name "lowMMQ.1" \
            --filter-expression 'MMQ.1 < ~{min_mapping_quality}' \
            --filter-name "lowMFRL.0" \
            --filter-expression 'MFRL.0 < ~{min_fragment_length}' \
            --filter-name "lowMFRL.1" \
            --filter-expression 'MFRL.1 < ~{min_fragment_length}' \
            --filter-name "lowMPOS" \
            --filter-expression 'MPOS.0 < ~{min_position_from_end_of_read}' \
            --filter-name "lowROQ" \
            --filter-expression 'ROQ < ~{min_read_orientation_quality}' \
            --filter-name "lowGERMQ" \
            --filter-expression 'GERMQ < ~{min_not_germline_quality}' \
            --filter-name "germline" \
            --filter-expression 'POPAF < ~{germline_min_population_af}' \
            ~{if (length(filter_names) > 0) then " --filter-name '" else ""}~{default="" sep="' --filter-name '" filter_names}~{if (length(filter_names) > 0) then "'" else ""} \
            ~{if (length(filter_expressions) > 0) then " --filter-expression '" else ""}~{default="" sep="' --filter-expression '" filter_expressions}~{if (length(filter_expressions) > 0) then "'" else ""} \
            --output '~{output_vcf}' \
            ~{variant_filtration_extra_args} \
            2> >(grep -v "WARN  JexlEngine - " >&2)

        rm -f 'filtered.left_aligned_and_trimmed.TAC.vcf.gz' 'filtered.left_aligned_and_trimmed.TAC.vcf.gz.tbi'
    >>>

    output {
        File filtered_vcf = output_vcf
        File filtered_vcf_idx = output_vcf_idx
        File filtering_stats = output_stats
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

    parameter_meta{
        interval_list: {localization_optional: true}
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
