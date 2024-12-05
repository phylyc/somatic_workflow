version development

import "patient.wdl" as p
import "patient.update_samples.wdl" as p_update_s
import "workflow_arguments.wdl" as wfargs
import "runtime_collection.wdl" as rtc
import "runtimes.wdl" as rt


workflow ModelSegments {
    input {
        Patient patient
        WorkflowArguments args
        RuntimeCollection runtime_collection
    }

    # The allelic counts were collected via GetPileupSummaries, which contains
    # more information. We first need to convert the pileup to allelic count format.

    scatter (sample in patient.samples) {
        if (defined(sample.snppanel_pileups)) {
            call PileupToAllelicCounts {
                input:
                    ref_dict = args.files.ref_dict,
                    pileup = sample.snppanel_pileups,
                    gvcf = patient.gvcf,
                    gvcf_idx = patient.gvcf_idx,
                    select_hets = true,
                    bam_name = sample.bam_name,
                    runtime_params = runtime_collection.pileup_to_allelic_counts
            }
        }
    }
    Array[File] sample_allelic_counts = select_all(PileupToAllelicCounts.allelic_counts)
    Array[Float] sample_error_probabilities = select_all(PileupToAllelicCounts.error_probability)

    # Save the allelic counts in each sample object (if they exist)

    if (length(sample_allelic_counts) > 0) {
        call p_update_s.UpdateSamples as AddAllelicCountsToSamples {
            input:
                patient = patient,
                snppanel_allelic_counts = sample_allelic_counts,
                genotype_error_probabilities = sample_error_probabilities
        }
    }
    Patient pat = select_first([AddAllelicCountsToSamples.updated_patient, patient])

    # Prepare ModelSegments input

    scatter (sample in pat.samples) {
        File? denoised_copy_ratios = sample.denoised_copy_ratios
        File? allelic_counts = sample.snppanel_allelic_counts
    }
    if (length(select_all(denoised_copy_ratios)) > 0) {
        Array[File] dcr = select_all(denoised_copy_ratios)
    }
    if (length(select_all(allelic_counts)) > 0) {
        Array[File] ac = select_all(allelic_counts)
    }
    if (defined(pat.matched_normal_sample)) {
        Sample matched_normal_sample = select_first([pat.matched_normal_sample])
        File? normal_allelic_counts = matched_normal_sample.snppanel_allelic_counts
    }

    # Now we can finally do the segmentation. First, we determine the patient-specific
    # segmentation, then we infer the sample-specific copy ratios and call amps/dels.
    # We set genotypung_homozygous_log_ratio_threshold = 10 since PileupToAllelicCounts
    # selects only hets, so the log odds ratio of hom/het should be very permissive.

    if (length(pat.samples) > 1) {
        call ModelSegmentsTask as MultiSampleModelSegments {
            input:
                denoised_copy_ratios = dcr,
                allelic_counts = ac,
                normal_allelic_counts = normal_allelic_counts,
                prefix = patient.name + ".segmentation",
                window_sizes = args.model_segments_window_sizes,
                genotyping_homozygous_log_ratio_threshold = 10,
                minimum_total_allele_count_normal = 0,  # all supplied loci are hets
                runtime_params = runtime_collection.model_segments
        }
    }

    scatter (sample in pat.samples) {
        if (defined(sample.denoised_copy_ratios))  {
            Array[File] dcr_list = select_all([sample.denoised_copy_ratios])
        }
        if (defined(sample.snppanel_allelic_counts)) {
            Array[File] ac_list = select_all([sample.snppanel_allelic_counts])
        }
        if (defined(sample.genotype_error_probabilities)) {
            Float error_probability = select_first([sample.genotype_error_probabilities])
        }
        call ModelSegmentsTask as SingleSampleInferCR {
            input:
                segments = MultiSampleModelSegments.multi_sample_segments,
                denoised_copy_ratios = dcr_list,
                allelic_counts = ac_list,
                normal_allelic_counts = normal_allelic_counts,
                genotyping_base_error_rate = error_probability,
                prefix = sample.name,
                window_sizes = args.model_segments_window_sizes,
                genotyping_homozygous_log_ratio_threshold = 10,
                minimum_total_allele_count_normal = 0,  # all supplied loci are hets
                runtime_params = runtime_collection.model_segments
        }

        call CallCopyRatioSegments {
            input:
                copy_ratio_segments = select_first([SingleSampleInferCR.cr_seg]),
                neutral_segment_copy_ratio_lower_bound = args.call_copy_ratios_neutral_segment_copy_ratio_lower_bound,
                neutral_segment_copy_ratio_upper_bound = args.call_copy_ratios_neutral_segment_copy_ratio_upper_bound,
                outlier_neutral_segment_copy_ratio_z_score_threshold = args.call_copy_ratios_outlier_neutral_segment_copy_ratio_z_score_threshold,
                calling_copy_ratio_z_score_threshold = args.call_copy_ratios_z_score_threshold,
                runtime_params = runtime_collection.call_copy_ratio_segments
        }

        call MergeCallsWithModeledSegments {
            input:
                copy_ratio_segments = select_first([SingleSampleInferCR.seg_final]),
                called_copy_ratio_segments = CallCopyRatioSegments.called_cr_seg,
                runtime_params = runtime_collection.merge_calls_with_modeled_segments
        }

        call PlotModeledSegments {
            input:
                ref_dict = args.files.ref_dict,
                sample_name = sample.name,
                segments = select_first([SingleSampleInferCR.seg_final]),
                denoised_copy_ratios = sample.denoised_copy_ratios,
                het_allelic_counts = SingleSampleInferCR.hets,
                runtime_params = runtime_collection.plot_modeled_segments
        }
    }

    call p_update_s.UpdateSamples as AddSegmentationResultsToSamples {
        input:
            patient = pat,
            af_model_parameters = select_all(SingleSampleInferCR.af_model_final_parameters),
            cr_model_parameters = select_all(SingleSampleInferCR.cr_model_final_parameters),
            called_copy_ratio_segmentations =  MergeCallsWithModeledSegments.merged_segments,
    }

    output {
        Patient updated_patient = AddSegmentationResultsToSamples.updated_patient

        Array[File] snppanel_allelic_counts = sample_allelic_counts

        File? modeled_segments = MultiSampleModelSegments.multi_sample_segments
        Array[File] hets = select_all(SingleSampleInferCR.hets)
        Array[File] af_model_begin_parameters = select_all(SingleSampleInferCR.af_model_begin_parameters)
        Array[File] cr_model_begin_parameters = select_all(SingleSampleInferCR.cr_model_begin_parameters)
        Array[File] af_model_final_parameters = select_all(SingleSampleInferCR.af_model_final_parameters)
        Array[File] cr_model_final_parameters = select_all(SingleSampleInferCR.cr_model_final_parameters)
        Array[File] igv_af = select_all(SingleSampleInferCR.igv_af)
        Array[File] igv_cr = select_all(SingleSampleInferCR.igv_cr)
        Array[File] seg_begin = select_all(SingleSampleInferCR.seg_begin)
        Array[File] called_copy_ratio_segmentations = MergeCallsWithModeledSegments.merged_segments
        Array[File] cr_plots = PlotModeledSegments.plot
    }
}

task PileupToAllelicCounts {
    input {
        File ref_dict
        String bam_name
        File? pileup
        File? gvcf
        File? gvcf_idx
        Int min_read_depth = 0
        Boolean select_hets = false

        Runtime runtime_params
    }

    String non_optional_pileup = select_first([pileup, "none"])
    String sample_name = basename(basename(basename(non_optional_pileup, ".gz"), ".pileup"), ".likelihoods")
    String output_file = sample_name + ".allelic_counts.tsv"

    command <<<
        set -euxo pipefail

        # Set default value
        printf "0.05" > error_probability.txt

        if [ "~{defined(pileup)}" == "true" ] ; then
            # HEADER
            cat '~{ref_dict}' > '~{output_file}'
            printf "@RG\tID:GATKCopyNumber\tSM:~{bam_name}\n" >> '~{output_file}'
            printf "CONTIG\tPOSITION\tREF_COUNT\tALT_COUNT\tREF_NUCLEOTIDE\tALT_NUCLEOTIDE\n" >> '~{output_file}'

            python <<EOF
import numpy as np
import pandas as pd

vcf_columns = ["contig", "position", "id", "ref", "alt", "qual", "filter", "info", "format", "genotype"]

pileup = pd.read_csv('~{pileup}', sep='\t', comment='#', low_memory=False).astype({"contig": str, "position": int, "ref_count": int, "alt_count": int, "other_alt_count": int, "allele_frequency": float})
gvcf = pd.read_csv('~{gvcf}', sep='\t', comment='#', header=None, low_memory=False, names=vcf_columns).astype({"contig": str, "position": int, "id": str, "ref": str, "alt": str, "info": str, "genotype": str})

if pileup.empty:
    print("No pileups found.")
else:
    df = pd.merge(pileup, gvcf, how='left', on=['contig', 'position'])
    df = df.loc[df[['ref_count', 'alt_count']].sum(axis=1) >= ~{min_read_depth}]
    if ~{if select_hets then "True" else "False"}:
        print("Selecting HETs")
        df = df.loc[df['genotype'] == '0/1']
    df[["contig", "position", "ref_count", "alt_count", "ref", "alt"]].to_csv('~{output_file}', sep='\t', index=False, header=False, mode='a')

    other_alt_counts = pileup["other_alt_count"].sum()
    total_counts = pileup[["ref_count", "alt_count", "other_alt_count"]].sum(axis=1).sum()
    error_probability = np.clip(1.5 * other_alt_counts / max(1, total_counts), a_min=0.001, a_max=0.05)
    np.savetxt("error_probability.txt", [error_probability])
EOF
        fi
    >>>

    output {
        File? allelic_counts = output_file
        Float error_probability = read_float("error_probability.txt")
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

task ModelSegmentsTask {
    input {
        File? segments
        Array[File]? denoised_copy_ratios
        Array[File]? allelic_counts
        File? normal_allelic_counts
        String prefix

        Float genotyping_base_error_rate = 0.05
        Float genotyping_homozygous_log_ratio_threshold = -10.0
        Int minimum_total_allele_count_normal = 30
        Array[Int] window_sizes = [8, 16, 32, 64, 128, 256]

        Runtime runtime_params
    }

    String output_dir = "."

    # todo: allow compressed input data

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{select_first([runtime_params.jar_override, "/root/gatk.jar"])}
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            ModelSegments \
            ~{"--segments '" + segments + "'"} \
            ~{true="--denoised-copy-ratios '" false="" defined(denoised_copy_ratios)}~{default="" sep="' --denoised-copy-ratios '" denoised_copy_ratios}~{true="'" false="" defined(denoised_copy_ratios)} \
            ~{true="--allelic-counts '" false="" defined(allelic_counts)}~{default="" sep="' --allelic-counts '" allelic_counts}~{true="'" false="" defined(allelic_counts)} \
            ~{"--normal-allelic-counts '" + normal_allelic_counts + "'"} \
            ~{"--output-prefix '" + prefix + "'"} \
            --genotyping-base-error-rate ~{genotyping_base_error_rate} \
            --genotyping-homozygous-log-ratio-threshold ~{genotyping_homozygous_log_ratio_threshold} \
            --minimum-total-allele-count-normal ~{minimum_total_allele_count_normal} \
            ~{sep=" " prefix("--window-size ", window_sizes)} \
            --output ~{output_dir}
    >>>

    output {
        File? multi_sample_segments = output_dir + "/" + prefix + ".interval_list"
        File? hets = output_dir + "/" + prefix + ".hets.tsv"
        File? af_model_begin_parameters = output_dir + "/" + prefix + ".modelBegin.af.param"
        File? cr_model_begin_parameters = output_dir + "/" + prefix + ".modelBegin.cr.param"
        File? af_model_final_parameters = output_dir + "/" + prefix + ".modelFinal.af.param"
        File? cr_model_final_parameters = output_dir + "/" + prefix + ".modelFinal.cr.param"
        File? seg_begin = output_dir + "/" + prefix + ".modelBegin.seg"
        File? seg_final = output_dir + "/" + prefix + ".modelFinal.seg"
        File? cr_seg = output_dir + "/" + prefix + ".cr.seg"
        File? igv_af = output_dir + "/" + prefix + ".af.igv.seg"
        File? igv_cr = output_dir + "/" + prefix + ".cr.igv.seg"
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
#        denoised_copy_ratios: {localization_optional: true}
#        allelic_counts: {localization_optional: true}
#        normal_allelic_counts: {localization_optional: true}
    }
}

task CallCopyRatioSegments {
    input {
        File copy_ratio_segments
        Float neutral_segment_copy_ratio_lower_bound = 0.9
        Float neutral_segment_copy_ratio_upper_bound = 1.1
        Float outlier_neutral_segment_copy_ratio_z_score_threshold = 2.0
        Float calling_copy_ratio_z_score_threshold = 2.0
        Runtime runtime_params
    }

    String called_segments = basename(copy_ratio_segments, ".seg") + ".called.seg"

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{select_first([runtime_params.jar_override, "/root/gatk.jar"])}
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            CallCopyRatioSegments \
            -I '~{copy_ratio_segments}' \
            -O '~{called_segments}' \
            --neutral-segment-copy-ratio-lower-bound ~{neutral_segment_copy_ratio_lower_bound} \
            --neutral-segment-copy-ratio-upper-bound ~{neutral_segment_copy_ratio_upper_bound} \
            --outlier-neutral-segment-copy-ratio-z-score-threshold ~{outlier_neutral_segment_copy_ratio_z_score_threshold} \
            --calling-copy-ratio-z-score-threshold ~{calling_copy_ratio_z_score_threshold}
    >>>

    output {
        File called_cr_seg = called_segments
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
#        copy_ratio_segments: {localization_optional: true}
    }
}

task MergeCallsWithModeledSegments {
    input {
        File copy_ratio_segments
        File called_copy_ratio_segments
        Runtime runtime_params
    }

    String dollar = "$"

    String output_file = basename(copy_ratio_segments, ".seg") + ".called.seg"

    command <<<
        # Extract headers (lines starting with "@") from the first file
        grep "^@" '~{copy_ratio_segments}' > '~{output_file}'

        # Merge the data (excluding headers) from both files
        paste <(grep -v "^@" '~{copy_ratio_segments}') <(grep -v "^@" '~{called_copy_ratio_segments}' | awk -F'\t' '{print $NF}') >> '~{output_file}'
    >>>

    output {
        File merged_segments = output_file
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

task PlotModeledSegments {
    input {
        File ref_dict

        String sample_name
        File segments
        File? denoised_copy_ratios
        File? het_allelic_counts

        Runtime runtime_params
    }

    String prefix = sample_name
    String output_dir = "."

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{select_first([runtime_params.jar_override, "/root/gatk.jar"])}
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            PlotModeledSegments \
            --sequence-dictionary '~{ref_dict}' \
            --segments '~{segments}' \
            ~{"--denoised-copy-ratios '" + denoised_copy_ratios + "'"} \
            ~{"--allelic-counts '" + het_allelic_counts + "'"} \
            --output-prefix '~{prefix}' \
            --output ~{output_dir}
    >>>

    output {
        File plot = output_dir + "/" + prefix + ".modeled.png"
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