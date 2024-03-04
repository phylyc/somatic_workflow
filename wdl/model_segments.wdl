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
        call PileupToAllelicCounts {
            input:
                ref_dict = args.files.ref_dict,
                pileup = sample.snp_array_pileups,
                bam_name = sample.bam_name,
                runtime_params = runtime_collection.pileup_to_allelic_counts
        }
    }
    Array[File] sample_allelic_counts = select_all(PileupToAllelicCounts.allelic_counts)

    # Save the allelic counts in each sample object (if they exist)

    if (length(sample_allelic_counts) > 1) {
        call p_update_s.UpdateSamples as AddAllelicCountsToSamples {
            input:
                patient = patient,
                snp_array_allelic_counts = sample_allelic_counts,
        }
    }
    Patient pat = select_first([AddAllelicCountsToSamples.updated_patient, patient])

    # Prepare ModelSegments input

    scatter (sample in pat.samples) {
        File? denoised_copy_ratios = sample.denoised_copy_ratios
        File? allelic_counts = sample.snp_array_allelic_counts
    }
    if (length(select_all(denoised_copy_ratios)) > 0) {
        Array[File] dcr = select_all(denoised_copy_ratios)
    }
    if (length(select_all(allelic_counts)) > 0) {
        Array[File] ac = select_all(allelic_counts)
    }
    if (defined(pat.matched_normal_sample)) {
        Sample matched_normal_sample = select_first([pat.matched_normal_sample])
        File? normal_allelic_counts = matched_normal_sample.snp_array_allelic_counts
    }

    # Now we can finally do the segmentation. First, we determine the patient-specific
    # segmentation, then we infer the sample-specific copy ratios and call amps/dels.

    if (length(pat.samples) > 1) {
        call ModelSegmentsTask as MultiSampleModelSegments {
            input:
                denoised_copy_ratios = dcr,
                allelic_counts = ac,
                normal_allelic_counts = normal_allelic_counts,
                prefix = patient.name + ".segmentation",
                runtime_params = runtime_collection.model_segments
        }
    }

    scatter (sample in pat.samples) {
        if (defined(sample.denoised_copy_ratios))  {
            Array[File] dcr_list = select_all([sample.denoised_copy_ratios])
        }
        if (defined(sample.snp_array_allelic_counts)) {
            Array[File] ac_list = select_all([sample.snp_array_allelic_counts])
        }
        call ModelSegmentsTask as SingleSampleInferCR {
            input:
                segments = MultiSampleModelSegments.multi_sample_segments,
                denoised_copy_ratios = dcr_list,
                allelic_counts = ac_list,
                normal_allelic_counts = normal_allelic_counts,
                prefix = sample.name,
                runtime_params = runtime_collection.model_segments
        }

        call CallCopyRatioSegments {
            input:
                copy_ratio_segments = select_first([SingleSampleInferCR.cr_seg]),
                runtime_params = runtime_collection.call_copy_ratio_segments
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

        call ModelSegmentsToACSConversion {
            input:
                seg_final = select_first([SingleSampleInferCR.seg_final]),
                af_model_parameters = select_first([SingleSampleInferCR.af_model_final_parameters]),
                runtime_params = runtime_collection.model_segments_to_acs_conversion
        }
    }

    call p_update_s.UpdateSamples as AddSegmentationResultsToSamples {
        input:
            patient = pat,
            copy_ratio_segmentations = select_all(SingleSampleInferCR.seg_final),
            af_model_parameters = select_all(SingleSampleInferCR.af_model_final_parameters),
            cr_model_parameters = select_all(SingleSampleInferCR.cr_model_final_parameters),
            called_copy_ratio_segmentations =  CallCopyRatioSegments.called_cr_seg,
            acs_copy_ratio_segmentations = ModelSegmentsToACSConversion.cnv_acs_conversion_seg,
            acs_copy_ratio_skews = ModelSegmentsToACSConversion.cnv_acs_conversion_skew
    }

    # todo: tag & filter germline events

    output {
        Patient updated_patient = AddSegmentationResultsToSamples.updated_patient

        Array[File] snp_array_allelic_counts = sample_allelic_counts

        File? modeled_segments = MultiSampleModelSegments.multi_sample_segments
        Array[File] hets = select_all(SingleSampleInferCR.hets)
        Array[File] af_model_begin_parameters = select_all(SingleSampleInferCR.af_model_begin_parameters)
        Array[File] cr_model_begin_parameters = select_all(SingleSampleInferCR.cr_model_begin_parameters)
        Array[File] af_model_final_parameters = select_all(SingleSampleInferCR.af_model_final_parameters)
        Array[File] cr_model_final_parameters = select_all(SingleSampleInferCR.cr_model_final_parameters)
        Array[File] igv_af = select_all(SingleSampleInferCR.igv_af)
        Array[File] igv_cr = select_all(SingleSampleInferCR.igv_cr)
        Array[File] seg_begin = select_all(SingleSampleInferCR.seg_begin)
        Array[File] seg_final = select_all(SingleSampleInferCR.seg_final)
        Array[File] cr_seg = select_all(SingleSampleInferCR.cr_seg)
        Array[File] called_cr_seg = CallCopyRatioSegments.called_cr_seg
        Array[File] cr_plots = PlotModeledSegments.plot
        Array[File] cr_converted_acs_segments = ModelSegmentsToACSConversion.cnv_acs_conversion_seg
        Array[File] cr_converted_acs_skews = ModelSegmentsToACSConversion.cnv_acs_conversion_skew
    }
}

task PileupToAllelicCounts {
    input {
        File ref_dict
        String bam_name
        File? pileup
        Runtime runtime_params
    }

    String non_optional_pileup = select_first([pileup, "none"])
    Boolean is_compressed = basename(non_optional_pileup, ".gz") != basename(non_optional_pileup)
    String uncompressed_pileup = basename(non_optional_pileup, ".gz")
    String sample_name = basename(basename(uncompressed_pileup, ".pileup"), ".likelihoods")
    String output_file = sample_name + ".allelic_counts.tsv"

    command <<<
        set -e
        if [ "~{defined(pileup)}" == "true" ] ; then
            if [ "~{is_compressed}" == "true" ] ; then
                gzip -cd '~{pileup}' > '~{uncompressed_pileup}'
            else
                cp '~{pileup}' '~{uncompressed_pileup}'
            fi

            # HEADER
            cat '~{ref_dict}' > '~{output_file}'
            printf "@RG\tID:GATKCopyNumber\tSM:~{bam_name}\n" >> '~{output_file}'
            printf "CONTIG\tPOSITION\tREF_COUNT\tALT_COUNT\tREF_NUCLEOTIDE\tALT_NUCLEOTIDE\n" >> '~{output_file}'
            # CONTENT
            # The ref and alt allele information is not available in the pileup file, so we set both to "N"
            tail -n +3 '~{uncompressed_pileup}' | awk -v OFS='\t' '{print $1, $2, $3, $4, "N", "N"}' >> '~{output_file}'

            # Calculate the error probability:
            # The error probability is calculated as a ratio of 'other_alt_count' to the
            # total counts, clipped within specified bounds.
#            ref_counts=$(awk '{sum+=$3} END {print sum}' '~{uncompressed_pileup}')
#            alt_counts=$(awk '{sum+=$4} END {print sum}' '~{uncompressed_pileup}')
#            other_alt_counts=$(awk '{sum+=$5} END {print sum}' '~{uncompressed_pileup}')
            # todo: finish this
        fi
    >>>

    output {
        File? allelic_counts = output_file
        Float error_probability = 0.05
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

        Float genotying_base_error_rate = 0.05
        Array[Int] window_sizes = [2, 4, 8, 16, 32, 64, 128, 256, 512]

        Runtime runtime_params
    }

    String output_dir = "."

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
            --genotyping-base-error-rate ~{genotying_base_error_rate} \
            --window-size [~{sep=", " window_sizes}] \
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

task ModelSegmentsToACSConversion {
    input {
        String script = "https://github.com/phylyc/genomics_workflows/raw/master/python/acs_conversion.py"

        File seg_final
        File af_model_parameters

        Int min_hets = 10
        Float maf90_threshold = 0.485

        Runtime runtime_params
    }

    String output_dir = "."

    String output_seg = output_dir + "/" + basename(basename(seg_final, ".seg"), ".modelFinal") + ".acs.seg"
    String output_skew = output_seg + ".skew"

    command <<<
        set -e
        wget -O acs_conversion.py ~{script}
        python acs_conversion.py \
            --output_dir ~{output_dir} \
            --seg ~{seg_final} \
            --af_parameters ~{af_model_parameters} \
            --min_hets ~{min_hets} \
            --maf90_threshold ~{maf90_threshold} \
            --verbose
    >>>

    output {
        File cnv_acs_conversion_seg = output_seg
        File cnv_acs_conversion_skew = output_skew
        Float cnv_acs_conversion_skew_float = read_float(output_skew)
        String cnv_acs_conversion_skew_string = read_string(output_skew)
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