version development

import "sample.wdl"
import "sequencing_run.wdl"
import "runtimes.wdl"
import "tasks.wdl"


struct Patient {
    String name
    Array[Sample] samples
    Array[Sample] tumor_samples
    Array[Sample] normal_samples

    Boolean has_tumor
    Boolean has_normal
    Sample? matched_normal_sample
}


workflow UpdatePatient {
    input {
        Patient patient
        String? name
        Array[Sample]? samples
        Array[Sample]? tumor_samples
        Array[Sample]? normal_samples
        Boolean? has_tumor
        Boolean? has_normal
        Sample? matched_normal_sample
    }

    Patient p = object {
        name: select_first([name, patient.name]),
        samples: select_first([samples, patient.samples]),
        tumor_samples: select_first([tumor_samples, patient.tumor_samples]),
        normal_samples: select_first([normal_samples, patient.normal_samples]),
        has_tumor: select_first([has_tumor, patient.has_tumor]),
        has_normal: select_first([has_normal, patient.has_normal]),
        # cannot use select_first for optional fields:
        matched_normal_sample: if defined(matched_normal_sample) then matched_normal_sample else patient.matched_normal_sample
    }

    output {
        Patient updated_patient = p
    }
}


workflow UpdateSamples {
    input {
        Patient patient
        Array[File]? read_counts
        Array[File]? denoised_copy_ratios
        Array[File]? standardized_copy_ratios
        Array[File]? snp_array_allelic_counts
        Array[File]? somatic_allelic_counts
        Array[File]? germline_allelic_counts
        Array[File]? contaminations
        Array[File]? af_segmentations
        Array[File]? copy_ratio_segmentations
    }

    # Iteratively update samples for each defined input

    if (defined(read_counts)) {
        scatter (pair in zip(patient.samples, select_first([read_counts, []]))) {
            if (pair.left.is_tumor) {
                call sample.UpdateSample as UpdateReadCountsTumor {
                    input:
                        sample = pair.left,
                        read_counts = pair.right,
                }
            }
            if (!pair.left.is_tumor) {
                call sample.UpdateSample as UpdateReadCountsNormal {
                    input:
                        sample = pair.left,
                        read_counts = pair.right,
                }
            }
        }
    }
    Array[Sample] tumor_samples_1 = select_all(select_first([UpdateReadCountsTumor.updated_sample, patient.tumor_samples]))
    Array[Sample] normal_samples_1 = select_all(select_first([UpdateReadCountsNormal.updated_sample, patient.normal_samples]))
    Array[Sample] samples_1 = flatten([tumor_samples_1, normal_samples_1])

    if (defined(denoised_copy_ratios)) {
        scatter (pair in zip(samples_1, select_first([denoised_copy_ratios, []]))) {
            if (pair.left.is_tumor) {
                call sample.UpdateSample as UpdateDenoisedCopyRatioTumor {
                    input:
                        sample = pair.left,
                        denoised_copy_ratios = pair.right,
                }
            }
            if (!pair.left.is_tumor) {
                call sample.UpdateSample as UpdateDenoisedCopyRatioNormal {
                    input:
                        sample = pair.left,
                        denoised_copy_ratios = pair.right,
                }
            }
        }
    }
    Array[Sample] tumor_samples_2 = select_all(select_first([UpdateDenoisedCopyRatioTumor.updated_sample, tumor_samples_1]))
    Array[Sample] normal_samples_2 = select_all(select_first([UpdateDenoisedCopyRatioNormal.updated_sample, normal_samples_1]))
    Array[Sample] samples_2 = flatten([tumor_samples_2, normal_samples_2])

    if (defined(standardized_copy_ratios)) {
        scatter (pair in zip(samples_2, select_first([standardized_copy_ratios, []]))) {
            if (pair.left.is_tumor) {
                call sample.UpdateSample as UpdateStandardizedCopyRatioTumor {
                    input:
                        sample = pair.left,
                        standardized_copy_ratios = pair.right,
                }
            }
            if (!pair.left.is_tumor) {
                call sample.UpdateSample as UpdateStandardizedCopyRatioNormal {
                    input:
                        sample = pair.left,
                        standardized_copy_ratios = pair.right,
                }
            }
        }
    }
    Array[Sample] tumor_samples_3 = select_all(select_first([UpdateStandardizedCopyRatioTumor.updated_sample, tumor_samples_2]))
    Array[Sample] normal_samples_3 = select_all(select_first([UpdateStandardizedCopyRatioNormal.updated_sample, normal_samples_2]))
    Array[Sample] samples_3 = flatten([tumor_samples_3, normal_samples_3])

    if (defined(snp_array_allelic_counts)) {
        scatter (pair in zip(samples_3, select_first([snp_array_allelic_counts, []]))) {
            if (pair.left.is_tumor) {
                call sample.UpdateSample as UpdateSnpArrayAllelicCountTumor {
                    input:
                        sample = pair.left,
                        snp_array_allelic_counts = pair.right,
                }
            }
            if (!pair.left.is_tumor) {
                call sample.UpdateSample as UpdateSnpArrayAllelicCountNormal {
                    input:
                        sample = pair.left,
                        snp_array_allelic_counts = pair.right,
                }
            }
        }
    }
    Array[Sample] tumor_samples_4 = select_all(select_first([UpdateSnpArrayAllelicCountTumor.updated_sample, tumor_samples_3]))
    Array[Sample] normal_samples_4 = select_all(select_first([UpdateSnpArrayAllelicCountNormal.updated_sample, normal_samples_3]))
    Array[Sample] samples_4 = flatten([tumor_samples_4, normal_samples_4])

    if (defined(somatic_allelic_counts)) {
        scatter (pair in zip(samples_4, select_first([somatic_allelic_counts, []]))) {
            if (pair.left.is_tumor) {
                call sample.UpdateSample as UpdateSomaticAllelicCountTumor {
                    input:
                        sample = pair.left,
                        somatic_allelic_counts = pair.right,
                }
            }
            if (!pair.left.is_tumor) {
                call sample.UpdateSample as UpdateSomaticAllelicCountNormal {
                    input:
                        sample = pair.left,
                        somatic_allelic_counts = pair.right,
                }
            }
        }
    }
    Array[Sample] tumor_samples_5 = select_all(select_first([UpdateSomaticAllelicCountTumor.updated_sample, tumor_samples_4]))
    Array[Sample] normal_samples_5 = select_all(select_first([UpdateSomaticAllelicCountNormal.updated_sample, normal_samples_4]))
    Array[Sample] samples_5 = flatten([tumor_samples_5, normal_samples_5])

    if (defined(germline_allelic_counts)) {
        scatter (pair in zip(samples_5, select_first([germline_allelic_counts, []]))) {
            if (pair.left.is_tumor) {
                call sample.UpdateSample as UpdateGermlineAllelicCountTumor {
                    input:
                        sample = pair.left,
                        germline_allelic_counts = pair.right,
                }
            }
            if (!pair.left.is_tumor) {
                call sample.UpdateSample as UpdateGermlineAllelicCountNormal {
                    input:
                        sample = pair.left,
                        germline_allelic_counts = pair.right,
                }
            }
        }
    }
    Array[Sample] tumor_samples_6 = select_all(select_first([UpdateGermlineAllelicCountTumor.updated_sample, tumor_samples_5]))
    Array[Sample] normal_samples_6 = select_all(select_first([UpdateGermlineAllelicCountNormal.updated_sample, normal_samples_5]))
    Array[Sample] samples_6 = flatten([tumor_samples_6, normal_samples_6])

    if (defined(contaminations)) {
        scatter (pair in zip(samples_6, select_first([contaminations, []]))) {
            if (pair.left.is_tumor) {
                call sample.UpdateSample as UpdateContaminationTumor {
                    input:
                        sample = pair.left,
                        contamination = pair.right,
                }
            }
            if (!pair.left.is_tumor) {
                call sample.UpdateSample as UpdateContaminationNormal {
                    input:
                        sample = pair.left,
                        contamination = pair.right,
                }
            }
        }
    }
    Array[Sample] tumor_samples_7 = select_all(select_first([UpdateContaminationTumor.updated_sample, tumor_samples_6]))
    Array[Sample] normal_samples_7 = select_all(select_first([UpdateContaminationNormal.updated_sample, normal_samples_6]))
    Array[Sample] samples_7 = flatten([tumor_samples_7, normal_samples_7])

    if (defined(af_segmentations)) {
        scatter (pair in zip(samples_7, select_first([af_segmentations, []]))) {
            if (pair.left.is_tumor) {
                call sample.UpdateSample as UpdateAfSegmentationTumor {
                    input:
                        sample = pair.left,
                        af_segmentation = pair.right,
                }
            }
            if (!pair.left.is_tumor) {
                call sample.UpdateSample as UpdateAfSegmentationNormal {
                    input:
                        sample = pair.left,
                        af_segmentation = pair.right,
                }
            }
        }
    }
    Array[Sample] tumor_samples_8 = select_all(select_first([UpdateAfSegmentationTumor.updated_sample, tumor_samples_7]))
    Array[Sample] normal_samples_8 = select_all(select_first([UpdateAfSegmentationNormal.updated_sample, normal_samples_7]))
    Array[Sample] samples_8 = flatten([tumor_samples_8, normal_samples_8])

    if (defined(copy_ratio_segmentations)) {
        scatter (pair in zip(samples_8, select_first([copy_ratio_segmentations, []]))) {
            if (pair.left.is_tumor) {
                call sample.UpdateSample as UpdateCopyRatioSegmentationTumor {
                    input:
                        sample = pair.left,
                        copy_ratio_segmentation = pair.right,
                }
            }
            if (!pair.left.is_tumor) {
                call sample.UpdateSample as UpdateCopyRatioSegmentationNormal {
                    input:
                        sample = pair.left,
                        copy_ratio_segmentation = pair.right,
                }
            }
        }
    }
    Array[Sample] tumor_samples_9 = select_all(select_first([UpdateCopyRatioSegmentationTumor.updated_sample, tumor_samples_8]))
    Array[Sample] normal_samples_9 = select_all(select_first([UpdateCopyRatioSegmentationNormal.updated_sample, normal_samples_8]))
    Array[Sample] samples_9 = flatten([tumor_samples_9, normal_samples_9])

    if (defined(patient.matched_normal_sample)) {
        Sample previous_matched_normal_sample = select_first([patient.matched_normal_sample])
        scatter (sample in samples_9) {
            if (sample.name == previous_matched_normal_sample.name) {
                Sample matched_normal_sample = sample
            }
        }
    }

    Patient p = object {
        name: patient.name,
        samples: samples_8,
        tumor_samples: tumor_samples_8,
        normal_samples: normal_samples_8,
        has_tumor: patient.has_tumor,
        has_normal: patient.has_normal,
        matched_normal_sample: if (patient.has_normal) then select_first(select_all(matched_normal_sample)) else patient.matched_normal_sample,
    }

    output {
        Patient updated_patient = p
    }
}


workflow DefinePatient {
    input {
        String individual_id

        Array[File]+ tumor_bams
        Array[File]+ tumor_bais
        Array[File]+ tumor_target_intervals
        Array[File]? tumor_annotated_target_intervals
        Array[File]? tumor_cnv_panel_of_normals
        Array[String]? tumor_sample_names

        Array[File]? normal_bams
        Array[File]? normal_bais
        Array[File]? normal_target_intervals
        Array[File]? normal_annotated_target_intervals
        Array[File]? normal_cnv_panel_of_normals
        Array[String]? normal_sample_names

        RuntimeCollection runtime_collection
    }

    Array[File] non_optional_normal_bams = select_first([normal_bams, []])
    Array[File] non_optional_normal_bais = select_first([normal_bais, []])
    Array[File] non_optional_normal_target_intervals = select_first([normal_target_intervals, []])

    Boolean has_tumor = length(tumor_bams) > 0
    Boolean has_normal = length(non_optional_normal_bams) > 0

    # We first define SequencingRuns for each tumor and normal bam, and then group them by sample name into Samples.

    scatter (s in transpose([tumor_bams, tumor_bais, tumor_target_intervals])) {
        call sequencing_run.DefineSequencingRun as DefineTumorSequencingRun {
            input:
                bam = s[0],
                bai = s[1],
                target_intervals = s[2],
                runtime_collection = runtime_collection
        }
        String tumor_bam_names = DefineTumorSequencingRun.sequencing_run.name
    }
    Array[SequencingRun] tumors_1 = DefineTumorSequencingRun.sequencing_run
    if (defined(tumor_annotated_target_intervals)) {
        scatter (pair in zip(tumors_1, select_first([tumor_annotated_target_intervals, []]))) {
            call sequencing_run.UpdateSequencingRun as UpdateAnnotatedTargetIntervalsTumorSeq {
                input:
                    sequencing_run = pair.left,
                    annotated_target_intervals = pair.right,
            }
        }
    }
    Array[SequencingRun] tumors_2 = select_first([UpdateAnnotatedTargetIntervalsTumorSeq.updated_sequencing_run, tumors_1])
    if (defined(tumor_cnv_panel_of_normals)) {
        scatter (pair in zip(tumors_2, select_first([tumor_cnv_panel_of_normals, []]))) {
            call sequencing_run.UpdateSequencingRun as UpdateCnvPanelOfNormalsTumorSeq {
                input:
                    sequencing_run = pair.left,
                    cnv_panel_of_normals = pair.right,
            }
        }
    }
    Array[SequencingRun] tumors_3 = select_first([UpdateCnvPanelOfNormalsTumorSeq.updated_sequencing_run, tumors_2])

    if (has_normal) {
        scatter (s in transpose([non_optional_normal_bams, non_optional_normal_bais, non_optional_normal_target_intervals])) {
            call sequencing_run.DefineSequencingRun as DefineNormalSequencingRun {
                input:
                    bam = s[0],
                    bai = s[1],
                    target_intervals = s[2],
                    runtime_collection = runtime_collection
            }
            String normal_bam_names = DefineNormalSequencingRun.sequencing_run.name
        }
    }
    Array[SequencingRun] normals_1 = select_first([DefineNormalSequencingRun.sequencing_run, []])
    if (defined(normal_annotated_target_intervals)) {
        scatter (pair in zip(normals_1, select_first([normal_annotated_target_intervals, []]))) {
            call sequencing_run.UpdateSequencingRun as UpdateAnnotatedTargetIntervalsNormalSeq {
                input:
                    sequencing_run = pair.left,
                    annotated_target_intervals = pair.right,
            }
        }
    }
    Array[SequencingRun] normals_2 = select_first([UpdateAnnotatedTargetIntervalsNormalSeq.updated_sequencing_run, normals_1])
    if (defined(normal_cnv_panel_of_normals)) {
        scatter (pair in zip(normals_2, select_first([normal_cnv_panel_of_normals, []]))) {
            call sequencing_run.UpdateSequencingRun as UpdateCnvPanelOfNormalsNormalSeq {
                input:
                    sequencing_run = pair.left,
                    cnv_panel_of_normals = pair.right,
            }
        }
    }
    Array[SequencingRun] normals_3 = select_first([UpdateCnvPanelOfNormalsNormalSeq.updated_sequencing_run, normals_2])

    # GroupBy sample name:
    # We assume that tumor_sample_names and tumor_bam_names share the same uniqueness,
    # that is if the supplied sample name is the same for two input bams, then the
    # bam names should also be the same, and vice versa.

    Array[String] tumor_names = select_first([tumor_sample_names, tumor_bam_names])
    scatter (pair in as_pairs(collect_by_key(zip(tumor_names, tumors_3)))) {
        Sample tumor_samples = object {
            name: pair.left,
            bam_name: pair.right[0].name,
            sequencing_runs: pair.right,
            is_tumor: true,
        }
    }
    Array[String] normal_names = select_first([normal_sample_names, normal_bam_names])
    scatter (pair in as_pairs(collect_by_key(zip(normal_names, normals_3)))) {
        Sample normal_samples = object {
            name: pair.left,
            bam_name: pair.right[0].name,
            sequencing_runs: pair.right,
            is_tumor: false,
        }
    }

    if (has_normal) {
        # We select the first normal sample to be the matched normal.
        # todo: select normal with greatest sequencing depth
        Sample best_matched_normal_sample = select_first(normal_samples)
    }

    Patient p = object {
        name: individual_id,
        samples: flatten([tumor_samples, normal_samples]),
        tumor_samples: tumor_samples,
        normal_samples: normal_samples,
        has_tumor: has_tumor,
        has_normal: has_normal,
        matched_normal_sample: best_matched_normal_sample,
    }

    output {
        Patient patient = p
    }
}