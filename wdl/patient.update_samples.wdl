version development

import "sample.wdl" as s
import "patient.wdl" as p


workflow UpdateSamples {
    input {
        Patient patient
        Array[File]? read_counts
        Array[File]? denoised_copy_ratios
        Array[File]? standardized_copy_ratios
        Array[File]? snp_array_pileups
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
                call s.UpdateSample as UpdateReadCountsTumor {
                    input:
                        sample = pair.left,
                        read_counts = pair.right,
                }
            }
            if (!pair.left.is_tumor) {
                call s.UpdateSample as UpdateReadCountsNormal {
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
                call s.UpdateSample as UpdateDenoisedCopyRatioTumor {
                    input:
                        sample = pair.left,
                        denoised_copy_ratios = pair.right,
                }
            }
            if (!pair.left.is_tumor) {
                call s.UpdateSample as UpdateDenoisedCopyRatioNormal {
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
                call s.UpdateSample as UpdateStandardizedCopyRatioTumor {
                    input:
                        sample = pair.left,
                        standardized_copy_ratios = pair.right,
                }
            }
            if (!pair.left.is_tumor) {
                call s.UpdateSample as UpdateStandardizedCopyRatioNormal {
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

    if (defined(snp_array_pileups)) {
        scatter (pair in zip(samples_3, select_first([snp_array_pileups, []]))) {
            if (pair.left.is_tumor) {
                call s.UpdateSample as UpdateSnpArrayPileupTumor {
                    input:
                        sample = pair.left,
                        snp_array_pileups = pair.right,
                }
            }
            if (!pair.left.is_tumor) {
                call s.UpdateSample as UpdateSnpArrayPileupNormal {
                    input:
                        sample = pair.left,
                        snp_array_pileups = pair.right,
                }
            }
        }
    }
    Array[Sample] tumor_samples_4 = select_all(select_first([UpdateSnpArrayPileupTumor.updated_sample, tumor_samples_3]))
    Array[Sample] normal_samples_4 = select_all(select_first([UpdateSnpArrayPileupNormal.updated_sample, normal_samples_3]))
    Array[Sample] samples_4 = flatten([tumor_samples_4, normal_samples_4])

    if (defined(snp_array_allelic_counts)) {
        scatter (pair in zip(samples_4, select_first([snp_array_allelic_counts, []]))) {
            if (pair.left.is_tumor) {
                call s.UpdateSample as UpdateSnpArrayAllelicCountTumor {
                    input:
                        sample = pair.left,
                        snp_array_allelic_counts = pair.right,
                }
            }
            if (!pair.left.is_tumor) {
                call s.UpdateSample as UpdateSnpArrayAllelicCountNormal {
                    input:
                        sample = pair.left,
                        snp_array_allelic_counts = pair.right,
                }
            }
        }
    }
    Array[Sample] tumor_samples_5 = select_all(select_first([UpdateSnpArrayAllelicCountTumor.updated_sample, tumor_samples_4]))
    Array[Sample] normal_samples_5 = select_all(select_first([UpdateSnpArrayAllelicCountNormal.updated_sample, normal_samples_4]))
    Array[Sample] samples_5 = flatten([tumor_samples_5, normal_samples_5])

    if (defined(somatic_allelic_counts)) {
        scatter (pair in zip(samples_5, select_first([somatic_allelic_counts, []]))) {
            if (pair.left.is_tumor) {
                call s.UpdateSample as UpdateSomaticAllelicCountTumor {
                    input:
                        sample = pair.left,
                        somatic_allelic_counts = pair.right,
                }
            }
            if (!pair.left.is_tumor) {
                call s.UpdateSample as UpdateSomaticAllelicCountNormal {
                    input:
                        sample = pair.left,
                        somatic_allelic_counts = pair.right,
                }
            }
        }
    }
    Array[Sample] tumor_samples_6 = select_all(select_first([UpdateSomaticAllelicCountTumor.updated_sample, tumor_samples_5]))
    Array[Sample] normal_samples_6 = select_all(select_first([UpdateSomaticAllelicCountNormal.updated_sample, normal_samples_5]))
    Array[Sample] samples_6 = flatten([tumor_samples_6, normal_samples_6])

    if (defined(germline_allelic_counts)) {
        scatter (pair in zip(samples_6, select_first([germline_allelic_counts, []]))) {
            if (pair.left.is_tumor) {
                call s.UpdateSample as UpdateGermlineAllelicCountTumor {
                    input:
                        sample = pair.left,
                        germline_allelic_counts = pair.right,
                }
            }
            if (!pair.left.is_tumor) {
                call s.UpdateSample as UpdateGermlineAllelicCountNormal {
                    input:
                        sample = pair.left,
                        germline_allelic_counts = pair.right,
                }
            }
        }
    }
    Array[Sample] tumor_samples_7 = select_all(select_first([UpdateGermlineAllelicCountTumor.updated_sample, tumor_samples_6]))
    Array[Sample] normal_samples_7 = select_all(select_first([UpdateGermlineAllelicCountNormal.updated_sample, normal_samples_6]))
    Array[Sample] samples_7 = flatten([tumor_samples_7, normal_samples_7])

    if (defined(contaminations)) {
        scatter (pair in zip(samples_7, select_first([contaminations, []]))) {
            if (pair.left.is_tumor) {
                call s.UpdateSample as UpdateContaminationTumor {
                    input:
                        sample = pair.left,
                        contamination = pair.right,
                }
            }
            if (!pair.left.is_tumor) {
                call s.UpdateSample as UpdateContaminationNormal {
                    input:
                        sample = pair.left,
                        contamination = pair.right,
                }
            }
        }
    }
    Array[Sample] tumor_samples_8 = select_all(select_first([UpdateContaminationTumor.updated_sample, tumor_samples_7]))
    Array[Sample] normal_samples_8 = select_all(select_first([UpdateContaminationNormal.updated_sample, normal_samples_7]))
    Array[Sample] samples_8 = flatten([tumor_samples_8, normal_samples_8])

    if (defined(af_segmentations)) {
        scatter (pair in zip(samples_8, select_first([af_segmentations, []]))) {
            if (pair.left.is_tumor) {
                call s.UpdateSample as UpdateAfSegmentationTumor {
                    input:
                        sample = pair.left,
                        af_segmentation = pair.right,
                }
            }
            if (!pair.left.is_tumor) {
                call s.UpdateSample as UpdateAfSegmentationNormal {
                    input:
                        sample = pair.left,
                        af_segmentation = pair.right,
                }
            }
        }
    }
    Array[Sample] tumor_samples_9 = select_all(select_first([UpdateAfSegmentationTumor.updated_sample, tumor_samples_8]))
    Array[Sample] normal_samples_9 = select_all(select_first([UpdateAfSegmentationNormal.updated_sample, normal_samples_8]))
    Array[Sample] samples_9 = flatten([tumor_samples_9, normal_samples_9])

    if (defined(copy_ratio_segmentations)) {
        scatter (pair in zip(samples_9, select_first([copy_ratio_segmentations, []]))) {
            if (pair.left.is_tumor) {
                call s.UpdateSample as UpdateCopyRatioSegmentationTumor {
                    input:
                        sample = pair.left,
                        copy_ratio_segmentation = pair.right,
                }
            }
            if (!pair.left.is_tumor) {
                call s.UpdateSample as UpdateCopyRatioSegmentationNormal {
                    input:
                        sample = pair.left,
                        copy_ratio_segmentation = pair.right,
                }
            }
        }
    }
    Array[Sample] tumor_samples_10 = select_all(select_first([UpdateCopyRatioSegmentationTumor.updated_sample, tumor_samples_9]))
    Array[Sample] normal_samples_10 = select_all(select_first([UpdateCopyRatioSegmentationNormal.updated_sample, normal_samples_9]))
    Array[Sample] samples_10 = flatten([tumor_samples_10, normal_samples_10])

    if (defined(patient.matched_normal_sample)) {
        Sample previous_matched_normal_sample = select_first([patient.matched_normal_sample])
        scatter (sample in samples_10) {
            if (sample.name == previous_matched_normal_sample.name) {
                Sample matched_normal_samples = sample
            }
        }
        Sample matched_normal_sample = select_all(matched_normal_samples)[0]
    }

    Patient pat = object {
        name: patient.name,
        samples: samples_10,
        tumor_samples: tumor_samples_10,
        normal_samples: normal_samples_10,
        has_tumor: patient.has_tumor,
        has_normal: patient.has_normal,
        matched_normal_sample: if defined(matched_normal_sample) then matched_normal_sample else patient.matched_normal_sample,
    }

    output {
        Patient updated_patient = pat
    }
}