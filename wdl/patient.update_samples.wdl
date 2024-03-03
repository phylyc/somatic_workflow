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
        Array[File]? called_copy_ratio_segmentations
    }

    # Split input arrays into tumor and normal arrays

    Int num_tumor_samples = length(patient.tumor_samples)
    Int num_normal_samples = length(patient.normal_samples)

    scatter (t in range(num_tumor_samples)) {
        File? t_read_counts = if defined(read_counts) then select_first([read_counts, []])[t] else None
        File? t_denoised_copy_ratios = if defined(denoised_copy_ratios) then select_first([denoised_copy_ratios, []])[t] else None
        File? t_standardized_copy_ratios = if defined(standardized_copy_ratios) then select_first([standardized_copy_ratios, []])[t] else None
        File? t_snp_array_pileups = if defined(snp_array_pileups) then select_first([snp_array_pileups, []])[t] else None
        File? t_snp_array_allelic_counts = if defined(snp_array_allelic_counts) then select_first([snp_array_allelic_counts, []])[t] else None
        File? t_somatic_allelic_counts = if defined(somatic_allelic_counts) then select_first([somatic_allelic_counts, []])[t] else None
        File? t_germline_allelic_counts = if defined(germline_allelic_counts) then select_first([germline_allelic_counts, []])[t] else None
        File? t_contaminations = if defined(contaminations) then select_first([contaminations, []])[t] else None
        File? t_af_segmentations = if defined(af_segmentations) then select_first([af_segmentations, []])[t] else None
        File? t_copy_ratio_segmentations = if defined(copy_ratio_segmentations) then select_first([copy_ratio_segmentations, []])[t] else None
        File? t_called_copy_ratio_segmentations = if defined(called_copy_ratio_segmentations) then select_first([called_copy_ratio_segmentations, []])[t] else None
    }
    Array[File] tumor_read_counts = select_all(t_read_counts)
    Array[File] tumor_denoised_copy_ratios = select_all(t_denoised_copy_ratios)
    Array[File] tumor_standardized_copy_ratios = select_all(t_standardized_copy_ratios)
    Array[File] tumor_snp_array_pileups = select_all(t_snp_array_pileups)
    Array[File] tumor_snp_array_allelic_counts = select_all(t_snp_array_allelic_counts)
    Array[File] tumor_somatic_allelic_counts = select_all(t_somatic_allelic_counts)
    Array[File] tumor_germline_allelic_counts = select_all(t_germline_allelic_counts)
    Array[File] tumor_contaminations = select_all(t_contaminations)
    Array[File] tumor_af_segmentations = select_all(t_af_segmentations)
    Array[File] tumor_copy_ratio_segmentations = select_all(t_copy_ratio_segmentations)
    Array[File] tumor_called_copy_ratio_segmentations = select_all(t_called_copy_ratio_segmentations)

    if (patient.has_normal) {
        scatter (n in range(num_normal_samples)) {
            Int m = num_tumor_samples + n
            File? n_read_counts = if defined(read_counts) then select_first([read_counts, []])[m] else None
            File? n_denoised_copy_ratios = if defined(denoised_copy_ratios) then select_first([denoised_copy_ratios, []])[m] else None
            File? n_standardized_copy_ratios = if defined(standardized_copy_ratios) then select_first([standardized_copy_ratios, []])[m] else None
            File? n_snp_array_pileups = if defined(snp_array_pileups) then select_first([snp_array_pileups, []])[m] else None
            File? n_snp_array_allelic_counts = if defined(snp_array_allelic_counts) then select_first([snp_array_allelic_counts, []])[m] else None
            File? n_somatic_allelic_counts = if defined(somatic_allelic_counts) then select_first([somatic_allelic_counts, []])[m] else None
            File? n_germline_allelic_counts = if defined(germline_allelic_counts) then select_first([germline_allelic_counts, []])[m] else None
            File? n_contaminations = if defined(contaminations) then select_first([contaminations, []])[m] else None
            File? n_af_segmentations = if defined(af_segmentations) then select_first([af_segmentations, []])[m] else None
            File? n_copy_ratio_segmentations = if defined(copy_ratio_segmentations) then select_first([copy_ratio_segmentations, []])[m] else None
            File? n_called_copy_ratio_segmentations = if defined(called_copy_ratio_segmentations) then select_first([called_copy_ratio_segmentations, []])[m] else None
        }
    }
    Array[File] normal_read_counts = select_all(select_first([n_read_counts, []]))
    Array[File] normal_denoised_copy_ratios = select_all(select_first([n_denoised_copy_ratios, []]))
    Array[File] normal_standardized_copy_ratios = select_all(select_first([n_standardized_copy_ratios, []]))
    Array[File] normal_snp_array_pileups = select_all(select_first([n_snp_array_pileups, []]))
    Array[File] normal_snp_array_allelic_counts = select_all(select_first([n_snp_array_allelic_counts, []]))
    Array[File] normal_somatic_allelic_counts = select_all(select_first([n_somatic_allelic_counts, []]))
    Array[File] normal_germline_allelic_counts = select_all(select_first([n_germline_allelic_counts, []]))
    Array[File] normal_contaminations = select_all(select_first([n_contaminations, []]))
    Array[File] normal_af_segmentations = select_all(select_first([n_af_segmentations, []]))
    Array[File] normal_copy_ratio_segmentations = select_all(select_first([n_copy_ratio_segmentations, []]))
    Array[File] normal_called_copy_ratio_segmentations = select_all(select_first([n_called_copy_ratio_segmentations, []]))

    # Update tumor samples:

    if (length(tumor_read_counts) > 0) {
        scatter (pair in zip(patient.tumor_samples, tumor_read_counts)) {
            call s.UpdateSample as UpdateReadCountsTumor {
                input:
                    sample = pair.left,
                    read_counts = pair.right,
            }
        }
    }
    Array[Sample] tumor_samples_1 = select_first([UpdateReadCountsTumor.updated_sample, patient.tumor_samples])

    if (length(tumor_denoised_copy_ratios) > 0) {
        scatter (pair in zip(tumor_samples_1, tumor_denoised_copy_ratios)) {
            call s.UpdateSample as UpdateDenoisedCopyRatioTumor {
                input:
                    sample = pair.left,
                    denoised_copy_ratios = pair.right,
            }
        }
    }
    Array[Sample] tumor_samples_2 = select_first([UpdateDenoisedCopyRatioTumor.updated_sample, tumor_samples_1])

    if (length(tumor_standardized_copy_ratios) > 0) {
        scatter (pair in zip(tumor_samples_2, tumor_standardized_copy_ratios)) {
            call s.UpdateSample as UpdateStandardizedCopyRatioTumor {
                input:
                    sample = pair.left,
                    standardized_copy_ratios = pair.right,
            }
        }
    }
    Array[Sample] tumor_samples_3 = select_first([UpdateStandardizedCopyRatioTumor.updated_sample, tumor_samples_2])

    if (length(tumor_snp_array_pileups) > 0) {
        scatter (pair in zip(tumor_samples_3, tumor_snp_array_pileups)) {
            call s.UpdateSample as UpdateSnpArrayPileupTumor {
                input:
                    sample = pair.left,
                    snp_array_pileups = pair.right,
            }
        }
    }
    Array[Sample] tumor_samples_4 = select_first([UpdateSnpArrayPileupTumor.updated_sample, tumor_samples_3])

    if (length(tumor_snp_array_allelic_counts) > 0) {
        scatter (pair in zip(tumor_samples_4, tumor_snp_array_allelic_counts)) {
            call s.UpdateSample as UpdateSnpArrayAllelicCountTumor {
                input:
                    sample = pair.left,
                    snp_array_allelic_counts = pair.right,
            }
        }
    }
    Array[Sample] tumor_samples_5 = select_first([UpdateSnpArrayAllelicCountTumor.updated_sample, tumor_samples_4])

    if (length(tumor_somatic_allelic_counts) > 0) {
        scatter (pair in zip(tumor_samples_5, tumor_somatic_allelic_counts)) {
            call s.UpdateSample as UpdateSomaticAllelicCountTumor {
                input:
                    sample = pair.left,
                    somatic_allelic_counts = pair.right,
            }
        }
    }
    Array[Sample] tumor_samples_6 = select_first([UpdateSomaticAllelicCountTumor.updated_sample, tumor_samples_5])

    if (length(tumor_germline_allelic_counts) > 0) {
        scatter (pair in zip(tumor_samples_6, tumor_germline_allelic_counts)) {
            call s.UpdateSample as UpdateGermlineAllelicCountTumor {
                input:
                    sample = pair.left,
                    germline_allelic_counts = pair.right,
            }
        }
    }
    Array[Sample] tumor_samples_7 = select_first([UpdateGermlineAllelicCountTumor.updated_sample, tumor_samples_6])

    if (length(tumor_contaminations) > 0) {
        scatter (pair in zip(tumor_samples_7, tumor_contaminations)) {
            call s.UpdateSample as UpdateContaminationTumor {
                input:
                    sample = pair.left,
                    contamination = pair.right,
            }
        }
    }
    Array[Sample] tumor_samples_8 = select_first([UpdateContaminationTumor.updated_sample, tumor_samples_7])

    if (length(tumor_af_segmentations) > 0) {
        scatter (pair in zip(tumor_samples_8, tumor_af_segmentations)) {
            call s.UpdateSample as UpdateAfSegmentationTumor {
                input:
                    sample = pair.left,
                    af_segmentation = pair.right,
            }
        }
    }
    Array[Sample] tumor_samples_9 = select_first([UpdateAfSegmentationTumor.updated_sample, tumor_samples_8])

    if (length(tumor_copy_ratio_segmentations) > 0) {
        scatter (pair in zip(tumor_samples_9, tumor_copy_ratio_segmentations)) {
            call s.UpdateSample as UpdateCopyRatioSegmentationTumor {
                input:
                    sample = pair.left,
                    copy_ratio_segmentation = pair.right,
            }
        }
    }
    Array[Sample] tumor_samples_10 = select_first([UpdateCopyRatioSegmentationTumor.updated_sample, tumor_samples_9])

    if (length(tumor_called_copy_ratio_segmentations) > 0) {
        scatter (pair in zip(tumor_samples_10, tumor_called_copy_ratio_segmentations)) {
            call s.UpdateSample as UpdateCalledCopyRatioSegmentationTumor {
                input:
                    sample = pair.left,
                    called_copy_ratio_segmentation = pair.right,
            }
        }
    }
    Array[Sample] tumor_samples_11 = select_first([UpdateCalledCopyRatioSegmentationTumor.updated_sample, tumor_samples_10])

    # Update normal samples:

    if (length(normal_read_counts) > 0) {
        scatter (pair in zip(patient.normal_samples, normal_read_counts)) {
            call s.UpdateSample as UpdateReadCountsNormal {
                input:
                    sample = pair.left,
                    read_counts = pair.right,
            }
        }
    }
    Array[Sample] normal_samples_1 = select_first([UpdateReadCountsNormal.updated_sample, patient.normal_samples])

    if (length(normal_denoised_copy_ratios) > 0) {
        scatter (pair in zip(normal_samples_1, normal_denoised_copy_ratios)) {
            call s.UpdateSample as UpdateDenoisedCopyRatioNormal {
                input:
                    sample = pair.left,
                    denoised_copy_ratios = pair.right,
            }
        }
    }
    Array[Sample] normal_samples_2 = select_first([UpdateDenoisedCopyRatioNormal.updated_sample, normal_samples_1])

    if (length(normal_standardized_copy_ratios) > 0) {
        scatter (pair in zip(normal_samples_2, normal_standardized_copy_ratios)) {
            call s.UpdateSample as UpdateStandardizedCopyRatioNormal {
                input:
                    sample = pair.left,
                    standardized_copy_ratios = pair.right,
            }
        }
    }
    Array[Sample] normal_samples_3 = select_first([UpdateStandardizedCopyRatioNormal.updated_sample, normal_samples_2])

    if (length(normal_snp_array_pileups) > 0) {
        scatter (pair in zip(normal_samples_3, normal_snp_array_pileups)) {
            call s.UpdateSample as UpdateSnpArrayPileupNormal {
                input:
                    sample = pair.left,
                    snp_array_pileups = pair.right,
            }
        }
    }
    Array[Sample] normal_samples_4 = select_first([UpdateSnpArrayPileupNormal.updated_sample, normal_samples_3])

    if (length(normal_snp_array_allelic_counts) > 0) {
        scatter (pair in zip(normal_samples_4, normal_snp_array_allelic_counts)) {
            call s.UpdateSample as UpdateSnpArrayAllelicCountNormal {
                input:
                    sample = pair.left,
                    snp_array_allelic_counts = pair.right,
            }
        }
    }
    Array[Sample] normal_samples_5 = select_first([UpdateSnpArrayAllelicCountNormal.updated_sample, normal_samples_4])

    if (length(normal_somatic_allelic_counts) > 0) {
        scatter (pair in zip(normal_samples_5, normal_somatic_allelic_counts)) {
            call s.UpdateSample as UpdateSomaticAllelicCountNormal {
                input:
                    sample = pair.left,
                    somatic_allelic_counts = pair.right,
            }
        }
    }
    Array[Sample] normal_samples_6 = select_first([UpdateSomaticAllelicCountNormal.updated_sample, normal_samples_5])

    if (length(normal_germline_allelic_counts) > 0) {
        scatter (pair in zip(normal_samples_6, normal_germline_allelic_counts)) {
            call s.UpdateSample as UpdateGermlineAllelicCountNormal {
                input:
                    sample = pair.left,
                    germline_allelic_counts = pair.right,
            }
        }
    }
    Array[Sample] normal_samples_7 = select_first([UpdateGermlineAllelicCountNormal.updated_sample, normal_samples_6])

    if (length(normal_contaminations) > 0) {
        scatter (pair in zip(normal_samples_7, normal_contaminations)) {
            call s.UpdateSample as UpdateContaminationNormal {
                input:
                    sample = pair.left,
                    contamination = pair.right,
            }
        }
    }
    Array[Sample] normal_samples_8 = select_first([UpdateContaminationNormal.updated_sample, normal_samples_7])

    if (length(normal_af_segmentations) > 0) {
        scatter (pair in zip(normal_samples_8, normal_af_segmentations)) {
            call s.UpdateSample as UpdateAfSegmentationNormal {
                input:
                    sample = pair.left,
                    af_segmentation = pair.right,
            }
        }
    }
    Array[Sample] normal_samples_9 = select_first([UpdateAfSegmentationNormal.updated_sample, normal_samples_8])

    if (length(normal_copy_ratio_segmentations) > 0) {
        scatter (pair in zip(normal_samples_9, normal_copy_ratio_segmentations)) {
            call s.UpdateSample as UpdateCopyRatioSegmentationNormal {
                input:
                    sample = pair.left,
                    copy_ratio_segmentation = pair.right,
            }
        }
    }
    Array[Sample] normal_samples_10 = select_first([UpdateCopyRatioSegmentationNormal.updated_sample, normal_samples_9])

    if (length(normal_called_copy_ratio_segmentations) > 0) {
        scatter (pair in zip(normal_samples_10, normal_called_copy_ratio_segmentations)) {
            call s.UpdateSample as UpdateCalledCopyRatioSegmentationNormal {
                input:
                    sample = pair.left,
                    called_copy_ratio_segmentation = pair.right,
            }
        }
    }
    Array[Sample] normal_samples_11 = select_first([UpdateCalledCopyRatioSegmentationNormal.updated_sample, normal_samples_10])

    # Update matched normal sample:

    if (defined(patient.matched_normal_sample)) {
        Sample previous_matched_normal_sample = select_first([patient.matched_normal_sample])
        scatter (sample in normal_samples_11) {
            if (sample.name == previous_matched_normal_sample.name) {
                Sample matched_normal_samples = sample
            }
        }
        Sample matched_normal_sample = select_all(matched_normal_samples)[0]
    }

    Patient pat = object {
        name: patient.name,
        samples: flatten([tumor_samples_11, normal_samples_11]),
        tumor_samples: tumor_samples_11,
        normal_samples: normal_samples_11,
        has_tumor: patient.has_tumor,
        has_normal: patient.has_normal,
        matched_normal_sample: if defined(matched_normal_sample) then matched_normal_sample else patient.matched_normal_sample,
    }

    output {
        Patient updated_patient = pat
    }
}