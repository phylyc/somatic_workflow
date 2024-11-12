version development

import "sample.wdl" as s
import "patient.wdl" as p


workflow UpdateSamples {
    input {
        Patient patient
        Array[File]? read_counts
        Array[File]? covered_regions
        Array[File]? denoised_copy_ratios
        Array[File]? standardized_copy_ratios
        Array[File]? snppanel_pileups
        Array[File]? snppanel_allelic_counts
        Array[Float]? genotype_error_probabilities
        Array[File]? somatic_allelic_counts
        Array[File]? germline_allelic_counts
        Array[File]? contaminations
        Array[File]? af_segmentations
        Array[File]? af_model_parameters
        Array[File]? cr_model_parameters
        Array[File]? called_copy_ratio_segmentations
        Array[File]? acs_copy_ratio_segmentations
        Array[Float]? acs_copy_ratio_skews
        Array[File]? annotated_variants
    }

    # Update samples:

    if (defined(read_counts)) {
        scatter (pair in zip(patient.samples, select_first([read_counts, []]))) {
            call s.UpdateSample as UpdateReadCounts {
                input:
                    sample = pair.left,
                    read_counts = pair.right,
            }
        }
    }
    Array[Sample] samples_1 = select_first([UpdateReadCounts.updated_sample, patient.samples])

    if (defined(covered_regions)) {
        scatter (pair in zip(samples_1, select_first([covered_regions, []]))) {
            call s.UpdateSample as UpdateCoveredRegions {
                input:
                    sample = pair.left,
                    covered_regions = pair.right,
            }
        }
    }
    Array[Sample] samples_2 = select_first([UpdateCoveredRegions.updated_sample, samples_1])

    if (defined(denoised_copy_ratios)) {
        scatter (pair in zip(samples_2, select_first([denoised_copy_ratios, []]))) {
            call s.UpdateSample as UpdateDenoisedCopyRatio {
                input:
                    sample = pair.left,
                    denoised_copy_ratios = pair.right,
            }
        }
    }
    Array[Sample] samples_3 = select_first([UpdateDenoisedCopyRatio.updated_sample, samples_2])

    if (defined(standardized_copy_ratios)) {
        scatter (pair in zip(samples_3, select_first([standardized_copy_ratios, []]))) {
            call s.UpdateSample as UpdateStandardizedCopyRatio {
                input:
                    sample = pair.left,
                    standardized_copy_ratios = pair.right,
            }
        }
    }
    Array[Sample] samples_4 = select_first([UpdateStandardizedCopyRatio.updated_sample, samples_3])

    if (defined(snppanel_pileups)) {
        scatter (pair in zip(samples_4, select_first([snppanel_pileups, []]))) {
            call s.UpdateSample as UpdateSnpPanelPileup {
                input:
                    sample = pair.left,
                    snppanel_pileups = pair.right,
            }
        }
    }
    Array[Sample] samples_5 = select_first([UpdateSnpPanelPileup.updated_sample, samples_4])

    if (defined(snppanel_allelic_counts)) {
        scatter (pair in zip(samples_5, select_first([snppanel_allelic_counts, []]))) {
            call s.UpdateSample as UpdateSnpPanelAllelicCounts {
                input:
                    sample = pair.left,
                    snppanel_allelic_counts = pair.right,
            }
        }
    }
    Array[Sample] samples_6 = select_first([UpdateSnpPanelAllelicCounts.updated_sample, samples_5])

    if (defined(genotype_error_probabilities)) {
        scatter (pair in zip(samples_6, select_first([genotype_error_probabilities, []]))) {
            call s.UpdateSample as UpdateGenotypeErrorProbabilities {
                input:
                    sample = pair.left,
                    genotype_error_probabilities = pair.right,
            }
        }
    }
    Array[Sample] samples_7 = select_first([UpdateGenotypeErrorProbabilities.updated_sample, samples_6])

    if (defined(somatic_allelic_counts)) {
        scatter (pair in zip(samples_7, select_first([somatic_allelic_counts, []]))) {
            call s.UpdateSample as UpdateSomaticAllelicCounts {
                input:
                    sample = pair.left,
                    somatic_allelic_counts = pair.right,
            }
        }
    }
    Array[Sample] samples_8 = select_first([UpdateSomaticAllelicCounts.updated_sample, samples_7])

    if (defined(germline_allelic_counts)) {
        scatter (pair in zip(samples_8, select_first([germline_allelic_counts, []]))) {
            call s.UpdateSample as UpdateGermlineAllelicCounts {
                input:
                    sample = pair.left,
                    germline_allelic_counts = pair.right,
            }
        }
    }
    Array[Sample] samples_9 = select_first([UpdateGermlineAllelicCounts.updated_sample, samples_8])

    if (defined(contaminations)) {
        scatter (pair in zip(samples_9, select_first([contaminations, []]))) {
            call s.UpdateSample as UpdateContamination {
                input:
                    sample = pair.left,
                    contamination = pair.right,
            }
        }
    }
    Array[Sample] samples_10 = select_first([UpdateContamination.updated_sample, samples_9])

    if (defined(af_segmentations)) {
        scatter (pair in zip(samples_10, select_first([af_segmentations, []]))) {
            call s.UpdateSample as UpdateAfSegmentation {
                input:
                    sample = pair.left,
                    af_segmentation = pair.right,
            }
        }
    }
    Array[Sample] samples_11 = select_first([UpdateAfSegmentation.updated_sample, samples_10])

    if (defined(af_model_parameters)) {
        scatter (pair in zip(samples_11, select_first([af_model_parameters, []]))) {
            call s.UpdateSample as UpdateAfModelParameters {
                input:
                    sample = pair.left,
                    af_model_parameters = pair.right,
            }
        }
    }
    Array[Sample] samples_12 = select_first([UpdateAfModelParameters.updated_sample, samples_11])

    if (defined(cr_model_parameters)) {
        scatter (pair in zip(samples_12, select_first([cr_model_parameters, []]))) {
            call s.UpdateSample as UpdateCrModelParameters {
                input:
                    sample = pair.left,
                    cr_model_parameters = pair.right,
            }
        }
    }
    Array[Sample] samples_13 = select_first([UpdateCrModelParameters.updated_sample, samples_12])

    if (defined(called_copy_ratio_segmentations)) {
        scatter (pair in zip(samples_13, select_first([called_copy_ratio_segmentations, []]))) {
            call s.UpdateSample as UpdateCalledCopyRatioSegmentation {
                input:
                    sample = pair.left,
                    called_copy_ratio_segmentation = pair.right,
            }
        }
    }
    Array[Sample] samples_14 = select_first([UpdateCalledCopyRatioSegmentation.updated_sample, samples_13])

    if (defined(acs_copy_ratio_segmentations)) {
        scatter (pair in zip(samples_14, select_first([acs_copy_ratio_segmentations, []]))) {
            call s.UpdateSample as UpdateAcsCopyRatioSegmentation {
                input:
                    sample = pair.left,
                    acs_copy_ratio_segmentation = pair.right,
            }
        }
    }
    Array[Sample] samples_15 = select_first([UpdateAcsCopyRatioSegmentation.updated_sample, samples_14])

    if (defined(acs_copy_ratio_skews)) {
        scatter (pair in zip(samples_15, select_first([acs_copy_ratio_skews, []]))) {
            call s.UpdateSample as UpdateAcsCopyRatioSkew {
                input:
                    sample = pair.left,
                    acs_copy_ratio_skew = pair.right,
            }
        }
    }
    Array[Sample] samples_16 = select_first([UpdateAcsCopyRatioSkew.updated_sample, samples_15])

    if (defined(annotated_variants)) {
        scatter (pair in zip(samples_16, select_first([annotated_variants, []]))) {
            call s.UpdateSample as UpdateAnnotatedVariants {
                input:
                    sample = pair.left,
                    annotated_variants = pair.right,
            }
        }
    }
    Array[Sample] samples_17 = select_first([UpdateAnnotatedVariants.updated_sample, samples_16])

    # Select tumor and normal samples:

    scatter (tumor_sample in patient.tumor_samples) {
        scatter (sample in samples_17) {
            if (sample.name == tumor_sample.name) {
                Sample selected_tumor_sample = sample
            }
        }
        Sample tumor_samples = select_all(selected_tumor_sample)[0]
    }

    if (patient.has_normal) {
        scatter (normal_sample in patient.normal_samples) {
            scatter (sample in samples_17) {
                if (sample.name == normal_sample.name) {
                    Sample selected_normal_sample = sample
                }
            }
            Sample normal_samples = select_all(selected_normal_sample)[0]
        }
    }
    Array[Sample] non_optional_normal_samples = select_first([normal_samples, []])

    if (defined(patient.matched_normal_sample)) {
        Sample previous_matched_normal_sample = select_first([patient.matched_normal_sample])
        scatter (sample in non_optional_normal_samples) {
            if (sample.name == previous_matched_normal_sample.name) {
                Sample matched_normal_samples = sample
            }
        }
        Sample matched_normal_sample = select_all(matched_normal_samples)[0]
    }

    # Update patient:

    Patient pat = object {
        name: patient.name,
        samples: flatten([tumor_samples, non_optional_normal_samples]),
        tumor_samples: tumor_samples,
        normal_samples: non_optional_normal_samples,
        has_tumor: patient.has_tumor,
        has_normal: patient.has_normal,
        matched_normal_sample: if defined(matched_normal_sample) then matched_normal_sample else patient.matched_normal_sample,
    }

    output {
        Patient updated_patient = pat
    }
}