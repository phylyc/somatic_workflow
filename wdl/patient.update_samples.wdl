version development

import "sample.wdl" as s
import "patient.wdl" as p


workflow UpdateSamples {
    input {
        Patient patient
        Array[Array[SequencingRun]]? sequencing_runs
        Array[File]? harmonized_callable_loci
        Array[File]? harmonized_denoised_total_copy_ratios
        Array[File]? harmonized_snppanel_allelic_pileup_summaries
        Array[File]? contamination_table
        Array[File]? af_segmentation_table
        Array[File]? allelic_pileup_summaries
        Array[File]? aggregated_allelic_read_counts
        Array[Float]? genotype_error_probabilities
        Array[File]? af_model_parameters
        Array[File]? cr_model_parameters
        Array[File]? called_copy_ratio_segmentation
        Array[File]? cr_plot
        Array[File]? acs_copy_ratio_segmentation
        Array[Float]? acs_copy_ratio_skew
        Array[File]? annotated_somatic_variants
        Array[File]? annotated_somatic_variants_idx
        Array[File]? absolute_acr_rdata
        Array[File]? absolute_acr_plot
        Array[File]? absolute_snv_maf
        Array[File]? absolute_indel_maf
        Array[Int]? absolute_solution
        Array[File]? absolute_maf
        Array[File]? absolute_segtab
        Array[File]? absolute_table
        Array[Float]? purity
        Array[Float]? ploidy
    }

    # Update samples:
    if (defined(sequencing_runs)) {
        scatter (pair in zip(patient.samples, select_first([sequencing_runs, []]))) {
            call s.UpdateSample as UpdateSequencingRuns {
                input:
                    sample = pair.left,
                    sequencing_runs = pair.right,
            }
        }
    }
    Array[Sample] samples_sr = select_first([UpdateSequencingRuns.updated_sample, patient.samples])

    if (defined(harmonized_callable_loci)) {
        scatter (pair in zip(samples_sr, select_first([harmonized_callable_loci, []]))) {
            call s.UpdateSample as UpdateHarmonizedCallableLoci {
                input:
                    sample = pair.left,
                    harmonized_callable_loci = pair.right,
            }
        }
    }
    Array[Sample] samples_hcl = select_first([UpdateHarmonizedCallableLoci.updated_sample, samples_sr])

    if (defined(harmonized_denoised_total_copy_ratios)) {
        scatter (pair in zip(samples_hcl, select_first([harmonized_denoised_total_copy_ratios, []]))) {
            call s.UpdateSample as UpdateHarmonizedDenoisedTotalCopyRatios {
                input:
                    sample = pair.left,
                    harmonized_denoised_total_copy_ratios = pair.right,
            }
        }
    }
    Array[Sample] samples_hdtcr = select_first([UpdateHarmonizedDenoisedTotalCopyRatios.updated_sample, samples_hcl])

    if (defined(harmonized_snppanel_allelic_pileup_summaries)) {
        scatter (pair in zip(samples_hdtcr, select_first([harmonized_snppanel_allelic_pileup_summaries, []]))) {
            call s.UpdateSample as UpdateHarmonizedSnpPanelAllelicPileupSummaries {
                input:
                    sample = pair.left,
                    harmonized_snppanel_allelic_pileup_summaries = pair.right,
            }
        }
    }
    Array[Sample] samples_hsap = select_first([UpdateHarmonizedSnpPanelAllelicPileupSummaries.updated_sample, samples_hdtcr])

    if (defined(contamination_table)) {
        scatter (pair in zip(samples_hsap, select_first([contamination_table, []]))) {
            call s.UpdateSample as UpdateContamination {
                input:
                    sample = pair.left,
                    contamination_table = pair.right,
            }
        }
    }
    Array[Sample] samples_ct = select_first([UpdateContamination.updated_sample, samples_hsap])

    if (defined(af_segmentation_table)) {
        scatter (pair in zip(samples_ct, select_first([af_segmentation_table, []]))) {
            call s.UpdateSample as UpdateAfSegmentation {
                input:
                    sample = pair.left,
                    af_segmentation_table = pair.right,
            }
        }
    }
    Array[Sample] samples_afst = select_first([UpdateAfSegmentation.updated_sample, samples_ct])

    if (defined(allelic_pileup_summaries)) {
        scatter (pair in zip(samples_afst, select_first([allelic_pileup_summaries, []]))) {
            call s.UpdateSample as UpdateAllelicPileupSummaries {
                input:
                    sample = pair.left,
                    allelic_pileup_summaries = pair.right,
            }
        }
    }
    Array[Sample] samples_aps = select_first([UpdateAllelicPileupSummaries.updated_sample, samples_afst])

    if (defined(aggregated_allelic_read_counts)) {
        scatter (pair in zip(samples_aps, select_first([aggregated_allelic_read_counts, []]))) {
            call s.UpdateSample as UpdateAggregatedAllelicReadCounts {
                input:
                    sample = pair.left,
                    aggregated_allelic_read_counts = pair.right,
            }
        }
    }
    Array[Sample] samples_aarc = select_first([UpdateAggregatedAllelicReadCounts.updated_sample, samples_aps])

    if (defined(genotype_error_probabilities)) {
        scatter (pair in zip(samples_aarc, select_first([genotype_error_probabilities, []]))) {
            call s.UpdateSample as UpdateGenotypeErrorProbabilities {
                input:
                    sample = pair.left,
                    genotype_error_probabilities = pair.right,
            }
        }
    }
    Array[Sample] samples_gep = select_first([UpdateGenotypeErrorProbabilities.updated_sample, samples_aarc])

    if (defined(af_model_parameters)) {
        scatter (pair in zip(samples_gep, select_first([af_model_parameters, []]))) {
            call s.UpdateSample as UpdateAfModelParameters {
                input:
                    sample = pair.left,
                    af_model_parameters = pair.right,
            }
        }
    }
    Array[Sample] samples_afmp = select_first([UpdateAfModelParameters.updated_sample, samples_gep])

    if (defined(cr_model_parameters)) {
        scatter (pair in zip(samples_afmp, select_first([cr_model_parameters, []]))) {
            call s.UpdateSample as UpdateCrModelParameters {
                input:
                    sample = pair.left,
                    cr_model_parameters = pair.right,
            }
        }
    }
    Array[Sample] samples_crmp = select_first([UpdateCrModelParameters.updated_sample, samples_afmp])

    if (defined(called_copy_ratio_segmentation)) {
        scatter (pair in zip(samples_crmp, select_first([called_copy_ratio_segmentation, []]))) {
            call s.UpdateSample as UpdateCalledCopyRatioSegmentation {
                input:
                    sample = pair.left,
                    called_copy_ratio_segmentation = pair.right,
            }
        }
    }
    Array[Sample] samples_ccrs = select_first([UpdateCalledCopyRatioSegmentation.updated_sample, samples_crmp])

    if (defined(cr_plot)) {
        scatter (pair in zip(samples_ccrs, select_first([cr_plot, []]))) {
            call s.UpdateSample as UpdateCrPlot {
                input:
                    sample = pair.left,
                    cr_plot = pair.right,
            }
        }
    }
    Array[Sample] samples_crp = select_first([UpdateCrPlot.updated_sample, samples_ccrs])

    if (defined(acs_copy_ratio_segmentation)) {
        scatter (pair in zip(samples_crp, select_first([acs_copy_ratio_segmentation, []]))) {
            call s.UpdateSample as UpdateAcsCopyRatioSegmentation {
                input:
                    sample = pair.left,
                    acs_copy_ratio_segmentation = pair.right,
            }
        }
    }
    Array[Sample] samples_acrs = select_first([UpdateAcsCopyRatioSegmentation.updated_sample, samples_crp])

    if (defined(acs_copy_ratio_skew)) {
        scatter (pair in zip(samples_acrs, select_first([acs_copy_ratio_skew, []]))) {
            call s.UpdateSample as UpdateAcsCopyRatioSkew {
                input:
                    sample = pair.left,
                    acs_copy_ratio_skew = pair.right,
            }
        }
    }
    Array[Sample] samples_acrsks = select_first([UpdateAcsCopyRatioSkew.updated_sample, samples_acrs])

    if (defined(annotated_somatic_variants)) {
        scatter (pair in zip(samples_acrsks, select_first([annotated_somatic_variants, []]))) {
            call s.UpdateSample as UpdateAnnotatedSomaticVariants {
                input:
                    sample = pair.left,
                    annotated_somatic_variants = pair.right,
            }
        }
    }
    Array[Sample] samples_asv = select_first([UpdateAnnotatedSomaticVariants.updated_sample, samples_acrsks])

    if (defined(annotated_somatic_variants_idx)) {
        scatter (pair in zip(samples_asv, select_first([annotated_somatic_variants_idx, []]))) {
            call s.UpdateSample as UpdateAnnotatedSomaticVariantsIdx {
                input:
                    sample = pair.left,
                    annotated_somatic_variants_idx = pair.right,
            }
        }
    }
    Array[Sample] samples_asvi = select_first([UpdateAnnotatedSomaticVariantsIdx.updated_sample, samples_asv])

    if (defined(absolute_acr_rdata)) {
        scatter (pair in zip(samples_asvi, select_first([absolute_acr_rdata, []]))) {
            call s.UpdateSample as UpdateAbsoluteRData {
                input:
                    sample = pair.left,
                    absolute_acr_rdata = pair.right,
            }
        }
    }
    Array[Sample] samples_ard = select_first([UpdateAbsoluteRData.updated_sample, samples_asvi])

    if (defined(absolute_acr_plot)) {
        scatter (pair in zip(samples_ard, select_first([absolute_acr_plot, []]))) {
            call s.UpdateSample as UpdateAbsolutePlot {
                input:
                    sample = pair.left,
                    absolute_acr_plot = pair.right,
            }
        }
    }
    Array[Sample] samples_acrp = select_first([UpdateAbsolutePlot.updated_sample, samples_ard])

    if (defined(absolute_snv_maf)) {
        scatter (pair in zip(samples_acrp, select_first([absolute_snv_maf, []]))) {
            call s.UpdateSample as UpdateAbsoluteSnvMaf {
                input:
                    sample = pair.left,
                    absolute_snv_maf = pair.right,
            }
        }
    }
    Array[Sample] samples_asnm = select_first([UpdateAbsoluteSnvMaf.updated_sample, samples_acrp])

    if (defined(absolute_indel_maf)) {
        scatter (pair in zip(samples_asnm, select_first([absolute_indel_maf, []]))) {
            call s.UpdateSample as UpdateAbsoluteIndelMaf {
                input:
                    sample = pair.left,
                    absolute_indel_maf = pair.right,
            }
        }
    }
    Array[Sample] samples_asim = select_first([UpdateAbsoluteIndelMaf.updated_sample, samples_asnm])

    if (defined(absolute_solution)) {
        scatter (pair in zip(samples_asim, select_first([absolute_solution, []]))) {
            call s.UpdateSample as UpdateAbsoluteSolution {
                input:
                    sample = pair.left,
                    absolute_solution = pair.right,
            }
        }
    }
    Array[Sample] samples_as = select_first([UpdateAbsoluteSolution.updated_sample, samples_asim])

    if (defined(absolute_maf)) {
        scatter (pair in zip(samples_as, select_first([absolute_maf, []]))) {
            call s.UpdateSample as UpdateAbsoluteMaf {
                input:
                    sample = pair.left,
                    absolute_maf = pair.right,
            }
        }
    }
    Array[Sample] samples_am = select_first([UpdateAbsoluteMaf.updated_sample, samples_as])

    if (defined(absolute_segtab)) {
        scatter (pair in zip(samples_am, select_first([absolute_segtab, []]))) {
            call s.UpdateSample as UpdateAbsoluteSegtab {
                input:
                    sample = pair.left,
                    absolute_segtab = pair.right,
            }
        }
    }
    Array[Sample] samples_asg = select_first([UpdateAbsoluteSegtab.updated_sample, samples_am])

    if (defined(absolute_table)) {
        scatter (pair in zip(samples_asg, select_first([absolute_table, []]))) {
            call s.UpdateSample as UpdateAbsoluteTable {
                input:
                    sample = pair.left,
                    absolute_table = pair.right,
            }
        }
    }
    Array[Sample] samples_at = select_first([UpdateAbsoluteTable.updated_sample, samples_asg])

    if (defined(purity)) {
        scatter (pair in zip(samples_at, select_first([purity, []]))) {
            call s.UpdateSample as UpdatePurity {
                input:
                    sample = pair.left,
                    purity = pair.right,
            }
        }
    }
    Array[Sample] samples_purity = select_first([UpdatePurity.updated_sample, samples_at])

    if (defined(ploidy)) {
        scatter (pair in zip(samples_purity, select_first([ploidy, []]))) {
            call s.UpdateSample as UpdatePloidy {
                input:
                    sample = pair.left,
                    ploidy = pair.right,
            }
        }
    }
    Array[Sample] samples_ploidy = select_first([UpdatePloidy.updated_sample, samples_purity])

    Array[Sample] samples = samples_ploidy

    # Select tumor and normal samples:

    scatter (tumor_sample in patient.tumor_samples) {
        scatter (sample in samples) {
            if (sample.name == tumor_sample.name) {
                Sample selected_tumor_sample = sample
            }
        }
        Sample tumor_samples = select_all(selected_tumor_sample)[0]
    }

    if (patient.has_normal) {
        scatter (normal_sample in patient.normal_samples) {
            scatter (sample in samples) {
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
    call p.UpdatePatient {
        input:
            patient = patient,
            samples = samples,
            tumor_samples = tumor_samples,
            normal_samples = non_optional_normal_samples,
            matched_normal_sample = matched_normal_sample
    }

    output {
        Patient updated_patient = UpdatePatient.updated_patient
    }
}