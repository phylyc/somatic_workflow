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
    Array[Array[SequencingRun]] sr = select_first([sequencing_runs, [[]]])
    if (length(sr) > 0) {
        scatter (pair in zip(patient.samples, sr)) {
            call s.UpdateSample as UpdateSequencingRuns {
                input:
                    sample = pair.left,
                    sequencing_runs = pair.right,
            }
        }
    }
    Array[Sample] samples_sr = select_first([UpdateSequencingRuns.updated_sample, patient.samples])

    Array[File] hcl = select_first([harmonized_callable_loci, []])
    if (length(hcl) > 0) {
        scatter (pair in zip(samples_sr, hcl)) {
            call s.UpdateSample as UpdateHarmonizedCallableLoci {
                input:
                    sample = pair.left,
                    harmonized_callable_loci = pair.right,
            }
        }
    }
    Array[Sample] samples_hcl = select_first([UpdateHarmonizedCallableLoci.updated_sample, samples_sr])

    Array[File] hdtcr = select_first([harmonized_denoised_total_copy_ratios, []])
    if (length(hdtcr) > 0) {
        scatter (pair in zip(samples_hcl, hdtcr)) {
            call s.UpdateSample as UpdateHarmonizedDenoisedTotalCopyRatios {
                input:
                    sample = pair.left,
                    harmonized_denoised_total_copy_ratios = pair.right,
            }
        }
    }
    Array[Sample] samples_hdtcr = select_first([UpdateHarmonizedDenoisedTotalCopyRatios.updated_sample, samples_hcl])

    Array[File] hsap = select_first([harmonized_snppanel_allelic_pileup_summaries, []])
    if (length(hsap) > 0) {
        scatter (pair in zip(samples_hdtcr, hsap)) {
            call s.UpdateSample as UpdateHarmonizedSnpPanelAllelicPileupSummaries {
                input:
                    sample = pair.left,
                    harmonized_snppanel_allelic_pileup_summaries = pair.right,
            }
        }
    }
    Array[Sample] samples_hsap = select_first([UpdateHarmonizedSnpPanelAllelicPileupSummaries.updated_sample, samples_hdtcr])

    Array[File] ct = select_first([harmonized_snppanel_allelic_pileup_summaries, []])
    if (length(ct) > 0) {
        scatter (pair in zip(samples_hsap, ct)) {
            call s.UpdateSample as UpdateContamination {
                input:
                    sample = pair.left,
                    contamination_table = pair.right,
            }
        }
    }
    Array[Sample] samples_ct = select_first([UpdateContamination.updated_sample, samples_hsap])

    Array[File] afst = select_first([af_segmentation_table, []])
    if (length(afst) > 0) {
        scatter (pair in zip(samples_ct, afst)) {
            call s.UpdateSample as UpdateAfSegmentation {
                input:
                    sample = pair.left,
                    af_segmentation_table = pair.right,
            }
        }
    }
    Array[Sample] samples_afst = select_first([UpdateAfSegmentation.updated_sample, samples_ct])

    Array[File] aps = select_first([allelic_pileup_summaries, []])
    if (length(aps) > 0) {
        scatter (pair in zip(samples_afst, aps)) {
            call s.UpdateSample as UpdateAllelicPileupSummaries {
                input:
                    sample = pair.left,
                    allelic_pileup_summaries = pair.right,
            }
        }
    }
    Array[Sample] samples_aps = select_first([UpdateAllelicPileupSummaries.updated_sample, samples_afst])

    Array[File] aacr = select_first([aggregated_allelic_read_counts, []])
    if (length(aacr) > 0) {
        scatter (pair in zip(samples_aps, aacr)) {
            call s.UpdateSample as UpdateAggregatedAllelicReadCounts {
                input:
                    sample = pair.left,
                    aggregated_allelic_read_counts = pair.right,
            }
        }
    }
    Array[Sample] samples_aarc = select_first([UpdateAggregatedAllelicReadCounts.updated_sample, samples_aps])

    Array[Float] gep = select_first([genotype_error_probabilities, []])
    if (length(gep) > 0) {
        scatter (pair in zip(samples_aarc, gep)) {
            call s.UpdateSample as UpdateGenotypeErrorProbabilities {
                input:
                    sample = pair.left,
                    genotype_error_probabilities = pair.right,
            }
        }
    }
    Array[Sample] samples_gep = select_first([UpdateGenotypeErrorProbabilities.updated_sample, samples_aarc])

    Array[File] afmp = select_first([af_model_parameters, []])
    if (length(afmp) > 0) {
        scatter (pair in zip(samples_gep, afmp)) {
            call s.UpdateSample as UpdateAfModelParameters {
                input:
                    sample = pair.left,
                    af_model_parameters = pair.right,
            }
        }
    }
    Array[Sample] samples_afmp = select_first([UpdateAfModelParameters.updated_sample, samples_gep])

    Array[File] crmp = select_first([cr_model_parameters, []])
    if (length(crmp) > 0) {
        scatter (pair in zip(samples_afmp, crmp)) {
            call s.UpdateSample as UpdateCrModelParameters {
                input:
                    sample = pair.left,
                    cr_model_parameters = pair.right,
            }
        }
    }
    Array[Sample] samples_crmp = select_first([UpdateCrModelParameters.updated_sample, samples_afmp])

    Array[File] ccrs = select_first([called_copy_ratio_segmentation, []])
    if (length(ccrs) > 0) {
        scatter (pair in zip(samples_crmp, ccrs)) {
            call s.UpdateSample as UpdateCalledCopyRatioSegmentation {
                input:
                    sample = pair.left,
                    called_copy_ratio_segmentation = pair.right,
            }
        }
    }
    Array[Sample] samples_ccrs = select_first([UpdateCalledCopyRatioSegmentation.updated_sample, samples_crmp])

    Array[File] crp = select_first([cr_plot, []])
    if (length(crp) > 0) {
        scatter (pair in zip(samples_ccrs, crp)) {
            call s.UpdateSample as UpdateCrPlot {
                input:
                    sample = pair.left,
                    cr_plot = pair.right,
            }
        }
    }
    Array[Sample] samples_crp = select_first([UpdateCrPlot.updated_sample, samples_ccrs])

    Array[File] acrs = select_first([acs_copy_ratio_segmentation, []])
    if (length(acrs) > 0) {
        scatter (pair in zip(samples_crp, acrs)) {
            call s.UpdateSample as UpdateAcsCopyRatioSegmentation {
                input:
                    sample = pair.left,
                    acs_copy_ratio_segmentation = pair.right,
            }
        }
    }
    Array[Sample] samples_acrs = select_first([UpdateAcsCopyRatioSegmentation.updated_sample, samples_crp])

    Array[Float] acrsks = select_first([acs_copy_ratio_skew, []])
    if (length(acrsks) > 0) {
        scatter (pair in zip(samples_acrs, acrsks)) {
            call s.UpdateSample as UpdateAcsCopyRatioSkew {
                input:
                    sample = pair.left,
                    acs_copy_ratio_skew = pair.right,
            }
        }
    }
    Array[Sample] samples_acrsks = select_first([UpdateAcsCopyRatioSkew.updated_sample, samples_acrs])

    Array[File] asv = select_first([annotated_somatic_variants, []])
    if (length(asv) > 0) {
        scatter (pair in zip(samples_acrsks, asv)) {
            call s.UpdateSample as UpdateAnnotatedSomaticVariants {
                input:
                    sample = pair.left,
                    annotated_somatic_variants = pair.right,
            }
        }
    }
    Array[Sample] samples_asv = select_first([UpdateAnnotatedSomaticVariants.updated_sample, samples_acrsks])

    Array[File] asvi = select_first([annotated_somatic_variants_idx, []])
    if (length(asvi) > 0) {
        scatter (pair in zip(samples_asv, asvi)) {
            call s.UpdateSample as UpdateAnnotatedSomaticVariantsIdx {
                input:
                    sample = pair.left,
                    annotated_somatic_variants_idx = pair.right,
            }
        }
    }
    Array[Sample] samples_asvi = select_first([UpdateAnnotatedSomaticVariantsIdx.updated_sample, samples_asv])

    Array[File] ard = select_first([aggregated_allelic_read_counts, []])
    if (length(ard) > 0) {
        scatter (pair in zip(samples_asvi, ard)) {
            call s.UpdateSample as UpdateAbsoluteRData {
                input:
                    sample = pair.left,
                    absolute_acr_rdata = pair.right,
            }
        }
    }
    Array[Sample] samples_ard = select_first([UpdateAbsoluteRData.updated_sample, samples_asvi])

    Array[File] acrp = select_first([absolute_acr_plot, []])
    if (length(acrp) > 0) {
        scatter (pair in zip(samples_ard, acrp)) {
            call s.UpdateSample as UpdateAbsolutePlot {
                input:
                    sample = pair.left,
                    absolute_acr_plot = pair.right,
            }
        }
    }
    Array[Sample] samples_acrp = select_first([UpdateAbsolutePlot.updated_sample, samples_ard])

    Array[File] asnm = select_first([absolute_snv_maf, []])
    if (length(asnm) > 0) {
        scatter (pair in zip(samples_acrp, asnm)) {
            call s.UpdateSample as UpdateAbsoluteSnvMaf {
                input:
                    sample = pair.left,
                    absolute_snv_maf = pair.right,
            }
        }
    }
    Array[Sample] samples_asnm = select_first([UpdateAbsoluteSnvMaf.updated_sample, samples_acrp])

    Array[File] asim = select_first([absolute_indel_maf, []])
    if (length(asim) > 0) {
        scatter (pair in zip(samples_asnm, asim)) {
            call s.UpdateSample as UpdateAbsoluteIndelMaf {
                input:
                    sample = pair.left,
                    absolute_indel_maf = pair.right,
            }
        }
    }
    Array[Sample] samples_asim = select_first([UpdateAbsoluteIndelMaf.updated_sample, samples_asnm])

    Array[Int] asol = select_first([absolute_solution, []])
    if (length(asol) > 0) {
        scatter (pair in zip(samples_asim, asol)) {
            call s.UpdateSample as UpdateAbsoluteSolution {
                input:
                    sample = pair.left,
                    absolute_solution = pair.right,
            }
        }
    }
    Array[Sample] samples_asol = select_first([UpdateAbsoluteSolution.updated_sample, samples_asim])

    Array[File] am = select_first([absolute_maf, []])
    if (length(am) > 0) {
        scatter (pair in zip(samples_asol, am)) {
            call s.UpdateSample as UpdateAbsoluteMaf {
                input:
                    sample = pair.left,
                    absolute_maf = pair.right,
            }
        }
    }
    Array[Sample] samples_am = select_first([UpdateAbsoluteMaf.updated_sample, samples_asol])

    Array[File] asg = select_first([absolute_segtab, []])
    if (length(asg) > 0) {
        scatter (pair in zip(samples_am, asg)) {
            call s.UpdateSample as UpdateAbsoluteSegtab {
                input:
                    sample = pair.left,
                    absolute_segtab = pair.right,
            }
        }
    }
    Array[Sample] samples_asg = select_first([UpdateAbsoluteSegtab.updated_sample, samples_am])

    Array[File] at = select_first([absolute_table, []])
    if (length(at) > 0) {
        scatter (pair in zip(samples_asg, at)) {
            call s.UpdateSample as UpdateAbsoluteTable {
                input:
                    sample = pair.left,
                    absolute_table = pair.right,
            }
        }
    }
    Array[Sample] samples_at = select_first([UpdateAbsoluteTable.updated_sample, samples_asg])

    Array[Float] purity_ = select_first([purity, []])
    if (length(purity_) > 0) {
        scatter (pair in zip(samples_at, purity_)) {
            call s.UpdateSample as UpdatePurity {
                input:
                    sample = pair.left,
                    purity = pair.right,
            }
        }
    }
    Array[Sample] samples_purity = select_first([UpdatePurity.updated_sample, samples_at])

    Array[Float] ploidy_ = select_first([ploidy, []])
    if (length(ploidy_) > 0) {
        scatter (pair in zip(samples_purity, ploidy_)) {
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
        if (defined(patient.matched_normal_sample)) {
            Sample previous_matched_normal_sample = select_first([patient.matched_normal_sample])
            scatter (sample in normal_samples) {
                if (sample.name == previous_matched_normal_sample.name) {
                    Sample matched_normal_samples = sample
                }
            }
            Sample matched_normal_sample = select_all(matched_normal_samples)[0]
        }
    }

    # Update patient:
    call p.UpdatePatient {
        input:
            patient = patient,
            samples = samples,
            tumor_samples = tumor_samples,
            normal_samples = normal_samples,
            matched_normal_sample = matched_normal_sample
    }

    output {
        Patient updated_patient = UpdatePatient.updated_patient
    }
}