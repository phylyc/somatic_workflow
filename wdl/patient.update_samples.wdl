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

        Array[File]? absolute_snv_maf
        Array[File]? absolute_indel_maf

        Array[File]? absolute_acr_rdata
        Array[File]? absolute_acr_plot
        Array[Int]? absolute_acr_solution
        Array[File]? absolute_acr_maf
        Array[File]? absolute_acr_segtab
        Array[File]? absolute_acr_segtab_igv
        Array[File]? absolute_acr_table
        Array[Float]? absolute_acr_purity
        Array[Float]? absolute_acr_ploidy

        Array[File]? absolute_tcr_rdata
        Array[File]? absolute_tcr_plot
        Array[Int]? absolute_tcr_solution
        Array[File]? absolute_tcr_maf
        Array[File]? absolute_tcr_segtab
        Array[File]? absolute_tcr_segtab_igv
        Array[File]? absolute_tcr_table
        Array[Float]? absolute_tcr_purity
        Array[Float]? absolute_tcr_ploidy

        Array[Int]? timepoint
    }

    # Build each updated Sample in a single pass. Every optional input array carries one
    # entry per sample (aligned by index) and overlays the matching Sample field; a field
    # whose input array is absent or empty is carried over from the existing sample
    # unchanged. Fields are overlaid independently because each input targets a distinct
    # field.
    #
    # Inputs are defaulted to empty arrays so that length() and indexing are always
    # well-defined. INVARIANT: any provided overlay array has length == length(patient.samples).
    Array[Array[SequencingRun]] sequencing_runs_arr = select_first([sequencing_runs, [[]]])
    Array[File] harmonized_callable_loci_arr = select_first([harmonized_callable_loci, []])
    Array[File] harmonized_denoised_total_copy_ratios_arr = select_first([harmonized_denoised_total_copy_ratios, []])
    Array[File] harmonized_snppanel_allelic_pileup_summaries_arr = select_first([harmonized_snppanel_allelic_pileup_summaries, []])
    Array[File] contamination_table_arr = select_first([contamination_table, []])
    Array[File] af_segmentation_table_arr = select_first([af_segmentation_table, []])
    Array[File] allelic_pileup_summaries_arr = select_first([allelic_pileup_summaries, []])
    Array[File] aggregated_allelic_read_counts_arr = select_first([aggregated_allelic_read_counts, []])
    Array[Float] genotype_error_probabilities_arr = select_first([genotype_error_probabilities, []])
    Array[File] af_model_parameters_arr = select_first([af_model_parameters, []])
    Array[File] cr_model_parameters_arr = select_first([cr_model_parameters, []])
    Array[File] called_copy_ratio_segmentation_arr = select_first([called_copy_ratio_segmentation, []])
    Array[File] cr_plot_arr = select_first([cr_plot, []])
    Array[File] acs_copy_ratio_segmentation_arr = select_first([acs_copy_ratio_segmentation, []])
    Array[Float] acs_copy_ratio_skew_arr = select_first([acs_copy_ratio_skew, []])
    Array[File] annotated_somatic_variants_arr = select_first([annotated_somatic_variants, []])
    Array[File] annotated_somatic_variants_idx_arr = select_first([annotated_somatic_variants_idx, []])
    Array[File] absolute_snv_maf_arr = select_first([absolute_snv_maf, []])
    Array[File] absolute_indel_maf_arr = select_first([absolute_indel_maf, []])
    Array[File] absolute_acr_rdata_arr = select_first([absolute_acr_rdata, []])
    Array[File] absolute_acr_plot_arr = select_first([absolute_acr_plot, []])
    Array[Int] absolute_acr_solution_arr = select_first([absolute_acr_solution, []])
    Array[File] absolute_acr_maf_arr = select_first([absolute_acr_maf, []])
    Array[File] absolute_acr_segtab_arr = select_first([absolute_acr_segtab, []])
    Array[File] absolute_acr_segtab_igv_arr = select_first([absolute_acr_segtab_igv, []])
    Array[File] absolute_acr_table_arr = select_first([absolute_acr_table, []])
    Array[Float] absolute_acr_purity_arr = select_first([absolute_acr_purity, []])
    Array[Float] absolute_acr_ploidy_arr = select_first([absolute_acr_ploidy, []])
    Array[File] absolute_tcr_rdata_arr = select_first([absolute_tcr_rdata, []])
    Array[File] absolute_tcr_plot_arr = select_first([absolute_tcr_plot, []])
    Array[Int] absolute_tcr_solution_arr = select_first([absolute_tcr_solution, []])
    Array[File] absolute_tcr_maf_arr = select_first([absolute_tcr_maf, []])
    Array[File] absolute_tcr_segtab_arr = select_first([absolute_tcr_segtab, []])
    Array[File] absolute_tcr_segtab_igv_arr = select_first([absolute_tcr_segtab_igv, []])
    Array[File] absolute_tcr_table_arr = select_first([absolute_tcr_table, []])
    Array[Float] absolute_tcr_purity_arr = select_first([absolute_tcr_purity, []])
    Array[Float] absolute_tcr_ploidy_arr = select_first([absolute_tcr_ploidy, []])
    Array[Int] timepoint_arr = select_first([timepoint, []])

    scatter (i in range(length(patient.samples))) {
        Sample sample_i = patient.samples[i]
        Sample updated_sample_i = object {
            name: sample_i.name,
            bam_name: sample_i.bam_name,
            is_tumor: sample_i.is_tumor,
            sequencing_runs: if length(flatten(sequencing_runs_arr)) > 0 then sequencing_runs_arr[i] else sample_i.sequencing_runs,

            harmonized_callable_loci: if length(harmonized_callable_loci_arr) > 0 then harmonized_callable_loci_arr[i] else sample_i.harmonized_callable_loci,
            harmonized_denoised_total_copy_ratios: if length(harmonized_denoised_total_copy_ratios_arr) > 0 then harmonized_denoised_total_copy_ratios_arr[i] else sample_i.harmonized_denoised_total_copy_ratios,
            harmonized_snppanel_allelic_pileup_summaries: if length(harmonized_snppanel_allelic_pileup_summaries_arr) > 0 then harmonized_snppanel_allelic_pileup_summaries_arr[i] else sample_i.harmonized_snppanel_allelic_pileup_summaries,
            contamination_table: if length(contamination_table_arr) > 0 then contamination_table_arr[i] else sample_i.contamination_table,
            af_segmentation_table: if length(af_segmentation_table_arr) > 0 then af_segmentation_table_arr[i] else sample_i.af_segmentation_table,
            allelic_pileup_summaries: if length(allelic_pileup_summaries_arr) > 0 then allelic_pileup_summaries_arr[i] else sample_i.allelic_pileup_summaries,
            aggregated_allelic_read_counts: if length(aggregated_allelic_read_counts_arr) > 0 then aggregated_allelic_read_counts_arr[i] else sample_i.aggregated_allelic_read_counts,
            genotype_error_probabilities: if length(genotype_error_probabilities_arr) > 0 then genotype_error_probabilities_arr[i] else sample_i.genotype_error_probabilities,
            af_model_parameters: if length(af_model_parameters_arr) > 0 then af_model_parameters_arr[i] else sample_i.af_model_parameters,
            cr_model_parameters: if length(cr_model_parameters_arr) > 0 then cr_model_parameters_arr[i] else sample_i.cr_model_parameters,
            called_copy_ratio_segmentation: if length(called_copy_ratio_segmentation_arr) > 0 then called_copy_ratio_segmentation_arr[i] else sample_i.called_copy_ratio_segmentation,
            cr_plot: if length(cr_plot_arr) > 0 then cr_plot_arr[i] else sample_i.cr_plot,
            acs_copy_ratio_segmentation: if length(acs_copy_ratio_segmentation_arr) > 0 then acs_copy_ratio_segmentation_arr[i] else sample_i.acs_copy_ratio_segmentation,
            acs_copy_ratio_skew: if length(acs_copy_ratio_skew_arr) > 0 then acs_copy_ratio_skew_arr[i] else sample_i.acs_copy_ratio_skew,
            annotated_somatic_variants: if length(annotated_somatic_variants_arr) > 0 then annotated_somatic_variants_arr[i] else sample_i.annotated_somatic_variants,
            annotated_somatic_variants_idx: if length(annotated_somatic_variants_idx_arr) > 0 then annotated_somatic_variants_idx_arr[i] else sample_i.annotated_somatic_variants_idx,

            absolute_snv_maf: if length(absolute_snv_maf_arr) > 0 then absolute_snv_maf_arr[i] else sample_i.absolute_snv_maf,
            absolute_indel_maf: if length(absolute_indel_maf_arr) > 0 then absolute_indel_maf_arr[i] else sample_i.absolute_indel_maf,

            absolute_acr_rdata: if length(absolute_acr_rdata_arr) > 0 then absolute_acr_rdata_arr[i] else sample_i.absolute_acr_rdata,
            absolute_acr_plot: if length(absolute_acr_plot_arr) > 0 then absolute_acr_plot_arr[i] else sample_i.absolute_acr_plot,
            absolute_acr_solution: if length(absolute_acr_solution_arr) > 0 then absolute_acr_solution_arr[i] else sample_i.absolute_acr_solution,
            absolute_acr_maf: if length(absolute_acr_maf_arr) > 0 then absolute_acr_maf_arr[i] else sample_i.absolute_acr_maf,
            absolute_acr_segtab: if length(absolute_acr_segtab_arr) > 0 then absolute_acr_segtab_arr[i] else sample_i.absolute_acr_segtab,
            absolute_acr_segtab_igv: if length(absolute_acr_segtab_igv_arr) > 0 then absolute_acr_segtab_igv_arr[i] else sample_i.absolute_acr_segtab_igv,
            absolute_acr_table: if length(absolute_acr_table_arr) > 0 then absolute_acr_table_arr[i] else sample_i.absolute_acr_table,
            absolute_acr_purity: if length(absolute_acr_purity_arr) > 0 then absolute_acr_purity_arr[i] else sample_i.absolute_acr_purity,
            absolute_acr_ploidy: if length(absolute_acr_ploidy_arr) > 0 then absolute_acr_ploidy_arr[i] else sample_i.absolute_acr_ploidy,

            absolute_tcr_rdata: if length(absolute_tcr_rdata_arr) > 0 then absolute_tcr_rdata_arr[i] else sample_i.absolute_tcr_rdata,
            absolute_tcr_plot: if length(absolute_tcr_plot_arr) > 0 then absolute_tcr_plot_arr[i] else sample_i.absolute_tcr_plot,
            absolute_tcr_solution: if length(absolute_tcr_solution_arr) > 0 then absolute_tcr_solution_arr[i] else sample_i.absolute_tcr_solution,
            absolute_tcr_maf: if length(absolute_tcr_maf_arr) > 0 then absolute_tcr_maf_arr[i] else sample_i.absolute_tcr_maf,
            absolute_tcr_segtab: if length(absolute_tcr_segtab_arr) > 0 then absolute_tcr_segtab_arr[i] else sample_i.absolute_tcr_segtab,
            absolute_tcr_segtab_igv: if length(absolute_tcr_segtab_igv_arr) > 0 then absolute_tcr_segtab_igv_arr[i] else sample_i.absolute_tcr_segtab_igv,
            absolute_tcr_table: if length(absolute_tcr_table_arr) > 0 then absolute_tcr_table_arr[i] else sample_i.absolute_tcr_table,
            absolute_tcr_purity: if length(absolute_tcr_purity_arr) > 0 then absolute_tcr_purity_arr[i] else sample_i.absolute_tcr_purity,
            absolute_tcr_ploidy: if length(absolute_tcr_ploidy_arr) > 0 then absolute_tcr_ploidy_arr[i] else sample_i.absolute_tcr_ploidy,

            timepoint: if length(timepoint_arr) > 0 then timepoint_arr[i] else sample_i.timepoint
        }
    }

    Array[Sample] samples = updated_sample_i

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

    # Assemble the updated patient inline (a plain struct construction, not a sub-workflow
    # call). Only the sample sets change here; every other field is carried over from the
    # input patient.
    Patient updated_patient_obj = object {
        name: patient.name,
        sex: patient.sex,
        samples: samples,
        tumor_samples: tumor_samples,
        normal_samples: select_first([normal_samples, patient.normal_samples]),
        has_tumor: patient.has_tumor,
        has_normal: patient.has_normal,
        matched_normal_sample: if defined(matched_normal_sample) then matched_normal_sample else patient.matched_normal_sample,
        shards: patient.shards,
        raw_snv_calls_vcf: patient.raw_snv_calls_vcf,
        raw_snv_calls_vcf_idx: patient.raw_snv_calls_vcf_idx,
        mutect2_stats: patient.mutect2_stats,
        orientation_bias: patient.orientation_bias,
        filtered_vcf: patient.filtered_vcf,
        filtered_vcf_idx: patient.filtered_vcf_idx,
        filtering_stats: patient.filtering_stats,
        somatic_vcf: patient.somatic_vcf,
        somatic_vcf_idx: patient.somatic_vcf_idx,
        num_somatic_variants: patient.num_somatic_variants,
        germline_vcf: patient.germline_vcf,
        germline_vcf_idx: patient.germline_vcf_idx,
        num_germline_variants: patient.num_germline_variants,
        somatic_calls_bam: patient.somatic_calls_bam,
        somatic_calls_bai: patient.somatic_calls_bai,
        rare_germline_alleles: patient.rare_germline_alleles,
        rare_germline_alleles_idx: patient.rare_germline_alleles_idx,
        gvcf: patient.gvcf,
        gvcf_idx: patient.gvcf_idx,
        snp_ref_counts: patient.snp_ref_counts,
        snp_alt_counts: patient.snp_alt_counts,
        snp_other_alt_counts: patient.snp_other_alt_counts,
        snp_sample_correlation: patient.snp_sample_correlation,
        snp_sample_correlation_min: patient.snp_sample_correlation_min,
        modeled_segments: patient.modeled_segments,
        mask_vcf: patient.mask_vcf,
        mask_vcf_idx: patient.mask_vcf_idx,
        mask_name: patient.mask_name
    }

    output {
        Patient updated_patient = updated_patient_obj
    }
}
