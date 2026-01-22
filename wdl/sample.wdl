version development

import "sequencing_run.wdl" as seqrun


struct Sample {
    String name
    String bam_name
    Array[SequencingRun] sequencing_runs
    Boolean is_tumor
                                                        ### ORIGIN
    File? harmonized_callable_loci                      # GATK CallableLoci | harmonization
    File? harmonized_denoised_total_copy_ratios         # GATK CollectReadCounts | DenoiseReadCounts | harmonization
    File? harmonized_snppanel_allelic_pileup_summaries  # GATK GetPileupSummaries | harmonization
    File? contamination_table                           # GATK CalculateContamination
    File? af_segmentation_table                         # GATK CalculateContamination / ModelSegments
    File? allelic_pileup_summaries                      # harmonized_snppanel_allelic_pileup_summaries + rare_germline_allelic_pileup_summaries
    File? aggregated_allelic_read_counts                # allelic_pileup_summaries | PileupToAllelicCounts
    Float? genotype_error_probabilities                 # PileupToAllelicCounts
    File? af_model_parameters                           # GATK ModelSegments
    File? cr_model_parameters                           # GATK ModelSegments
    File? called_copy_ratio_segmentation                # GATK ModelSegments + CallCopyRatioSegments + merge
    File? cr_plot                                       # GATK PlotModeledSegments
    File? acs_copy_ratio_segmentation                   # ModelSegmentsToACSConversion
    Float? acs_copy_ratio_skew                          # ModelSegmentsToACSConversion
    File? annotated_somatic_variants                    # GATK Funcotator
    File? annotated_somatic_variants_idx                # GATK Funcotator
    File? absolute_acr_rdata                            # ABSOLUTE
    File? absolute_acr_plot                             # ABSOLUTE
    File? absolute_snv_maf                              # ABSOLUTE
    File? absolute_indel_maf                            # ABSOLUTE
    Int? absolute_solution                              # manual
    File? absolute_maf                                  # ABSOLUTE
    File? absolute_segtab                               # ABSOLUTE
    File? absolute_maf_postprocessed                    # ABSOLUTE + Postprocess
    File? absolute_segtab_postprocessed                 # ABSOLUTE + Postprocess
    File? absolute_segtab_igv_postprocessed             # ABSOLUTE + Postprocess
    File? absolute_table                                # ABSOLUTE
    Float? purity                                       # ABSOLUTE
    Float? ploidy                                       # ABSOLUTE
    Int? timepoint                                      # manual
}


workflow UpdateSample {
    input {
        Sample sample
        String? name
        String? bam_name
        Array[SequencingRun]? sequencing_runs
        Boolean? is_tumor

        File? harmonized_callable_loci
        File? harmonized_denoised_total_copy_ratios
        File? harmonized_snppanel_allelic_pileup_summaries
        File? contamination_table
        File? af_segmentation_table
        File? allelic_pileup_summaries
        File? aggregated_allelic_read_counts
        Float? genotype_error_probabilities
        File? af_model_parameters
        File? cr_model_parameters
        File? called_copy_ratio_segmentation
        File? cr_plot
        File? acs_copy_ratio_segmentation
        Float? acs_copy_ratio_skew
        File? annotated_somatic_variants
        File? annotated_somatic_variants_idx
        File? absolute_acr_rdata
        File? absolute_acr_plot
        File? absolute_snv_maf
        File? absolute_indel_maf
        Int? absolute_solution
        File? absolute_maf
        File? absolute_segtab
        File? absolute_maf_postprocessed
        File? absolute_segtab_postprocessed
        File? absolute_segtab_igv_postprocessed
        File? absolute_table
        Float? purity
        Float? ploidy
        Int? timepoint
    }

    Sample s = object {
        name: select_first([name, sample.name]),
        bam_name: select_first([bam_name, sample.bam_name]),
        sequencing_runs: select_first([sequencing_runs, sample.sequencing_runs]),
        is_tumor: select_first([is_tumor, sample.is_tumor]),
        # cannot use select_first for optional fields:
        harmonized_callable_loci: if defined(harmonized_callable_loci) then harmonized_callable_loci else sample.harmonized_callable_loci,
        harmonized_denoised_total_copy_ratios: if defined(harmonized_denoised_total_copy_ratios) then harmonized_denoised_total_copy_ratios else sample.harmonized_denoised_total_copy_ratios,
        harmonized_snppanel_allelic_pileup_summaries: if defined(harmonized_snppanel_allelic_pileup_summaries) then harmonized_snppanel_allelic_pileup_summaries else sample.harmonized_snppanel_allelic_pileup_summaries,
        contamination_table: if defined(contamination_table) then contamination_table else sample.contamination_table,
        af_segmentation_table: if defined(af_segmentation_table) then af_segmentation_table else sample.af_segmentation_table,
        allelic_pileup_summaries: if defined(allelic_pileup_summaries) then allelic_pileup_summaries else sample.allelic_pileup_summaries,
        aggregated_allelic_read_counts: if defined(aggregated_allelic_read_counts) then aggregated_allelic_read_counts else sample.aggregated_allelic_read_counts,
        genotype_error_probabilities: if defined(genotype_error_probabilities) then genotype_error_probabilities else sample.genotype_error_probabilities,
        af_model_parameters: if defined(af_model_parameters) then af_model_parameters else sample.af_model_parameters,
        cr_model_parameters: if defined(cr_model_parameters) then cr_model_parameters else sample.cr_model_parameters,
        called_copy_ratio_segmentation: if defined(called_copy_ratio_segmentation) then called_copy_ratio_segmentation else sample.called_copy_ratio_segmentation,
        cr_plot: if defined(cr_plot) then cr_plot else sample.cr_plot,
        acs_copy_ratio_segmentation: if defined(acs_copy_ratio_segmentation) then acs_copy_ratio_segmentation else sample.acs_copy_ratio_segmentation,
        acs_copy_ratio_skew: if defined(acs_copy_ratio_skew) then acs_copy_ratio_skew else sample.acs_copy_ratio_skew,
        annotated_somatic_variants: if defined(annotated_somatic_variants) then annotated_somatic_variants else sample.annotated_somatic_variants,
        annotated_somatic_variants_idx: if defined(annotated_somatic_variants_idx) then annotated_somatic_variants_idx else sample.annotated_somatic_variants_idx,
        absolute_acr_rdata: if defined(absolute_acr_rdata) then absolute_acr_rdata else sample.absolute_acr_rdata,
        absolute_acr_plot: if defined(absolute_acr_plot) then absolute_acr_plot else sample.absolute_acr_plot,
        absolute_snv_maf: if defined(absolute_snv_maf) then absolute_snv_maf else sample.absolute_snv_maf,
        absolute_indel_maf: if defined(absolute_indel_maf) then absolute_indel_maf else sample.absolute_indel_maf,
        absolute_solution: if defined(absolute_solution) then absolute_solution else sample.absolute_solution,
        absolute_maf: if defined(absolute_maf) then absolute_maf else sample.absolute_maf,
        absolute_segtab: if defined(absolute_segtab) then absolute_segtab else sample.absolute_segtab,
        absolute_maf_postprocessed: if defined(absolute_maf_postprocessed) then absolute_maf_postprocessed else sample.absolute_maf_postprocessed,
        absolute_segtab_postprocessed: if defined(absolute_segtab_postprocessed) then absolute_segtab_postprocessed else sample.absolute_segtab_postprocessed,
        absolute_segtab_igv_postprocessed: if defined(absolute_segtab_igv_postprocessed) then absolute_segtab_igv_postprocessed else sample.absolute_segtab_igv_postprocessed,
        absolute_table: if defined(absolute_table) then absolute_table else sample.absolute_table,
        purity: if defined(purity) then purity else sample.purity,
        ploidy: if defined(ploidy) then ploidy else sample.ploidy,
        timepoint: if defined(timepoint) then timepoint else sample.timepoint
    }

    output {
        Sample updated_sample = s
    }
}
