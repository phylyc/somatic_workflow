version development

import "sequencing_run.wdl" as seqrun


struct Sample {
    String name
    String bam_name
    Array[SequencingRun] sequencing_runs
    Boolean is_tumor

    File? read_counts                       # GATK CollectReadCounts
    File? denoised_copy_ratios              # GATK DenoiseReadCounts
    File? standardized_copy_ratios          # GATK DenoiseReadCounts
    File? snp_array_pileups                 # GATK GetPileupSummaries
    File? snp_array_allelic_counts          # PileupToAllelicCounts / GATK CollectAllelicCounts
    Float? genotype_error_probabilities     # PileupToAllelicCounts
    File? somatic_allelic_counts            # GATK GetPileupSummaries
    File? germline_allelic_counts           # GATK GetPileupSummaries
    File? contamination                     # GATK CalculateContamination
    File? af_segmentation                   # GATK CalculateContamination
    File? copy_ratio_segmentation           # GATK ModelSegments
    File? af_model_parameters               # GATK ModelSegments
    File? cr_model_parameters               # GATK ModelSegments
    File? called_copy_ratio_segmentation    # GATK CallCopyRatioSegments
    File? acs_copy_ratio_segmentation       # ModelSegmentsToACSConversion
    Float? acs_copy_ratio_skew              # ModelSegmentsToACSConversion
    File? annotated_variants                # GATK Funcotator
}


workflow UpdateSample {
    input {
        Sample sample
        String? name
        String? bam_name
        Array[SequencingRun]? sequencing_runs
        Boolean? is_tumor

        File? read_counts
        File? denoised_copy_ratios
        File? standardized_copy_ratios
        File? snp_array_pileups
        File? snp_array_allelic_counts
        Float? genotype_error_probabilities
        File? somatic_allelic_counts
        File? germline_allelic_counts
        File? contamination
        File? af_segmentation
        File? copy_ratio_segmentation
        File? af_model_parameters
        File? cr_model_parameters
        File? called_copy_ratio_segmentation
        File? acs_copy_ratio_segmentation
        Float? acs_copy_ratio_skew
        File? annotated_variants
    }

    Sample s = object {
        name: select_first([name, sample.name]),
        bam_name: select_first([bam_name, sample.bam_name]),
        sequencing_runs: select_first([sequencing_runs, sample.sequencing_runs]),
        is_tumor: select_first([is_tumor, sample.is_tumor]),
        # cannot use select_first for optional fields:
        read_counts: if defined(read_counts) then read_counts else sample.read_counts,
        denoised_copy_ratios: if defined(denoised_copy_ratios) then denoised_copy_ratios else sample.denoised_copy_ratios,
        standardized_copy_ratios: if defined(standardized_copy_ratios) then standardized_copy_ratios else sample.standardized_copy_ratios,
        snp_array_pileups: if defined(snp_array_pileups) then snp_array_pileups else sample.snp_array_pileups,
        snp_array_allelic_counts: if defined(snp_array_allelic_counts) then snp_array_allelic_counts else sample.snp_array_allelic_counts,
        genotype_error_probabilities: if defined(genotype_error_probabilities) then genotype_error_probabilities else sample.genotype_error_probabilities,
        somatic_allelic_counts: if defined(somatic_allelic_counts) then somatic_allelic_counts else sample.somatic_allelic_counts,
        germline_allelic_counts: if defined(germline_allelic_counts) then germline_allelic_counts else sample.germline_allelic_counts,
        contamination: if defined(contamination) then contamination else sample.contamination,
        af_segmentation: if defined(af_segmentation) then af_segmentation else sample.af_segmentation,
        copy_ratio_segmentation: if defined(copy_ratio_segmentation) then copy_ratio_segmentation else sample.copy_ratio_segmentation,
        af_model_parameters: if defined(af_model_parameters) then af_model_parameters else sample.af_model_parameters,
        cr_model_parameters: if defined(cr_model_parameters) then cr_model_parameters else sample.cr_model_parameters,
        called_copy_ratio_segmentation: if defined(called_copy_ratio_segmentation) then called_copy_ratio_segmentation else sample.called_copy_ratio_segmentation,
        acs_copy_ratio_segmentation: if defined(acs_copy_ratio_segmentation) then acs_copy_ratio_segmentation else sample.acs_copy_ratio_segmentation,
        acs_copy_ratio_skew: if defined(acs_copy_ratio_skew) then acs_copy_ratio_skew else sample.acs_copy_ratio_skew,
        annotated_variants: if defined(annotated_variants) then annotated_variants else sample.annotated_variants
    }

    output {
        Sample updated_sample = s
    }
}
