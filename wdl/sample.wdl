version development

import "sequencing_run.wdl"
import "runtimes.wdl"


struct Sample {
    String name
    String bam_name
    Array[SequencingRun] sequencing_runs
    Boolean is_tumor

    File? read_counts
    File? denoised_copy_ratios
    File? standardized_copy_ratios
    File? snp_array_allelic_counts
    File? somatic_allelic_counts
    File? germline_allelic_counts
    File? contamination
    File? af_segmentation
    File? copy_ratio_segmentation
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
        File? snp_array_allelic_counts
        File? somatic_allelic_counts
        File? germline_allelic_counts
        File? contamination
        File? af_segmentation
        File? copy_ratio_segmentation
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
        snp_array_allelic_counts: if defined(snp_array_allelic_counts) then snp_array_allelic_counts else sample.snp_array_allelic_counts,
        somatic_allelic_counts: if defined(somatic_allelic_counts) then somatic_allelic_counts else sample.somatic_allelic_counts,
        germline_allelic_counts: if defined(germline_allelic_counts) then germline_allelic_counts else sample.germline_allelic_counts,
        contamination: if defined(contamination) then contamination else sample.contamination,
        af_segmentation: if defined(af_segmentation) then af_segmentation else sample.af_segmentation,
        copy_ratio_segmentation: if defined(copy_ratio_segmentation) then copy_ratio_segmentation else sample.copy_ratio_segmentation
    }

    output {
        Sample updated_sample = s
    }
}
