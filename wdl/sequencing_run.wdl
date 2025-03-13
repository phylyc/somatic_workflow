version development


struct SequencingRun {
    String name
    String sample_name
    File bam
    File bai
    File target_intervals
    File? annotated_target_intervals
    File? cnv_panel_of_normals
    Boolean? is_paired_end
    Boolean use_for_tCR
    Boolean use_for_aCR

    # CACHE
    File? callable_loci
    File? total_read_counts
    File? denoised_total_copy_ratios
    File? snppanel_allelic_pileup_summaries
}


workflow UpdateSequencingRun {
    input {
        SequencingRun sequencing_run
        String? name
        String? sample_name
        File? bam
        File? bai
        File? target_intervals
        File? annotated_target_intervals
        File? cnv_panel_of_normals
        Boolean? is_paired_end
        Boolean? use_for_tCR
        Boolean? use_for_aCR

        File? callable_loci
        File? total_read_counts
        File? denoised_total_copy_ratios
        File? snppanel_allelic_pileup_summaries
    }

    SequencingRun seq_run = object {
        name: select_first([name, sequencing_run.name]),
        sample_name: select_first([sample_name, sequencing_run.sample_name]),
        bam: select_first([bam, sequencing_run.bam]),
        bai: select_first([bai, sequencing_run.bai]),
        target_intervals: select_first([target_intervals, sequencing_run.target_intervals]),
        # cannot use select_first for optional fields:
        annotated_target_intervals: if defined(annotated_target_intervals) then annotated_target_intervals else sequencing_run.annotated_target_intervals,
        cnv_panel_of_normals: if defined(cnv_panel_of_normals) then cnv_panel_of_normals else sequencing_run.cnv_panel_of_normals,
        is_paired_end: if defined(is_paired_end) then is_paired_end else sequencing_run.is_paired_end,
        use_for_tCR: select_first([use_for_tCR, sequencing_run.use_for_tCR]),
        use_for_aCR: select_first([use_for_aCR, sequencing_run.use_for_aCR]),
        callable_loci: if defined(callable_loci) then callable_loci else sequencing_run.callable_loci,
        total_read_counts: if defined(total_read_counts) then total_read_counts else sequencing_run.total_read_counts,
        denoised_total_copy_ratios: if defined(denoised_total_copy_ratios) then denoised_total_copy_ratios else sequencing_run.denoised_total_copy_ratios,
        snppanel_allelic_pileup_summaries: if defined(snppanel_allelic_pileup_summaries) then snppanel_allelic_pileup_summaries else sequencing_run.snppanel_allelic_pileup_summaries,
    }

    output {
        SequencingRun updated_sequencing_run = seq_run
    }
}