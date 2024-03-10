version development


struct SequencingRun {
    String name
    File bam
    File bai
    File target_intervals
    File? annotated_target_intervals
    File? cnv_panel_of_normals
    Boolean? is_paired_end
    Boolean use_for_dCR
    Boolean use_for_aCR
}


workflow UpdateSequencingRun {
    input {
        SequencingRun sequencing_run
        String? name
        File? bam
        File? bai
        File? target_intervals
        File? annotated_target_intervals
        File? cnv_panel_of_normals
        Boolean? is_paired_end
        Boolean? use_for_dCR
        Boolean? use_for_aCR
    }

    SequencingRun seq_run = object {
        name: select_first([name, sequencing_run.name]),
        bam: select_first([bam, sequencing_run.bam]),
        bai: select_first([bai, sequencing_run.bai]),
        target_intervals: select_first([target_intervals, sequencing_run.target_intervals]),
        # cannot use select_first for optional fields:
        annotated_target_intervals: if defined(annotated_target_intervals) then annotated_target_intervals else sequencing_run.annotated_target_intervals,
        cnv_panel_of_normals: if defined(cnv_panel_of_normals) then cnv_panel_of_normals else sequencing_run.cnv_panel_of_normals,
        is_paired_end: if defined(is_paired_end) then is_paired_end else sequencing_run.is_paired_end,
        use_for_dCR: select_first([use_for_dCR, sequencing_run.use_for_dCR]),
        use_for_aCR: select_first([use_for_aCR, sequencing_run.use_for_aCR])
    }

    output {
        SequencingRun updated_sequencing_run = seq_run
    }
}