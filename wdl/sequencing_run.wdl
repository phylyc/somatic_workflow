version development

import "runtimes.wdl"
import "tasks.wdl"


struct SequencingRun {
    String name
    File bam
    File bai
    File target_intervals
    File? annotated_target_intervals
    File? cnv_panel_of_normals
    Boolean? is_paired_end
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
    }

    SequencingRun seq_run = object {
        name: select_first([name, sequencing_run.name]),
        bam: select_first([bam, sequencing_run.bam]),
        bai: select_first([bai, sequencing_run.bai]),
        target_intervals: select_first([target_intervals, sequencing_run.target_intervals]),
        # cannot use select_first for optional fields:
        annotated_target_intervals: if defined(annotated_target_intervals) then annotated_target_intervals else sequencing_run.annotated_target_intervals,
        cnv_panel_of_normals: if defined(cnv_panel_of_normals) then cnv_panel_of_normals else sequencing_run.cnv_panel_of_normals,
        is_paired_end: if defined(is_paired_end) then is_paired_end else sequencing_run.is_paired_end
    }

    output {
        SequencingRun updated_sequencing_run = seq_run
    }
}


workflow DefineSequencingRun {
    input {
        String? name
        File bam
        File bai
        File target_intervals
        File? annotated_target_intervals
        File? cnv_panel_of_normals
        Boolean? is_paired_end

        RuntimeCollection runtime_collection
    }

    if (!defined(name)) {
        call tasks.GetSampleName {
            input:
                bam = bam,
                runtime_params = runtime_collection.get_sample_name,
        }
    }

    # The object syntax is likely going to be deprecated in newer WDL versions!
    # Definition via {"bam": sample.left.left, ..., "is_tumor": true} leads to
    # Error(s): No coercion defined from wom value(s) '"true"' of type 'String' to 'Boolean'.
    # :(
    SequencingRun seq_run = object {
        name: select_first([name, GetSampleName.sample_name]),
        bam: bam,
        bai: bai,
        target_intervals: target_intervals,
        annotated_target_intervals: annotated_target_intervals,
        cnv_panel_of_normals: cnv_panel_of_normals,
        is_paired_end: is_paired_end
    }

    output {
        SequencingRun sequencing_run = seq_run
    }
}