version development

import "runtime_collection.wdl" as rtc
import "sequencing_run.wdl" as seqrun
import "tasks.wdl"


workflow DefineSequencingRun {
    input {
        String? name
        String? sample_name
        File bam
        File bai
        File target_intervals
        File? annotated_target_intervals
        File? cnv_panel_of_normals
        Boolean? is_paired_end
        Boolean use_for_dCR = true
        Boolean use_for_aCR = true

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
    # Replace "object" with name of the struct in the future.
    SequencingRun seq_run = object {
        name: select_first([name, GetSampleName.sample_name]),
        sample_name: select_first([sample_name, GetSampleName.sample_name]),
        bam: bam,
        bai: bai,
        target_intervals: target_intervals,
        annotated_target_intervals: annotated_target_intervals,
        cnv_panel_of_normals: cnv_panel_of_normals,
        is_paired_end: is_paired_end,
        use_for_dCR: use_for_dCR,
        use_for_aCR: use_for_aCR,
    }

    output {
        SequencingRun sequencing_run = seq_run
    }
}