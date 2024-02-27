version development

import "runtime_collection.wdl"
import "sequencing_run.wdl"
import "tasks.wdl"


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