version development

import "runtime_collection.wdl" as rtc
import "sequencing_run.wdl" as seqrun
import "tasks.wdl"


workflow DefineSequencingRun {
    input {
        String? name
        String? sample_name
        Int timepoint = 0
        File bam
        File bai
        File target_intervals
        File? annotated_target_intervals
        File? cnv_panel_of_normals
        Boolean? is_paired_end
        Boolean use_for_tCR = true
        Boolean use_for_aCR = true

        # CACHE
        File? callable_loci
        File? total_read_counts
        File? denoised_total_copy_ratios
        File? snppanel_allelic_pileup_summaries
        File? rare_germline_allelic_pileup_summaries

        RuntimeCollection runtime_collection
    }

    if (!defined(name)) {
        call tasks.GetSampleName {
            input:
                bam = bam,
                bai = bai,
                runtime_params = runtime_collection.get_sample_name,
        }
    }

    # The object syntax is likely going to be deprecated in newer WDL versions!
    # Replace "object" with name of the struct in the future.
    SequencingRun seq_run = object {
        name: select_first([name, GetSampleName.sample_name]),
        sample_name: select_first([sample_name, GetSampleName.sample_name]),
        timepoint: timepoint,
        bam: bam,
        bai: bai,
        target_intervals: target_intervals,
        annotated_target_intervals: annotated_target_intervals,
        cnv_panel_of_normals: cnv_panel_of_normals,
        is_paired_end: is_paired_end,
        use_for_tCR: use_for_tCR,
        use_for_aCR: use_for_aCR,
        callable_loci: callable_loci,
        total_read_counts: total_read_counts,
        denoised_total_copy_ratios: denoised_total_copy_ratios,
        snppanel_allelic_pileup_summaries: snppanel_allelic_pileup_summaries,
        rare_germline_allelic_pileup_summaries: rare_germline_allelic_pileup_summaries,
    }

    output {
        SequencingRun sequencing_run = seq_run
    }
}