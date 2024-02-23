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


#workflow DefineSample {
#    input {
#        String? name
#        String? bam_name
#        Array[File] bams
#        Array[File] bais
#        Array[File] target_intervals
#        Array[File]? annotated_target_intervals
#        Array[File]? cnv_panel_of_normals
#        Array[Boolean]? is_paired_end
#        Boolean is_tumor
#
#        RuntimeCollection runtime_collection
#    }
#
#    # Define the sequencing runs with minimal required information:
#
#    scatter (sequencing_run in transpose([bams, bais, target_intervals])) {
#        call sequencing_run.DefineSequencingRun {
#            input:
#                name = bam_name,
#                bam = sequencing_run[0],
#                bai = sequencing_run[1],
#                target_intervals = sequencing_run[2],
#                runtime_collection = runtime_collection,
#        }
#    }
#    Array[SequencingRun] seq_runs_1 = DefineSequencingRun.sequencing_run
#
#    String this_bam_name = seq_runs_1[0].name  # This is the same as bam_name if defined.
#
#    # Update the sequencing runs with optional information:
#
#    if (defined(annotated_target_intervals)) {
#        scatter (pair in zip(seq_runs_1, select_first([annotated_target_intervals, []]))) {
#            call sequencing_run.UpdateSequencingRun as AddAnnotatedTargetIntervals {
#                input:
#                    sequencing_run = pair.left,
#                    annotated_target_intervals = pair.right
#            }
#        }
#    }
#    Array[SequencingRun] seq_runs_2 = select_first([AddAnnotatedTargetIntervals.updated_sequencing_run, seq_runs_1])
#
#    if (defined(cnv_panel_of_normals)) {
#        scatter (pair in zip(seq_runs_2, select_first([cnv_panel_of_normals, []]))) {
#            call sequencing_run.UpdateSequencingRun as AddCnvPanelOfNormals {
#                input:
#                    sequencing_run = pair.left,
#                    cnv_panel_of_normals = pair.right
#            }
#        }
#    }
#    Array[SequencingRun] seq_runs_3 = select_first([AddCnvPanelOfNormals.updated_sequencing_run, seq_runs_2])
#
#    if (defined(is_paired_end)) {
#        scatter (pair in zip(DefineSequencingRun.sequencing_run, select_first([is_paired_end, []]))) {
#            call sequencing_run.UpdateSequencingRun as AddPairedEnd {
#                input:
#                    sequencing_run = pair.left,
#                    is_paired_end = pair.right
#            }
#        }
#    }
#    Array[SequencingRun] seq_runs_4 = select_first([AddPairedEnd.updated_sequencing_run, seq_runs_3])
#
#    Sample s = object {
#        name: select_first([name, this_bam_name]),
#        bam_name: this_bam_name,
#        sequencing_runs: seq_runs_4,
#        is_tumor: is_tumor
#    }
#
#    output {
#        Sample sample = s
#    }
#}
