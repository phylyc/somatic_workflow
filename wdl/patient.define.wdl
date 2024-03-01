version development

import "sequencing_run.wdl" as seq_run
import "sequencing_run.define.wdl" as seq_run_def
import "sample.wdl" as s
import "patient.wdl" as p
import "runtime_collection.wdl" as rtc


workflow DefinePatient {
    input {
        String individual_id

        Array[File]+ tumor_bams
        Array[File]+ tumor_bais
        Array[File]+ tumor_target_intervals
        Array[File]? tumor_annotated_target_intervals
        Array[File]? tumor_cnv_panel_of_normals
        Array[String]? tumor_sample_names

        Array[File]? normal_bams
        Array[File]? normal_bais
        Array[File]? normal_target_intervals
        Array[File]? normal_annotated_target_intervals
        Array[File]? normal_cnv_panel_of_normals
        Array[String]? normal_sample_names

        RuntimeCollection runtime_collection
    }

    Array[File] non_optional_normal_bams = select_first([normal_bams, []])
    Array[File] non_optional_normal_bais = select_first([normal_bais, []])
    Array[File] non_optional_normal_target_intervals = select_first([normal_target_intervals, []])

    Boolean has_tumor = length(tumor_bams) > 0
    Boolean has_normal = length(non_optional_normal_bams) > 0

    # We first define SequencingRuns for each tumor and normal bam, and then group them by sample name into Samples.

    scatter (tuple in transpose([tumor_bams, tumor_bais, tumor_target_intervals])) {
        call seq_run_def.DefineSequencingRun as DefineTumorSequencingRun {
            input:
                bam = tuple[0],
                bai = tuple[1],
                target_intervals = tuple[2],
                runtime_collection = runtime_collection
        }
        String tumor_bam_names = DefineTumorSequencingRun.sequencing_run.name
    }
    Array[SequencingRun] tumors_1 = DefineTumorSequencingRun.sequencing_run
    if (defined(tumor_annotated_target_intervals)) {
        scatter (pair in zip(tumors_1, select_first([tumor_annotated_target_intervals, []]))) {
            call seq_run.UpdateSequencingRun as UpdateAnnotatedTargetIntervalsTumorSeq {
                input:
                    sequencing_run = pair.left,
                    annotated_target_intervals = pair.right,
            }
        }
    }
    Array[SequencingRun] tumors_2 = select_first([UpdateAnnotatedTargetIntervalsTumorSeq.updated_sequencing_run, tumors_1])
    if (defined(tumor_cnv_panel_of_normals)) {
        scatter (pair in zip(tumors_2, select_first([tumor_cnv_panel_of_normals, []]))) {
            if (size(pair.right) > 0) {
                # For some sequencing platforms a panel of normals may not be available.
                # The denoise read counts task will then just use the anntated target
                # intervals to do GC correction.
                File t_cnv_panel_of_normals = pair.right
            }
            call seq_run.UpdateSequencingRun as UpdateCnvPanelOfNormalsTumorSeq {
                input:
                    sequencing_run = pair.left,
                    cnv_panel_of_normals = t_cnv_panel_of_normals,
            }
        }
    }
    Array[SequencingRun] tumors_3 = select_first([UpdateCnvPanelOfNormalsTumorSeq.updated_sequencing_run, tumors_2])

    if (has_normal) {
        scatter (tuple in transpose([non_optional_normal_bams, non_optional_normal_bais, non_optional_normal_target_intervals])) {
            call seq_run_def.DefineSequencingRun as DefineNormalSequencingRun {
                input:
                    bam = tuple[0],
                    bai = tuple[1],
                    target_intervals = tuple[2],
                    runtime_collection = runtime_collection
            }
            String normal_bam_names = DefineNormalSequencingRun.sequencing_run.name
        }
    }
    Array[SequencingRun] normals_1 = select_first([DefineNormalSequencingRun.sequencing_run, []])
    if (defined(normal_annotated_target_intervals)) {
        scatter (pair in zip(normals_1, select_first([normal_annotated_target_intervals, []]))) {
            call seq_run.UpdateSequencingRun as UpdateAnnotatedTargetIntervalsNormalSeq {
                input:
                    sequencing_run = pair.left,
                    annotated_target_intervals = pair.right,
            }
        }
    }
    Array[SequencingRun] normals_2 = select_first([UpdateAnnotatedTargetIntervalsNormalSeq.updated_sequencing_run, normals_1])
    if (defined(normal_cnv_panel_of_normals)) {
        scatter (pair in zip(normals_2, select_first([normal_cnv_panel_of_normals, []]))) {
            if (size(pair.right) > 0) {
                # For some sequencing platforms a panel of normals may not be available.
                # The denoise read counts task will then just use the anntated target
                # intervals to do GC correction.
                File n_cnv_panel_of_normals = pair.right
            }
            call seq_run.UpdateSequencingRun as UpdateCnvPanelOfNormalsNormalSeq {
                input:
                    sequencing_run = pair.left,
                    cnv_panel_of_normals = n_cnv_panel_of_normals,
            }
        }
    }
    Array[SequencingRun] normals_3 = select_first([UpdateCnvPanelOfNormalsNormalSeq.updated_sequencing_run, normals_2])

    # GroupBy sample name:
    # We assume that tumor_sample_names and tumor_bam_names share the same uniqueness,
    # that is if the supplied sample name is the same for two input bams, then the
    # bam names should also be the same, and vice versa.

    Array[String] tumor_names = select_first([tumor_sample_names, tumor_bam_names])
    scatter (pair in as_pairs(collect_by_key(zip(tumor_names, tumors_3)))) {
        Sample tumor_samples = object {
            name: pair.left,
            bam_name: pair.right[0].name,
            sequencing_runs: pair.right,
            is_tumor: true,
        }
    }
    Array[String] normal_names = select_first([normal_sample_names, normal_bam_names])
    scatter (pair in as_pairs(collect_by_key(zip(normal_names, normals_3)))) {
        Sample normal_samples = object {
            name: pair.left,
            bam_name: pair.right[0].name,
            sequencing_runs: pair.right,
            is_tumor: false,
        }
    }

    if (has_normal) {
        # We select the first normal sample to be the matched normal.
        # todo: select normal with greatest sequencing depth
        Sample best_matched_normal_sample = select_first(normal_samples)
    }

    Patient pat = object {
        name: individual_id,
        samples: flatten([tumor_samples, normal_samples]),
        tumor_samples: tumor_samples,
        normal_samples: normal_samples,
        has_tumor: has_tumor,
        has_normal: has_normal,
        matched_normal_sample: best_matched_normal_sample,
    }

    output {
        Patient patient = pat
    }
}