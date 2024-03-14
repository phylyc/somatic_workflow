version development

import "patient.wdl" as p
import "patient.update_samples.wdl" as p_update_s
import "workflow_arguments.wdl" as wfargs
import "runtime_collection.wdl" as rtc
import "runtimes.wdl" as rt


workflow FilterSegments {
    input {
        Patient patient
        WorkflowArguments args
        RuntimeCollection runtime_collection
    }

    if (defined(patient.matched_normal_sample)) {
        scatter (sample in patient.tumor_samples) {
            Sample matched_normal_sample = select_first([patient.matched_normal_sample])
            if (defined(sample.called_copy_ratio_segmentation) && defined(matched_normal_sample.called_copy_ratio_segmentation)) {
                call FilterGermlineCNVs {
                    input:
                        script = args.filter_germline_cnvs_script,
                        tumor_called_copy_ratio_segmentation = select_first([sample.called_copy_ratio_segmentation]),
                        normal_called_copy_ratio_segmentation = select_first([matched_normal_sample.called_copy_ratio_segmentation]),
                        runtime_params = runtime_collection.filter_germline_cnvs
                }
            }
        }

        scatter (sample in patient.normal_samples) {
            File? normal_called_copy_ratio_segmentation = sample.called_copy_ratio_segmentation
        }

        call p_update_s.UpdateSamples as UpdateSegmentations {
            input:
                patient = patient,
                called_copy_ratio_segmentations = flatten([
                    select_all(FilterGermlineCNVs.filtered_called_copy_ratio_segmentation),
                    select_all(normal_called_copy_ratio_segmentation)
                ])
        }
    }

    output {
        Patient updated_patient = select_first([UpdateSegmentations.updated_patient, patient])
    }
}

task FilterGermlineCNVs {
    input {
        String script = "https://github.com/phylyc/somatic_workflow/raw/master/python/filter_germline_cnvs.py"

        File tumor_called_copy_ratio_segmentation
        File normal_called_copy_ratio_segmentation

        Int min_segment_length = 100

        Boolean verbose = true

        Runtime runtime_params
    }

    String filtered_called_copy_ratio_segmentation = basename(tumor_called_copy_ratio_segmentation, ".seg") + ".germline_filtered.seg"

    command <<<
        set -euxo pipefail
        wget -O filter_germline_cnvs.py ~{script}
        python filter_germline_cnvs.py \
            --tumor_called_copy_ratio_segmentation '~{tumor_called_copy_ratio_segmentation}' \
            --normal_called_copy_ratio_segmentation '~{normal_called_copy_ratio_segmentation}' \
            --output '~{filtered_called_copy_ratio_segmentation}' \
            --min_segment_length ~{min_segment_length} \
            ~{if verbose then "--verbose" else ""}
    >>>

    output {
        File filtered_called_copy_ratio_segmentation = filtered_called_copy_ratio_segmentation
    }

    runtime {
        docker: runtime_params.docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: runtime_params.machine_mem + " MB"
        runtime_minutes: runtime_params.runtime_minutes
        disks: "local-disk " + runtime_params.disk + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }
}