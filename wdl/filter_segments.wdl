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

    scatter (sample in patient.samples) {
        String sample_name = sample.name
        File? called_copy_ratio_segmentation = sample.called_copy_ratio_segmentation
    }
    Array[File] sample_names = select_all(sample_name)
    Array[File] copy_ratios = select_all(called_copy_ratio_segmentation)

    if (defined(patient.matched_normal_sample)) {
        Sample matched_normal_sample = select_first([patient.matched_normal_sample])
        call FilterGermlineCNVs {
            input:
                script = args.harmonize_copy_ratios_script,
                sample_names = sample_names,
                normal_sample = matched_normal_sample.name,
                copy_ratios = copy_ratios,
                compress_output = args.compress_output,
                min_segment_length = args.filter_germline_cnvs_min_segment_length,
                suffix = ".modelFinal.called.germline_filtered.seg",
                runtime_params = runtime_collection.harmonize_copy_ratios
        }

        # sort output to match order of sample_names since glob doesn't guarantee order
        scatter (sample in patient.samples) {
            scatter (fccr in FilterGermlineCNVs.filtered_called_copy_ratio) {
                String this_sample_name = basename(fccr, ".modelFinal.called.germline_filtered.seg")
                if (sample.name == this_sample_name) {
                    File this_cr = fccr
                }
            }
            Array[File] this_sample_cr = select_all(this_cr)
        }
        Array[File] sorted_filtered_copy_ratios = flatten(this_sample_cr)

        call p_update_s.UpdateSamples as UpdateSegmentations {
            input:
                patient = patient,
                called_copy_ratio_segmentations = sorted_filtered_copy_ratios
        }
    }

    output {
        Patient updated_patient = select_first([UpdateSegmentations.updated_patient, patient])
        Array[File]? filtered_called_copy_ratio_segmentations = sorted_filtered_copy_ratios
    }
}

task FilterGermlineCNVs {
    input {
        String script = "https://github.com/phylyc/somatic_workflow/raw/master/python/harmonize_copy_ratios.py"

        Array[String]+ sample_names
        Array[File]+ copy_ratios
        String? normal_sample
        Int min_target_length = 100
        Boolean compress_output = false
        Boolean verbose = true

        Int min_segment_length = 100

        Boolean compress_output = false
        Boolean verbose = true

        String suffix = ".germline_filtered.seg"

        Runtime runtime_params
    }

    String output_dir = "."

    command <<<
        set -euxo pipefail
        n_threads=$(nproc)
        wget -O harmonize_copy_ratios.py ~{script}
        python harmonize_copy_ratios.py \
            --output_dir '~{output_dir}' \
            ~{sep="' " prefix("--sample '", sample_names)}' \
            ~{sep="' " prefix("--copy_ratio '", copy_ratios)}' \
            ~{"--normal_sample '" + normal_sample + "'"} \
            --suffix ".germline_filtered.seg" \
            --column_names \
                CONTIG START END \
                NUM_POINTS_COPY_RATIO NUM_POINTS_ALLELE_FRACTION \
                LOG2_COPY_RATIO_POSTERIOR_10 LOG2_COPY_RATIO_POSTERIOR_50 LOG2_COPY_RATIO_POSTERIOR_90 \
                MINOR_ALLELE_FRACTION_POSTERIOR_10 MINOR_ALLELE_FRACTION_POSTERIOR_50 MINOR_ALLELE_FRACTION_POSTERIOR_90 \
                CALL \
            --column_types \
                str int int \
                int int \
                float float float \
                float float float \
                str \
            --agg_col NUM_POINTS_COPY_RATIO \
            --agg_func mean \
            --agg_col NUM_POINTS_ALLELE_FRACTION \
            --agg_func mean \
            --agg_col LOG2_COPY_RATIO_POSTERIOR_10 \
            --agg_func mean \
            --agg_col LOG2_COPY_RATIO_POSTERIOR_50 \
            --agg_func mean \
            --agg_col LOG2_COPY_RATIO_POSTERIOR_90 \
            --agg_func mean \
            --agg_col MINOR_ALLELE_FRACTION_POSTERIOR_10 \
            --agg_func mean \
            --agg_col MINOR_ALLELE_FRACTION_POSTERIOR_50 \
            --agg_func mean \
            --agg_col MINOR_ALLELE_FRACTION_POSTERIOR_90 \
            --agg_func mean \
            --agg_col CALL \
            --agg_func first \
            --filter_germline_calls \
            --threads $n_threads \
            --min_target_length ~{min_segment_length} \
            ~{if compress_output then "--compress_output" else ""} \
            ~{if verbose then "--verbose" else ""}
    >>>

    output {
        Array[File] filtered_called_copy_ratio = glob("*" + suffix)
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