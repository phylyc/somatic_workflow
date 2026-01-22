version development

import "sample.wdl" as s
import "patient.wdl" as p
import "patient.update_samples.wdl" as p_update_s
import "harmonize_samples.wdl" as hs
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
        if (defined(sample.harmonized_denoised_total_copy_ratios) && defined(sample.called_copy_ratio_segmentation)) {
            call FilterCopyRatios {
                input:
                    copy_ratios = select_first([sample.harmonized_denoised_total_copy_ratios]),
                    segmentations = select_first([sample.called_copy_ratio_segmentation]),
                    min_probes = args.filter_segments_min_probes,
                    runtime_params = runtime_collection.filter_copy_ratios
            }
        }

        call s.UpdateSample as AddFilteredCopyRatios {
            input:
                sample = sample,
                harmonized_denoised_total_copy_ratios = FilterCopyRatios.filtered_denoised_copy_ratios
        }

        String sample_dcr_names = sample.name
        File? sample_dcr = FilterCopyRatios.filtered_denoised_copy_ratios
    }

    if (length(select_all(sample_dcr)) > 0) {
        call hs.HarmonizeCopyRatios {
            input:
                script = args.script_harmonize_copy_ratios,
                ref_dict = args.files.ref_dict,
                sample_names = sample_dcr_names,
                denoised_copy_ratios = select_all(sample_dcr),
                min_target_length = args.harmonize_min_target_length,
                compress_output = false,
                runtime_params = runtime_collection.harmonize_copy_ratios
        }

        # sort output to match order of sample_names since glob doesn't guarantee order
        scatter (sample in AddFilteredCopyRatios.updated_sample) {
            scatter (h_dcr in HarmonizeCopyRatios.harmonized_denoised_copy_ratios) {
                String this_dcr_sample_name = basename(basename(basename(basename(h_dcr, ".hdf5"), ".tsv"), ".denoised_CR"), ".harmonized")
                if (sample.name == this_dcr_sample_name) {
                    File this_dcr = h_dcr
                }
            }
            Array[File] this_sample_dcr = select_all(this_dcr)
        }
        Array[File] sorted_harmonized_denoised_copy_ratios = flatten(this_sample_dcr)

        call p_update_s.UpdateSamples as UpdateCopyRatios {
            input:
                patient = patient,
                harmonized_denoised_total_copy_ratios = sorted_harmonized_denoised_copy_ratios
        }
    }

    output {
        Patient updated_patient = select_first([UpdateCopyRatios.updated_patient, patient])
    }
}

task FilterCopyRatios {
    input {
        File copy_ratios
        File segmentations
        Int min_probes = 0
        Int max_num_segments = 1000

        Runtime runtime_params
    }

    String filtered_copy_ratios = basename(copy_ratios, ".tsv") + ".filtered.tsv"

    command <<<
        set -euxo pipefail

        num_segments=$(( $(grep -v "^@" "~{segmentations}" | wc -l ) - 1 ))

        if [[ $num_segments -gt ~{max_num_segments} ]]; then
            cp ~{copy_ratios} ~{filtered_copy_ratios}
            echo "This sample has $num_segments segments and is likely oversegmented (> ~{max_num_segments}). No filtering performed."
            exit 0
        fi

        echo ""

        awk -F'\t' -v OFS='\t' -v min_probes="~{min_probes}" '
          /^@/ || $1=="CONTIG" { next }
          ($4+0) < min_probes { print $1, $2-1, $3 }
        ' "~{segmentations}" > "bad_segments.bed"

        awk -F'\t' -v OFS='\t' '
          /^@/ || $1=="CONTIG" { next }
          { print $1, $2-1, $3, $0 }    # bed start=end-1; payload=$0
        ' "~{copy_ratios}" > "copy_ratios.bed"

        num_bad_segments=$(wc -l < "bad_segments.bed")
        echo "Removing denoised read count observations in $num_bad_segments segments with < ~{min_probes} markers."

        # copy header
        awk '/^@/ || /^CONTIG/{print}' "~{copy_ratios}" \
          > "~{filtered_copy_ratios}"

        # filter copy ratios
        if [ -s "bad_segments.bed" ]; then
          bedtools intersect \
            -a "copy_ratios.bed" \
            -b "bad_segments.bed" \
            -v -f 1.0 | cut -f4- \
            >> "~{filtered_copy_ratios}"
        else
          awk '!/^@/ && $1!="CONTIG"{print}' "~{copy_ratios}" \
            >> "~{filtered_copy_ratios}"
        fi

        echo ""
    >>>

    output {
        File filtered_denoised_copy_ratios = filtered_copy_ratios
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
