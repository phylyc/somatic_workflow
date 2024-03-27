version development

import "sample.wdl" as s
import "sequencing_run.wdl" as seqrun
import "runtime_collection.wdl" as rtc
import "runtimes.wdl" as rt


workflow HarmonizeSamples {
    input {
        String harmonize_copy_ratios_script
        String merge_pileups_script
        Array[Sample] samples
        Array[Array[File?]] denoised_copy_ratios
        Array[Array[File?]] allelic_counts

        Int harmonize_min_target_length = 100
        Boolean compress_output = false

        RuntimeCollection runtime_collection
    }

    scatter (sample in samples) {
        scatter (sequencing_run in sample.sequencing_runs) {
            SequencingRun this_sequencing_run = sequencing_run
        }
    }
    Array[SequencingRun] sequencing_runs = flatten(this_sequencing_run)

    Array[File] all_dCRs = select_all(flatten(denoised_copy_ratios))
    Array[File] all_ACs = select_all(flatten(allelic_counts))

    Boolean has_dCR = length(all_dCRs) > 0
    Boolean has_AC = length(all_ACs) > 0

    if (has_dCR) {
        scatter (pair in zip(sequencing_runs, all_dCRs)) {
            if (pair.left.use_for_dCR) {
                String chosen_dcr_name = pair.left.sample_name
                File chosen_dcr = pair.right
            }
        }
        Array[String] seq_dcr_sample_names = select_all(chosen_dcr_name)
        Array[File] seq_denoised_copy_ratios = select_all(chosen_dcr)

        call HarmonizeCopyRatios {
            input:
                script = harmonize_copy_ratios_script,
                sample_names = seq_dcr_sample_names,
                denoised_copy_ratios = seq_denoised_copy_ratios,
                min_target_length = harmonize_min_target_length,
                compress_output = compress_output,
                runtime_params = runtime_collection.harmonize_copy_ratios
        }

        # sort output to match order of sample_names since glob doesn't guarantee order
        scatter (sample in samples) {
            scatter (h_dcr in HarmonizeCopyRatios.harmonized_denoised_copy_ratios) {
                String this_dcr_sample_name = basename(basename(basename(basename(h_dcr, ".hdf5"), ".tsv"), ".denoised_CR"), ".harmonized")
                if (sample.name == this_dcr_sample_name) {
                    File this_dcr = h_dcr
                }
            }
            Array[File] this_sample_dcr = select_all(this_dcr)
        }
        Array[File] sorted_harmonized_denoised_copy_ratios = flatten(this_sample_dcr)
    }

    if (has_AC) {
        scatter (pair in zip(sequencing_runs, all_ACs)) {
            if (pair.left.use_for_aCR) {
                String chosen_ac_name = pair.left.sample_name
                File chosen_ac = pair.right
            }
        }
        Array[String] seq_acr_sample_names = select_all(chosen_ac_name)
        Array[File] seq_allelic_counts = select_all(chosen_ac)

        call MergeAllelicCounts {
            input:
                script = merge_pileups_script,
                sample_names = seq_acr_sample_names,
                allelic_counts = seq_allelic_counts,
                compress_output = compress_output,
                runtime_params = runtime_collection.merge_allelic_counts
        }

        # sort output to match order of sample_names since glob doesn't guarantee order
        scatter (sample in samples) {
            scatter (allelic_count in MergeAllelicCounts.merged_allelic_counts) {
                String this_sample_name = basename(basename(allelic_count, ".gz"), ".pileup")
                if (sample.name == this_sample_name) {
                    File this_allelic_counts = allelic_count
                }
            }
            Array[File] this_sample_allelic_counts = select_all(this_allelic_counts)
        }
        Array[File] sorted_allelic_counts = flatten(this_sample_allelic_counts)
    }

    output {
        Array[File]? harmonized_denoised_copy_ratios = sorted_harmonized_denoised_copy_ratios
        Array[File]? merged_allelic_counts = sorted_allelic_counts
    }
}


task HarmonizeCopyRatios {
    input {
        String script = "https://github.com/phylyc/genomics_workflows/raw/master/python/harmonize_copy_ratios.py"

        Array[String]+ sample_names
        Array[File]+ denoised_copy_ratios
        Int min_target_length = 100
        Boolean compress_output = false
        Boolean verbose = true

        Runtime runtime_params
    }

    String output_dir = "."

#    Array[String] harmonized_denoised_cr = suffix(sample_names, ".denoised_CR.tsv")

    command <<<
        set -euxo pipefail
        wget -O harmonize_copy_ratios.py ~{script}
        python harmonize_copy_ratios.py \
            --output_dir '~{output_dir}' \
            ~{sep="' " prefix("--sample '", sample_names)}' \
            ~{sep="' " prefix("--copy_ratio '", denoised_copy_ratios)}' \
            --suffix ".harmonized.denoised_CR" \
            --threads ~{runtime_params.cpu} \
            --min_target_length ~{min_target_length} \
            ~{if compress_output then "--compress_output" else ""} \
            ~{if verbose then "--verbose" else ""}
    >>>

    output {
        Array[File] harmonized_denoised_copy_ratios = glob("*.harmonized.denoised_CR.tsv")
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

task MergeAllelicCounts {
    input {
        String script = "https://github.com/phylyc/genomics_workflows/raw/master/python/merge_pileups.py"

        Array[String]+ sample_names
        Array[File]+ allelic_counts
        Boolean compress_output = false
        Boolean verbose = true

        Runtime runtime_params
    }

    String output_dir = "."

#    Array[String] merged_pileups = suffix(sample_names, ".pileup")

    command <<<
        set -euxo pipefail
        wget -O merge_pileups.py ~{script}
        python merge_pileups.py \
            --output_dir '~{output_dir}' \
            ~{sep="' " prefix("--sample '", sample_names)}' \
            ~{sep="' " prefix("-P '", allelic_counts)}' \
            ~{if compress_output then "--compress_output" else ""} \
            ~{if verbose then "--verbose" else ""}
    >>>

    output {
        Array[File] merged_allelic_counts = glob(output_dir + "/*.pileup" + (if compress_output then ".gz" else ""))
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