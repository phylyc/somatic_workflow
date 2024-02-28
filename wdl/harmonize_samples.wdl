version development

import "sample.wdl" as s
import "runtime_collection.wdl" as rtc
import "runtimes.wdl" as rt


workflow HarmonizeSamples {
    input {
        Array[Sample] samples
        Array[Array[File?]] denoised_copy_ratios
        Array[Array[File?]] standardized_copy_ratios
        Array[Array[File?]] allelic_counts

        Boolean compress_output = false

        RuntimeCollection runtime_collection
    }

    # The file inputs are per sequencing run!
    scatter (sample in samples) {
        scatter (sequencing_run in sample.sequencing_runs) {
            String this_seq_sample_names = sample.name
        }
        String sample_names = sample.name
    }
    Array[String] seq_sample_names = flatten(this_seq_sample_names)

    Array[File] non_optional_denoised_copy_ratios = select_all(flatten(denoised_copy_ratios))
    Array[File] non_optional_standardized_copy_ratios = select_all(flatten(standardized_copy_ratios))
    Array[File] non_optional_allelic_counts = select_all(flatten(allelic_counts))

    Boolean has_dCR = length(non_optional_denoised_copy_ratios) > 0
    Boolean has_sCR = length(non_optional_standardized_copy_ratios) > 0
    Boolean has_AC = length(non_optional_allelic_counts) > 0

    if (has_dCR && has_sCR) {
        call HarmonizeCopyRatios {
            input:
                sample_names = seq_sample_names,
                denoised_copy_ratios = non_optional_denoised_copy_ratios,
                standardized_copy_ratios = non_optional_standardized_copy_ratios,
                compress_output = compress_output,
                runtime_params = runtime_collection.harmonize_copy_ratios
        }

        # sort output to match order of sample_names since glob doesn't guarantee order
        scatter (sample in samples) {
            scatter (pair in zip(HarmonizeCopyRatios.harmonized_denoised_copy_ratios, HarmonizeCopyRatios.harmonized_standardized_copy_ratios)) {
                String this_dcr_sample_name = basename(basename(basename(pair.left, ".hdf5"), ".tsv"), ".denoised_CR")
                String this_scr_sample_name = basename(basename(basename(pair.right, ".hdf5"), ".tsv"), ".standardized_CR")
                if (sample.name == this_dcr_sample_name) {
                    File? this_dcr = pair.left
                }
                if (sample.name == this_scr_sample_name) {
                    File? this_scr = pair.right
                }
            }
            Array[File] this_sample_dcr = select_all(this_dcr)
            Array[File] this_sample_scr = select_all(this_scr)
        }
        Array[File]? sorted_harmonized_denoised_copy_ratios = flatten(this_sample_dcr)
        Array[File]? sorted_harmonized_standardized_copy_ratios = flatten(this_sample_scr)
    }

    if (has_AC) {
        call MergeAllelicCounts {
            input:
                sample_names = seq_sample_names,
                allelic_counts = non_optional_allelic_counts,
                compress_output = compress_output,
                runtime_params = runtime_collection.merge_allelic_counts
        }

        # sort output to match order of sample_names since glob doesn't guarantee order
        scatter (sample in samples) {
            scatter (allelic_count in MergeAllelicCounts.merged_allelic_counts) {
                String this_sample_name = basename(basename(allelic_count, ".gz"), ".pileup")
                if (sample.name == this_sample_name) {
                    File? this_allelic_counts = allelic_count
                }
            }
            Array[File] this_sample_allelic_counts = select_all(this_allelic_counts)
        }
        Array[File]? sorted_allelic_counts = flatten(this_sample_allelic_counts)
    }

    output {
        Array[File]? harmonized_denoised_copy_ratios = sorted_harmonized_denoised_copy_ratios
        Array[File]? harmonized_standardized_copy_ratios = sorted_harmonized_standardized_copy_ratios
        Array[File]? merged_allelic_counts = sorted_allelic_counts
    }
}


task HarmonizeCopyRatios {
    input {
        String script = "https://github.com/phylyc/genomics_workflows/raw/master/python/harmonize_copy_ratios.py"

        Array[String]+ sample_names
        Array[File]+ denoised_copy_ratios
        Array[File]+ standardized_copy_ratios
        Boolean compress_output = false
        Boolean verbose = true

        Runtime runtime_params
    }

    String output_dir = "."

#    Array[String] harmonized_denoised_cr = suffix(sample_names, ".denoised_CR.tsv")
#    Array[String] harmonized_standardized_cr = suffix(sample_names, ".standardized_CR.tsv")

    command <<<
        set -e
        wget -O harmonize_copy_ratios.py ~{script}
        python harmonize_copy_ratios.py \
            --output_dir '~{output_dir}' \
            ~{sep="' " prefix("--sample '", sample_names)}' \
            ~{sep="' " prefix("--denoised_cr '", denoised_copy_ratios)}' \
            ~{sep="' " prefix("--standardized_cr '", standardized_copy_ratios)}' \
            ~{if compress_output then "--compress_output" else ""} \
            ~{if verbose then "--verbose" else ""}
    >>>

    output {
        Array[File] harmonized_denoised_copy_ratios = glob("*.denoised_CR.tsv")
        Array[File] harmonized_standardized_copy_ratios = glob("*.standardized_CR.tsv")
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
        set -e
        wget -O merge_pileups.py ~{script}
        python merge_pileups.py \
            --output_dir '~{output_dir}' \
            ~{sep="' " prefix("--sample '", sample_names)}' \
            ~{sep="' " prefix("-P '", allelic_counts)}' \
            ~{if compress_output then "--compress_output" else ""} \
            ~{if verbose then "--verbose" else ""}
    >>>

    output {
        Array[File] merged_allelic_counts = glob("*.pileup" + (if compress_output then ".gz" else ""))
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