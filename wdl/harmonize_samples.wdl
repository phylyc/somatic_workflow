version development

import "sample.wdl" as s
import "sequencing_run.wdl" as seqrun
import "runtime_collection.wdl" as rtc
import "runtimes.wdl" as rt


workflow HarmonizeSamples {
    input {
        File ref_dict
        String harmonize_copy_ratios_script
        String merge_pileups_script
        Array[Sample] samples

        Int harmonize_min_target_length = 0
        Int pileups_min_read_depth = 0
        Boolean compress_output = false

        RuntimeCollection runtime_collection
    }

    scatter (sample in samples) {
        scatter (sequencing_run in sample.sequencing_runs) {
            File? cl = sequencing_run.callable_loci
            if (sequencing_run.use_for_tCR && (size(sequencing_run.denoised_total_copy_ratios) > 0)) {
                String dcr_sample_name = sequencing_run.sample_name
                File? dcr = sequencing_run.denoised_total_copy_ratios
            }
            if (sequencing_run.use_for_aCR && (size(sequencing_run.snppanel_allelic_pileup_summaries) > 0)) {
                String aps_sample_name = sequencing_run.sample_name
                File? aps = sequencing_run.snppanel_allelic_pileup_summaries
            }
        }
        Array[File] sample_cl = select_all(cl)
        Array[String] sample_dcr_names = select_all(dcr_sample_name)
        Array[File] sample_dcr = select_all(dcr)
        Array[String] sample_aps_name = select_all(aps_sample_name)
        Array[File] sample_aps = select_all(aps)
    }

    if (length(flatten(sample_cl)) > 0) {
        Array[File] sorted_harmonized_callable_loci = flatten(sample_cl)
        # TODO: harmonize bed files
    }

    if (length(flatten(sample_dcr_names)) > 0) {
        call HarmonizeCopyRatios {
            input:
                script = harmonize_copy_ratios_script,
                ref_dict = ref_dict,
                sample_names = flatten(sample_dcr_names),
                denoised_copy_ratios = flatten(sample_dcr),
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

    if (length(flatten(sample_aps_name)) > 0) {
        call MergeAllelicCounts {
            input:
                script = merge_pileups_script,
                ref_dict = ref_dict,
                sample_names = flatten(sample_aps_name),
                allelic_counts = flatten(sample_aps),
                min_read_depth = pileups_min_read_depth,
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
        Array[File]? harmonized_callable_loci = sorted_harmonized_callable_loci
        Array[File]? harmonized_denoised_copy_ratios = sorted_harmonized_denoised_copy_ratios
        Array[File]? merged_allelic_counts = sorted_allelic_counts
    }
}


task HarmonizeCopyRatios {
    input {
        String script = "https://github.com/phylyc/somatic_workflow/raw/master/python/harmonize_copy_ratios.py"

        File ref_dict
        Array[String]+ sample_names
        Array[File]+ denoised_copy_ratios
        Int min_target_length = 0
        Boolean compress_output = false
        Boolean verbose = true

        Runtime runtime_params
    }

    String output_dir = "."
    String suffix = ".harmonized.denoised_CR.tsv"

    command <<<
        set -euxo pipefail
        n_threads=$(nproc)
        wget -O harmonize_copy_ratios.py ~{script}
        python harmonize_copy_ratios.py \
            --output_dir '~{output_dir}' \
            --ref_dict '~{ref_dict}' \
            ~{sep="' " prefix("--sample '", sample_names)}' \
            ~{sep="' " prefix("--copy_ratio '", denoised_copy_ratios)}' \
            --suffix ~{suffix} \
            --threads $n_threads \
            --min_target_length ~{min_target_length} \
            --column_names CONTIG START END LOG2_COPY_RATIO \
            --column_types str int int float \
            --agg_col LOG2_COPY_RATIO \
            --agg_func strongest_signal \
            ~{if compress_output then "--compress_output" else ""} \
            ~{if verbose then "--verbose" else ""}
    >>>

    output {
        Array[File] harmonized_denoised_copy_ratios = glob(output_dir + "/*" + suffix)
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
        String script = "https://github.com/phylyc/somatic_workflow/raw/master/python/merge_pileups.py"

        File ref_dict
        Array[String]+ sample_names
        Array[File]+ allelic_counts
        Int min_read_depth = 0
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
            -D ~{ref_dict} \
            --min_read_depth ~{min_read_depth} \
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