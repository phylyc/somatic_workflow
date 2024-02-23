version development

import "sample.wdl"
import "sequencing_run.wdl"
import "runtimes.wdl"


workflow MergeSequencingRuns {
    input {
        Sample sample
        Array[File]? read_counts
        Array[File]? denoised_copy_ratios
        Array[File]? standardized_copy_ratios
        Array[File]? snp_array_allelic_counts
        Array[File]? somatic_allelic_counts
        Array[File]? germline_allelic_counts

        RuntimeCollection runtime_collection
    }

    call MergeReadCountData {
        input:
            sample_name = sample.name,
            read_counts = read_counts,
            denoised_copy_ratios = denoised_copy_ratios,
            standardized_copy_ratios = standardized_copy_ratios,
            runtime_params = runtime_collection.merge_sample
    }

    call MergeAllelicCounts {
        input:
            sample_name = sample.name,
            snp_array_allelic_counts = snp_array_allelic_counts,
            somatic_allelic_counts = somatic_allelic_counts,
            germline_allelic_counts = germline_allelic_counts,
            runtime_params = runtime_collection.merge_sample
    }

    call sample.UpdateSample {
        input:
            sample = sample,
            read_counts = MergeReadCountData.merged_read_counts,
            denoised_copy_ratios = MergeReadCountData.merged_denoised_copy_ratios,
            standardized_copy_ratios = MergeReadCountData.merged_standardized_copy_ratios,
            snp_array_allelic_counts = MergeAllelicCounts.merged_snp_array_allelic_counts,
            somatic_allelic_counts = MergeAllelicCounts.merged_somatic_allelic_counts,
            germline_allelic_counts = MergeAllelicCounts.merged_germline_allelic_counts
    }

    output {
        Sample updated_sample = UpdateSample.updated_sample
        File? merged_read_counts = MergeReadCountData.merged_read_counts
        File? merged_denoised_copy_ratios = MergeReadCountData.merged_denoised_copy_ratios
        File? merged_standardized_copy_ratios = MergeReadCountData.merged_standardized_copy_ratios
        File? merged_snp_array_allelic_counts = MergeAllelicCounts.merged_snp_array_allelic_counts
        File? merged_somatic_allelic_counts = MergeAllelicCounts.merged_somatic_allelic_counts
        File? merged_germline_allelic_counts = MergeAllelicCounts.merged_germline_allelic_counts
    }
}

task MergeReadCountData {
    input {
        String sample_name
        Array[File]? read_counts
        Array[File]? denoised_copy_ratios
        Array[File]? standardized_copy_ratios

        Boolean compress_output = true
        Runtime runtime_params
    }

    String output_read_counts = sample_name + ".read_counts.tsv"
    String output_denoised_copy_ratios = sample_name + ".denoised_CR.tsv"
    String output_standardized_copy_ratios = sample_name + ".standardized_CR.tsv"

    # For now, just select the first file in the array
    # TODO: Select the file with highest signal to noise ratio. How do we measure this? Variance across COUNT column?

    String command_read_count = if (defined(read_counts)) then "cp " + select_first([read_counts])[0] + " " + output_read_counts else ""
    String command_denoised_copy_ratios = if (defined(denoised_copy_ratios)) then "cp " + select_first([denoised_copy_ratios])[0] + " " + output_denoised_copy_ratios else ""
    String command_standardized_copy_ratios = if (defined(standardized_copy_ratios)) then "cp " + select_first([standardized_copy_ratios])[0] + " " + output_standardized_copy_ratios else ""

    command <<<
        set -e

        ~{command_read_count}
        ~{command_denoised_copy_ratios}
        ~{command_standardized_copy_ratios}
    >>>

    output {
        File? merged_read_counts = output_read_counts
        File? merged_denoised_copy_ratios = output_denoised_copy_ratios
        File? merged_standardized_copy_ratios = output_standardized_copy_ratios
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
        String sample_name
        Array[File]? snp_array_allelic_counts
        Array[File]? somatic_allelic_counts
        Array[File]? germline_allelic_counts

        Boolean compress_output = true
        Runtime runtime_params
    }

    String output_snp_array_allelic_counts = sample_name + ".snp_array.pileup"
    String output_somatic_allelic_counts = sample_name + ".somatic.pileup"
    String output_germline_allelic_counts = sample_name + ".germline.pileup"

    command <<<
        set -e
        python <<CODE
import argparse
import gzip
import pandas as pd
import warnings


def get_header_and_df(file_path: str, columns: list[str] = None):
    open_func = gzip.open if file_path.endswith(".gz") else open
    try:
        with open_func(file_path, "rt") as file:
            header = file.readline()  # assumes header is first line
            df = pd.read_csv(file, sep="\t", comment="#", header=0, names=columns, low_memory=False)
    except Exception as e:
        warnings.warn(f"Exception reading file {file_path}: {e}")
        warnings.warn(f"Setting header to None and df to empty DataFrame.")
        header = None
        df = pd.DataFrame(columns=columns)
    return header, df


def write_header_and_df(header: str, df: pd.DataFrame, file_path: str):
    open_func = gzip.open if file_path.endswith(".gz") else open
    with open_func(file_path, "wt") as file:
        file.write(header) if header is not None else None
        df.to_csv(file, sep="\t", index=False)


def merge_pileups(file_paths: list[str], output_file: str):
    columns = ["contig", "position", "ref_count", "alt_count", "other_alt_count", "allele_frequency"]
    headers = []
    pileups = []
    for file_path in file_paths:
        pileup_header, df = get_header_and_df(file_path=file_path, columns=columns)
        headers.append(pileup_header)
        pileups.append(df)
    headers = [h for h in headers if h is not None]
    merged_pileup = pd.concat(
        pileups, ignore_index=True
    ).groupby(["contig", "position"]).agg(
        {
            "ref_count": "sum",
            "alt_count": "sum",
            "other_alt_count": "sum",
            "allele_frequency": "max"
        }
    ).reset_index()
    new_header = headers[0] if len(headers) > 0 else None
    write_header_and_df(header=new_header, df=merged_pileup, file_path=output_file)


if ~{defined(snp_array_allelic_counts)}:
    merge_pileups(~{snp_array_allelic_counts}, "~{output_snp_array_allelic_counts}")
if ~{defined(somatic_allelic_counts)}:
    merge_pileups(~{somatic_allelic_counts}, "~{output_somatic_allelic_counts}")
if ~{defined(germline_allelic_counts)}:
    merge_pileups(~{germline_allelic_counts}, "~{output_germline_allelic_counts}")

        CODE
    >>>

    output {
        File? merged_snp_array_allelic_counts = output_snp_array_allelic_counts
        File? merged_somatic_allelic_counts = output_somatic_allelic_counts
        File? merged_germline_allelic_counts = output_germline_allelic_counts
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
