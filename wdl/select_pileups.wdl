version development

import "runtime_collection.wdl" as rtc
import "runtimes.wdl" as rt


workflow SelectPileupsWorkflow {
    input {
        File pileup_summaries
        String sample_name = basename(basename(pileup_summaries, ".gz"), ".pileup")

        Float minimum_population_allele_frequency = 0.01
        Float maximum_population_allele_frequency = 0.2
        Int minimum_read_depth = 0

        Boolean compress_output = false

        RuntimeCollection runtime_collection = RuntimeParameters.rtc
    }

    call rtc.DefineRuntimeCollection as RuntimeParameters

    call SelectPileups {
        input:
            pileup_summaries = pileup_summaries,
            sample_name = sample_name,
            minimum_population_allele_frequency = minimum_population_allele_frequency,
            maximum_population_allele_frequency = maximum_population_allele_frequency,
            minimum_read_depth = minimum_read_depth,
            compress_output = compress_output,
            runtime_params = runtime_collection.select_pileup_summaries,
    }

    output {
        File selected_pileup_summaries = SelectPileups.selected_pileup_summaries
    }
}

task SelectPileups {
    input {
        File pileup_summaries
        String sample_name

        Float minimum_population_allele_frequency = 0.01
        Float maximum_population_allele_frequency = 0.2
        Int minimum_read_depth = 0

        Boolean compress_output

        Runtime runtime_params
    }

    String uncompressed_pileup_summaries = basename(pileup_summaries, ".gz")
    Boolean is_compressed = uncompressed_pileup_summaries != basename(pileup_summaries)

    String pileup_file = sample_name + ".pileup"
    String output_file = pileup_file + if compress_output then ".gz" else ""
    String tmp_pileup_file = "tmp." + pileup_file

    command <<<
        set -uxo pipefail

        if [ "~{is_compressed}" == "true" ] ; then
            bgzip -cd '~{pileup_summaries}' > '~{uncompressed_pileup_summaries}'
        else
            mv '~{pileup_summaries}' '~{uncompressed_pileup_summaries}'
        fi

        # Extract leading comment lines
        grep '^#' '~{uncompressed_pileup_summaries}' > '~{tmp_pileup_file}'

        # Extract column headers
        grep -v '^#' '~{uncompressed_pileup_summaries}' | head -n 1 >> '~{tmp_pileup_file}'

        # Count the number of lines that are not comments (headers) or column headers
        num_variants=$(( $(grep -vc '^#' '~{uncompressed_pileup_summaries}') - 1 ))

        if [ "$num_variants" -gt 0 ]; then
            # Extract table and select lines with read depth >= min_read_depth
            grep -v '^#' '~{uncompressed_pileup_summaries}' | tail -n +2 \
                | awk -F"\t" '$3 + $4 + $5 >= ~{minimum_read_depth}' \
                | awk -F"\t" '$6 >= ~{minimum_population_allele_frequency}' \
                | awk -F"\t" '$6 <= ~{maximum_population_allele_frequency}' \
                >> '~{tmp_pileup_file}'
        fi

        mv '~{tmp_pileup_file}' '~{pileup_file}'

        # Count the number of lines that are not comments (headers) or column headers
        num_selected_variants=$(( $(grep -vc '^#' '~{pileup_file}') - 1 ))

        echo ">> Selected $num_selected_variants loci out of $num_variants."

        if [ "~{compress_output}" == "true" ] ; then
            bgzip -c '~{pileup_file}' > '~{output_file}'
        fi
    >>>

    output {
        File selected_pileup_summaries = output_file
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