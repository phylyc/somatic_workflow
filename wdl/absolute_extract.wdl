version development

import "runtime_collection.wdl" as rtc


workflow absolute_extract_workflow {
    input {
        String sample_name
        File absolute_summary_data_gatk_acnv
        Int absolute_called_solution
        String analyst_id

        RuntimeCollection runtime_collection = RuntimeParameters.rtc
    }

    call rtc.DefineRuntimeCollection as RuntimeParameters

    call absolute_extract {
        input:
            sample_name = sample_name,
            absolute_summary_data_gatk_acnv = absolute_summary_data_gatk_acnv,
            absolute_called_solution = absolute_called_solution,
            analyst_id = analyst_id,
            runtime_params = runtime_collection.absolute_extract

    }
}


task absolute_extract{
    input {
        String sample_name
        File absolute_summary_data_gatk_acnv
        Int absolute_called_solution
        String analyst_id

        Runtime runtime_params
    }

    command <<<
        set -euxo pipefail

        Rscript /usr/local/bin/ABSOLUTE_extract_cli_start.R \
            --solution_num $absolute_called_solution \
            --analyst_id "${analyst_id}" \
            --rdata_modes_fn "${absolute_summary_data_gatk_acnv}" \
            --sample_name "${sample_name}" \
            --results_dir . \
            --abs_lib_dir /xchip/tcga/Tools/absolute/releases/v1.5/

        cp "reviewed/SEG_MAF/${sample_name}_ABS_MAF.txt" .
        cp "reviewed/SEG_MAF/${sample_name}.segtab.txt" .
        cp "reviewed/samples/${sample_name}.ABSOLUTE.${analyst_id}.called.RData" .
        cp "reviewed/${sample_name}.${analyst_id}.ABSOLUTE.table.txt" .

        cut -f4 "${sample_name}.${analyst_id}.ABSOLUTE.table.txt" | tail -n 1 > purity
        cut -f5 "${sample_name}.${analyst_id}.ABSOLUTE.table.txt" | tail -n 1 > ploidy
    >>>

    output {
        File absolute_annotated_maf_capture = sample_name + "_ABS_MAF.txt"
        File absolute_seg_file = sample_name + ".segtab.txt"
        File absolute_segdat_file = sample_name + ".ABSOLUTE." + analyst_id + ".called.RData"
        File absolute_table = sample_name + "." + analyst_id + ".ABSOLUTE.table.txt"
        String purity = read_string("purity")
        String ploidy = read_string("ploidy")
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