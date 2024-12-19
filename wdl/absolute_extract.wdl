version development

import "runtime_collection.wdl" as rtc


workflow AbsoluteExtract {
    input {
        File rdata
        Int called_solution
        String analyst_id
        String? copy_ratio_type

        RuntimeCollection runtime_collection = RuntimeParameters.rtc
    }

    call rtc.DefineRuntimeCollection as RuntimeParameters

    call AbsoluteExtractTask {
        input:
            rdata = rdata,
            called_solution = called_solution,
            analyst_id = analyst_id,
#            copy_ratio_type = copy_ratio_type,
            runtime_params = runtime_collection.absolute_extract

    }

    output {
        File? absolute_maf = AbsoluteExtractTask.abs_maf
        File? absolute_segtab = AbsoluteExtractTask.segtab
        File? absolute_called_rdata = AbsoluteExtractTask.called_rdata
        File? absolute_table = AbsoluteExtractTask.table
        File? absolute_gene_corrected_cn = AbsoluteExtractTask.gene_corrected_cn
        File? absolute_rescaled_total_cn = AbsoluteExtractTask.rescaled_total_cn
        String? absolute_purity = AbsoluteExtractTask.purity
        String? absolute_ploidy = AbsoluteExtractTask.ploidy
    }
}


task AbsoluteExtractTask {
    input {
        File rdata
        Int called_solution
        String analyst_id
        String copy_ratio_type = "allelic"

        Runtime runtime_params
    }

    String sample_name = basename(rdata, "." + copy_ratio_type + ".ABSOLUTE.RData")
    String output_dir = "."
    String output_table = output_dir + "/reviewed/" + sample_name + "." + analyst_id + ".ABSOLUTE.table.txt"

    command <<<
        set -euxo pipefail

        if [[ "~{called_solution}" =~ ^[0-9]+$ ]] && (( ~{called_solution} > 0 )); then
            Rscript /library/scripts/extract_solution.R \
                --solution_num ~{called_solution} \
                --results_dir ~{output_dir} \
                --analyst_id ~{analyst_id} \
                --pkg_dir "/" \
                --sample ~{sample_name} \
                --rdata ~{rdata} \
                --copy_num_type ~{copy_ratio_type}

            cut -f4 "~{output_table}" | tail -n 1 > purity
            cut -f5 "~{output_table}" | tail -n 1 > ploidy
        else
            echo "Called solution needs to be a positive integer but is: ~{called_solution}"
            printf "-1" > purity
            printf "-1" > ploidy
        fi
    >>>

    output {
        File? abs_maf = output_dir + "/reviewed/SEG_MAF/" + sample_name + "_ABS_MAF.txt"
        File? segtab = output_dir + "/reviewed/SEG_MAF/" + sample_name + ".segtab.txt"
        File? called_rdata = output_dir + "/reviewed/samples/" + sample_name + ".ABSOLUTE." + analyst_id + ".called.RData"
        File? table = output_table
        File? gene_corrected_cn = output_dir + "/reviewed/" + sample_name + "_gene_corrected_CN.txt"
        File? rescaled_total_cn = output_dir + "/reviewed/" + sample_name + "_rescaled_total_cn.IGV.seg.txt"
        String? purity = read_string("purity")
        String? ploidy = read_string("ploidy")
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