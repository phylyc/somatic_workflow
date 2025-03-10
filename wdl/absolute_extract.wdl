version development

import "runtime_collection.wdl" as rtc


workflow AbsoluteExtract {
    input {
        String sample_name
        String? sex

        File rdata
        Int called_solution
        String analyst_id
        String? copy_ratio_type

        File? copy_ratio_segmentation
        File? annotated_variants
        File? gvcf

        RuntimeCollection runtime_collection = RuntimeParameters.rtc
    }

    call rtc.DefineRuntimeCollection as RuntimeParameters

    call AbsoluteExtractTask {
        input:
            rdata = rdata,
            called_solution = called_solution,
            analyst_id = analyst_id,
            copy_ratio_type = copy_ratio_type,
            runtime_params = runtime_collection.absolute_extract
    }

    call Postprocess {
        input:
            sample_name = sample_name,
            sex = sex,
            maf = AbsoluteExtractTask.abs_maf,
            seg = AbsoluteExtractTask.segtab,
            copy_ratio_segmentation = copy_ratio_segmentation,
            annotated_variants = annotated_variants,
            gvcf = gvcf,
            purity = AbsoluteExtractTask.purity,
            ploidy = AbsoluteExtractTask.ploidy,
            runtime_params = runtime_collection.absolute_extract_postprocess
    }

    output {
        File absolute_maf = Postprocess.abs_maf
        File absolute_segtab = Postprocess.segtab
        File absolute_table = AbsoluteExtractTask.table
        Float absolute_purity = AbsoluteExtractTask.purity
        Float absolute_ploidy = AbsoluteExtractTask.ploidy
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
    >>>

    output {
        File abs_maf = output_dir + "/reviewed/SEG_MAF/" + sample_name + ".ABS_MAF.txt"
        File segtab = output_dir + "/reviewed/SEG_MAF/" + sample_name + ".segtab.txt"
        File segtab_igv = output_dir + "/reviewed/SEG_MAF/" + sample_name + ".IGV.seg.txt"
        File called_rdata = output_dir + "/reviewed/samples/" + sample_name + ".ABSOLUTE." + analyst_id + ".called.RData"
        File table = output_table
        File gene_corrected_cn = output_dir + "/reviewed/" + sample_name + ".gene_corrected_CN.txt"
        Float purity = read_float("purity")
        Float ploidy = read_float("ploidy")
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

task Postprocess {
    input {
        String script = "https://github.com/phylyc/somatic_workflow/raw/master/python/map_to_absolute_copy_number.py"

        String sample_name
        String? sex
        File maf
        File seg
        File? copy_ratio_segmentation
        File? annotated_variants
        File? gvcf
        Float? purity
        Float? ploidy

        Runtime runtime_params
    }

    command <<<
        set -euxo pipefail
        wget -O map_to_absolute_copy_number.py ~{script}
        python map_to_absolute_copy_number.py \
            --sample '~{sample_name}' \
            ~{"--sex  " + sex} \
            --absolute_seg '~{seg}' \
            --cr_seg '~{copy_ratio_segmentation}' \
            --gvcf '~{gvcf}' \
            --purity ~{purity} \
            --ploidy ~{ploidy} \
            --outdir "."
    >>>

    output {
        File abs_maf = maf
        File segtab = sample_name + ".segtab.completed.txt"
        File segtab_igv = sample_name + ".IGV.seg.completed.txt"
        File rescued_intervals = sample_name + ".rescued_intervals.txt"
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