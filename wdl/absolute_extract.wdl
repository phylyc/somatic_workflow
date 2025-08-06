version development

import "runtime_collection.wdl" as rtc


workflow AbsoluteExtract {
    input {
        String? sample_name
        String? sex

        File rdata
        Int called_solution
        String analyst_id
        String? copy_ratio_type

        File acs_copy_ratio_segmentation
        Float acs_copy_ratio_skew
        File? snv_maf
        File? indel_maf
        File? gvcf
        String? genome_build
        String map_to_absolute_copy_number_script = "https://github.com/phylyc/somatic_workflow/raw/master/python/map_to_absolute_copy_number.py"

        RuntimeCollection runtime_collection = RuntimeParameters.rtc
    }

    call rtc.DefineRuntimeCollection as RuntimeParameters

    call AbsoluteExtractTask {
        input:
            rdata = rdata,
            called_solution = called_solution,
            analyst_id = analyst_id,
            copy_ratio_type = copy_ratio_type,
            acs_copy_ratio_segmentation = acs_copy_ratio_segmentation,
            acs_copy_ratio_skew = acs_copy_ratio_skew,
            snv_maf = snv_maf,
            indel_maf = indel_maf,
            sex = sex,
            genome_build = genome_build,
            runtime_params = runtime_collection.absolute_extract
    }

    call Postprocess {
        input:
            script = map_to_absolute_copy_number_script,
            sample_name = sample_name,
            sex = sex,
            maf = AbsoluteExtractTask.abs_maf,
            seg = AbsoluteExtractTask.segtab,
            seg_igv = AbsoluteExtractTask.segtab_igv,
            copy_ratio_segmentation = acs_copy_ratio_segmentation,
#            snv_maf = snv_maf,
#            indel_maf = indel_maf,
#            gvcf = gvcf,
            purity = AbsoluteExtractTask.purity,
            ploidy = AbsoluteExtractTask.ploidy,
            runtime_params = runtime_collection.absolute_extract_postprocess
    }

    output {
        File absolute_table = AbsoluteExtractTask.table
        Float absolute_purity = AbsoluteExtractTask.purity
        Float absolute_ploidy = AbsoluteExtractTask.ploidy
        File absolute_maf = Postprocess.abs_maf
        File absolute_segtab = Postprocess.segtab
    }
}


task AbsoluteExtractTask {
    input {
        File rdata
        Int called_solution
        String analyst_id
        String copy_ratio_type = "allelic"

        File? acs_copy_ratio_segmentation
        Float? acs_copy_ratio_skew
        File? snv_maf
        File? indel_maf
        String? sex
        String? platform
        String? genome_build

        Runtime runtime_params
    }

    String sample_name = basename(rdata, "." + copy_ratio_type + ".ABSOLUTE.RData")
    String output_dir = "."
    String output_table = output_dir + "/reviewed/" + sample_name + "." + analyst_id + ".ABSOLUTE.table.txt"
    String output_abs_maf = output_dir + "/reviewed/SEG_MAF/" + sample_name + ".ABS_MAF.txt"
    String output_segtab = output_dir + "/reviewed/SEG_MAF/" + sample_name + ".segtab.txt"
    String output_segtab_igv = output_dir + "/reviewed/SEG_MAF/" + sample_name + ".IGV.seg.txt"
    String output_called_rdata = output_dir + "/reviewed/samples/" + sample_name + ".ABSOLUTE." + analyst_id + ".called.RData"
    String output_gene_corrected_cn = output_dir + "/reviewed/" + sample_name + ".gene_corrected_CN.txt"

    command <<<
        set -euxo pipefail

        mkdir -p "~{output_dir}/reviewed/"
        mkdir -p "~{output_dir}/reviewed/SEG_MAF/"

        # Create dummy files to guarantee output files
        touch "~{output_table}"
        touch "~{output_abs_maf}"
        touch "~{output_segtab}"
        touch "~{output_segtab_igv}"
        touch "~{output_called_rdata}"
        touch "~{output_gene_corrected_cn}"

        if [[ "~{called_solution}" == "-1" ]] ; then
            echo -e "array\tsample\tcall status\tpurity\tploidy\tGenome doublings\tdelta\tCoverage for 80% power\tCancer DNA fraction\tSubclonal genome fraction\ttau\tE_CR\n" \
                > "~{output_table}"
            echo -e "~{sample_name}\t~{sample_name}\tfailed\t\t\t\t\t\t\t\t\t\n" \
                >> "~{output_table}"
            echo -1 > purity
            echo -1 > ploidy
            exit 0

        elif [[ "~{called_solution}" == "0" ]] ; then
            Rscript /library/scripts/run_absolute.R \
                --results_dir "~{output_dir}/~{sample_name}.force-call" \
                --sample "~{sample_name}" \
                --seg_dat_fn "~{acs_copy_ratio_segmentation}" \
                ~{"--maf '" + snv_maf + "'"} \
                ~{"--indel_maf '" + indel_maf + "'"} \
                --alpha 1 \
                --tau 2 \
                ~{"--gender " + sex} \
                ~{"--platform " + platform} \
                --ssnv_skew ~{acs_copy_ratio_skew} \
                --copy_num_type ~{copy_ratio_type} \
                ~{"--genome_build '" + genome_build + "'"} \
                --pkg_dir "/"

            this_rdata="~{output_dir}/~{sample_name}.force-call/~{sample_name}.allelic.ABSOLUTE.RData"
            this_called_solution=1

        else
            this_rdata='~{rdata}'
            this_called_solution=~{called_solution}
        fi

        Rscript /library/scripts/extract_solution.R \
            --solution_num $this_called_solution \
            --rdata "$this_rdata" \
            --results_dir ~{output_dir} \
            --analyst_id '~{analyst_id}' \
            --sample '~{sample_name}' \
            --copy_num_type ~{copy_ratio_type} \
            --pkg_dir "/"

        if [[ "~{called_solution}" == "0" ]] ; then
            sed -i -E "s/\tcalled\t/\tlow purity\t/" "~{output_table}"
        fi

        cut -f4 "~{output_table}" | tail -n 1 > purity
        cut -f5 "~{output_table}" | tail -n 1 > ploidy
    >>>

    output {
        File table = output_table
        Float purity = read_float("purity")
        Float ploidy = read_float("ploidy")
        File abs_maf = output_abs_maf
        File segtab = output_segtab
        File segtab_igv = output_segtab_igv
        File called_rdata = output_called_rdata
        File gene_corrected_cn = output_gene_corrected_cn
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

        String? sample_name
        String? sex
        File maf
        File seg
        File seg_igv
        File copy_ratio_segmentation
#        File? snv_maf
#        File? indel_maf
#        File? gvcf
        Float purity
        Float ploidy

        Runtime runtime_params
    }

    String this_sample_name = if defined(sample_name) then sample_name else basename(seg, ".segtab.txt")
    String output_maf = this_sample_name + ".ABS_MAF.completed.txt"
    String output_segtab = this_sample_name + ".segtab.completed.txt"
    String output_segtab_igv = this_sample_name + ".IGV.seg.completed.txt"
    String output_rescued_intervals = this_sample_name + ".rescued_intervals.txt"

    command <<<
        set -euxo pipefail

        # Create default output files
        cp "~{maf}" "~{output_maf}"
        cp "~{seg}" "~{output_segtab}"
        cp "~{seg_igv}" "~{output_segtab_igv}"
        echo "Chromosome\tStart.bp\tEnd.bp\n" > "~{output_rescued_intervals}"

        wget -O map_to_absolute_copy_number.py ~{script}
        python map_to_absolute_copy_number.py \
            --sample '~{this_sample_name}' \
            ~{"--sex  " + sex} \
            --absolute_seg '~{seg}' \
            --cr_seg '~{copy_ratio_segmentation}' \
            --purity ~{purity} \
            --ploidy ~{ploidy} \
            --outdir "."
    >>>

    output {
        File abs_maf = output_maf
        File segtab = output_segtab
        File segtab_igv = output_segtab_igv
        File rescued_intervals = output_rescued_intervals
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