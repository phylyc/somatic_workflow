version development

import "runtime_collection.wdl" as rtc


workflow AbsoluteExtract {
    input {
        String? sample_name
        String? sex

        File rdata
        Int called_solution
        String analyst_id
        String copy_ratio_type = "allelic"

        File acs_copy_ratio_segmentation
        Float acs_copy_ratio_skew
        File? snv_maf
        File? indel_maf
        File? gvcf
        String? genome_build
        String map_to_absolute_copy_number_script = "https://github.com/phylyc/somatic_workflow/raw/master/python/map_to_absolute_copy_number.py"
        String calculate_cancer_cell_fraction_script = "https://github.com/phylyc/somatic_workflow/raw/master/python/calculate_cancer_cell_fraction.py"

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

    # Postprocessing (segment rescue + allelic re-split + CCF recomputation) only
    # applies to allelic copy-ratio output. In total copy-ratio mode ABSOLUTE drops no
    # segments and the segtab already carries total copy number and CCF, and the
    # allelic-split assumptions the rescue makes do not hold genome-wide — so skip it
    # and pass ABSOLUTE's output through unchanged.
    Boolean is_allelic = copy_ratio_type == "allelic"

    if (is_allelic) {
        call Postprocess {
            input:
                cnv_script = map_to_absolute_copy_number_script,
                snv_script = calculate_cancer_cell_fraction_script,
                sample_name = sample_name,
                sex = sex,
                maf = AbsoluteExtractTask.abs_maf,
                seg = AbsoluteExtractTask.segtab,
                seg_igv = AbsoluteExtractTask.segtab_igv,
                copy_ratio_segmentation = acs_copy_ratio_segmentation,
                acs_copy_ratio_skew = acs_copy_ratio_skew,
                snv_maf = snv_maf,
                indel_maf = indel_maf,
                purity = AbsoluteExtractTask.purity,
                ploidy = AbsoluteExtractTask.ploidy,
                copy_ratio_type = copy_ratio_type,
                runtime_params = runtime_collection.absolute_extract_postprocess
        }
    }

    output {
        File absolute_table = AbsoluteExtractTask.table
        Float absolute_purity = AbsoluteExtractTask.purity
        Float absolute_ploidy = AbsoluteExtractTask.ploidy
        File absolute_maf = select_first([Postprocess.abs_maf, AbsoluteExtractTask.abs_maf])
        File absolute_segtab = select_first([Postprocess.segtab, AbsoluteExtractTask.segtab])
        File absolute_segtab_igv = select_first([Postprocess.segtab_igv, AbsoluteExtractTask.segtab_igv])
        File? absolute_rescued_intervals = Postprocess.rescued_intervals
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

        Int organism_normal_ploidy = 2

        Runtime runtime_params
    }

    String sample_name = basename(rdata, ".ABSOLUTE." + copy_ratio_type + ".RData")
    String output_dir = "."
    String output_table = output_dir + "/reviewed/" + sample_name + "." + analyst_id + ".ABSOLUTE.table." + copy_ratio_type + ".txt"
    String output_abs_maf = output_dir + "/reviewed/SEG_MAF/" + sample_name + ".ABS_MAF." + copy_ratio_type + ".txt"
    String output_segtab = output_dir + "/reviewed/SEG_MAF/" + sample_name + ".segtab." + copy_ratio_type + ".txt"
    String output_segtab_igv = output_dir + "/reviewed/SEG_MAF/" + sample_name + ".IGV.seg." + copy_ratio_type + ".txt"
    String output_called_rdata = output_dir + "/reviewed/samples/" + sample_name + ".ABSOLUTE." + copy_ratio_type + "." + analyst_id + ".called.RData"
    String output_gene_corrected_cn = output_dir + "/reviewed/" + sample_name + ".gene_corrected_CN." + copy_ratio_type + ".txt"

    command <<<
        set -euxo pipefail

        mkdir -p "~{output_dir}/reviewed/samples/"
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
            # This is a low-purity sample. We force purity to 1 for ABSOLUTE to
            # still get the same type of output files as for all other samples.
            # We set purity to 0 below.
            Rscript /opt/absolute/library/scripts/run_absolute.R \
                --results_dir "~{output_dir}/~{sample_name}.force-call" \
                --sample "~{sample_name}" \
                --seg_dat_fn "~{acs_copy_ratio_segmentation}" \
                ~{"--maf '" + snv_maf + "'"} \
                ~{"--indel_maf '" + indel_maf + "'"} \
                --alpha 1 \
                --tau ~{organism_normal_ploidy} \
                ~{"--gender " + sex} \
                ~{"--platform " + platform} \
                ~{if (defined(acs_copy_ratio_skew) && (acs_copy_ratio_skew > 0)) then "--ssnv_skew " + acs_copy_ratio_skew else ""} \
                --copy_num_type ~{copy_ratio_type} \
                ~{"--genome_build '" + genome_build + "'"} \
                --pkg_dir "/opt/absolute"

            this_rdata="~{output_dir}/~{sample_name}.force-call/~{sample_name}.ABSOLUTE.~{copy_ratio_type}.RData"
            this_called_solution=1

        else
            this_rdata='~{rdata}'
            this_called_solution=~{called_solution}
        fi

        Rscript /opt/absolute/library/scripts/extract_solution.R \
            --solution_num $this_called_solution \
            --rdata "$this_rdata" \
            --results_dir ~{output_dir} \
            --analyst_id '~{analyst_id}' \
            --sample '~{sample_name}' \
            --copy_num_type ~{copy_ratio_type} \
            --pkg_dir "/opt/absolute"

        if [[ "~{called_solution}" == "0" ]] ; then
            # Set call status to "low purity" and purity to 0:
            sed -i -E "s/\tcalled\t1/\tlow purity\t0/" "~{output_table}"
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
        String cnv_script = "https://github.com/phylyc/somatic_workflow/raw/master/python/map_to_absolute_copy_number.py"
        String snv_script = "https://github.com/phylyc/somatic_workflow/raw/master/python/calculate_cancer_cell_fraction.py"

        String? sample_name
        String? sex
        File maf
        File seg
        File seg_igv
        File copy_ratio_segmentation
        Float? acs_copy_ratio_skew
        File? snv_maf
        File? indel_maf
        Float purity
        Float ploidy
        String copy_ratio_type = "allelic"

        Int organism_normal_ploidy = 2

        Runtime runtime_params
    }

    String this_sample_name = if defined(sample_name) then sample_name else basename(seg, ".segtab.txt")
    String output_maf = this_sample_name + ".ABS_MAF." + copy_ratio_type + ".completed.txt"
    String output_segtab = this_sample_name + ".segtab." + copy_ratio_type + ".completed.txt"
    String output_segtab_igv = this_sample_name + ".IGV.seg." + copy_ratio_type + ".completed.txt"
    String output_rescued_intervals = this_sample_name + ".rescued_intervals." + copy_ratio_type + ".txt"

    command <<<
        set -euxo pipefail

        # Create default output files
        cp "~{maf}" "~{output_maf}"
        cp "~{seg}" "~{output_segtab}"
        cp "~{seg_igv}" "~{output_segtab_igv}"
        echo "Chromosome\tStart.bp\tEnd.bp\n" > "~{output_rescued_intervals}"

        wget -O map_to_absolute_copy_number.py ~{cnv_script}
        # If purity == 0, we rescue intervals assuming purity 1,
        python map_to_absolute_copy_number.py \
            --sample '~{this_sample_name}' \
            ~{"--sex  " + sex} \
            --absolute_segtab '~{seg}' \
            --acs_cr_seg '~{copy_ratio_segmentation}' \
            --purity ~{purity} \
            --ploidy ~{ploidy} \
            --normal_ploidy ~{organism_normal_ploidy} \
            --copy_num_type ~{copy_ratio_type} \
            --allelic_resplit_focals \
            --outdir "."

        if [[ "~{defined(snv_maf)}" == "true" || "~{defined(indel_maf)}" == "true" ]] ; then
            wget -O calculate_cancer_cell_fraction.py ~{snv_script}
            python calculate_cancer_cell_fraction.py \
                --sample '~{this_sample_name}' \
                ~{"--sex " + sex} \
                --absolute_maf '~{maf}' \
                --absolute_segtab '~{output_segtab}' \
                ~{if (defined(acs_copy_ratio_skew) && (acs_copy_ratio_skew > 0)) then "--ssnv_skew " + acs_copy_ratio_skew else ""} \
                ~{"--snv_maf '" + snv_maf + "'"} \
                ~{"--indel_maf '" + indel_maf + "'"} \
                --purity ~{purity} \
                --ploidy ~{ploidy} \
                --normal_ploidy ~{organism_normal_ploidy} \
                --copy_num_type ~{copy_ratio_type} \
                --outdir '.'
        fi
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