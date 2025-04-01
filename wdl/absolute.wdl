version development

import "runtime_collection.wdl" as rtc


workflow Absolute {
    input {
        String sample_name
        File copy_ratio_segmentation
        File af_model_parameters
        File? annotated_variants
        Float? purity
        Float? ploidy
        String? sex

        String acs_conversion_script = "https://github.com/phylyc/somatic_workflow/raw/master/python/acs_conversion.py"
        Int min_hets = 0
        Int min_probes = 2
        Float maf90_threshold = 0.485

        RuntimeCollection runtime_collection = RuntimeParameters.rtc
    }

    call rtc.DefineRuntimeCollection as RuntimeParameters

    call ModelSegmentsToACSConversion {
        input:
            script = acs_conversion_script,
            seg_final = copy_ratio_segmentation,
            af_model_parameters = af_model_parameters,
            min_hets = min_hets,
            min_probes = min_probes,
            maf90_threshold = maf90_threshold,
            runtime_params = runtime_collection.model_segments_to_acs_conversion
    }

    if (defined(annotated_variants)) {
        call ProcessMAFforAbsolute {
            input:
                sample_name = sample_name,
                maf = select_first([annotated_variants]),
                runtime_params = runtime_collection.process_maf_for_absolute
        }
    }

    call AbsoluteTask as AbsoluteACRTask {
        input:
            sample_name = sample_name,
            seg_file = ModelSegmentsToACSConversion.acs_converted_seg,
            skew = ModelSegmentsToACSConversion.skew,
            snv_maf = ProcessMAFforAbsolute.snv_maf,
            indel_maf = ProcessMAFforAbsolute.indel_maf,
            copy_ratio_type = "allelic",
            purity = purity,
            ploidy = ploidy,
            sex = sex,
            runtime_params = runtime_collection.absolute
    }

    # TODO: fix Absolute package for tCR
#    call AbsoluteTask as AbsoluteTCRTask {
#        input:
#            sample_name = sample_name,
#            seg_file = ModelSegmentsToACSConversion.acs_converted_seg,
#            skew = ModelSegmentsToACSConversion.skew,
#            snv_maf = ProcessMAFforAbsolute.snv_maf,
#            indel_maf = ProcessMAFforAbsolute.indel_maf,
#            copy_number_type = "total",
#            sex = sex,
#            runtime_params = runtime_collection.absolute
#    }

    output {
        File acs_copy_ratio_segmentation = ModelSegmentsToACSConversion.acs_converted_seg
        Float acs_copy_ratio_skew = ModelSegmentsToACSConversion.skew
        File? snv_maf = ProcessMAFforAbsolute.snv_maf
        File? indel_maf = ProcessMAFforAbsolute.indel_maf
        File acr_plot = AbsoluteACRTask.plot
        File acr_rdata = AbsoluteACRTask.rdata
#        File? tcr_plot = AbsoluteTCRTask.plot
#        File? tcr_rdata = AbsoluteTCRTask.rdata
    }
}

task ModelSegmentsToACSConversion {
    input {
        String script = "https://github.com/phylyc/somatic_workflow/raw/master/python/acs_conversion.py"

        File seg_final
        File af_model_parameters

        Int min_hets = 10
        Int min_probes = 4
        Float maf90_threshold = 0.485

        Runtime runtime_params
    }

    String output_dir = "."

    String output_seg = basename(seg_final, ".seg") + ".acs.seg"
    String output_skew = output_seg + ".skew"

    command <<<
        set -e
        wget -O acs_conversion.py ~{script}
        python acs_conversion.py \
            --output_dir ~{output_dir} \
            --seg '~{seg_final}' \
            --af_parameters '~{af_model_parameters}' \
            --min_hets ~{min_hets} \
            --min_probes ~{min_probes} \
            --maf90_threshold ~{maf90_threshold} \
            --verbose
    >>>

    output {
        File acs_converted_seg = output_seg
        Float skew = read_float(output_skew)
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

task ProcessMAFforAbsolute {
    input {
        String sample_name
        File maf

        Runtime runtime_params
    }

    String uncompressed_maf = basename(maf, ".gz")
    Boolean is_compressed = uncompressed_maf != basename(maf)

    String output_snv_maf = sample_name + ".snv.maf"
    String output_indel_maf = sample_name + ".indel.maf"

    command <<<
        set -euxo pipefail

        if [ "~{is_compressed}" == "true" ] ; then
            gzip -cd '~{maf}' > '~{uncompressed_maf}'
        else
            mv '~{maf}' '~{uncompressed_maf}'
        fi

        grep "^#" '~{uncompressed_maf}' > '~{output_snv_maf}'
        grep "^#" '~{uncompressed_maf}' > '~{output_indel_maf}'

        python <<EOF
import pandas as pd

maf = pd.read_csv('~{uncompressed_maf}', sep='\t', comment='#')
if maf.empty:
    print("No variants found in the input MAF file.")
else:
    cols_to_keep = [
        col
        for col in maf.columns
        if maf[col].astype(str).map(len).max() < 1000 | maf[col].isna().all()
    ]
    print("Removing columns:", sorted(set(maf.columns) - set(cols_to_keep)))
    print("Keeping columns:", cols_to_keep)

    maf = maf[cols_to_keep].rename(columns={"Start_Position": "Start_position", "End_Position": "End_position"})
    maf.loc[maf["Variant_Type"].isin(["SNP", "DNP", "TNP", "MNP", "ONP"])].to_csv('~{output_snv_maf}', sep='\t', index=False, mode='a')
    maf.loc[maf["Variant_Type"].isin(["INS", "DEL"])].to_csv('~{output_indel_maf}', sep='\t', index=False, mode='a')
EOF
    >>>

    output {
        File snv_maf = output_snv_maf
        File indel_maf = output_indel_maf
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

task AbsoluteTask {
    input {
        String sample_name
        File seg_file
        Float skew
        File? snv_maf
        File? indel_maf
        Float? purity
        Float? ploidy
        String? sex
        String? platform
        String copy_ratio_type = "allelic"

        Runtime runtime_params
    }

    String output_dir = "."
    String output_plot = output_dir + "/" + sample_name + "." + copy_ratio_type + ".ABSOLUTE_plot.pdf"
    String output_rdata = output_dir + "/" + sample_name + "." + copy_ratio_type + ".ABSOLUTE.RData"
    String output_mode_res = output_dir + "/" + sample_name + "." + copy_ratio_type + ".ABSOLUTE_mode.res.Rds"
    String output_mode_tab = output_dir + "/" + sample_name + "." + copy_ratio_type + ".ABSOLUTE_mode.tab.Rds"
    String output_ssnv_mode_tab = output_dir + "/" + sample_name + "." + copy_ratio_type + ".ABSOLUTE_SSNV.mode.res.Rds"

    command <<<
        # ABSOLUTE may fail for various edge cases, but we still want to capture
        # the output files for all other samples if run within the somatic workflow.
        set +e
        set -uxo pipefail

        touch '~{output_plot}'
        touch '~{output_rdata}'

        num_segments=$(( $(wc -l < '~{seg_file}') - 1 ))

        if [ $num_segments -gt 0 ] ; then
            Rscript /library/scripts/run_absolute.R \
                --results_dir '~{output_dir}' \
                --sample '~{sample_name}' \
                --seg_dat_fn '~{seg_file}' \
                ~{"--maf '" + snv_maf + "'"} \
                ~{"--indel_maf '" + indel_maf + "'"} \
                ~{"--alpha " + purity} \
                ~{"--tau " + ploidy} \
                ~{"--gender  " + sex} \
                ~{"--platform " + platform} \
                --ssnv_skew ~{skew} \
                --copy_num_type ~{copy_ratio_type} \
                --pkg_dir "/"
        else
            echo "No segments found in the input segmentation file. Exiting." >&2
        fi

        exit 0
    >>>

    output {
        File plot = output_plot
        File rdata = output_rdata
        File? mode_res = output_mode_res
        File? mode_tab = output_mode_tab
        File? ssnv_mode_res = output_ssnv_mode_tab
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