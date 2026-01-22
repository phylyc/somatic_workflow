version development

import "sample.wdl" as s
import "patient.wdl" as p
import "patient.update_samples.wdl" as p_update_s
import "workflow_arguments.wdl" as wfargs
import "runtime_collection.wdl" as rtc
import "absolute.wdl" as abs
import "absolute_extract.wdl" as abs_extract


workflow ClonalAnalysisWorkflow {
    input {
        Patient patient
        WorkflowArguments args
        RuntimeCollection runtime_collection
    }

    if (args.run_clonal_decomposition) {
        scatter (sample in patient.samples) {
            if (defined(sample.called_copy_ratio_segmentation) && defined(sample.af_model_parameters) && !defined(sample.absolute_acr_rdata) && !defined(sample.absolute_acr_plot)) {
                call abs.Absolute {
                    input:
                        acs_conversion_script = args.script_acs_conversion,
                        sample_name = sample.name,
                        copy_ratio_segmentation = select_first([sample.called_copy_ratio_segmentation]),
                        af_model_parameters = select_first([sample.af_model_parameters]),
                        annotated_variants = sample.annotated_somatic_variants,
                        purity = sample.purity,
                        ploidy = sample.ploidy,
                        sex = patient.sex,
                        min_hets = args.absolute_min_hets,
                        min_probes = args.absolute_min_probes,
                        maf90_threshold = args.absolute_maf90_threshold,
                        genome_build = args.genome_build,
                        runtime_collection = runtime_collection
                }
            }

            File acs_copy_ratio_segmentation = select_first([Absolute.acs_copy_ratio_segmentation, sample.acs_copy_ratio_segmentation])
            Float acs_copy_ratio_skew = select_first([Absolute.acs_copy_ratio_skew, sample.acs_copy_ratio_skew])
            File? snv_maf = if defined(Absolute.snv_maf) then Absolute.snv_maf else sample.absolute_snv_maf
            File? indel_maf = if defined(Absolute.indel_maf) then Absolute.indel_maf else sample.absolute_indel_maf
            File acr_rdata = select_first([Absolute.acr_rdata, sample.absolute_acr_rdata])
            File acr_plot = select_first([Absolute.acr_plot, sample.absolute_acr_plot])

            if (defined(sample.absolute_solution)) {
                call abs_extract.AbsoluteExtract {
                    input:
                        map_to_absolute_copy_number_script = args.script_map_to_absolute_copy_number,
                        sample_name = sample.name,
                        sex = patient.sex,
                        rdata = acr_rdata,
                        called_solution = select_first([sample.absolute_solution]),
                        analyst_id = args.analyst_id,
                        copy_ratio_type = "allelic",
                        acs_copy_ratio_segmentation = acs_copy_ratio_segmentation,
                        acs_copy_ratio_skew = acs_copy_ratio_skew,
                        snv_maf = snv_maf,
                        indel_maf = indel_maf,
                        gvcf = patient.gvcf,
                        genome_build = args.absolute_genome_build,
                }
            }
        }

        if (length(select_all(snv_maf)) > 0) {
            Array[File] absolute_snv_maf = select_all(snv_maf)
        }
        if (length(select_all(indel_maf)) > 0) {
            Array[File] absolute_indel_maf = select_all(indel_maf)
        }
        if (length(select_all(AbsoluteExtract.absolute_maf)) > 0) {
            Array[File] absolute_maf = select_all(AbsoluteExtract.absolute_maf)
        }
        if (length(select_all(AbsoluteExtract.absolute_segtab)) > 0) {
            Array[File] absolute_segtab = select_all(AbsoluteExtract.absolute_segtab)
        }
        if (length(select_all(AbsoluteExtract.absolute_table)) > 0) {
            Array[File] absolute_table = select_all(AbsoluteExtract.absolute_table)
        }
        if (length(select_all(AbsoluteExtract.absolute_purity)) > 0) {
            Array[Float] purity = select_all(AbsoluteExtract.absolute_purity)
        }
        if (length(select_all(AbsoluteExtract.absolute_ploidy)) > 0) {
            Array[Float] ploidy = select_all(AbsoluteExtract.absolute_ploidy)
        }

        call p_update_s.UpdateSamples {
            input:
                patient = patient,
                acs_copy_ratio_segmentation = acs_copy_ratio_segmentation,
                acs_copy_ratio_skew = acs_copy_ratio_skew,
                absolute_acr_rdata = acr_rdata,
                absolute_acr_plot = acr_plot,
                absolute_snv_maf = absolute_snv_maf,
                absolute_indel_maf = absolute_indel_maf,
                absolute_maf = absolute_maf,
                absolute_segtab = absolute_segtab,
                absolute_table = absolute_table,
                purity = purity,
                ploidy = ploidy,
        }
    }

    # phylogicNDT

    output {
        Patient updated_patient = select_first([UpdateSamples.updated_patient, patient])
    }
}