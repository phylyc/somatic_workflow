version development

import "sample.wdl" as s
import "patient.wdl" as p
import "patient.update_samples.wdl" as p_update_s
import "workflow_arguments.wdl" as wfargs
import "runtime_collection.wdl" as rtc
import "absolute.wdl" as abs
import "absolute_extract.wdl" as abs_extract
import "phylogicndt.wdl" as phylogicndt


workflow ClonalAnalysisWorkflow {
    input {
        Patient patient
        WorkflowArguments args
        RuntimeCollection runtime_collection
    }

    if (args.run_clonal_decomposition) {
        scatter (sample in patient.samples) {
            # Force ABSOLUTE's purity (alpha) and ploidy (tau) only when the user supplied
            # both; a lone value is ignored so ABSOLUTE infers both (alpha and tau are
            # forced together or not at all).
            if (defined(sample.user_purity) && defined(sample.user_ploidy)) {
                Float forced_purity = select_first([sample.user_purity])
                Float forced_ploidy = select_first([sample.user_ploidy])
            }
            if (defined(sample.called_copy_ratio_segmentation) && defined(sample.af_model_parameters) && !defined(sample.absolute_acr_rdata) && !defined(sample.absolute_acr_plot)) {
                call abs.Absolute {
                    input:
                        acs_conversion_script = args.script_acs_conversion,
                        sample_name = sample.name,
                        copy_ratio_segmentation = select_first([sample.called_copy_ratio_segmentation]),
                        af_model_parameters = sample.af_model_parameters,
                        annotated_variants = sample.annotated_somatic_variants,
                        purity = forced_purity,
                        ploidy = forced_ploidy,
                        sex = patient.sex,
                        min_hets = args.absolute_min_hets,
                        min_probes = args.absolute_min_probes,
                        maf90_threshold = args.absolute_maf90_threshold,
                        genome_build = args.genome_build,
                        runtime_collection = runtime_collection
                }
            }

            File acs_cr_segmentation = select_first([Absolute.acs_copy_ratio_segmentation, sample.acs_copy_ratio_segmentation])
            Float acs_cr_skew = select_first([Absolute.acs_copy_ratio_skew, sample.acs_copy_ratio_skew])
            File? snv_maf = if defined(Absolute.snv_maf) then Absolute.snv_maf else sample.absolute_snv_maf
            File? indel_maf = if defined(Absolute.indel_maf) then Absolute.indel_maf else sample.absolute_indel_maf
            File acr_rdata = select_first([Absolute.acr_rdata, sample.absolute_acr_rdata])
            File acr_plot = select_first([Absolute.acr_plot, sample.absolute_acr_plot])
            File tcr_rdata = select_first([Absolute.tcr_rdata, sample.absolute_tcr_rdata])
            File tcr_plot = select_first([Absolute.tcr_plot, sample.absolute_tcr_plot])

            if (defined(sample.absolute_acr_solution)) {
                call abs_extract.AbsoluteExtract as AbsoluteACRExtract {
                    input:
                        map_to_absolute_copy_number_script = args.script_map_to_absolute_copy_number,
                        calculate_cancer_cell_fraction_script = args.script_calculate_cancer_cell_fraction,
                        sample_name = sample.name,
                        sex = patient.sex,
                        rdata = acr_rdata,
                        called_solution = select_first([sample.absolute_acr_solution]),
                        analyst_id = args.analyst_id,
                        copy_ratio_type = "allelic",
                        acs_copy_ratio_segmentation = acs_cr_segmentation,
                        acs_copy_ratio_skew = acs_cr_skew,
                        snv_maf = snv_maf,
                        indel_maf = indel_maf,
                        gvcf = patient.gvcf,
                        genome_build = args.genome_build,
                        runtime_collection = runtime_collection
                }
            }

            if (defined(sample.absolute_tcr_solution)) {
                call abs_extract.AbsoluteExtract as AbsoluteTCRExtract {
                    input:
                        map_to_absolute_copy_number_script = args.script_map_to_absolute_copy_number,
                        calculate_cancer_cell_fraction_script = args.script_calculate_cancer_cell_fraction,
                        sample_name = sample.name,
                        sex = patient.sex,
                        rdata = tcr_rdata,
                        called_solution = select_first([sample.absolute_tcr_solution]),
                        analyst_id = args.analyst_id,
                        copy_ratio_type = "total",
                        acs_copy_ratio_segmentation = acs_cr_segmentation,
                        acs_copy_ratio_skew = acs_cr_skew,
                        snv_maf = snv_maf,
                        indel_maf = indel_maf,
                        gvcf = patient.gvcf,
                        genome_build = args.genome_build,
                        runtime_collection = runtime_collection
                }
            }
        }

        # ABSOLUTE input
        if (length(select_all(snv_maf)) > 0) {
            Array[File] abs_snv_maf = select_all(snv_maf)
        }
        if (length(select_all(indel_maf)) > 0) {
            Array[File] abs_indel_maf = select_all(indel_maf)
        }

        # ACR
        if (length(select_all(AbsoluteACRExtract.absolute_maf)) > 0) {
            Array[File] abs_acr_maf = select_all(AbsoluteACRExtract.absolute_maf)
        }
        if (length(select_all(AbsoluteACRExtract.absolute_segtab)) > 0) {
            Array[File] abs_acr_segtab = select_all(AbsoluteACRExtract.absolute_segtab)
        }
        if (length(select_all(AbsoluteACRExtract.absolute_segtab_igv)) > 0) {
            Array[File] abs_acr_segtab_igv = select_all(AbsoluteACRExtract.absolute_segtab_igv)
        }
        if (length(select_all(AbsoluteACRExtract.absolute_table)) > 0) {
            Array[File] abs_acr_table = select_all(AbsoluteACRExtract.absolute_table)
        }
        if (length(select_all(AbsoluteACRExtract.absolute_purity)) > 0) {
            Array[Float] abs_acr_purity = select_all(AbsoluteACRExtract.absolute_purity)
        }
        if (length(select_all(AbsoluteACRExtract.absolute_ploidy)) > 0) {
            Array[Float] abs_acr_ploidy = select_all(AbsoluteACRExtract.absolute_ploidy)
        }

        # TCR
        if (length(select_all(AbsoluteTCRExtract.absolute_maf)) > 0) {
            Array[File] abs_tcr_maf = select_all(AbsoluteTCRExtract.absolute_maf)
        }
        if (length(select_all(AbsoluteTCRExtract.absolute_segtab)) > 0) {
            Array[File] abs_tcr_segtab = select_all(AbsoluteTCRExtract.absolute_segtab)
        }
        if (length(select_all(AbsoluteTCRExtract.absolute_segtab_igv)) > 0) {
            Array[File] abs_tcr_segtab_igv = select_all(AbsoluteTCRExtract.absolute_segtab_igv)
        }
        if (length(select_all(AbsoluteTCRExtract.absolute_table)) > 0) {
            Array[File] abs_tcr_table = select_all(AbsoluteTCRExtract.absolute_table)
        }
        if (length(select_all(AbsoluteTCRExtract.absolute_purity)) > 0) {
            Array[Float] abs_tcr_purity = select_all(AbsoluteTCRExtract.absolute_purity)
        }
        if (length(select_all(AbsoluteTCRExtract.absolute_ploidy)) > 0) {
            Array[Float] abs_tcr_ploidy = select_all(AbsoluteTCRExtract.absolute_ploidy)
        }

        call p_update_s.UpdateSamples as AddAbsoluteResultsToSamples {
            input:
                patient = patient,
                acs_copy_ratio_segmentation = acs_cr_segmentation,
                acs_copy_ratio_skew = acs_cr_skew,
                absolute_snv_maf = abs_snv_maf,
                absolute_indel_maf = abs_indel_maf,

                absolute_acr_rdata = acr_rdata,
                absolute_acr_plot = acr_plot,
                absolute_acr_maf = abs_acr_maf,
                absolute_acr_segtab = abs_acr_segtab,
                absolute_acr_segtab_igv = abs_acr_segtab_igv,
                absolute_acr_table = abs_acr_table,
                absolute_acr_purity = abs_acr_purity,
                absolute_acr_ploidy = abs_acr_ploidy,

                absolute_tcr_rdata = tcr_rdata,
                absolute_tcr_plot = tcr_plot,
                absolute_tcr_maf = abs_tcr_maf,
                absolute_tcr_segtab = abs_tcr_segtab,
                absolute_tcr_segtab_igv = abs_tcr_segtab_igv,
                absolute_tcr_table = abs_tcr_table,
                absolute_tcr_purity = abs_tcr_purity,
                absolute_tcr_ploidy = abs_tcr_ploidy,
        }

        # Only run PhylogicNDT if there are MAFs with ccf annotation
        if (length(select_first([abs_acr_maf, []])) > 0) {
            scatter (sample in AddAbsoluteResultsToSamples.updated_patient.samples) {
                if (defined(sample.absolute_acr_maf) && defined(sample.absolute_acr_purity) && (sample.absolute_acr_purity > 0)) {
                    String? phylogic_sample_name = sample.name
                    File? sample_absolute_maf = sample.absolute_acr_maf
                    File? sample_absolute_segtab = sample.absolute_acr_segtab
                    Float? sample_purity = sample.absolute_acr_purity
                    Int? sample_timepoint = sample.timepoint
                }
            }
            if (args.phylogic_use_segtab && length(select_all(sample_absolute_segtab)) > 0) {
                Array[File]? phylogic_absolute_segtabs = select_all(sample_absolute_segtab)
            }
            if (length(select_all(sample_timepoint)) > 0) {
                Array[Int]? phylogic_timepoints = select_all(sample_timepoint)
            }

            if (length(select_all(phylogic_sample_name)) > 0) {
                call phylogicndt.PhylogicNDT {
                    input:
                        patient_id = AddAbsoluteResultsToSamples.updated_patient.name,
                        sample_names = select_all(phylogic_sample_name),
                        absolute_mafs = select_all(sample_absolute_maf),
                        absolute_segtabs = phylogic_absolute_segtabs,
                        absolute_purities = select_all(sample_purity),
                        timepoints = phylogic_timepoints,
                        use_indels = args.phylogic_use_indels,
                        impute_missing_snvs = args.phylogic_impute_missing_snvs,
                        min_coverage = args.phylogic_min_coverage,
                        driver_genes_file = args.files.phylogic_driver_genes_file,
                        focal_cnv_intervals = args.files.phylogic_focal_cnv_intervals,
                        genome_build = args.genome_build,
                        runtime_collection = runtime_collection
                }
            }
        }
    }

    output {
        Patient updated_patient = select_first([AddAbsoluteResultsToSamples.updated_patient, patient])
    }
}