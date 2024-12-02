version development

import "sample.wdl" as s
import "patient.wdl" as p
import "patient.update_samples.wdl" as p_update_s
import "workflow_arguments.wdl" as wfargs
import "runtime_collection.wdl" as rtc
import "collect_read_counts.wdl" as crc
import "collect_allelic_counts.wdl" as cac
import "harmonize_samples.wdl" as hs
import "genotype_snp_panel.wdl" as gsp
import "model_segments.wdl" as ms
import "filter_segments.wdl" as fs


workflow CNVWorkflow {
    input {
        Patient patient
        WorkflowArguments args
        RuntimeCollection runtime_collection
    }

    # ModelSegments requires the allelic counts to be pulled down at the same
    # set of loci for all samples. GetPileupSummaries does not guarantee this,
    # however, GenotypeSNPArray enforces this. Hence, if we want to run the
    # copy-ratio segmentation workflow, then we also need to run the contamination
    # workflow. Estimating the contamination is helpful for genotyping. Only HETs
    # are supplied to ModelSegements, with high genotyping_homozygous_log_ratio_threshold
    # to ensure that all loci are taken as HETs.

    # Prepare input for the contamination model / SNP array genotying

    scatter (sample in patient.samples) {
        String sample_names = sample.name
    }
    scatter (tumor_sample in patient.tumor_samples) {
        File? t_pileups = tumor_sample.snppanel_pileups
    }
    Array[File] tumor_pileups = select_all(t_pileups)
    if (patient.has_normal) {
        scatter (normal_sample in patient.normal_samples) {
            File? n_pileups = normal_sample.snppanel_pileups
            String? n_sample_name = normal_sample.name
        }
        Array[File] normal_pileups = select_all(n_pileups)
        Array[String] normal_sample_names = select_all(n_sample_name)
    }

    # Run the contamination model

    if (length(tumor_pileups) > 0) {
        call gsp.GenotypeSNPPanel {
            input:
                genotype_variants_script = args.genotype_variants_script,
                scattered_interval_list = args.scattered_interval_list,
                ref_dict = args.files.ref_dict,
                individual_id = patient.name,
                sample_names = sample_names,
                normal_sample_names = normal_sample_names,
                tumor_pileups = tumor_pileups,
                normal_pileups = normal_pileups,
                common_germline_alleles = args.files.common_germline_alleles,
                common_germline_alleles_idx = args.files.common_germline_alleles_idx,
                rare_germline_alleles = patient.rare_germline_alleles,
                rare_germline_alleles_idx = patient.rare_germline_alleles_idx,
                genotype_variants_min_read_depth = args.min_snppanel_read_depth,
                genotype_variants_min_genotype_likelihood = args.genotype_variants_min_genotype_likelihood,
                genotype_variants_overdispersion = args.genotype_variants_overdispersion,
                genotype_variants_ref_bias = args.genotype_variants_ref_bias,
                genotype_variants_save_sample_genotype_likelihoods = true,
                compress_output = args.compress_output,
                runtime_collection = runtime_collection,
        }

        call p_update_s.UpdateSamples as AddPileupsAndContaminationToSamples {
            input:
                patient = patient,
                snppanel_pileups = GenotypeSNPPanel.sample_genotype_likelihoods,  # Careful: This is not technically in a pileup format!
                contaminations = GenotypeSNPPanel.contaminations,
                af_segmentations = GenotypeSNPPanel.segmentations,
        }

        call p.UpdatePatient as AddGVCFtoPatient {
            input:
                patient = AddPileupsAndContaminationToSamples.updated_patient,
                gvcf = GenotypeSNPPanel.genotyped_vcf,
                gvcf_idx = GenotypeSNPPanel.genotyped_vcf_idx,
        }
    }

    if (args.run_model_segments) {
        call ms.ModelSegments {
            input:
                patient = select_first([AddGVCFtoPatient.updated_patient, patient]),
                args = args,
                runtime_collection = runtime_collection,
        }

        # todo: fix germline filter
#        call fs.FilterSegments {
#            input:
#                patient = ModelSegments.updated_patient,
#                args = args,
#                runtime_collection = runtime_collection,
#        }
    }

    # todo: FuncotateSegments

    output {
        Patient updated_patient = select_first([ModelSegments.updated_patient, AddGVCFtoPatient.updated_patient, patient])
#        Patient updated_patient = select_first([FilterSegments.updated_patient, ModelSegments.updated_patient, AddGVCFtoPatient.updated_patient, patient])

        File? snppanel_genotyped_vcf = GenotypeSNPPanel.genotyped_vcf
        File? snppanel_genotyped_vcf_idx = GenotypeSNPPanel.genotyped_vcf_idx
        File? snppanel_ref_counts = GenotypeSNPPanel.ref_counts
        File? snppanel_alt_counts = GenotypeSNPPanel.alt_counts
        File? snppanel_other_alt_counts = GenotypeSNPPanel.other_alt_counts
        File? snppanel_sample_correlation = GenotypeSNPPanel.sample_correlation
        Array[File]? snppanel_sample_genotype_likelihoods = GenotypeSNPPanel.sample_genotype_likelihoods
        Array[File]? snppanel_pileups = GenotypeSNPPanel.pileups
        Array[File]? contamination_tables = GenotypeSNPPanel.contaminations
        Array[File]? segmentation_tables = GenotypeSNPPanel.segmentations

        Array[File]? snppanel_allelic_counts = ModelSegments.snppanel_allelic_counts

        File? modeled_segments = ModelSegments.modeled_segments
        Array[File]? cr_segmentations = ModelSegments.called_copy_ratio_segmentations
#        Array[File]? cr_segmentations = FilterSegments.filtered_called_copy_ratio_segmentations
        Array[File]? cr_plots = ModelSegments.cr_plots
        Array[File]? af_model_parameters = ModelSegments.af_model_final_parameters
        Array[File]? cr_model_parameters = ModelSegments.cr_model_final_parameters
    }
}