version development

import "patient.wdl" as p
import "patient.update_samples.wdl" as p_update_s
import "workflow_arguments.wdl" as wfargs
import "runtime_collection.wdl" as rtc
import "genotype_variants.wdl" as gv
import "model_segments.wdl" as ms
#import "filter_segments.wdl" as fs


workflow CNVWorkflow {
    input {
        Patient patient
        WorkflowArguments args
        RuntimeCollection runtime_collection
    }

    # ModelSegments requires the allelic counts to be pulled down at the same
    # set of loci for all samples. GetPileupSummaries does not guarantee this,
    # however, GenotypeVariants enforces this. Estimating the contamination is
    # helpful for genotyping.

    scatter (sample in patient.samples) {
        String sample_names = sample.name
        File? pileups = sample.allelic_pileup_summaries
        File? contaminations = sample.contamination_table
        File? af_segmentations = select_first([sample.called_copy_ratio_segmentation, sample.af_segmentation_table])
        File? af_model_params = sample.af_model_parameters
    }
    Array[File] gt_pileups = select_all(pileups)
    if (length(select_all(contaminations)) > 0) {
        Array[File] contamination_tables = select_all(contaminations)
    }
    if (length(select_all(af_segmentations)) > 0) {
        Array[File] segmentation_tables = select_all(af_segmentations)
    }
    if (length(select_all(af_model_params)) > 0) {
        Array[File] af_pre_model_parameters = select_all(af_model_params)
    }

    if (patient.has_normal) {
        scatter (normal_sample in patient.normal_samples) {
            String? n_sample_name = normal_sample.name
        }
        Array[String] normal_sample_names = select_all(n_sample_name)
    }

    if (length(gt_pileups) > 0) {
        call gv.GenotypeVariants {
            input:
                script = args.genotype_variants_script,
                patient_id = patient.name,
                sex = patient.sex,
                sample_names = sample_names,
                normal_sample_names = normal_sample_names,
                pileups = gt_pileups,
                contamination_tables = contamination_tables,
                segmentation_tables = segmentation_tables,
                af_model_parameters = af_pre_model_parameters,
                common_germline_alleles = args.files.common_germline_alleles,
                common_germline_alleles_idx = args.files.common_germline_alleles_idx,
                rare_germline_alleles = patient.rare_germline_alleles,
                rare_germline_alleles_idx = patient.rare_germline_alleles_idx,
                compress_output = args.compress_output,
                min_read_depth = args.min_snppanel_read_depth,
                min_genotype_likelihood = args.genotype_variants_min_genotype_likelihood,
                outlier_prior = args.genotype_variants_outlier_prior,
                overdispersion = args.genotype_variants_overdispersion,
                ref_bias = args.genotype_variants_ref_bias,
                select_hets = false,
                save_sample_genotype_likelihoods = true,
                verbose = true,
                runtime_collection = runtime_collection,
        }

        # todo: phase gvcf

        call p_update_s.UpdateSamples as AddPileupsToSamples {
            input:
                patient = patient,
                allelic_pileup_summaries = GenotypeVariants.sample_genotype_likelihoods,  # Careful: This is not technically in a pileup format!
        }

        call p.UpdatePatient as AddGVCFtoPatient {
            input:
                patient = AddPileupsToSamples.updated_patient,
                gvcf = GenotypeVariants.vcf,
                gvcf_idx = GenotypeVariants.vcf_idx,
                snp_ref_counts = GenotypeVariants.ref_counts,
                snp_alt_counts = GenotypeVariants.alt_counts,
                snp_other_alt_counts = GenotypeVariants.other_alt_counts,
                snp_sample_correlation = GenotypeVariants.sample_correlation,
        }
    }

    if (args.run_model_segments) {
        call ms.ModelSegments {
            input:
                patient = select_first([AddGVCFtoPatient.updated_patient, patient]),
                args = args,
                runtime_collection = runtime_collection,
        }

#        if (args.run_filter_segments) {
#            call fs.FilterSegments {
#                input:
#                    patient = ModelSegments.updated_patient,
#                    args = args,
#                    runtime_collection = runtime_collection,
#            }
#        }
    }

    # todo: FuncotateSegments

    output {
        Patient updated_patient = select_first([ModelSegments.updated_patient, AddGVCFtoPatient.updated_patient, patient])
    }
}