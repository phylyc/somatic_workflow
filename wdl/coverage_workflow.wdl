version development

import "sequencing_run.wdl" as seqrun
import "sample.wdl" as s
import "patient.wdl" as p
import "patient.update_samples.wdl" as p_update_s
import "workflow_arguments.wdl" as wfargs
import "runtime_collection.wdl" as rtc
import "collect_callable_loci.wdl" as ccl
import "collect_read_counts.wdl" as crc
import "collect_allelic_counts.wdl" as cac
import "harmonize_samples.wdl" as hs
import "calculate_contamination.wdl" as cc
import "genotype_variants.wdl" as gv
import "model_segments.wdl" as ms


workflow CoverageWorkflow {
    input {
        Patient patient
        WorkflowArguments args
        RuntimeCollection runtime_collection
    }

    scatter (sample in patient.samples) {
        scatter (sequencing_run in sample.sequencing_runs) {
            if (args.run_collect_callable_loci && (size(sequencing_run.callable_loci) == 0)) {
                call ccl.CollectCallableLoci {
                    input:
                        ref_fasta = args.files.ref_fasta,
                        ref_fasta_index = args.files.ref_fasta_index,
                        ref_dict = args.files.ref_dict,
                        sample_name = sample.name,
                        bam = sequencing_run.bam,
                        bai = sequencing_run.bai,
                        is_paired_end = sequencing_run.is_paired_end,
                        runtime_collection = runtime_collection,
                }
            }

            if (args.run_collect_total_read_counts && (size(sequencing_run.total_read_counts) == 0)) {
                call crc.CollectReadCounts {
                    input:
                        ref_fasta = args.files.ref_fasta,
                        ref_fasta_index = args.files.ref_fasta_index,
                        ref_dict = args.files.ref_dict,
                        sample_name = sample.name,
                        bam = sequencing_run.bam,
                        bai = sequencing_run.bai,
                        interval_list = sequencing_run.target_intervals,
                        annotated_interval_list = sequencing_run.annotated_target_intervals,
                        read_count_panel_of_normals = sequencing_run.cnv_panel_of_normals,
                        is_paired_end = sequencing_run.is_paired_end,
                        max_soft_clipped_bases = args.collect_read_counts_max_soft_clipped_bases,
                        runtime_collection = runtime_collection,
                }
            }

            if (args.run_collect_allelic_read_counts && (size(sequencing_run.snppanel_allelic_pileup_summaries) == 0)) {
                call cac.CollectAllelicCounts {
                    input:
                        ref_dict = args.files.ref_dict,
                        bam = sequencing_run.bam,
                        bai = sequencing_run.bai,
                        is_paired_end = sequencing_run.is_paired_end,
                        sample_name = sample.name + ".snppanel",
                        interval_list = sequencing_run.target_intervals,
                        scattered_interval_list = args.files.scattered_intervals,
                        variants = args.files.common_germline_alleles,
                        variants_idx = args.files.common_germline_alleles_idx,
                        getpileupsummaries_extra_args = args.getpileupsummaries_extra_args,
                        minimum_population_allele_frequency = args.min_snppanel_pop_af,
                        maximum_population_allele_frequency = args.max_snppanel_pop_af,
                        minimum_read_depth = args.min_snppanel_read_depth,
                        runtime_collection = runtime_collection,
                }
            }

            call seqrun.UpdateSequencingRun as SeqAddCoverage {
                input:
                    sequencing_run = sequencing_run,
                    callable_loci = CollectCallableLoci.bed,
                    total_read_counts = CollectReadCounts.read_counts,
                    denoised_total_copy_ratios = CollectReadCounts.denoised_copy_ratios,
                    snppanel_allelic_pileup_summaries = CollectAllelicCounts.pileup_summaries,
            }
        }
    }

    call p_update_s.UpdateSamples as PatientAddCoverage {
        input:
            patient = patient,
            sequencing_runs = SeqAddCoverage.updated_sequencing_run,
    }

    # todo: FilterIntervals

    call hs.HarmonizeSamples {
        input:
            ref_dict = args.files.ref_dict,
            harmonize_copy_ratios_script = args.harmonize_copy_ratios_script,
            merge_pileups_script = args.merge_pileups_script,
            samples = PatientAddCoverage.updated_patient.samples,
            harmonize_min_target_length = args.harmonize_min_target_length,
            pileups_min_read_depth = args.min_snppanel_read_depth,
            compress_output = false,
            runtime_collection = runtime_collection,
    }

    call p_update_s.UpdateSamples as ConsensusPatient {
        input:
            patient = patient,
            harmonized_callable_loci = HarmonizeSamples.harmonized_callable_loci,
            harmonized_denoised_total_copy_ratios = HarmonizeSamples.harmonized_denoised_copy_ratios,
            harmonized_snppanel_allelic_pileup_summaries = HarmonizeSamples.merged_allelic_counts,
            allelic_pileup_summaries = HarmonizeSamples.merged_allelic_counts,  # Will be overwritten later
    }

    if (defined(HarmonizeSamples.merged_allelic_counts) && args.run_contamination_model) {
        scatter (sample in ConsensusPatient.updated_patient.samples) {
            if (size(sample.contamination_table) == 0) {
                if (sample.is_tumor && defined(patient.matched_normal_sample)) {
                    Sample matched_normal_sample = select_first([patient.matched_normal_sample])
                    File? matched_normal_pileups = matched_normal_sample.harmonized_snppanel_allelic_pileup_summaries
                }

                call cc.CalculateContamination {
                    input:
                        tumor_pileups = sample.harmonized_snppanel_allelic_pileup_summaries,
                        normal_pileups = matched_normal_pileups,
                        runtime_collection = runtime_collection,
                }
            }

            File contamination_table = select_first([CalculateContamination.contamination_table, sample.contamination_table])
        }

        call p_update_s.UpdateSamples as AddContaminationToSamples {
            input:
                patient = ConsensusPatient.updated_patient,
                contamination_table = contamination_table,
        }

        # Perform a first-pass single-sample segmentation to get prior allelic
        # copy ratio segmentations for genotyping.
        call ms.ModelSegments as FirstPassSegmentation {
            input:
                patient = AddContaminationToSamples.updated_patient,
                args = args,
                runtime_collection = runtime_collection,
                pre_select_hets = false,
                gvcf = args.files.common_germline_alleles,
                gvcf_idx = args.files.common_germline_alleles_idx,
        }

        scatter (sample in FirstPassSegmentation.updated_patient.samples) {
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

        if (FirstPassSegmentation.updated_patient.has_normal) {
            scatter (normal_sample in FirstPassSegmentation.updated_patient.normal_samples) {
                String? n_sample_name = normal_sample.name
            }
            Array[String] normal_sample_names = select_all(n_sample_name)
        }

        call gv.GenotypeVariants {
            input:
                script = args.genotype_variants_script,
                patient_id = FirstPassSegmentation.updated_patient.name,
                sex = FirstPassSegmentation.updated_patient.sex,
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
                patient = FirstPassSegmentation.updated_patient,
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

        # Perform a second-pass multi-sample segmentation to get better priors for
        # allelic copy ratio segmentations: (FilterMutectCalls & GenotypeVariants)
        call ms.ModelSegments as SecondPassSegmentation {
            input:
                patient = AddGVCFtoPatient.updated_patient,
                args = args,
                runtime_collection = runtime_collection,
        }
    }

    output {
        Patient updated_patient = select_first([SecondPassSegmentation.updated_patient, ConsensusPatient.updated_patient])

        Array[File]? first_pass_cr_segmentations = FirstPassSegmentation.called_copy_ratio_segmentations
        Array[File]? first_pass_cr_plots = FirstPassSegmentation.cr_plots
        Array[File]? first_pass_af_model_parameters = FirstPassSegmentation.af_model_final_parameters
        Array[File]? first_pass_cr_model_parameters = FirstPassSegmentation.cr_model_final_parameters

        Array[File]? second_pass_cr_segmentations = SecondPassSegmentation.called_copy_ratio_segmentations
        Array[File]? second_pass_cr_plots = SecondPassSegmentation.cr_plots
        Array[File]? second_pass_af_model_parameters = SecondPassSegmentation.af_model_final_parameters
        Array[File]? second_pass_cr_model_parameters = SecondPassSegmentation.cr_model_final_parameters
    }
}