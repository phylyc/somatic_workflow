version development

import "sample.wdl" as s
import "patient.wdl" as p
import "patient.update_samples.wdl" as p_update_s
import "workflow_arguments.wdl" as wfargs
import "runtime_collection.wdl" as rtc
import "collect_read_counts.wdl" as crc
import "collect_allelic_counts.wdl" as cac
import "harmonize_samples.wdl" as hs
import "genotype_snp_array.wdl" as gsa
import "model_segments.wdl" as ms


workflow CNVWorkflow {
    input {
        Patient patient
        WorkflowArguments args
        RuntimeCollection runtime_collection
    }

    scatter (sample in patient.samples) {
        scatter (sequencing_run in sample.sequencing_runs) {
            if (args.run_collect_target_coverage) {
                call crc.CollectReadCounts {
                    input:
                        ref_fasta = args.ref_fasta,
                        ref_fasta_index = args.ref_fasta_index,
                        ref_dict = args.ref_dict,
                        sample_name = sample.name,
                        sequencing_run = sequencing_run,
                        runtime_collection = runtime_collection,
                }
            }

            if (args.run_collect_allelic_coverage) {
                call cac.CollectAllelicCounts {
                    input:
                        ref_dict = args.ref_dict,
                        bam = sequencing_run.bam,
                        bai = sequencing_run.bai,
                        sample_name = sample.name + ".snp_array",
                        interval_list = sequencing_run.target_intervals,
                        scattered_interval_list = args.scattered_interval_list,
                        common_germline_alleles = args.common_germline_alleles,
                        common_germline_alleles_idx = args.common_germline_alleles_idx,
                        getpileupsummaries_extra_args = args.getpileupsummaries_extra_args,
                        minimum_population_allele_frequency = args.min_snp_array_pop_af,
                        maximum_population_allele_frequency = args.max_snp_array_pop_af,
                        minimum_read_depth = args.min_snp_array_read_depth,
                        runtime_collection = runtime_collection,
                }
            }
        }
    }

    if (args.run_collect_target_coverage) {
        # There is no way to harmonize the total read count data. This needs to be done
        # via a hierarchical model for the segmentation / copy ratio inference.
        # Careful: The there is one read count file per sequencing run, not per sample!
        Array[File]? read_counts = select_all(flatten(CollectReadCounts.read_counts))
    }

    call hs.HarmonizeSamples {
        input:
            samples = patient.samples,
            denoised_copy_ratios = CollectReadCounts.denoised_copy_ratios,
            allelic_counts = CollectAllelicCounts.pileup_summaries,
            compress_output = false,
            runtime_collection = runtime_collection,
    }

    call p_update_s.UpdateSamples as ConsensusPatient {
        input:
            patient = patient,
            denoised_copy_ratios = HarmonizeSamples.harmonized_denoised_copy_ratios,
            snp_array_pileups = HarmonizeSamples.merged_allelic_counts,
    }

    if (args.run_contamination_model) {
        scatter (sample in ConsensusPatient.updated_patient.samples) {
            String sample_names = sample.name
        }
        scatter (tumor_sample in ConsensusPatient.updated_patient.tumor_samples) {
            File? t_pileups = tumor_sample.snp_array_pileups
        }
        Array[File] tumor_pileups = select_all(t_pileups)
        if (patient.has_normal) {
            scatter (normal_sample in ConsensusPatient.updated_patient.normal_samples) {
                File? n_pileups = normal_sample.snp_array_pileups
            }
            Array[File]? normal_pileups = select_all(n_pileups)
        }
        call gsa.GenotypeSNPArray {
            input:
                scattered_interval_list = args.scattered_interval_list,
                ref_dict = args.ref_dict,
                individual_id = ConsensusPatient.updated_patient.name,
                sample_names = sample_names,
                tumor_pileups = tumor_pileups,
                normal_pileups = normal_pileups,
                common_germline_alleles = args.common_germline_alleles,
                common_germline_alleles_idx = args.common_germline_alleles_idx,
                getpileupsummaries_extra_args = args.getpileupsummaries_extra_args,
                genotype_variants_script = args.genotype_variants_script,
                genotype_variants_save_sample_genotype_likelihoods = args.genotype_variants_save_sample_genotype_likelihoods,
                compress_output = args.compress_output,
                runtime_collection = runtime_collection,
        }

        call p_update_s.UpdateSamples as AddPileupsAndContaminationToSamples {
            input:
                patient = ConsensusPatient.updated_patient,
                snp_array_pileups = GenotypeSNPArray.pileups,
                contaminations = GenotypeSNPArray.contaminations,
                af_segmentations = GenotypeSNPArray.segmentations,
        }
    }

    Patient updated_patient_ = select_first([AddPileupsAndContaminationToSamples.updated_patient, ConsensusPatient.updated_patient, patient])

    if (args.run_model_segments) {
        call ms.ModelSegments {
            input:
                patient = updated_patient_,
                args = args,
                runtime_collection = runtime_collection,
        }
    }

    output {
        Patient updated_patient = select_first([ModelSegments.updated_patient, updated_patient_])

        File? genotyped_snparray_vcf = GenotypeSNPArray.genotyped_vcf
        File? genotyped_snparray_vcf_idx = GenotypeSNPArray.genotyped_vcf_idx
        File? snparray_ref_counts = GenotypeSNPArray.ref_counts
        File? snparray_alt_counts = GenotypeSNPArray.alt_counts
        File? snparray_other_alt_counts = GenotypeSNPArray.other_alt_counts
        File? sample_snp_correlation = GenotypeSNPArray.sample_correlation
        Array[File]? sample_snparray_genotype_likelihoods = GenotypeSNPArray.sample_genotype_likelihoods
        Array[File]? snparray_pileups = GenotypeSNPArray.pileups
        Array[File]? snparray_allelic_counts = ModelSegments.snp_array_allelic_counts
        Array[File]? contamination_tables = GenotypeSNPArray.contaminations
        Array[File]? segmentation_tables = GenotypeSNPArray.segmentations

        File? modeled_segments = ModelSegments.modeled_segments
        Array[File]? cr_segments = ModelSegments.seg_final
        Array[File]? called_cr_segments = ModelSegments.called_cr_seg
        Array[File]? af_model_parameters = ModelSegments.af_model_final_parameters
        Array[File]? cr_model_parameters = ModelSegments.cr_model_final_parameters

        Array[File]? target_read_counts = read_counts
        Array[File]? denoised_copy_ratios = HarmonizeSamples.harmonized_denoised_copy_ratios
    }
}