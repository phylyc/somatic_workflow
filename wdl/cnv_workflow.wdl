version development

import "sample.wdl"
import "patient.wdl"
import "workflow_arguments.wdl"
import "runtimes.wdl"
import "collect_read_counts.wdl" as crc
import "collect_allelic_counts.wdl" as cac
import "merge_sequencing_runs.wdl" as msr
import "genotype_snp_array.wdl" as gsa
#import "model_segments.wdl" as ms


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
                        bam = sequencing_run.bam,
                        bai = sequencing_run.bai,
                        sample_name = sample.name,
                        interval_list = sequencing_run.target_intervals,
                        annotated_interval_list = sequencing_run.annotated_target_intervals,
                        read_count_panel_of_normals = sequencing_run.cnv_panel_of_normals,
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

        call msr.MergeSequencingRuns {
            input:
                sample = sample,
                read_counts = select_all(CollectReadCounts.read_counts),
                denoised_copy_ratios = select_all(CollectReadCounts.denoised_copy_ratios),
                standardized_copy_ratios = select_all(CollectReadCounts.standardized_copy_ratios),
                snp_array_allelic_counts = select_all(CollectAllelicCounts.pileup_summaries),
                runtime_collection = runtime_collection,
        }

        # Define samples to update patient with
        Sample samples = MergeSequencingRuns.updated_sample
        if (sample.is_tumor) {
            Sample tumor_samples = MergeSequencingRuns.updated_sample
        }
        if (!sample.is_tumor) {
            Sample normal_samples = MergeSequencingRuns.updated_sample
        }
        if (defined(patient.matched_normal_sample)) {
            Sample previous_matched_normal_sample = select_first([patient.matched_normal_sample])
            if (sample.name == previous_matched_normal_sample.name) {
                Sample matched_normal_sample = MergeSequencingRuns.updated_sample
            }
        }
    }

    call patient.UpdatePatient as ConsensusPatient {
        input:
            patient = patient,
            samples = select_all(samples),
            tumor_samples = select_all(tumor_samples),
            normal_samples = select_all(normal_samples),
            matched_normal_sample = if (patient.has_normal) then select_first(select_all(matched_normal_sample)) else patient.matched_normal_sample,
    }

    if (args.run_contamination_model) {
        scatter (sample in ConsensusPatient.updated_patient.samples) {
            String sample_names = sample.name
        }
        scatter (tumor_sample in ConsensusPatient.updated_patient.tumor_samples) {
            File? t_pileups = tumor_sample.snp_array_allelic_counts
        }
        Array[File] tumor_pileups = select_all(t_pileups)
        if (patient.has_normal) {
            scatter (normal_sample in ConsensusPatient.updated_patient.normal_samples) {
                File? n_pileups = normal_sample.snp_array_allelic_counts
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

        # Attach pileups, contamination, and segmentation tables to samples
        call patient.UpdateSamples as AddPileupsAndContaminationToSamples {
            input:
                patient = ConsensusPatient.updated_patient,
                snp_array_allelic_counts = GenotypeSNPArray.pileups,
                contaminations = GenotypeSNPArray.contaminations,
                af_segmentations = GenotypeSNPArray.segmentations,
        }
    }

#    call ms.ModelSegments {
#        input:
#            patient = AddPileupsAndContaminationToSamples.updated_patient,
#    }

#    call tasks.CallCopyRatioSegments {
#        input:
#            copy_ratios_segments = ModelSegments.copy_ratios_segments,
#    }

    output {
        Patient updated_patient = select_first([AddPileupsAndContaminationToSamples.updated_patient, ConsensusPatient.updated_patient, patient])

        File? genotyped_snparray_vcf = GenotypeSNPArray.genotyped_vcf
        File? genotyped_snparray_vcf_idx = GenotypeSNPArray.genotyped_vcf_idx
        File? snparray_ref_counts = GenotypeSNPArray.ref_counts
        File? snparray_alt_counts = GenotypeSNPArray.alt_counts
        File? snparray_other_alt_counts = GenotypeSNPArray.other_alt_counts
        File? sample_snp_correlation = GenotypeSNPArray.sample_correlation
        Array[File]? sample_snparray_genotype_likelihoods = GenotypeSNPArray.sample_genotype_likelihoods
        Array[File]? snparray_allelic_counts = GenotypeSNPArray.pileups
        Array[File]? contamination_tables = GenotypeSNPArray.contaminations
        Array[File]? segmentation_tables = GenotypeSNPArray.segmentations

        Array[File?]? target_read_counts = MergeSequencingRuns.merged_read_counts
        Array[File?]? denoised_copy_ratios = MergeSequencingRuns.merged_denoised_copy_ratios
        Array[File?]? standardized_copy_ratios = MergeSequencingRuns.merged_standardized_copy_ratios
    }
}