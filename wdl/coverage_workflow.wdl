version development

import "sample.wdl" as s
import "patient.wdl" as p
import "patient.update_samples.wdl" as p_update_s
import "workflow_arguments.wdl" as wfargs
import "runtime_collection.wdl" as rtc
import "collect_read_counts.wdl" as crc
import "collect_allelic_counts.wdl" as cac
import "collect_covered_regions.wdl" as ccr
import "harmonize_samples.wdl" as hs
import "calculate_contamination.wdl" as cc


workflow CoverageWorkflow {
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

            if (args.run_collect_allelic_coverage) {
                call cac.CollectAllelicCounts {
                    input:
                        ref_dict = args.files.ref_dict,
                        bam = sequencing_run.bam,
                        bai = sequencing_run.bai,
                        is_paired_end = sequencing_run.is_paired_end,
                        sample_name = sample.name + ".snp_array",
                        interval_list = sequencing_run.target_intervals,
                        scattered_interval_list = args.scattered_interval_list,
                        variants = args.files.common_germline_alleles,
                        variants_idx = args.files.common_germline_alleles_idx,
                        getpileupsummaries_extra_args = args.getpileupsummaries_extra_args,
                        minimum_population_allele_frequency = args.min_snp_array_pop_af,
                        maximum_population_allele_frequency = args.max_snp_array_pop_af,
                        minimum_read_depth = args.min_snp_array_read_depth,
                        runtime_collection = runtime_collection,
                }
            }

            File sample_bams = sequencing_run.bam
            File sample_bais = sequencing_run.bai

            if (select_first([sequencing_run.is_paired_end, false])) {
                Boolean run_is_paired_end = true
            }
        }

        # hacky way of implementing "ANY: Array[Boolean] -> Boolean"
        Boolean sample_is_paired_end = length(select_all(run_is_paired_end)) > 0

        if (args.run_collect_covered_regions) {
            call ccr.CollectCoveredRegions {
                input:
                    ref_fasta = args.files.ref_fasta,
                    ref_fasta_index = args.files.ref_fasta_index,
                    ref_dict = args.files.ref_dict,
                    sample_name = sample.name,
                    bams = sample_bams,
                    bais = sample_bais,
                    min_read_depth_threshold = args.min_read_depth,
                    is_paired_end = sample_is_paired_end,
                    output_format = "interval_list",
                    runtime_collection = runtime_collection,
            }
        }
    }

    if (args.run_collect_target_coverage) {
        # There is no way to harmonize the total read count data. This needs to be done
        # via a hierarchical model for the segmentation / copy ratio inference.
        Array[File]? read_counts = select_all(flatten(CollectReadCounts.read_counts))
    }

    if (args.run_collect_covered_regions) {
        Array[File]? covered_regions = select_all(CollectCoveredRegions.regions_interval_list)
    }

    # todo: FilterIntervals

    call hs.HarmonizeSamples {
        input:
            ref_dict = args.files.ref_dict,
            harmonize_copy_ratios_script = args.harmonize_copy_ratios_script,
            merge_pileups_script = args.merge_pileups_script,
            samples = patient.samples,
            denoised_copy_ratios = CollectReadCounts.denoised_copy_ratios,
            allelic_counts = CollectAllelicCounts.pileup_summaries,
            harmonize_min_target_length = args.harmonize_min_target_length,
            compress_output = false,
            runtime_collection = runtime_collection,
    }

    call p_update_s.UpdateSamples as ConsensusPatient {
        input:
            patient = patient,
            covered_regions = covered_regions,
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
            Array[File] normal_pileups = select_all(n_pileups)
        }

        if (defined(normal_pileups)) {
            scatter (normal_pileup in select_first([normal_pileups])) {
                call cc.CalculateContamination as CalculateNormalContamination {
                    input:
                        tumor_pileups = normal_pileup,
                        runtime_collection = runtime_collection,
                }
            }

            # todo: Choose the normal with the greatest sequencing depth.
            File matched_normal_pileup = select_first(select_first([normal_pileups, []]))
        }

        scatter (tumor_pileup in tumor_pileups) {
            # The only reason for supplying the matched normal pileup is to select
            # sites that have been confidently genotyped as homozygous SNPs in the normal.
            call cc.CalculateContamination as CalculateTumorContamination {
                input:
                    tumor_pileups = tumor_pileup,
                    normal_pileups = matched_normal_pileup,
                    runtime_collection = runtime_collection,
            }
        }

        Array[File] contaminations = flatten([
            CalculateTumorContamination.contamination_table,
            select_first([CalculateNormalContamination.contamination_table, []])
        ])
        Array[File] af_segmentations = flatten([
            CalculateTumorContamination.segmentation,
            select_first([CalculateNormalContamination.segmentation, []])
        ])

        call p_update_s.UpdateSamples as AddContaminationToSamples {
            input:
                patient = ConsensusPatient.updated_patient,
                contaminations = contaminations,
                af_segmentations = af_segmentations,
        }
    }

    output {
        Patient updated_patient = select_first([AddContaminationToSamples.updated_patient, ConsensusPatient.updated_patient, patient])

        Array[File]? target_read_counts = read_counts
        Array[File]? denoised_copy_ratios = HarmonizeSamples.harmonized_denoised_copy_ratios
        Array[File]? snp_array_pileups = HarmonizeSamples.merged_allelic_counts
        Array[File?]? covered_regions_interval_list = CollectCoveredRegions.regions_interval_list
        Array[File]? contamination_tables = contaminations
        Array[File]? segmentation_tables = af_segmentations
    }
}