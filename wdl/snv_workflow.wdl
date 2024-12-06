version development

import "patient.wdl" as p
import "patient.update_samples.wdl" as p_update_s
import "workflow_arguments.wdl" as wfargs
import "runtime_collection.wdl" as rtc
import "call_variants.wdl" as cv
import "filter_variants.wdl" as fv
import "collect_allelic_counts.wdl" as cac
import "harmonize_samples.wdl" as hs
import "annotate_variants.wdl" as av
import "tasks.wdl"
#import "calculate_tumor_mutation_burden.wdl" as tmb


workflow SNVWorkflow {
    input {
        Patient patient
        WorkflowArguments args
        RuntimeCollection runtime_collection
    }

    if (args.run_variant_calling) {
        call cv.CallVariants {
            input:
                patient = patient,
                args = args,
                runtime_collection = runtime_collection,
        }

        if (args.run_variant_filter) {
            call fv.FilterVariants {
                input:
                    patient = patient,
                    args = args,
                    vcf = CallVariants.vcf,
                    vcf_idx = CallVariants.vcf_idx,
                    mutect_stats = CallVariants.mutect_stats,
                    orientation_bias = CallVariants.orientation_bias,
                    runtime_collection = runtime_collection,
            }

            # subset CallVariants.bam to reads covering the FilteredVariants.somatic_vcf only
            if (length(CallVariants.bams) > 0) {
                call tasks.PrintReads as SomaticBam {
                    input:
                        ref_fasta = args.files.ref_fasta,
                        ref_fasta_index = args.files.ref_fasta_index,
                        ref_dict = args.files.ref_dict,
                        patient_name = patient.name,
                        bams = select_all(CallVariants.bams),
                        bais = select_all(CallVariants.bais),
                        vcf = FilterVariants.somatic_vcf,
                        vcf_idx = FilterVariants.somatic_vcf_idx,
                        runtime_params = runtime_collection.print_reads
                }
            }

            if (args.keep_germline && defined(FilterVariants.germline_vcf)) {
                # Collect allelic pileups for all putative germline sites that were
                # not yet collected via the coverage workflow, then merge them.
                # This allows for more sensitive aCR segmentation.
                call tasks.SelectVariants as SelectGermlineNotInResource {
                    input:
                        ref_fasta = args.files.ref_fasta,
                        ref_fasta_index = args.files.ref_fasta_index,
                        ref_dict = args.files.ref_dict,
                        vcf = select_first([FilterVariants.germline_vcf]),
                        vcf_idx = select_first([FilterVariants.germline_vcf_idx]),
                        interval_blacklist = args.files.common_germline_alleles,
                        interval_blacklist_idx = args.files.common_germline_alleles_idx,
                        compress_output = args.compress_output,
                        runtime_params = runtime_collection.select_variants
                }

                if (SelectGermlineNotInResource.num_selected_variants > 0) {
                    call cac.VcfToPileupVariants as GermlineVariantsNotInResource {
                        input:
                            vcf = SelectGermlineNotInResource.selected_vcf,
                            vcf_idx = SelectGermlineNotInResource.selected_vcf_idx,
                            runtime_params = runtime_collection.vcf_to_pileup_variants,
                    }

                    scatter (sample in patient.samples) {
                        scatter (sequencing_run in sample.sequencing_runs) {
                            call cac.CollectAllelicCounts as GermlineAllelicCounts {
                                input:
                                    scattered_interval_list = args.scattered_interval_list,
                                    bam = sequencing_run.bam,
                                    bai = sequencing_run.bai,
                                    is_paired_end = sequencing_run.is_paired_end,
                                    ref_dict = args.files.ref_dict,
                                    variants = GermlineVariantsNotInResource.variants,
                                    variants_idx = GermlineVariantsNotInResource.variants_idx,
                                    getpileupsummaries_extra_args = args.getpileupsummaries_extra_args,
                                    sample_name = sequencing_run.name,
                                    runtime_collection = runtime_collection,
                            }
                            String seq_sample_names = sample.name
                        }

                        Array[String] mgac_task_sample_names = (
                            if defined(sample.snppanel_pileups)
                            then flatten([seq_sample_names, [sample.name]])
                            else seq_sample_names
                        )
                        Array[File] mgac_task_allelic_counts = (
                            if defined(sample.snppanel_pileups)
                            then flatten([GermlineAllelicCounts.pileup_summaries, [select_first([sample.snppanel_pileups])]])
                            else GermlineAllelicCounts.pileup_summaries
                        )
                        call hs.MergeAllelicCounts as MergeGermlineAllelicCounts {
                            input:
                                ref_dict = args.files.ref_dict,
                                script = args.merge_pileups_script,
                                sample_names = mgac_task_sample_names,
                                allelic_counts = mgac_task_allelic_counts,
                                min_read_depth = args.min_snppanel_read_depth,
                                compress_output = args.compress_output,
                                runtime_params = runtime_collection.merge_allelic_counts,
                        }
                        # We select the first file since we only supplied one unique sample name, so all counts were merged into the same file.
                        File germline_allelic_counts_ = select_first(MergeGermlineAllelicCounts.merged_allelic_counts)
                    }

                    call p_update_s.UpdateSamples as ExtendPileupsForSamples {
                        input:
                            patient = patient,
                            snppanel_pileups = germline_allelic_counts_,
                    }
                    call p.UpdatePatient as AddGermlineAlleles {
                        input:
                            patient = ExtendPileupsForSamples.updated_patient,
                            rare_germline_alleles = GermlineVariantsNotInResource.variants,
                            rare_germline_alleles_idx = GermlineVariantsNotInResource.variants_idx
                    }
                }
            }

            if (FilterVariants.num_somatic_variants > 0) {
                scatter (sample in patient.samples) {
                    scatter (sequencing_run in sample.sequencing_runs) {
                        call cac.CollectAllelicCounts as SomaticAllelicCounts {
                            input:
                                scattered_interval_list = args.scattered_interval_list,
                                bam = sequencing_run.bam,
                                bai = sequencing_run.bai,
                                is_paired_end = sequencing_run.is_paired_end,
                                ref_dict = args.files.ref_dict,
                                vcf = FilterVariants.somatic_vcf,
                                vcf_idx = FilterVariants.somatic_vcf_idx,
                                getpileupsummaries_extra_args = args.getpileupsummaries_extra_args,
                                sample_name = sequencing_run.name + ".somatic",
                                runtime_collection = runtime_collection,
                        }
                        String somatic_seq_sample_names = sample.name + ".somatic"
                    }

                    call hs.MergeAllelicCounts as MergeSomaticAllelicCounts {
                        input:
                            ref_dict = args.files.ref_dict,
                            script = args.merge_pileups_script,
                            sample_names = somatic_seq_sample_names,
                            allelic_counts = SomaticAllelicCounts.pileup_summaries,
                            runtime_params = runtime_collection.merge_allelic_counts,
                    }
                    # We select the first file since we only supplied one unique sample name, so all counts were merged.
                    File? somatic_allelic_counts_ = select_first(MergeSomaticAllelicCounts.merged_allelic_counts)
                }
            }

            if (args.run_variant_annotation) {
                # The sample scatter needs to be outside of the call to AnnotateVariants
                # since cromwell shits the bed for piping optional inputs into a nested scatter.
                scatter (sample in patient.samples) {
                    if (sample.is_tumor && defined(patient.matched_normal_sample)) {
                        Sample matched_normal_sample = select_first([patient.matched_normal_sample])
                        String? matched_normal_sample_name = matched_normal_sample.name
                        String? matched_normal_bam_name = matched_normal_sample.bam_name
                    }

                    call av.AnnotateVariants {
                        input:
                            vcf = FilterVariants.somatic_vcf,
                            vcf_idx = FilterVariants.somatic_vcf_idx,
                            num_variants = FilterVariants.num_somatic_variants,
                            individual_id = patient.name,
                            tumor_sample_name = sample.name,
                            tumor_bam_name = sample.bam_name,
                            normal_sample_name = matched_normal_sample_name,
                            normal_bam_name = matched_normal_bam_name,
                            args = args,
                            runtime_collection = runtime_collection,
                    }
                }

                call p_update_s.UpdateSamples as AddAnnotatedVariantsToSamples {
                    input:
                        patient = select_first([AddGermlineAlleles.updated_patient, patient]),
                        annotated_variants = AnnotateVariants.annotated_variants,
                }
            }
        }

#        call tmb.CalculateTumorMutationBurden as TMB {
#
#        }
    }

    output {
        Patient updated_patient = select_first([AddAnnotatedVariantsToSamples.updated_patient, AddGermlineAlleles.updated_patient, patient])

        File? unfiltered_vcf = CallVariants.vcf
        File? unfiltered_vcf_idx = CallVariants.vcf_idx
        File? mutect_stats = CallVariants.mutect_stats
        File? somatic_calls_bam = SomaticBam.output_bam
        File? somatic_calls_bai = SomaticBam.output_bai

        File? orientation_bias = CallVariants.orientation_bias
        File? filtered_vcf = FilterVariants.filtered_vcf
        File? filtered_vcf_idx = FilterVariants.filtered_vcf_idx
        File? somatic_vcf = FilterVariants.somatic_vcf
        File? somatic_vcf_idx = FilterVariants.somatic_vcf_idx
        File? germline_vcf = FilterVariants.germline_vcf
        File? germline_vcf_idx = FilterVariants.germline_vcf_idx
        File? filtering_stats = FilterVariants.filtering_stats

        Array[File?]? germline_allelic_counts = germline_allelic_counts_
        Array[File?]? somatic_allelic_counts = somatic_allelic_counts_

        Array[File]? annotated_variants = AnnotateVariants.annotated_variants
        Array[File?]? annotated_variants_idx = AnnotateVariants.annotated_variants_idx

#        Array[File?]? tmb = TMB.tmb
    }
}
