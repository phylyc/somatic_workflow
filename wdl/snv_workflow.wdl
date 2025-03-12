version development

import "sequencing_run.wdl" as seqrun
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
        if (!defined(patient.raw_snv_calls_vcf)) {
            call cv.CallVariants {
                input:
                    patient = patient,
                    args = args,
                    runtime_collection = runtime_collection,
            }
        }

        if (args.run_variant_filter) {
            call fv.FilterVariants {
                input:
                    patient = select_first([CallVariants.updated_patient, patient]),
                    args = args,
                    runtime_collection = runtime_collection,
            }

            if (args.keep_germline && defined(FilterVariants.updated_patient.germline_vcf)) {
                # Collect allelic pileups for all putative germline sites that were
                # not yet collected via the coverage workflow, then merge them.
                # This allows for more sensitive aCR segmentation.
                call tasks.SelectVariants as SelectGermlineNotInResource {
                    input:
                        ref_fasta = args.files.ref_fasta,
                        ref_fasta_index = args.files.ref_fasta_index,
                        ref_dict = args.files.ref_dict,
                        vcf = select_first([FilterVariants.updated_patient.germline_vcf]),
                        vcf_idx = select_first([FilterVariants.updated_patient.germline_vcf_idx]),
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

                            call seqrun.UpdateSequencingRun as SeqAddRareGermlinePileups {
                                input:
                                    sequencing_run = sequencing_run,
                                    rare_germline_allelic_pileup_summaries = GermlineAllelicCounts.pileup_summaries,
                            }
                        }

                        Array[String] mgac_task_sample_names = (
                            if defined(sample.allelic_pileup_summaries)
                            then flatten([seq_sample_names, [sample.name]])
                            else seq_sample_names
                        )
                        Array[File] mgac_task_allelic_counts = (
                            if defined(sample.allelic_pileup_summaries)
                            then flatten([GermlineAllelicCounts.pileup_summaries, [select_first([sample.allelic_pileup_summaries])]])
                            else GermlineAllelicCounts.pileup_summaries
                        )
                        call hs.MergeAllelicCounts as MergeGermlineAllelicCounts {
                            input:
                                ref_dict = args.files.ref_dict,
                                script = args.merge_pileups_script,
                                sample_names = mgac_task_sample_names,
                                allelic_counts = mgac_task_allelic_counts,
                                compress_output = args.compress_output,
                                runtime_params = runtime_collection.merge_allelic_counts,
                        }
                        # We select the first file since we only supplied one unique sample name, so all counts were merged into the same file.
                        File extended_allelic_pileup_summaries = select_first(MergeGermlineAllelicCounts.merged_allelic_counts)
                    }
                    Array[File] rg_acs = select_all(flatten(GermlineAllelicCounts.pileup_summaries))

                    call p_update_s.UpdateSamples as ExtendAllelicPileups {
                        input:
                            patient = FilterVariants.updated_patient,
                            sequencing_runs = SeqAddRareGermlinePileups.updated_sequencing_run,
                            allelic_pileup_summaries = extended_allelic_pileup_summaries,
                    }
                }
            }

            call p.UpdatePatient as AddGermlineAlleles {
                input:
                    patient = select_first([ExtendAllelicPileups.updated_patient, FilterVariants.updated_patient]),
                    rare_germline_alleles = GermlineVariantsNotInResource.variants,
                    rare_germline_alleles_idx = GermlineVariantsNotInResource.variants_idx
            }

            if (args.run_variant_annotation) {
                # The sample scatter needs to be outside of the call to AnnotateVariants
                # since cromwell shits the bed for piping optional inputs into a nested scatter.
                Patient pat = select_first([AddGermlineAlleles.updated_patient, FilterVariants.updated_patient])
                scatter (sample in pat.samples) {
                    if (sample.is_tumor && defined(pat.matched_normal_sample)) {
                        Sample matched_normal_sample = select_first([pat.matched_normal_sample])
                        String? matched_normal_sample_name = matched_normal_sample.name
                        String? matched_normal_bam_name = matched_normal_sample.bam_name
                    }

                    call av.AnnotateVariants {
                        input:
                            vcf = select_first([pat.somatic_vcf]),
                            vcf_idx = select_first([pat.somatic_vcf_idx]),
                            num_variants = pat.num_somatic_variants,
                            individual_id = patient.name,
                            tumor_sample_name = sample.name,
                            tumor_bam_name = sample.bam_name,
                            normal_sample_name = matched_normal_sample_name,
                            normal_bam_name = matched_normal_bam_name,
                            args = args,
                            runtime_collection = runtime_collection,
                    }
                }
                if (length(select_all(AnnotateVariants.annotated_variants_idx)) > 0) {
                    Array[File] annotated_variants_idx = select_all(AnnotateVariants.annotated_variants_idx)
                }

                call p_update_s.UpdateSamples as AddAnnotatedVariantsToSamples {
                    input:
                        patient = select_first([AddGermlineAlleles.updated_patient, patient]),
                        annotated_somatic_variants = AnnotateVariants.annotated_variants,
                        annotated_somatic_variants_idx = annotated_variants_idx,
                }
            }
        }

#        call tmb.CalculateTumorMutationBurden as TMB {
#
#        }
    }

    output {
        Patient updated_patient = select_first([AddAnnotatedVariantsToSamples.updated_patient, AddGermlineAlleles.updated_patient, CallVariants.updated_patient, patient])
    }
}
