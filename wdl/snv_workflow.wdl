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
            # Only collect SNPs since Indels or MNVs are too likely misclassified.
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
                    select_variants_extra_args = "--select-type-to-include SNP",
                    runtime_params = runtime_collection.select_variants
            }

            if (SelectGermlineNotInResource.num_selected_variants > 0) {
                scatter (sample in FilterVariants.updated_patient.samples) {
                    String bam_names = sample.bam_name
                    String sample_names = sample.name
                    File? allelic_pileup_summaries = sample.allelic_pileup_summaries
                }
                call cac.VcfToPileupVariants as GermlineVariantsNotInResource {
                    input:
                        vcf = SelectGermlineNotInResource.selected_vcf,
                        vcf_idx = SelectGermlineNotInResource.selected_vcf_idx,
                        sample_names = bam_names,
                        compress_output = args.compress_output,
                        runtime_params = runtime_collection.vcf_to_pileup_variants,
                }

                call hs.MergeAllelicCounts as MergeGermlineAllelicCounts {
                    input:
                        ref_dict = args.files.ref_dict,
                        script = args.script_merge_pileups,
                        sample_names = flatten([sample_names, sample_names]),
                        allelic_counts = select_all(flatten([allelic_pileup_summaries, GermlineVariantsNotInResource.pileups])),
                        compress_output = args.compress_output,
                        runtime_params = runtime_collection.merge_allelic_counts,
                }

                # sort output to match order of sample_names since glob doesn't guarantee order
                scatter (sample in FilterVariants.updated_patient.samples) {
                    scatter (allelic_count in MergeGermlineAllelicCounts.merged_allelic_counts) {
                        String this_sample_name = basename(basename(allelic_count, ".gz"), ".pileup")
                        if (sample.name == this_sample_name) {
                            File this_allelic_counts = allelic_count
                        }
                    }
                    Array[File] this_sample_allelic_counts = select_all(this_allelic_counts)
                }
                Array[File] sorted_allelic_counts = flatten(this_sample_allelic_counts)

                call p_update_s.UpdateSamples as ExtendAllelicPileups {
                    input:
                        patient = FilterVariants.updated_patient,
                        allelic_pileup_summaries = sorted_allelic_counts,
                }

                call p.UpdatePatient as AddGermlineAlleles {
                    input:
                        patient = ExtendAllelicPileups.updated_patient,
                        rare_germline_alleles = GermlineVariantsNotInResource.variants,
                        rare_germline_alleles_idx = GermlineVariantsNotInResource.variants_idx
                }
            }
        }

        if (args.run_variant_annotation) {
            # The sample scatter needs to be outside of the call to AnnotateVariants
            # since cromwell shits the bed for piping optional inputs into a nested scatter.
            Patient pat = select_first([AddGermlineAlleles.updated_patient, FilterVariants.updated_patient])
            scatter (sample in pat.samples) {
                if (size(sample.annotated_somatic_variants) == 0) {
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
                File annotated_somatic_variants = select_first([AnnotateVariants.annotated_variants, sample.annotated_somatic_variants])
                File? annotated_somatic_variants_idx = if defined(AnnotateVariants.annotated_variants_idx) then AnnotateVariants.annotated_variants_idx else sample.annotated_somatic_variants_idx
            }
            if (length(select_all(annotated_somatic_variants_idx)) > 0) {
                Array[File] annotated_variants_idx = select_all(annotated_somatic_variants_idx)
            }

            call p_update_s.UpdateSamples as AddAnnotatedVariantsToSamples {
                input:
                    patient = pat,
                    annotated_somatic_variants = annotated_somatic_variants,
                    annotated_somatic_variants_idx = annotated_variants_idx,
            }
        }
    }

#        call tmb.CalculateTumorMutationBurden as TMB {
#
#        }

    output {
        Patient updated_patient = select_first([
            AddAnnotatedVariantsToSamples.updated_patient,
            AddGermlineAlleles.updated_patient,
            FilterVariants.updated_patient,
            CallVariants.updated_patient,
            patient
        ])
    }
}
