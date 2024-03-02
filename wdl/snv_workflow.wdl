version development

import "patient.wdl" as p
import "workflow_arguments.wdl" as wfargs
import "runtime_collection.wdl" as rtc
import "call_variants.wdl" as cv
import "filter_variants.wdl" as fv
import "collect_allelic_counts.wdl" as cac
import "harmonize_samples.wdl" as hs
import "annotate_variants.wdl" as av


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

            if (args.run_collect_called_variants_allelic_coverage) {
                scatter (sample in patient.samples) {
                    if (args.keep_germline) {
                        scatter (sequencing_run in sample.sequencing_runs) {
                            call cac.CollectAllelicCounts as GermlineAllelicCounts {
                                input:
                                    scattered_interval_list = args.scattered_interval_list,
                                    bam = sequencing_run.bam,
                                    bai = sequencing_run.bai,
                                    ref_dict = args.ref_dict,
                                    vcf = FilterVariants.germline_vcf,
                                    vcf_idx = FilterVariants.germline_vcf_idx,
                                    getpileupsummaries_extra_args = args.getpileupsummaries_extra_args,
                                    sample_name = sequencing_run.name + ".germline",
                                    runtime_collection = runtime_collection,
                            }
                            String germline_seq_sample_names = sample.name + ".germline"
                        }

                        call hs.MergeAllelicCounts as MergeGermlineAllelicCounts {
                            input:
                                sample_names = germline_seq_sample_names,
                                allelic_counts = GermlineAllelicCounts.pileup_summaries,
                                compress_output = args.compress_output,
                                runtime_params = runtime_collection.merge_allelic_counts,
                        }
                        # We select the first file since we only supplied one unique sample name, so all counts were merged.
                        File? germline_allelic_counts = select_first(MergeGermlineAllelicCounts.merged_allelic_counts)
                    }

                    scatter (sequencing_run in sample.sequencing_runs) {
                        call cac.CollectAllelicCounts as SomaticAllelicCounts {
                            input:
                                scattered_interval_list = args.scattered_interval_list,
                                bam = sequencing_run.bam,
                                bai = sequencing_run.bai,
                                ref_dict = args.ref_dict,
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
                            sample_names = somatic_seq_sample_names,
                            allelic_counts = SomaticAllelicCounts.pileup_summaries,
                            compress_output = args.compress_output,
                            runtime_params = runtime_collection.merge_allelic_counts,
                    }
                    # We select the first file since we only supplied one unique sample name, so all counts were merged.
                    File? somatic_allelic_counts = select_first(MergeSomaticAllelicCounts.merged_allelic_counts)
                }
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
                        vcf = select_first([FilterVariants.somatic_vcf, CallVariants.vcf]),
                        vcf_idx = select_first([FilterVariants.somatic_vcf_idx, CallVariants.vcf_idx]),
                        individual_id = patient.name,
                        tumor_sample_name = sample.name,
                        tumor_bam_name = sample.bam_name,
                        normal_sample_name = matched_normal_sample_name,
                        normal_bam_name = matched_normal_bam_name,
                        args = args,
                        runtime_collection = runtime_collection,
                }
            }
        }
    }

    output {
        File? unfiltered_vcf = CallVariants.vcf
        File? unfiltered_vcf_idx = CallVariants.vcf_idx
        File? mutect_stats = CallVariants.mutect_stats
        File? bam = CallVariants.bam
        File? bai = CallVariants.bai

        File? orientation_bias = CallVariants.orientation_bias
        File? filtered_vcf = FilterVariants.filtered_vcf
        File? filtered_vcf_idx = FilterVariants.filtered_vcf_idx
        File? somatic_vcf = FilterVariants.somatic_vcf
        File? somatic_vcf_idx = FilterVariants.somatic_vcf_idx
        File? germline_vcf = FilterVariants.germline_vcf
        File? germline_vcf_idx = FilterVariants.germline_vcf_idx
        File? filtering_stats = FilterVariants.filtering_stats

        Array[File?]? called_germline_allelic_counts = germline_allelic_counts
        Array[File?]? called_somatic_allelic_counts = somatic_allelic_counts

        Array[File]? annotated_variants = AnnotateVariants.annotated_variants
        Array[File?]? annotated_variants_idx = AnnotateVariants.annotated_variants_idx

    }
}
