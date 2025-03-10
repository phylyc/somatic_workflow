version development

import "sequencing_run.wdl" as seq_run
import "sequencing_run.define.wdl" as seq_run_def
import "sample.wdl" as s
import "shard.wdl" as sh
import "patient.wdl" as p
import "patient.update_samples.wdl" as p_update_s
import "patient.update_shards.wdl" as p_update_sh
import "runtime_collection.wdl" as rtc


workflow DefinePatient {
    input {
        String name
        String? sex

        Array[String]? normal_sample_names
        # for each sequencing run:
        Array[File]+ bams
        Array[File]+ bais
        Array[File]+ target_intervals
        Array[File]? annotated_target_intervals
        Array[File]? cnv_panel_of_normals
        Array[Boolean]? is_paired_end
        Array[Boolean]? use_for_tCR
        Array[Boolean]? use_for_aCR
        Array[String]? sample_names
        # CACHE
        Array[File]? callable_loci
        Array[File]? total_read_counts
        Array[File]? denoised_total_copy_ratios
        Array[File]? snppanel_allelic_pileup_summaries
        Array[File]? rare_germline_allelic_pileup_summaries

        # for each sample:
        # CACHE (as returned by the workflow)
        Array[File]? harmonized_callable_loci
        Array[File]? harmonized_denoised_total_copy_ratios
        Array[File]? harmonized_snppanel_allelic_pileup_summaries
        Array[File]? contamination_table
        Array[File]? af_segmentation_table
        Array[File]? allelic_pileup_summaries
        Array[File]? aggregated_allelic_read_counts
        Array[Float]? genotype_error_probabilities
        Array[File]? af_model_parameters
        Array[File]? cr_model_parameters
        Array[File]? called_copy_ratio_segmentation
        Array[File]? acs_copy_ratio_segmentation
        Array[Float]? acs_copy_ratio_skew
        Array[File]? annotated_somatic_variants
        Array[File]? annotated_somatic_variants_idx
        Array[File]? absolute_acr_rdata
        Array[File]? absolute_acr_plot
        Array[Int]? absolute_solution
        Array[File]? absolute_maf
        Array[File]? absolute_segtab
        Array[File]? absolute_table
        Array[Float]? purity
        Array[Float]? ploidy

        # for the patient-level shards:
        # CACHE
        Array[File]? scattered_intervals
        Array[File]? raw_calls_mutect2_vcf_scattered
        Array[File]? raw_calls_mutect2_vcf_idx_scattered
        Array[File]? raw_mutect2_stats_scattered
        Array[File]? raw_mutect2_bam_out_scattered
        Array[File]? raw_mutect2_bai_out_scattered
        Array[File]? raw_mutect2_artifact_priors_scattered

        # for the patient:
        # CACHE
        File? raw_snv_calls_vcf
        File? raw_snv_calls_vcf_idx
        File? mutect2_stats
        File? orientation_bias
        File? filtered_vcf
        File? filtered_vcf_idx
        File? filtering_stats
        File? somatic_vcf
        File? somatic_vcf_idx
        Int? num_somatic_variants
        File? germline_vcf
        File? germline_vcf_idx
        Int? num_germline_variants
        File? somatic_calls_bam
        File? somatic_calls_bai
        File? rare_germline_alleles
        File? rare_germline_alleles_idx
        File? gvcf
        File? gvcf_idx
        File? snp_ref_counts
        File? snp_alt_counts
        File? snp_other_alt_counts
        File? snp_sample_correlation

        RuntimeCollection runtime_collection
    }

    Array[String] non_optional_normal_sample_names = select_first([normal_sample_names, []])
    Boolean has_normal = length(non_optional_normal_sample_names) > 0

    # We first define SequencingRuns for each bam, and then group them by sample name into Samples.

    scatter (tuple in transpose([bams, bais, target_intervals])) {
        call seq_run_def.DefineSequencingRun {
            input:
                bam = tuple[0],
                bai = tuple[1],
                target_intervals = tuple[2],
                runtime_collection = runtime_collection
        }
        String bam_names = DefineSequencingRun.sequencing_run.name
    }
    Array[SequencingRun] seqruns_dsr = DefineSequencingRun.sequencing_run

    if (defined(annotated_target_intervals)) {
        scatter (pair in zip(seqruns_dsr, select_first([annotated_target_intervals, []]))) {
            call seq_run.UpdateSequencingRun as UpdateAnnotatedTargetIntervals {
                input:
                    sequencing_run = pair.left,
                    annotated_target_intervals = pair.right,
            }
        }
    }
    Array[SequencingRun] seqruns_ati = select_first([UpdateAnnotatedTargetIntervals.updated_sequencing_run, seqruns_dsr])

    if (defined(cnv_panel_of_normals)) {
        scatter (pair in zip(seqruns_ati, select_first([cnv_panel_of_normals, []]))) {
            if (size(pair.right) > 0) {
                # For some sequencing platforms a panel of normals may not be available.
                # The denoise read counts task will then just use the anntated target
                # intervals to do GC correction. The convention here is to pass
                # a file of size 0B in place of the cnv PoN.
                File this_cnv_panel_of_normals = pair.right
            }
            call seq_run.UpdateSequencingRun as UpdateCnvPanelOfNormals {
                input:
                    sequencing_run = pair.left,
                    cnv_panel_of_normals = this_cnv_panel_of_normals,
            }
        }
    }
    Array[SequencingRun] seqruns_cpn = select_first([UpdateCnvPanelOfNormals.updated_sequencing_run, seqruns_ati])

    if (defined(is_paired_end)) {
        scatter (pair in zip(seqruns_cpn, select_first([is_paired_end, []]))) {
            call seq_run.UpdateSequencingRun as UpdateIsPairedEnd {
                input:
                    sequencing_run = pair.left,
                    is_paired_end = pair.right,
            }
        }
    }
    Array[SequencingRun] seqruns_ipe = select_first([UpdateIsPairedEnd.updated_sequencing_run, seqruns_cpn])

    if (defined(use_for_tCR)) {
        scatter (pair in zip(seqruns_ipe, select_first([use_for_tCR, []]))) {
            call seq_run.UpdateSequencingRun as UpdateUseForDCR {
                input:
                    sequencing_run = pair.left,
                    use_for_tCR = pair.right,
            }
        }
    }
    Array[SequencingRun] seqruns_udcr = select_first([UpdateUseForDCR.updated_sequencing_run, seqruns_ipe])

    if (defined(use_for_aCR)) {
        scatter (pair in zip(seqruns_udcr, select_first([use_for_aCR, []]))) {
            call seq_run.UpdateSequencingRun as UpdateUseForACR {
                input:
                    sequencing_run = pair.left,
                    use_for_aCR = pair.right,
            }
        }
    }
    Array[SequencingRun] seqruns_uacr = select_first([UpdateUseForACR.updated_sequencing_run, seqruns_udcr])

    if (defined(sample_names)) {
        scatter (pair in zip(seqruns_uacr, select_first([sample_names, []]))) {
            call seq_run.UpdateSequencingRun as UpdateSampleName {
                input:
                    sequencing_run = pair.left,
                    sample_name = pair.right,
            }
        }
    }
    Array[SequencingRun] seqruns_sn = select_first([UpdateSampleName.updated_sequencing_run, seqruns_uacr])

    if (defined(callable_loci)) {
        scatter (pair in zip(seqruns_sn, select_first([callable_loci, []]))) {
            call seq_run.UpdateSequencingRun as UpdateCallableLoci {
                input:
                    sequencing_run = pair.left,
                    callable_loci = pair.right,
            }
        }
    }
    Array[SequencingRun] seqruns_cl = select_first([UpdateCallableLoci.updated_sequencing_run, seqruns_sn])

    if (defined(total_read_counts)) {
        scatter (pair in zip(seqruns_cl, select_first([total_read_counts, []]))) {
            call seq_run.UpdateSequencingRun as UpdateTotalReadCounts {
                input:
                    sequencing_run = pair.left,
                    total_read_counts = pair.right,
            }
        }
    }
    Array[SequencingRun] seqruns_trc = select_first([UpdateTotalReadCounts.updated_sequencing_run, seqruns_cl])

    if (defined(denoised_total_copy_ratios)) {
        scatter (pair in zip(seqruns_trc, select_first([denoised_total_copy_ratios, []]))) {
            call seq_run.UpdateSequencingRun as UpdateDenoisedTotalCopyRatios {
                input:
                    sequencing_run = pair.left,
                    denoised_total_copy_ratios = pair.right,
            }
        }
    }
    Array[SequencingRun] seqruns_dtcr = select_first([UpdateDenoisedTotalCopyRatios.updated_sequencing_run, seqruns_trc])

    if (defined(snppanel_allelic_pileup_summaries)) {
        scatter (pair in zip(seqruns_dtcr, select_first([snppanel_allelic_pileup_summaries, []]))) {
            call seq_run.UpdateSequencingRun as UpdateAllelicPileupSummaries {
                input:
                    sequencing_run = pair.left,
                    snppanel_allelic_pileup_summaries = pair.right,
            }
        }
    }
    Array[SequencingRun] seqruns_saps = select_first([UpdateAllelicPileupSummaries.updated_sequencing_run, seqruns_dtcr])

    if (defined(rare_germline_allelic_pileup_summaries)) {
        scatter (pair in zip(seqruns_saps, select_first([rare_germline_allelic_pileup_summaries, []]))) {
            call seq_run.UpdateSequencingRun as UpdateRareGermlineAllelicPileupSummaries {
                input:
                    sequencing_run = pair.left,
                    rare_germline_allelic_pileup_summaries = pair.right,
            }
        }
    }
    Array[SequencingRun] seqruns_rgaps = select_first([UpdateRareGermlineAllelicPileupSummaries.updated_sequencing_run, seqruns_saps])

    # GroupBy sample name:
    # We assume that sample_names and bam_names share the same uniqueness,
    # that is if the supplied sample name is the same for two input bams, then the
    # bam names should also be the same, and vice versa.

    Array[String] theses_sample_names = select_first([sample_names, bam_names])
    Array[Pair[String, Array[SequencingRun]]] sample_dict = as_pairs(collect_by_key(zip(theses_sample_names, seqruns_rgaps)))

    # Pick tumor and normal samples apart:

    call GetUniqueSampleNameSets {
        input:
            sample_names = theses_sample_names,
            normal_sample_names = non_optional_normal_sample_names,
            runtime_params = runtime_collection.get_tumor_sample_names
    }

    scatter (tumor_sample_name in GetUniqueSampleNameSets.unique_tumor_sample_names) {
        scatter (pair in sample_dict) {
            if (pair.left == tumor_sample_name) {
                Sample selected_tumor_sample = object {
                    name: pair.left,
                    bam_name: pair.right[0].name,
                    sequencing_runs: pair.right,
                    is_tumor: true,
                }
            }
        }
        Sample tumor_samples = select_all(selected_tumor_sample)[0]
    }

    if (has_normal) {
        scatter (normal_sample_name in GetUniqueSampleNameSets.unique_normal_sample_names) {
            scatter (pair in sample_dict) {
                if (pair.left == normal_sample_name) {
                    Sample selected_normal_sample = object {
                        name: pair.left,
                        bam_name: pair.right[0].name,
                        sequencing_runs: pair.right,
                        is_tumor: false,
                    }
                }
            }
            Sample this_normal_samples = select_all(selected_normal_sample)[0]
        }
        Sample best_matched_normal_sample = select_first(this_normal_samples)
    }
    Array[Sample] normal_samples = select_first([this_normal_samples, []])

    if (defined(scattered_intervals)) {
        scatter (intervals in select_first([scattered_intervals, []])) {
            Shard shards = object {
                intervals: intervals,
            }
        }
    }

    Patient pat = object {
        name: name,
        sex: sex,
        samples: flatten([tumor_samples, normal_samples]),
        tumor_samples: tumor_samples,
        normal_samples: normal_samples,
        has_tumor: length(tumor_samples) > 0,
        has_normal: has_normal,
        matched_normal_sample: best_matched_normal_sample,
        # CACHE
        shards: shards,
        raw_snv_calls_vcf: raw_snv_calls_vcf,
        raw_snv_calls_vcf_idx: raw_snv_calls_vcf_idx,
        mutect2_stats: mutect2_stats,
        orientation_bias: orientation_bias,
        filtered_vcf: filtered_vcf,
        filtered_vcf_idx: filtered_vcf_idx,
        filtering_stats: filtering_stats,
        somatic_vcf: somatic_vcf,
        somatic_vcf_idx: somatic_vcf_idx,
        num_somatic_variants: num_somatic_variants,
        germline_vcf: germline_vcf,
        germline_vcf_idx: germline_vcf_idx,
        num_germline_variants: num_germline_variants,
        somatic_calls_bam: somatic_calls_bam,
        somatic_calls_bai: somatic_calls_bai,
        rare_germline_alleles: rare_germline_alleles,
        rare_germline_alleles_idx: rare_germline_alleles_idx,
        gvcf: gvcf,
        gvcf_idx: gvcf_idx,
        snp_ref_counts: snp_ref_counts,
        snp_alt_counts: snp_alt_counts,
        snp_other_alt_counts: snp_other_alt_counts,
    }

    call p_update_sh.UpdateShards {
        input:
            patient = pat,
            raw_calls_mutect2_vcf_scattered = raw_calls_mutect2_vcf_scattered,
            raw_calls_mutect2_vcf_idx_scattered = raw_calls_mutect2_vcf_idx_scattered,
            raw_mutect2_stats_scattered = raw_mutect2_stats_scattered,
            raw_mutect2_bam_out_scattered = raw_mutect2_bam_out_scattered,
            raw_mutect2_bai_out_scattered = raw_mutect2_bai_out_scattered,
            raw_mutect2_artifact_priors_scattered = raw_mutect2_artifact_priors_scattered,
    }

    call p_update_s.UpdateSamples {
        input:
            patient = UpdateShards.updated_patient,
            harmonized_callable_loci = harmonized_callable_loci,
            harmonized_denoised_total_copy_ratios = harmonized_denoised_total_copy_ratios,
            harmonized_snppanel_allelic_pileup_summaries = harmonized_snppanel_allelic_pileup_summaries,
            contamination_table = contamination_table,
            af_segmentation_table = af_segmentation_table,
            allelic_pileup_summaries = allelic_pileup_summaries,
            aggregated_allelic_read_counts = aggregated_allelic_read_counts,
            genotype_error_probabilities = genotype_error_probabilities,
            af_model_parameters = af_model_parameters,
            cr_model_parameters = cr_model_parameters,
            called_copy_ratio_segmentation = called_copy_ratio_segmentation,
            acs_copy_ratio_segmentation = acs_copy_ratio_segmentation,
            acs_copy_ratio_skew = acs_copy_ratio_skew,
            annotated_somatic_variants = annotated_somatic_variants,
            annotated_somatic_variants_idx = annotated_somatic_variants_idx,
            absolute_acr_rdata = absolute_acr_rdata,
            absolute_acr_plot = absolute_acr_plot,
            absolute_solution = absolute_solution,
            absolute_maf = absolute_maf,
            absolute_segtab = absolute_segtab,
            absolute_table = absolute_table,
            purity = purity,
            ploidy = ploidy
    }

    output {
        Patient patient = UpdateSamples.updated_patient
    }
}

task GetUniqueSampleNameSets {
    input {
        Array[String] sample_names
        Array[String] normal_sample_names
        Runtime runtime_params
    }

    String dollar = "$"

    command <<<
        set -euxo pipefail

        # Convert comma-separated strings to arrays
        IFS=',' read -r -a all_names <<< "~{sep="," sample_names}"
        IFS=',' read -r -a normal_names <<< "~{sep="," normal_sample_names}"

        # Create an associative array to hold unique names
        declare -A unique_tumor_names
        declare -A unique_normal_names

        # Loop through all names
        for name in "~{dollar}{all_names[@]}"; do
            # Check if name is in the normal names list
            if [[ " ~{dollar}{normal_names[*]} " =~ " ${name} " ]]; then
                unique_normal_names["$name"]=1
            else
                unique_tumor_names["$name"]=1
            fi
        done

        # Write unique tumor sample names to file
        touch "tumor_sample_names.txt"
        for tumor_name in "~{dollar}{!unique_tumor_names[@]}"; do
            echo "$tumor_name" >> "tumor_sample_names.txt"
        done

        # Write unique normal sample names to file
        touch "normal_sample_names.txt"
        for normal_name in "~{dollar}{!unique_normal_names[@]}"; do
            echo "$normal_name" >> "normal_sample_names.txt"
        done
    >>>

    output {
        Array[String] unique_tumor_sample_names = read_lines("tumor_sample_names.txt")
        Array[String] unique_normal_sample_names = read_lines("normal_sample_names.txt")
    }

    runtime {
        docker: runtime_params.docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: runtime_params.machine_mem + " MB"
        runtime_minutes: runtime_params.runtime_minutes
        disks: "local-disk " + runtime_params.disk + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }
}