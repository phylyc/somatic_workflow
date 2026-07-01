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
        Array[Int]? timepoints

        # CACHE
        Array[Array[File]]? callable_loci
        Array[Array[File]]? total_read_counts
        Array[Array[File]]? denoised_total_copy_ratios
        Array[Array[File]]? snppanel_allelic_pileup_summaries

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
        Array[File]? cr_plot
        Array[File]? acs_copy_ratio_segmentation
        Array[Float]? acs_copy_ratio_skew
        Array[File]? annotated_somatic_variants
        Array[File]? annotated_somatic_variants_idx
        Array[File]? absolute_snv_maf
        Array[File]? absolute_indel_maf

        Array[File]? absolute_acr_rdata
        Array[File]? absolute_acr_plot
        Array[Int]? absolute_acr_solution
        Array[File]? absolute_acr_maf
        Array[File]? absolute_acr_segtab
        Array[File]? absolute_acr_table
        Array[Float]? absolute_acr_purity
        Array[Float]? absolute_acr_ploidy

        Array[File]? absolute_tcr_rdata
        Array[File]? absolute_tcr_plot
        Array[Int]? absolute_tcr_solution
        Array[File]? absolute_tcr_maf
        Array[File]? absolute_tcr_segtab
        Array[File]? absolute_tcr_table
        Array[Float]? absolute_tcr_purity
        Array[Float]? absolute_tcr_ploidy

        # user-supplied purity/ploidy that force ABSOLUTE's alpha/tau (both or neither):
        Array[Float]? user_purity
        Array[Float]? user_ploidy

        # for the patient-level shards:
        # CACHE
        Array[File]? scattered_intervals_for_variant_calling
        Array[Int]? skip_shards
        Array[Int]? high_mem_shards
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
        # TODO: ADD peddy file names here
        # DONE
        String? ancestry_pred
        Float? ancestry_prob
        File? background_pca_table
        File? pca_plot

        File? mask_vcf
        File? mask_vcf_idx
        String? mask_name

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

    # Overlay the per-sequencing-run inputs in a single pass. Each provided input array
    # carries one entry per run (aligned by index) and sets the matching field; an absent
    # or empty input leaves the field unchanged. cnv_panel_of_normals follows the 0-byte
    # placeholder convention: a size-0 file means "no PoN", so the field is left undefined.
    # INVARIANT: any provided input array has length == length(seqruns_dsr).
    Array[File] annotated_target_intervals_arr = select_first([annotated_target_intervals, []])
    Array[File] cnv_panel_of_normals_arr = select_first([cnv_panel_of_normals, []])
    Array[Boolean] is_paired_end_arr = select_first([is_paired_end, []])
    Array[Boolean] use_for_tCR_arr = select_first([use_for_tCR, []])
    Array[Boolean] use_for_aCR_arr = select_first([use_for_aCR, []])
    Array[String] sample_names_arr = select_first([sample_names, []])
    Array[Int] timepoints_arr = select_first([timepoints, []])

    scatter (i in range(length(seqruns_dsr))) {
        SequencingRun seq_run_i = seqruns_dsr[i]
        if (length(cnv_panel_of_normals_arr) > 0) {
            if (size(cnv_panel_of_normals_arr[i]) > 0) {
                File this_cnv_panel_of_normals = cnv_panel_of_normals_arr[i]
            }
        }
        SequencingRun seq_run_with_metadata = object {
            name: seq_run_i.name,
            sample_name: if length(sample_names_arr) > 0 then sample_names_arr[i] else seq_run_i.sample_name,
            timepoint: if length(timepoints_arr) > 0 then timepoints_arr[i] else seq_run_i.timepoint,
            bam: seq_run_i.bam,
            bai: seq_run_i.bai,
            target_intervals: seq_run_i.target_intervals,
            annotated_target_intervals: if length(annotated_target_intervals_arr) > 0 then annotated_target_intervals_arr[i] else seq_run_i.annotated_target_intervals,
            cnv_panel_of_normals: if defined(this_cnv_panel_of_normals) then this_cnv_panel_of_normals else seq_run_i.cnv_panel_of_normals,
            is_paired_end: if length(is_paired_end_arr) > 0 then is_paired_end_arr[i] else seq_run_i.is_paired_end,
            use_for_tCR: if length(use_for_tCR_arr) > 0 then use_for_tCR_arr[i] else seq_run_i.use_for_tCR,
            use_for_aCR: if length(use_for_aCR_arr) > 0 then use_for_aCR_arr[i] else seq_run_i.use_for_aCR,
            callable_loci: seq_run_i.callable_loci,
            total_read_counts: seq_run_i.total_read_counts,
            denoised_total_copy_ratios: seq_run_i.denoised_total_copy_ratios,
            snppanel_allelic_pileup_summaries: seq_run_i.snppanel_allelic_pileup_summaries,
        }
    }
    Array[SequencingRun] seqruns_with_metadata = seq_run_with_metadata

    # GroupBy sample name:
    # We assume that sample_names and bam_names share the same uniqueness,
    # that is if the supplied sample name is the same for two input bams, then the
    # bam names should also be the same, and vice versa.

    Array[String] theses_sample_names = select_first([sample_names, bam_names])
    Array[Pair[String, Array[SequencingRun]]] sample_dict = as_pairs(collect_by_key(zip(theses_sample_names, seqruns_with_metadata)))

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
                    timepoint: pair.right[0].timepoint,
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
                        timepoint: pair.right[0].timepoint,
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

    # TODO: This technically has the potential to fail since the Patient object
    # requires Array[Shard] shards to be non-optional! In this workflow this will
    # always be defined though.
    if (defined(scattered_intervals_for_variant_calling)) {
        Array[File] this_scattered_intervals = select_first([scattered_intervals_for_variant_calling, []])
        scatter (pair in zip(range(length(this_scattered_intervals)), this_scattered_intervals)) {
            # replace by "contain(skip_shards, pair.left)" in future WDL versions
            Array[Int] this_skip_shards = select_first([skip_shards, []])
            if (length(this_skip_shards) > 0) {
                scatter (skip in this_skip_shards) {
                    if (pair.left == skip) {
                        Boolean this_skip_shard = true
                    }
                }
            }
            # replace by "contain(high_mem_shards, pair.left)" in future WDL versions
            Array[Int] this_high_mem_shards = select_first([high_mem_shards, []])
            if (length(this_high_mem_shards) > 0) {
                scatter (high_mem_shard in this_high_mem_shards) {
                    if (pair.left == high_mem_shard) {
                        Boolean this_high_mem_shard = true
                    }
                }
            }
            Shard shards = object {
                id: pair.left,
                intervals: pair.right,
                skip: length(select_all(select_first([this_skip_shard, []]))) > 0,
                is_high_mem: length(select_all(select_first([this_high_mem_shard, []]))) > 0,
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
        mask_vcf: mask_vcf,
        mask_vcf_idx: mask_vcf_idx,
        mask_name: mask_name,
        ancestry_pred: ancestry_pred,
        ancestry_prob: ancestry_prob,
        background_pca_table: background_pca_table,
        pca_plot: pca_plot
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

    scatter (sample in UpdateShards.updated_patient.samples) {
        Array[SequencingRun] seqruns = sample.sequencing_runs
    }

    # Overlay the per-run CACHE coverage fields (callable loci, read counts, denoised
    # copy ratios, allelic pileups) in a single nested pass over samples and their runs.
    # Each provided input is an Array[Array[File]] aligned by [sample][run]; an absent or
    # empty input leaves the field unchanged.
    Array[Array[File]] callable_loci_arr = select_first([callable_loci, [[]]])
    Array[Array[File]] total_read_counts_arr = select_first([total_read_counts, [[]]])
    Array[Array[File]] denoised_total_copy_ratios_arr = select_first([denoised_total_copy_ratios, [[]]])
    Array[Array[File]] snppanel_allelic_pileup_summaries_arr = select_first([snppanel_allelic_pileup_summaries, [[]]])

    scatter (sample_idx in range(length(seqruns))) {
        scatter (run_idx in range(length(seqruns[sample_idx]))) {
            SequencingRun existing_seq_run = seqruns[sample_idx][run_idx]
            SequencingRun seq_run_with_coverage = object {
                name: existing_seq_run.name,
                sample_name: existing_seq_run.sample_name,
                timepoint: existing_seq_run.timepoint,
                bam: existing_seq_run.bam,
                bai: existing_seq_run.bai,
                target_intervals: existing_seq_run.target_intervals,
                annotated_target_intervals: existing_seq_run.annotated_target_intervals,
                cnv_panel_of_normals: existing_seq_run.cnv_panel_of_normals,
                is_paired_end: existing_seq_run.is_paired_end,
                use_for_tCR: existing_seq_run.use_for_tCR,
                use_for_aCR: existing_seq_run.use_for_aCR,
                callable_loci: if length(flatten(callable_loci_arr)) > 0 then callable_loci_arr[sample_idx][run_idx] else existing_seq_run.callable_loci,
                total_read_counts: if length(flatten(total_read_counts_arr)) > 0 then total_read_counts_arr[sample_idx][run_idx] else existing_seq_run.total_read_counts,
                denoised_total_copy_ratios: if length(flatten(denoised_total_copy_ratios_arr)) > 0 then denoised_total_copy_ratios_arr[sample_idx][run_idx] else existing_seq_run.denoised_total_copy_ratios,
                snppanel_allelic_pileup_summaries: if length(flatten(snppanel_allelic_pileup_summaries_arr)) > 0 then snppanel_allelic_pileup_summaries_arr[sample_idx][run_idx] else existing_seq_run.snppanel_allelic_pileup_summaries,
            }
        }
    }
    Array[Array[SequencingRun]] seqruns_with_coverage = seq_run_with_coverage

    call p_update_s.UpdateSamples {
        input:
            patient = UpdateShards.updated_patient,
            sequencing_runs = seqruns_with_coverage,
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
            cr_plot = cr_plot,
            acs_copy_ratio_segmentation = acs_copy_ratio_segmentation,
            acs_copy_ratio_skew = acs_copy_ratio_skew,
            annotated_somatic_variants = annotated_somatic_variants,
            annotated_somatic_variants_idx = annotated_somatic_variants_idx,

            absolute_snv_maf = absolute_snv_maf,
            absolute_indel_maf = absolute_indel_maf,

            absolute_acr_rdata = absolute_acr_rdata,
            absolute_acr_plot = absolute_acr_plot,
            absolute_acr_solution = absolute_acr_solution,
            absolute_acr_maf = absolute_acr_maf,
            absolute_acr_segtab = absolute_acr_segtab,
            absolute_acr_table = absolute_acr_table,
            absolute_acr_purity = absolute_acr_purity,
            absolute_acr_ploidy = absolute_acr_ploidy,

            absolute_tcr_rdata = absolute_tcr_rdata,
            absolute_tcr_plot = absolute_tcr_plot,
            absolute_tcr_solution = absolute_tcr_solution,
            absolute_tcr_maf = absolute_tcr_maf,
            absolute_tcr_segtab = absolute_tcr_segtab,
            absolute_tcr_table = absolute_tcr_table,
            absolute_tcr_purity = absolute_tcr_purity,
            absolute_tcr_ploidy = absolute_tcr_ploidy,

            user_purity = user_purity,
            user_ploidy = user_ploidy
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

        # Use WDL's built-in write_lines() to safely dump arrays to temp files.
        # This creates 0-byte files if the input arrays are empty.
        # We then sort those files directly. This is equivalent to python's sorted().
        LC_ALL=C sort -u ~{write_lines(sample_names)} > all_sorted.txt
        LC_ALL=C sort -u ~{write_lines(normal_sample_names)} > normal_sorted.txt

        # Extract Tumor Samples
        comm -23 all_sorted.txt normal_sorted.txt > tumor_sample_names.txt

        # Extract Normal Samples
        comm -12 all_sorted.txt normal_sorted.txt > normal_sample_names.txt
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