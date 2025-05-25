version development

import "patient.wdl" as p


workflow Output {
    input {
        Patient patient
    }

    scatter (sample in patient.samples) {
        scatter (seqrun in sample.sequencing_runs) {
            File? cl = seqrun.callable_loci
            File? trc = seqrun.total_read_counts
            File? dtcr = seqrun.denoised_total_copy_ratios
            File? saps = seqrun.snppanel_allelic_pileup_summaries
        }
        Array[File] out_callable_loci = select_all(cl)
        Array[File] out_total_read_counts = select_all(trc)
        Array[File] out_denoised_total_copy_ratios = select_all(dtcr)
        Array[File] out_snppanel_allelic_pileup_summaries = select_all(saps)
        String sample_name = sample.name
        File? out_harmonized_callable_loci = sample.harmonized_callable_loci
        File? out_harmonized_denoised_total_copy_ratios = sample.harmonized_denoised_total_copy_ratios
        File? out_harmonized_snppanel_allelic_pileup_summaries = sample.harmonized_snppanel_allelic_pileup_summaries
        File? out_contamination_table = sample.contamination_table
        File? out_af_segmentation_table = sample.af_segmentation_table
        File? out_allelic_pileup_summaries = sample.allelic_pileup_summaries
        File? out_aggregated_allelic_read_counts = sample.aggregated_allelic_read_counts
        Float? out_genotype_error_probabilities = sample.genotype_error_probabilities
        File? out_af_model_parameters = sample.af_model_parameters
        File? out_cr_model_parameters = sample.cr_model_parameters
        File? out_called_copy_ratio_segmentation = sample.called_copy_ratio_segmentation
        File? out_cr_plot = sample.cr_plot
        File? out_acs_copy_ratio_segmentation = sample.acs_copy_ratio_segmentation
        Float? out_acs_copy_ratio_skew = sample.acs_copy_ratio_skew
        File? out_annotated_somatic_variants = sample.annotated_somatic_variants
        File? out_annotated_somatic_variants_idx = sample.annotated_somatic_variants_idx
        File? out_absolute_acr_rdata = sample.absolute_acr_rdata
        File? out_absolute_acr_plot = sample.absolute_acr_plot
        File? out_absolute_snv_maf = sample.absolute_snv_maf
        File? out_absolute_indel_maf = sample.absolute_indel_maf
        Int? out_absolute_solution = sample.absolute_solution
        File? out_absolute_maf = sample.absolute_maf
        File? out_absolute_segtab = sample.absolute_segtab
        File? out_absolute_table = sample.absolute_table
        Float? out_purity = sample.purity
        Float? out_ploidy = sample.ploidy
    }

    if (length(flatten(out_callable_loci)) > 0) {
        Array[Array[File]] cl_out = out_callable_loci
    }
    if (length(flatten(out_total_read_counts)) > 0) {
        Array[Array[File]] trc_out = out_total_read_counts
    }
    if (length(flatten(out_denoised_total_copy_ratios)) > 0) {
        Array[Array[File]] dtcr_out = out_denoised_total_copy_ratios
    }
    if (length(flatten(out_snppanel_allelic_pileup_summaries)) > 0) {
        Array[Array[File]] snpaps_out = out_snppanel_allelic_pileup_summaries
    }

    if (length(select_all(out_harmonized_callable_loci)) > 0) {
        Array[File] hcl_out = select_all(out_harmonized_callable_loci)
    }
    if (length(select_all(out_harmonized_denoised_total_copy_ratios)) > 0) {
        Array[File] hdtrc_out = select_all(out_harmonized_denoised_total_copy_ratios)
    }
    if (length(select_all(out_harmonized_snppanel_allelic_pileup_summaries)) > 0) {
        Array[File] haps_out = select_all(out_harmonized_snppanel_allelic_pileup_summaries)
    }
    if (length(select_all(out_contamination_table)) > 0) {
        Array[File] ct_out = select_all(out_contamination_table)
    }
    if (length(select_all(out_af_segmentation_table)) > 0) {
        Array[File] afst_out = select_all(out_af_segmentation_table)
    }
    if (length(select_all(out_allelic_pileup_summaries)) > 0) {
        Array[File] aps_out = select_all(out_allelic_pileup_summaries)
    }
    if (length(select_all(out_aggregated_allelic_read_counts)) > 0) {
        Array[File] aarc_out = select_all(out_aggregated_allelic_read_counts)
    }
    if (length(select_all(out_genotype_error_probabilities)) > 0) {
        Array[Float] gep_out = select_all(out_genotype_error_probabilities)
    }
    if (length(select_all(out_af_model_parameters)) > 0) {
        Array[File] afmp_out = select_all(out_af_model_parameters)
    }
    if (length(select_all(out_cr_model_parameters)) > 0) {
        Array[File] crmp_out = select_all(out_cr_model_parameters)
    }
    if (length(select_all(out_called_copy_ratio_segmentation)) > 0) {
        Array[File] ccrs_out = select_all(out_called_copy_ratio_segmentation)
    }
    if (length(select_all(out_cr_plot)) > 0) {
        Array[File] crp_out = select_all(out_cr_plot)
    }
    if (length(select_all(out_acs_copy_ratio_segmentation)) > 0) {
        Array[File] acrs_out = select_all(out_acs_copy_ratio_segmentation)
    }
    if (length(select_all(out_acs_copy_ratio_skew)) > 0) {
        Array[Float] acrskew_out = select_all(out_acs_copy_ratio_skew)
    }
    if (length(select_all(out_annotated_somatic_variants)) > 0) {
        Array[File] asv_out = select_all(out_annotated_somatic_variants)
    }
    if (length(select_all(out_annotated_somatic_variants_idx)) > 0) {
        Array[File] asvi_out = select_all(out_annotated_somatic_variants_idx)
    }
    if (length(select_all(out_absolute_acr_rdata)) > 0) {
        Array[File] aar_out = select_all(out_absolute_acr_rdata)
    }
    if (length(select_all(out_absolute_acr_plot)) > 0) {
        Array[File] acr_plot_out = select_all(out_absolute_acr_plot)
    }
    if (length(select_all(out_absolute_snv_maf)) > 0) {
        Array[File] asm_out = select_all(out_absolute_snv_maf)
    }
    if (length(select_all(out_absolute_indel_maf)) > 0) {
        Array[File] aim_out = select_all(out_absolute_indel_maf)
    }
    if (length(select_all(out_absolute_solution)) > 0) {
        Array[Int] as_out = select_all(out_absolute_solution)
    }
    if (length(select_all(out_absolute_maf)) > 0) {
        Array[File] am_out = select_all(out_absolute_maf)
    }
    if (length(select_all(out_absolute_segtab)) > 0) {
        Array[File] ast_out = select_all(out_absolute_segtab)
    }
    if (length(select_all(out_absolute_table)) > 0) {
        Array[File] at_out = select_all(out_absolute_table)
    }
    if (length(select_all(out_purity)) > 0) {
        Array[Float] purity_out = select_all(out_purity)
    }
    if (length(select_all(out_ploidy)) > 0) {
        Array[Float] ploidy_out = select_all(out_ploidy)
    }

    scatter (shard in patient.shards) {
        File? raw_calls_mutect2_vcf = shard.raw_calls_mutect2_vcf
        File? raw_calls_mutect2_vcf_idx = shard.raw_calls_mutect2_vcf_idx
        File? raw_mutect2_stats = shard.raw_mutect2_stats
        File? raw_mutect2_bam_out = shard.raw_mutect2_bam_out
        File? raw_mutect2_bai_out = shard.raw_mutect2_bai_out
        File? raw_mutect2_artifact_priors = shard.raw_mutect2_artifact_priors
    }

    if (length(select_all(raw_calls_mutect2_vcf)) > 0) {
        Array[File] rcm2vcf_out = select_all(raw_calls_mutect2_vcf)
    }
    if (length(select_all(raw_calls_mutect2_vcf_idx)) > 0) {
        Array[File] rcm2vcf_idx_out = select_all(raw_calls_mutect2_vcf_idx)
    }
    if (length(select_all(raw_mutect2_stats)) > 0) {
        Array[File] rm2_stats_out = select_all(raw_mutect2_stats)
    }
    if (length(select_all(raw_mutect2_bam_out)) > 0) {
        Array[File] rm2bam_out = select_all(raw_mutect2_bam_out)
    }
    if (length(select_all(raw_mutect2_bai_out)) > 0) {
        Array[File] rm2bai_out = select_all(raw_mutect2_bai_out)
    }
    if (length(select_all(raw_mutect2_artifact_priors)) > 0) {
        Array[File] rm2ap_out = select_all(raw_mutect2_artifact_priors)
    }

    output {
        # for each sequencing run:
        # CACHE (as returned by the workflow)
        Array[Array[File]]? callable_loci = cl_out
        Array[Array[File]]? total_read_counts = trc_out
        Array[Array[File]]? denoised_total_copy_ratios = dtcr_out
        Array[Array[File]]? snppanel_allelic_pileup_summaries = snpaps_out

        # for each sample:
        # CACHE (as returned by the workflow)
        Array[String]? sample_names_ordered = sample_name
        Array[File]? harmonized_callable_loci = hcl_out
        Array[File]? harmonized_denoised_total_copy_ratios = hdtrc_out
        Array[File]? harmonized_snppanel_allelic_pileup_summaries = haps_out
        Array[File]? contamination_table = ct_out
        Array[File]? af_segmentation_table = afst_out
        Array[File]? allelic_pileup_summaries = aps_out
        Array[File]? aggregated_allelic_read_counts = aarc_out
        Array[Float]? genotype_error_probabilities = gep_out
        Array[File]? af_model_parameters = afmp_out
        Array[File]? cr_model_parameters = crmp_out
        Array[File]? called_copy_ratio_segmentation = ccrs_out
        Array[File]? cr_plot = crp_out
        Array[File]? acs_copy_ratio_segmentation = acrs_out
        Array[Float]? acs_copy_ratio_skew = acrskew_out
        Array[File]? annotated_somatic_variants = asv_out
        Array[File?]? annotated_somatic_variants_idx = asvi_out
        Array[File]? absolute_acr_rdata = aar_out
        Array[File]? absolute_acr_plot = acr_plot_out
        Array[File]? absolute_snv_maf = asm_out
        Array[File]? absolute_indel_maf = aim_out
        Array[Int]? absolute_solution = as_out
        Array[File]? absolute_maf = am_out
        Array[File]? absolute_segtab = ast_out
        Array[File]? absolute_table = at_out
        Array[Float]? purity = purity_out
        Array[Float]? ploidy = ploidy_out

        # for each shard:
        # CACHE (as returned by the workflow)
        Array[File]? raw_calls_mutect2_vcf_scattered = rcm2vcf_out
        Array[File]? raw_calls_mutect2_vcf_idx_scattered = rcm2vcf_idx_out
        Array[File]? raw_mutect2_stats_scattered = rm2_stats_out
        Array[File]? raw_mutect2_bam_out_scattered = rm2bam_out
        Array[File]? raw_mutect2_bai_out_scattered = rm2bai_out
        Array[File]? raw_mutect2_artifact_priors_scattered = rm2ap_out
    }
}