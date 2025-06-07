version development

import "patient.wdl" as p


workflow MergePatients {
    input {
        Patient patient
        Patient other
    }

    Patient pat = object {
        name: patient.name,
        samples: patient.samples,
        tumor_samples: patient.tumor_samples,
        normal_samples: patient.normal_samples,
        has_tumor: patient.has_tumor,
        has_normal: patient.has_normal,
        # cannot use select_first for optional fields:
        sex: if defined(other.sex) then other.sex else patient.sex,
        matched_normal_sample: if defined(other.matched_normal_sample) then other.matched_normal_sample else patient.matched_normal_sample,
        shards: if defined(other.shards) then other.shards else patient.shards,
        raw_snv_calls_vcf: if defined(other.raw_snv_calls_vcf) then other.raw_snv_calls_vcf else patient.raw_snv_calls_vcf,
        raw_snv_calls_vcf_idx: if defined(other.raw_snv_calls_vcf_idx) then other.raw_snv_calls_vcf_idx else patient.raw_snv_calls_vcf_idx,
        mutect2_stats: if defined(other.mutect2_stats) then other.mutect2_stats else patient.mutect2_stats,
        orientation_bias: if defined(other.orientation_bias) then other.orientation_bias else patient.orientation_bias,
        filtered_vcf: if defined(other.filtered_vcf) then other.filtered_vcf else patient.filtered_vcf,
        filtered_vcf_idx: if defined(other.filtered_vcf_idx) then other.filtered_vcf_idx else patient.filtered_vcf_idx,
        filtering_stats: if defined(other.filtering_stats) then other.filtering_stats else patient.filtering_stats,
        somatic_vcf: if defined(other.somatic_vcf) then other.somatic_vcf else patient.somatic_vcf,
        somatic_vcf_idx: if defined(other.somatic_vcf_idx) then other.somatic_vcf_idx else patient.somatic_vcf_idx,
        num_somatic_variants: if defined(other.num_somatic_variants) then other.num_somatic_variants else patient.num_somatic_variants,
        germline_vcf: if defined(other.germline_vcf) then other.germline_vcf else patient.germline_vcf,
        germline_vcf_idx: if defined(other.germline_vcf_idx) then other.germline_vcf_idx else patient.germline_vcf_idx,
        num_germline_variants: if defined(other.num_germline_variants) then other.num_germline_variants else patient.num_germline_variants,
        somatic_calls_bam: if defined(other.somatic_calls_bam) then other.somatic_calls_bam else patient.somatic_calls_bam,
        somatic_calls_bai: if defined(other.somatic_calls_bai) then other.somatic_calls_bai else patient.somatic_calls_bai,
        rare_germline_alleles: if defined(other.rare_germline_alleles) then other.rare_germline_alleles else patient.rare_germline_alleles,
        rare_germline_alleles_idx: if defined(other.rare_germline_alleles_idx) then other.rare_germline_alleles_idx else patient.rare_germline_alleles_idx,
        gvcf: if defined(other.gvcf) then other.gvcf else patient.gvcf,
        gvcf_idx: if defined(other.gvcf_idx) then other.gvcf_idx else patient.gvcf_idx,
        snp_ref_counts: if defined(other.snp_ref_counts) then other.snp_ref_counts else patient.snp_ref_counts,
        snp_alt_counts: if defined(other.snp_alt_counts) then other.snp_alt_counts else patient.snp_alt_counts,
        snp_other_alt_counts: if defined(other.snp_other_alt_counts) then other.snp_other_alt_counts else patient.snp_other_alt_counts,
        snp_sample_correlation: if defined(other.snp_sample_correlation) then other.snp_sample_correlation else patient.snp_sample_correlation,
        snp_sample_correlation_min: if defined(other.snp_sample_correlation_min) then other.snp_sample_correlation_min else patient.snp_sample_correlation_min,
        modeled_segments: if defined(other.modeled_segments) then other.modeled_segments else patient.modeled_segments
    }

    output {
        Patient updated_patient = pat
    }
}