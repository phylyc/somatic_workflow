version development

import "sample.wdl" as s
import "shard.wdl" as sh


struct Patient {
    String name
    String? sex
    Array[Sample] samples
    Array[Sample] tumor_samples
    Array[Sample] normal_samples

    Boolean has_tumor
    Boolean has_normal
    Sample? matched_normal_sample

    # CACHE
    Array[Shard] shards
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
    File? modeled_segments
}


workflow UpdatePatient {
    input {
        Patient patient
        String? name
        Array[Sample]? samples
        Array[Sample]? tumor_samples
        Array[Sample]? normal_samples
        String? sex
        Boolean? has_tumor
        Boolean? has_normal
        Sample? matched_normal_sample

        # CACHE
        Array[Shard]? shards
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
        File? modeled_segments
    }

    Patient pat = object {
        name: select_first([name, patient.name]),
        samples: select_first([samples, patient.samples]),
        tumor_samples: select_first([tumor_samples, patient.tumor_samples]),
        normal_samples: select_first([normal_samples, patient.normal_samples]),
        has_tumor: select_first([has_tumor, patient.has_tumor]),
        has_normal: select_first([has_normal, patient.has_normal]),
        # cannot use select_first for optional fields:
        sex: if defined(sex) then sex else patient.sex,
        matched_normal_sample: if defined(matched_normal_sample) then matched_normal_sample else patient.matched_normal_sample,
        shards: if defined(shards) then shards else patient.shards,
        raw_snv_calls_vcf: if defined(raw_snv_calls_vcf) then raw_snv_calls_vcf else patient.raw_snv_calls_vcf,
        raw_snv_calls_vcf_idx: if defined(raw_snv_calls_vcf_idx) then raw_snv_calls_vcf_idx else patient.raw_snv_calls_vcf_idx,
        mutect2_stats: if defined(mutect2_stats) then mutect2_stats else patient.mutect2_stats,
        orientation_bias: if defined(orientation_bias) then orientation_bias else patient.orientation_bias,
        filtered_vcf: if defined(filtered_vcf) then filtered_vcf else patient.filtered_vcf,
        filtered_vcf_idx: if defined(filtered_vcf_idx) then filtered_vcf_idx else patient.filtered_vcf_idx,
        filtering_stats: if defined(filtering_stats) then filtering_stats else patient.filtering_stats,
        somatic_vcf: if defined(somatic_vcf) then somatic_vcf else patient.somatic_vcf,
        somatic_vcf_idx: if defined(somatic_vcf_idx) then somatic_vcf_idx else patient.somatic_vcf_idx,
        num_somatic_variants: if defined(num_somatic_variants) then num_somatic_variants else patient.num_somatic_variants,
        germline_vcf: if defined(germline_vcf) then germline_vcf else patient.germline_vcf,
        germline_vcf_idx: if defined(germline_vcf_idx) then germline_vcf_idx else patient.germline_vcf_idx,
        num_germline_variants: if defined(num_germline_variants) then num_germline_variants else patient.num_germline_variants,
        somatic_calls_bam: if defined(somatic_calls_bam) then somatic_calls_bam else patient.somatic_calls_bam,
        somatic_calls_bai: if defined(somatic_calls_bai) then somatic_calls_bai else patient.somatic_calls_bai,
        rare_germline_alleles: if defined(rare_germline_alleles) then rare_germline_alleles else patient.rare_germline_alleles,
        rare_germline_alleles_idx: if defined(rare_germline_alleles_idx) then rare_germline_alleles_idx else patient.rare_germline_alleles_idx,
        gvcf: if defined(gvcf) then gvcf else patient.gvcf,
        gvcf_idx: if defined(gvcf_idx) then gvcf_idx else patient.gvcf_idx,
        snp_ref_counts: if defined(snp_ref_counts) then snp_ref_counts else patient.snp_ref_counts,
        snp_alt_counts: if defined(snp_alt_counts) then snp_alt_counts else patient.snp_alt_counts,
        snp_other_alt_counts: if defined(snp_other_alt_counts) then snp_other_alt_counts else patient.snp_other_alt_counts,
        snp_sample_correlation: if defined(snp_sample_correlation) then snp_sample_correlation else patient.snp_sample_correlation,
        modeled_segments: if defined(modeled_segments) then modeled_segments else patient.modeled_segments
    }

    output {
        Patient updated_patient = pat
    }
}