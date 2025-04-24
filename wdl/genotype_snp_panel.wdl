version development

import "runtime_collection.wdl" as rtc
import "calculate_contamination.wdl" as cc
import "genotype_variants.wdl" as gv


workflow GenotypeSNPPanel {
    input {
        Array[File]? scattered_interval_list

        File? ref_dict

        String patient_id
        String? sex
        Array[String] sample_names
        Array[File] tumor_pileups
        Array[File]? normal_pileups
        Array[String]? normal_sample_names

        File? common_germline_alleles  # SNP panel
        File? common_germline_alleles_idx
        File? rare_germline_alleles  # From variant calling
        File? rare_germline_alleles_idx

        String genotype_variants_script = "https://github.com/phylyc/genomics_workflows/raw/master/python/genotype.py"
        Boolean genotype_variants_save_sample_genotype_likelihoods = false

        Int genotype_variants_min_read_depth = 10
        Float genotype_variants_min_genotype_likelihood = 0.995
        Float genotype_variants_outlier_prior = 0.00001
        Int genotype_variants_overdispersion = 10
        Float genotype_variants_ref_bias = 1.05

        Boolean compress_output = true

        RuntimeCollection runtime_collection
    }

    if (defined(normal_pileups)) {
        scatter (normal_pileup in select_first([normal_pileups])) {
            call cc.CalculateContamination as CalculateNormalContamination {
                input:
                    scattered_interval_list = scattered_interval_list,
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
                scattered_interval_list = scattered_interval_list,
                tumor_pileups = tumor_pileup,
                normal_pileups = matched_normal_pileup,
                runtime_collection = runtime_collection,
        }
    }

    Array[File] common_germline_allele_pileups = flatten([
        tumor_pileups,
        select_first([normal_pileups, []])
    ])

    Array[File] contamination_tables = flatten([
        CalculateTumorContamination.contamination_table,
        select_first([CalculateNormalContamination.contamination_table, []])
    ])
    Array[File] segmentation_tables = flatten([
        CalculateTumorContamination.af_segmentation_table,
        select_first([CalculateNormalContamination.af_segmentation_table, []])
    ])

    call gv.GenotypeVariants as GenotypeSNPPanelVariants {
        input:
            script = genotype_variants_script,
            patient_id = patient_id,
            sex = sex,
            sample_names = sample_names,
            normal_sample_names = normal_sample_names,
            pileups = common_germline_allele_pileups,
            contamination_tables = contamination_tables,
            segmentation_tables = segmentation_tables,
            common_germline_alleles = common_germline_alleles,
            common_germline_alleles_idx = common_germline_alleles_idx,
            rare_germline_alleles = rare_germline_alleles,
            rare_germline_alleles_idx = rare_germline_alleles_idx,
            compress_output = compress_output,
            min_read_depth = genotype_variants_min_read_depth,
            min_genotype_likelihood = genotype_variants_min_genotype_likelihood,
            outlier_prior = genotype_variants_outlier_prior,
            overdispersion = genotype_variants_overdispersion,
            ref_bias = genotype_variants_ref_bias,
            select_hets = false,
            save_sample_genotype_likelihoods = genotype_variants_save_sample_genotype_likelihoods,
            runtime_collection = runtime_collection,
    }

    output {
        Array[File] pileups = common_germline_allele_pileups
        Array[File] contaminations = contamination_tables
        Array[File] segmentations = segmentation_tables

        File genotyped_vcf = GenotypeSNPPanelVariants.vcf
        File genotyped_vcf_idx = GenotypeSNPPanelVariants.vcf_idx
        File ref_counts = GenotypeSNPPanelVariants.ref_counts
        File alt_counts = GenotypeSNPPanelVariants.alt_counts
        File other_alt_counts = GenotypeSNPPanelVariants.other_alt_counts
        File sample_correlation = GenotypeSNPPanelVariants.sample_correlation
        Boolean samples_are_from_same_patient = GenotypeSNPPanelVariants.samples_are_from_same_patient
        Array[File]? sample_genotype_likelihoods = GenotypeSNPPanelVariants.sample_genotype_likelihoods
    }
}
