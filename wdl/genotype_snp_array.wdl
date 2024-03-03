version development

import "runtime_collection.wdl" as rtc
import "calculate_contamination.wdl" as cc
import "genotype_variants.wdl" as gv


workflow GenotypeSNPArray {
    input {
        Array[File]? scattered_interval_list

        File? ref_dict

        String individual_id
        Array[String] sample_names
        Array[File] tumor_pileups
        Array[File]? normal_pileups

        File? common_germline_alleles  # SNP array
        File? common_germline_alleles_idx

        String genotype_variants_script = "https://github.com/phylyc/genomics_workflows/raw/master/python/genotype.py"
        Boolean genotype_variants_save_sample_genotype_likelihoods = false

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
        CalculateTumorContamination.segmentation,
        select_first([CalculateNormalContamination.segmentation, []])
    ])

    call gv.GenotypeVariants as GenotypeSNPArrayVariants {
        input:
            script = genotype_variants_script,
            individual_id = individual_id,
            sample_names = sample_names,
            pileups = common_germline_allele_pileups,
            contamination_tables = contamination_tables,
            segmentation_tables = segmentation_tables,
            common_germline_alleles = select_first([common_germline_alleles]),
            common_germline_alleles_idx = select_first([common_germline_alleles_idx]),
            compress_output = compress_output,
            select_hets = false,
            save_sample_genotype_likelihoods = genotype_variants_save_sample_genotype_likelihoods,
            verbose = true,
            runtime_collection = runtime_collection,
    }

    output {
        Array[File] pileups = common_germline_allele_pileups
        Array[File] contaminations = contamination_tables
        Array[File] segmentations = segmentation_tables

        File genotyped_vcf = GenotypeSNPArrayVariants.vcf
        File genotyped_vcf_idx = GenotypeSNPArrayVariants.vcf_idx
        File ref_counts = GenotypeSNPArrayVariants.ref_counts
        File alt_counts = GenotypeSNPArrayVariants.alt_counts
        File other_alt_counts = GenotypeSNPArrayVariants.other_alt_counts
        File sample_correlation = GenotypeSNPArrayVariants.sample_correlation
        Array[File]? sample_genotype_likelihoods = GenotypeSNPArrayVariants.sample_genotype_likelihoods
    }
}
