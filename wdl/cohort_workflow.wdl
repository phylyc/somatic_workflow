workflow SomaticWorkflow {
    input {
        File sif
        File snp_panel
    }

    call SIF {
        input:
            sif = sif
    }

    scatter (platform in SIF.platforms) {
        call SplitIntervals as PaddedScatteredIntervals {
            input:
                interval_list = platform.interval_list,
                padding = platform.padding
        }
        call SplitIntervals as UnpaddedScatteredIntervals {
            input:
                interval_list = platform.interval_list,
                padding = 0
        }

        scatter (sample in platform.samples) {
            call CollectReadCounts {
                input:
                    scattered_interval_list = PaddedScatteredIntervals.interval_files,
                    bam = sample.bam,
                    bai = sample.bai,
                    name = sample.name
            }
            call CollectAllelicCounts {
                input:
                    scattered_interval_list = UnpaddedScatteredIntervals.interval_files,
                    snp_panel = snp_panel,
                    bam = sample.bam,
                    bai = sample.bai,
                    name = sample.name
            }
        }
        call PhaseAlleles {
            input:
                allelic_counts = CollectAllelicCounts.allelic_counts,
                snp_panel = snp_panel
        }
        call Somix as MSACS {
            input:
                recipe = "segmentation",
                interval_list = PaddedIntervals.preprocessed_interval_list,
                read_counts = CollectReadCounts.read_counts,
                phased_allelic_counts = PhaseAlleles.phased_allelic_counts,
        }
    }

    scatter (patient in SIF.patients) {
        call MergeSegmentations {
            input:
                patient = patient,
                segmentations = MSACS.segmentations
        }

        # todo: circular: Need ICA_factorization for ploidy (variant call) and need somatic variants for ICA_factorization!
        call VariantCall {
            input:
                tumors = patient.tumor_samples,
                normals = patient.normal_samples,
                snp_panel = snp_panel,
                cr_segmentation = MergeSegmentations.cr_segmentation
        }

        call Somix as MSABS {
            input:
                recipe = "ICA_factorization",
                cr_segmentation = MergeSegmentations.cr_segmentation,
                somatic_variants = VariantCall.variants
        }

        call Phylogic {
            input:
                somatic_variants = VariantCall.variants,
                germline_variants = MergeSegmentations.cr_segmentation.hets,
                cn_segmentation = MSABS.segmentations,
        }
    }
}
