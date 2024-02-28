version development

import "patient.wdl" as p
import "runtimes.wdl" as rt
import "tasks.wdl"
import "multi-sample_somatic_workflow.wdl" as mssw
import "phase_variants.wdl" as phv


workflow SomaticWorkflow {
    input {
        File sif
        File snp_panel
    }

    call SIF {
        input:
            sif = sif
    }

    scatter (patient in SIF.patients) {
        call mssw.MultiSampleSomaticWorkflow as MSSW {
            input:
                patient = patient
        }
    }

    call MergeVCFs {
        input:
            vcfs = MSSW.genotyped_snparray_vcf,
            vcfs_idx = MSSW.genotyped_snparray_vcf_idx
    }

    call phv.PhaseVariants {
        input:
            vcf = MergeVCFs.vcf,
            shapeit_reference = snp_panel
    }

    scatter (patient in SIF.patients) {
        call MSACS {
            input:
                patient = patient,
                vcf = PhaseVariants.vcf,
                ref_counts = MSSW.snparray_ref_counts,
                alt_counts = MSSW.snparray_alt_counts,
                other_alt_counts = MSSW.snparray_other_alt_counts
        }
    }

    call Somix {
        input:
            recipe = "tCR_segmentation",
            interval_list = SplitIntervals.preprocessed_interval_list,
            read_counts = MSSW.target_read_counts,
            allelic_segmentation = MSACS.allelic_segmentation,
    }

    scatter (patient in SIF.patients) {
        call MSABS {
            input:
                recipe = "ICA_factorization",
                tcr_segmentation = Somix.segmentation,
                somatic_variants = VariantCall.variants
        }

        call Phylogic {
            input:
                somatic_variants = MSSW.called_somatic_allelic_counts,
                germline_variants = PhaseAlleles.variants,
                cn_segmentation = MSABS.segmentations,
        }
    }
}
