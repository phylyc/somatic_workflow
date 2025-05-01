version development

import "create_AF_only_germline_resource.wdl" as singleAFonly
import "../runtime_collection.wdl" as rtc
import "../tasks.wdl"


workflow CreateAFonlyVCFfromFile {
    input {
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        File vcfs_file
        File vcfs_idx_file
        File? minAF_file

        Float minAF = 0.05

        Boolean compress_output = true

        String af_only_output_name = "af_only.filtered"
        String af_only_common_output_name = "af_only.filtered.common.biallelic"

        RuntimeCollection runtime_collection = RuntimeParameters.rtc
    }

    call rtc.DefineRuntimeCollection as RuntimeParameters

    Array[File]+ vcfs = read_lines(vcfs_file)
    Array[File]+ vcfs_idx = read_lines(vcfs_idx_file)

    if (defined(minAF_file)) {
        Array[String]+ minAFs_from_file = read_lines(select_first([minAF_file]))
    }
    scatter (vcf in vcfs) {
        String default_minAFs = minAF
    }
    Array[String] minAFs = select_first([minAFs_from_file, default_minAFs])

    scatter (pair in zip(zip(vcfs, vcfs_idx), minAFs)) {
        call singleAFonly.CreateAFonlyVCF {
            input:
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                ref_dict = ref_dict,
                vcf = pair.left.left,
                vcf_idx = pair.left.right,
                minAF = pair.right,
                compress_output = compress_output,
                runtime_collection = runtime_collection
        }
    }

    call tasks.GatherVCFs as GatherAFonlyVCFs {
        input:
            vcfs = CreateAFonlyVCF.af_only_vcf,
            vcfs_idx = CreateAFonlyVCF.af_only_vcf_idx,
            output_name = af_only_output_name,
            compress_output = compress_output,
            runtime_params = runtime_collection.gather_vcfs
    }

    call tasks.GatherVCFs as GatherCommonAFonlyVCFs {
        input:
            vcfs = CreateAFonlyVCF.common_af_only_vcf,
            vcfs_idx = CreateAFonlyVCF.common_af_only_vcf_idx,
            output_name = af_only_common_output_name,
            compress_output = compress_output,
            drop_duplicate_sites = true,
            runtime_params = runtime_collection.gather_vcfs
    }

    output {
        File af_only_vcf = GatherAFonlyVCFs.merged_vcf
        File af_only_vcf_idx = GatherAFonlyVCFs.merged_vcf_idx
        File common_af_only_vcf = GatherCommonAFonlyVCFs.merged_vcf
        File common_af_only_vcf_idx = GatherCommonAFonlyVCFs.merged_vcf_idx
    }
}