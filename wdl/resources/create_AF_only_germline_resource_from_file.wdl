version development

import "create_AF_only_germline_resource.wdl" as singleAFonly
import "runtime_collection.wdl" as rtc
import "tasks.wdl"


workflow CreateAFonlyVCFfromFile {
    input {
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        File vcfs_file
        File vcfs_idx_file

        Boolean compress_output = true
        Boolean create_biallelic = false
        Boolean create_multiallelic = false

        RuntimeCollection runtime_collection = RuntimeParameters.rtc
    }

    call rtc.DefineRuntimeCollection as RuntimeParameters

    Array[File]+ vcfs = read_lines(vcfs_file)
    Array[File]+ vcfs_idx = read_lines(vcfs_idx_file)

    scatter (vcf_pair in zip(vcfs, vcfs_idx)) {
        call singleAFonly.CreateAFonlyVCF {
            input:
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                ref_dict = ref_dict,
                vcf = vcf_pair.left,
                vcf_idx = vcf_pair.right,
                compress_output = compress_output,
                create_biallelic = create_biallelic,
                create_multiallelic = create_multiallelic,
                runtime_collection = runtime_collection
        }
    }

    call tasks.MergeVCFs as MergeAFonlyVCFs {
        input:
            vcfs = CreateAFonlyVCF.af_only_vcf,
            vcfs_idx = CreateAFonlyVCF.af_only_vcf_idx,
            output_name = "af_only.filtered",
            compress_output = compress_output,
            runtime_params = runtime_collection.merge_vcfs
    }

    if (create_biallelic) {
        call tasks.MergeVCFs as MergeBiallelicAFonlyVCFs {
            input:
                vcfs = select_all(CreateAFonlyVCF.biallelic_af_only_vcf),
                vcfs_idx = select_all(CreateAFonlyVCF.biallelic_af_only_vcf_idx),
                output_name = "af_only.filtered.biallelic",
                compress_output = compress_output,
                runtime_params = runtime_collection.merge_vcfs
        }
    }

    if (create_multiallelic) {
        call tasks.MergeVCFs as MergeMultiallelicAFonlyVCFs {
            input:
                vcfs = select_all(CreateAFonlyVCF.multiallelic_af_only_vcf),
                vcfs_idx = select_all(CreateAFonlyVCF.multiallelic_af_only_vcf_idx),
                output_name = "af_only.filtered.multiallelic",
                compress_output = compress_output,
                runtime_params = runtime_collection.merge_vcfs
        }
    }

    output {
        File af_only_vcf = MergeAFonlyVCFs.merged_vcf
        File af_only_vcf_idx = MergeAFonlyVCFs.merged_vcf_idx
  		File? biallelic_af_only_vcf = MergeBiallelicAFonlyVCFs.merged_vcf
  		File? biallelic_af_only_vcf_idx = MergeBiallelicAFonlyVCFs.merged_vcf
        File? multiallelic_af_only_vcf = MergeMultiallelicAFonlyVCFs.merged_vcf
        File? multiallelic_af_only_vcf_idx = MergeMultiallelicAFonlyVCFs.merged_vcf
    }
}
