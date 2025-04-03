version development

import "absolute_extract.wdl" as ae
import "runtime_collection.wdl" as rtc


workflow AbsoluteExtractPatient {
    input {
        Array[File] rdata
        Array[Int] called_solutions
        String analyst_id
        Array[String]? copy_ratio_types

        Array[File] acs_copy_ratio_segmentation
        Array[Float] acs_copy_ratio_skew
        Array[File]? snv_maf
        Array[File]? indel_maf
        File? gvcf

        RuntimeCollection runtime_collection = RuntimeParameters.rtc
    }

    call rtc.DefineRuntimeCollection as RuntimeParameters

    scatter (tuple in transpose(select_all([rdata, called_solutions, acs_copy_ratio_segmentation, acs_copy_ratio_skew, snv_maf, indel_maf]))) {
        if (defined(snv_maf)) {
            File this_snv_maf = tuple[4]
        }
        if (defined(indel_maf)) {
            File this_indel_maf = tuple[5]
        }
        call ae.AbsoluteExtract {
            input:
                rdata = tuple[0],
                called_solution = tuple[1],
                acs_copy_ratio_segmentation = tuple[2],
                acs_copy_ratio_skew = tuple[3],
                snv_maf = this_snv_maf,
                indel_maf = this_indel_maf,
                analyst_id = analyst_id,
                # copy_ratio_type = copy_ratio_type,
                # sample_name = sample_name,
                runtime_collection = runtime_collection
        }
    }

    output {
        Array[File?] absolute_maf = AbsoluteExtract.absolute_maf
        Array[File?] absolute_segtab = AbsoluteExtract.absolute_segtab
        Array[File] absolute_table = AbsoluteExtract.absolute_table
        Array[Float] absolute_purity = AbsoluteExtract.absolute_purity
        Array[Float] absolute_ploidy = AbsoluteExtract.absolute_ploidy
    }
}