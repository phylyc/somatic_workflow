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
        File? snv_maf
        File? indel_maf
        File? gvcf

        RuntimeCollection runtime_collection = RuntimeParameters.rtc
    }

    call rtc.DefineRuntimeCollection as RuntimeParameters

    scatter (tuple in transpose([rdata, called_solutions, acs_copy_ratio_segmentation, acs_copy_ratio_skew])) {
        call ae.AbsoluteExtract {
            input:
                rdata = tuple[0],
                called_solution = tuple[1],
                acs_copy_ratio_segmentation = tuple[2],
                acs_copy_ratio_skew = tuple[3],
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