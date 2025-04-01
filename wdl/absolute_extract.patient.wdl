version development

import "absolute_extract.wdl" as ae
import "runtime_collection.wdl" as rtc


workflow AbsoluteExtractPatient {
    input {
        Array[File] rdata
        Array[Int] called_solutions
        String analyst_id
        Array[String]? copy_ratio_types

        RuntimeCollection runtime_collection = RuntimeParameters.rtc
    }

    call rtc.DefineRuntimeCollection as RuntimeParameters

    scatter (pair in zip(rdata, called_solutions)) {
        if (pair.right > 0) {
            call ae.AbsoluteExtract {
                input:
                    rdata = pair.left,
                    called_solution = pair.right,
                    analyst_id = analyst_id,
    #                copy_ratio_type = copy_ratio_type,
                    runtime_collection = runtime_collection
            }
        }
    }

    output {
        Array[File] absolute_maf = select_all(AbsoluteExtract.absolute_maf)
        Array[File] absolute_segtab = select_all(AbsoluteExtract.absolute_segtab)
        Array[File] absolute_table = select_all(AbsoluteExtract.absolute_table)
        Array[String] absolute_purity = select_all(AbsoluteExtract.absolute_purity)
        Array[String] absolute_ploidy = select_all(AbsoluteExtract.absolute_ploidy)
    }
}