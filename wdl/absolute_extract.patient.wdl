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
        call ae.AbsoluteExtract {
            input:
                rdata = pair.left,
                called_solution = pair.right,
                analyst_id = analyst_id,
#                copy_ratio_type = copy_ratio_type,
                runtime_collection = runtime_collection
        }
    }

    output {
        Array[File?] absolute_maf = AbsoluteExtract.absolute_maf
        Array[File?] absolute_segtab = AbsoluteExtract.absolute_segtab
        Array[File?] absolute_called_rdata = AbsoluteExtract.absolute_called_rdata
        Array[File?] absolute_table = AbsoluteExtract.absolute_table
        Array[File?] absolute_gene_corrected_cn = AbsoluteExtract.absolute_gene_corrected_cn
        Array[File?] absolute_rescaled_total_cn = AbsoluteExtract.absolute_rescaled_total_cn
        Array[String?] absolute_purity = AbsoluteExtract.absolute_purity
        Array[String?] absolute_ploidy = AbsoluteExtract.absolute_ploidy
    }
}