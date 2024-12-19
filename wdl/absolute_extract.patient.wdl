version development

import "absolute_extract.wdl" as ae
import "runtime_collection.wdl" as rtc


workflow AbsoluteExtractPatient {
    input {
        Array[String]? sample_name
        Array[File] rdata
        Array[Int] called_solution
        String analyst_id
        Array[String]? copy_ratio_type

        RuntimeCollection runtime_collection = RuntimeParameters.rtc
    }

    call rtc.DefineRuntimeCollection as RuntimeParameters

    scatter (pair in zip(rdata, called_solution)) {
        call ae.AbsoluteExtract {
            input:
                rdata = pair.left,
                called_solution = pair.right,
                analyst_id = analyst_id,
                runtime_collection = runtime_collection

        }
    }

    output {
        Array[File?] abs_maf = AbsoluteExtract.abs_maf
        Array[File] segtab = AbsoluteExtract.segtab
        Array[File] called_rdata = AbsoluteExtract.called_rdata
        Array[File] table = AbsoluteExtract.table
        Array[File] gene_corrected_cn = AbsoluteExtract.gene_corrected_cn
        Array[File] rescaled_total_cn = AbsoluteExtract.rescaled_total_cn
        Array[String] purity = AbsoluteExtract.purity
        Array[String] ploidy = AbsoluteExtract.ploidy
    }
}