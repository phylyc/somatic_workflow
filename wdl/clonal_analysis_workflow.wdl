version development

import "sample.wdl" as s
import "patient.wdl" as p
import "patient.update_samples.wdl" as p_update_s
import "workflow_arguments.wdl" as wfargs
import "runtime_collection.wdl" as rtc
import "absolute.wdl" as abs


workflow ClonalAnalysisWorkflow {
    input {
        Patient patient
        WorkflowArguments args
        RuntimeCollection runtime_collection
    }

    if (args.run_clonal_decomposition) {
        scatter (sample in patient.tumor_samples) {
            if (defined(sample.called_copy_ratio_segmentation) && defined(sample.af_model_parameters) && defined(sample.annotated_variants)) {
                call abs.Absolute {
                    input:
                        acs_conversion_script = args.acs_conversion_script,
                        sample_name = sample.name,
                        copy_ratio_segmentation = select_first([sample.called_copy_ratio_segmentation]),
                        af_model_parameters = select_first([sample.af_model_parameters]),
                        annotated_variants = sample.annotated_variants,
                        sex = patient.sex,
                        min_hets = args.absolute_min_hets,
                        min_probes = args.absolute_min_probes,
                        maf90_threshold = args.absolute_maf90_threshold,
                        runtime_collection = runtime_collection
                }
            }
        }
        Array[File] acr_plots = select_all(Absolute.acr_plot)
        Array[File] acr_rdata = select_all(Absolute.acr_rdata)
#        Array[File] tcr_plots = select_all(Absolute.tcr_plot)
#        Array[File] tcr_rdata = select_all(Absolute.tcr_rdata)
    }

    # phylogicNDT

    output {
        Array[File]? absolute_acr_plots = acr_plots
        Array[File]? absolute_acr_rdata = acr_rdata
#        Array[File]? absolute_tcr_plots = tcr_plots
#        Array[File]? absolute_tcr_rdata = tcr_rdata
    }
}