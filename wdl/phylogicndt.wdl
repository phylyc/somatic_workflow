version development

import "runtime_collection.wdl" as rtc

workflow PhylogicNDT {
    input {
        String patient_id

        Array[String]? sample_names
        Array[File] absolute_mafs
        Array[File]? absolute_segtabs
        Array[Float] absolute_purities
        Array[Int]? timepoints
        Array[Float]? tumor_mutation_burdens
        Boolean run_with_BuildTree = true
        Boolean use_indels = false
        Boolean impute_missing_snvs = false
        Int min_coverage = 8
        File? driver_genes_file

        RuntimeCollection runtime_collection = RuntimeParameters.rtc
    }

    call rtc.DefineRuntimeCollection as RuntimeParameters

    # Create patient SIF and runs PhylogicNDT
    call PhylogicNDTTask {
        input:
            patient_id = patient_id,
            sample_names = sample_names,
            absolute_mafs = absolute_mafs,
            absolute_segtabs = absolute_segtabs,
            absolute_purities = absolute_purities,
            timepoints = timepoints,
            tumor_mutation_burdens = tumor_mutation_burdens,
            run_with_BuildTree = run_with_BuildTree,
            use_indels = use_indels,
            impute_missing_snvs = impute_missing_snvs,
            min_coverage = min_coverage,
            driver_genes_file = driver_genes_file,
            runtime_params = runtime_collection.phylogicndt_task
    }

    output {
        File sif_file = PhylogicNDTTask.sif_file
        
        # Outputs from PhylogicNDT Cluster
        File cnvs = PhylogicNDTTask.cnvs
        File mut_ccfs = PhylogicNDTTask.mut_ccfs
        File unclustered = PhylogicNDTTask.unclustered
        File cluster_ccfs = PhylogicNDTTask.cluster_ccfs
        # Array[File] pie_plots = PhylogicNDTTask.pie_plots
        Array[File] mutation_plots = PhylogicNDTTask.one_d_mutation_plots
        Array[File] cluster_plots = PhylogicNDTTask.one_d_cluster_plots

        # Outputs from PhylogicNDT BuildTree
        # Note: these are optional and may not be present if run_with_BuildTree is false
        File? report = PhylogicNDTTask.report
        File? cell_population_abundances = PhylogicNDTTask.cell_population_abundances
        File? cell_population_mcmc_trace = PhylogicNDTTask.cell_population_mcmc_trace
        File? constrained_ccf = PhylogicNDTTask.constrained_ccf
        File? build_tree_posteriors = PhylogicNDTTask.build_tree_posteriors

        File? growth_rates = PhylogicNDTTask.growth_rates
        File? growth_rate_plot = PhylogicNDTTask.growth_rate_plot

        File? timing_comparison = PhylogicNDTTask.timing_comparison
        File? timing_table = PhylogicNDTTask.timing_table
        File? timing_graph = PhylogicNDTTask.timing_graph
        File? timing_report = PhylogicNDTTask.timing_report
        File? timing_wgd_supporting_events = PhylogicNDTTask.timing_wgd_supporting_events
    }
}

task PhylogicNDTTask {
    input {
        String patient_id

        Array[String]? sample_names
        Array[File] absolute_mafs
        Array[File]? absolute_segtabs
        Array[Float] absolute_purities
        Array[Int]? timepoints
        Array[Float]? tumor_mutation_burdens
        Boolean run_with_BuildTree = true
        Boolean use_indels = false
        Boolean impute_missing_snvs = false
        Int min_coverage = 8
        Float Pk_k_r = 3.0
        Float Pk_k_mu = 3.0
        File? driver_genes_file
        Runtime runtime_params
    }

    String sif = patient_id + ".sif"
    String timing_sif = patient_id + ".sif.timing.txt"

    command <<<
        set -e
        python /build/PhylogicNDT/create_patient_sif.py \
            --patient_id '~{patient_id}' \
            ~{if defined(sample_names) then "--sample_names '" else ""}~{default="" sep="' '" sample_names}~{if defined(sample_names) then "'" else ""} \
            --absolute_mafs '~{sep="' '" absolute_mafs}' \
            ~{if defined(absolute_segtabs) then "--absolute_segtabs '" else ""}~{default="" sep="' '" absolute_segtabs}~{if defined(absolute_segtabs) then "'" else ""} \
            --absolute_purities ~{sep=" " absolute_purities} \
            ~{if defined(timepoints) then "--timepoints " else ""}~{default="" sep=" " timepoints} \
            ~{if defined(tumor_mutation_burdens) then "--tumor_mutation_burdens " else ""}~{default="" sep=" " tumor_mutation_burdens} \
            --outfile '~{sif}'

        python /build/PhylogicNDT/PhylogicNDT.py Cluster \
            -i '~{patient_id}' \
            -sif '~{sif}' \
            ~{"--driver_genes_file " + driver_genes_file} \
            ~{if use_indels then "--use_indels" else ""} \
            ~{if impute_missing_snvs then "--impute" else ""} \
            ~{if run_with_BuildTree then "--run_with_BuildTree" else ""} \
            --Pk_k_r ~{Pk_k_r} \
            --Pk_k_mu ~{Pk_k_mu} \
            --min_coverage ~{min_coverage}

        # Cell populations are already being inferred.

        if [ -f "~{patient_id}_cell_population_mcmc_trace.tsv" ]; then
            # May not yield any result for single samples
            python /build/PhylogicNDT/PhylogicNDT.py GrowthKinetics \
                -i "~{patient_id}" \
                -ab "~{patient_id}_cell_population_mcmc_trace.tsv" \
                -sif "~{sif}"
        fi

        if [ "~{defined(absolute_segtabs)}" = "true" ]; then
            python /build/PhylogicNDT/PhylogicNDT.py Timing \
                -i "~{patient_id}" \
                -sif "~{timing_sif}"
        fi
    >>>

    output {
        File sif_file = sif

        # PhylogicNDT Cluster outputs
        File cnvs = '~{patient_id}.cnvs.txt'
        File mut_ccfs = '~{patient_id}.mut_ccfs.txt'
        File unclustered = '~{patient_id}.unclustered.txt'
        File cluster_ccfs = '~{patient_id}.cluster_ccfs.txt'
        # Array[File] pie_plots = glob("~{patient_id}_pie_plots/*.pieplot.svg") # Pie plots already in HTML report
        Array[File] one_d_mutation_plots = glob("~{patient_id}_1d_mutation_plots/*.mutations_ccfs.svg")
        Array[File] one_d_cluster_plots = glob("~{patient_id}_1d_cluster_plots/*.cluster_ccfs.svg")

        # PhylogicNDT BuildTree outputs
        File? report = "~{patient_id}.phylogic_report.html"
        File? cell_population_abundances = "~{patient_id}_cell_population_abundances.tsv"
        File? cell_population_mcmc_trace = "~{patient_id}_cell_population_mcmc_trace.tsv"
        File? constrained_ccf = "~{patient_id}_constrained_ccf.tsv"
        File? build_tree_posteriors = "~{patient_id}_build_tree_posteriors.tsv"

        File? growth_rates = "~{patient_id}_growth_rate.tsv"
        File? growth_rate_plot = "~{patient_id}.growth_rate.pdf"

        File? timing_comparison = "~{patient_id}.comp.tsv"
        File? timing_table = "~{patient_id}.timing.tsv"
        File? timing_graph = "~{patient_id}.comp_graph.pdf"
        File? timing_report = "~{patient_id}.phylogic_timing_report.html"
        File? timing_wgd_supporting_events = "~{patient_id}.WGD_supporting_events.timing.tsv"
    }

    runtime {
        docker: runtime_params.docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: runtime_params.machine_mem + " MB"
        runtime_minutes: runtime_params.runtime_minutes
        disks: "local-disk " + runtime_params.disk + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }
}