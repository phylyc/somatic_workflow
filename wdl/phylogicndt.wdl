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
        Boolean run_with_BuildTree = true
        String phylogicndt_create_sif_script = "https://github.com/phylyc/somatic_workflow/raw/master/python/create_patient_sif.py"

        RuntimeCollection runtime_collection = RuntimeParameters.rtc
    }

    call rtc.DefineRuntimeCollection as RuntimeParameters

    # Create patient SIF and runs PhylogicNDT
    call PhylogicNDTTask {
        input:
            script = phylogicndt_create_sif_script,
            patient_id = patient_id,
            sample_names = sample_names,
            absolute_mafs = absolute_mafs,
            absolute_segtabs = absolute_segtabs,
            absolute_purities = absolute_purities,
            timepoints = timepoints,
            run_with_BuildTree = run_with_BuildTree,
            runtime_params = runtime_collection.phylogicndt_task
    }

    output {
        File phylogic_sif_file = PhylogicNDTTask.sif_file
        
        # Outputs from PhylogicNDT Cluster
        Array[File] phylogic_pie_plots = PhylogicNDTTask.pie_plots
        Array[File] phylogic_mutation_plots = PhylogicNDTTask.one_d_mutation_plots
        Array[File] phylogic_cluster_plots = PhylogicNDTTask.one_d_cluster_plots
        File phylogic_cnvs = PhylogicNDTTask.cnvs
        File phylogic_mut_ccfs = PhylogicNDTTask.mut_ccfs
        File phylogic_unclustered = PhylogicNDTTask.unclustered
        File phylogic_cluster_ccfs = PhylogicNDTTask.cluster_ccfs

        # Outputs from PhylogicNDT BuildTree
        # Note: these are optional and may not be present if run_with_BuildTree is false
        File? phylogic_report = PhylogicNDTTask.phylogic_report
        File? phylogic_cell_population_abundances = PhylogicNDTTask.cell_population_abundances
        File? phylogic_cell_population_mcmc_trace = PhylogicNDTTask.cell_population_mcmc_trace
        File? phylogic_constrained_ccf = PhylogicNDTTask.constrained_ccf
        File? phylogic_build_tree_posteriors = PhylogicNDTTask.build_tree_posteriors
    }
}

task PhylogicNDTTask {
    input {
        String script = "https://github.com/phylyc/somatic_workflow/raw/master/python/create_patient_sif.py"
        String patient_id

        Array[String]? sample_names
        Array[File] absolute_mafs
        Array[File]? absolute_segtabs
        Array[Float] absolute_purities
        Array[Int]? timepoints
        Boolean run_with_BuildTree = true
        Runtime runtime_params
    }

    String sif = patient_id + ".sif"

    command <<<
        set -e
        wget -O create_patient_sif.py ~{script}

        # phylogicndt_docker uses python for Python 3.8.5 
        python create_patient_sif.py \
            --patient_id '~{patient_id}' \
            ~{if defined(sample_names) then "--sample_names '" else ""}~{default="" sep="' '" sample_names}~{if defined(sample_names) then "'" else ""} \
            --absolute_mafs '~{sep="' '" absolute_mafs}' \
            ~{if defined(absolute_segtabs) then "--absolute_segtabs '" else ""}~{default="" sep="' '" absolute_segtabs}~{if defined(absolute_segtabs) then "'" else ""} \
            --absolute_purities ~{sep=" " absolute_purities} \
            ~{if defined(timepoints) then "--timepoints " else ""}~{default="" sep=" " timepoints} \
            --outfile '~{sif}'

        # phylogicndt_docker uses python2 for Python 2.7.18 
        python2 PhylogicNDT.py Cluster \
            -i '~{patient_id}' \
            -sif '~{sif}' \
            ~{if run_with_BuildTree then "--run_with_BuildTree" else ""}
    >>>

    output {
        File sif_file = sif

        # PhylogicNDT Cluster outputs
        Array[File] pie_plots = glob("~{patient_id}_pie_plots/*.pieplot.svg")
        Array[File] one_d_mutation_plots = glob("~{patient_id}_1d_mutation_plots/*.mutations_ccfs.svg")
        Array[File] one_d_cluster_plots = glob("~{patient_id}_1d_cluster_plots/*.cluster_ccfs.svg")
        File cnvs = '~{patient_id}.cnvs.txt'
        File mut_ccfs = '~{patient_id}.mut_ccfs.txt'
        File unclustered = '~{patient_id}.unclustered.txt'
        File cluster_ccfs = '~{patient_id}.cluster_ccfs.txt'

        # PhylogicNDT BuildTree outputs
        File? phylogic_report = select_first([glob("~{patient_id}.phylogic_report.html"), []])
        File? cell_population_abundances = select_first([glob("~{patient_id}_cell_population_abundances.tsv"), []])
        File? cell_population_mcmc_trace = select_first([glob("~{patient_id}_cell_population_mcmc_trace.tsv"), []])
        File? constrained_ccf = select_first([glob("~{patient_id}_constrained_ccf.tsv"), []])
        File? build_tree_posteriors = select_first([glob("~{patient_id}_build_tree_posteriors.tsv"), []])
    }

    # Need docker with PhylogicNDT !!
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