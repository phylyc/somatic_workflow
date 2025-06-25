version development

import "runtime_collection.wdl" as rtc

workflow PhylogicNDT {
    input {
        String patient_id

        Array[File] absolute_maf
        Array[File]? absolute_segtab
        Array[Float] absolute_purity
        Array[Int]? timepoints
        Boolean run_with_BuildTree = true

        RuntimeCollection runtime_collection = RuntimeParameters.rtc
    }

    call rtc.DefineRuntimeCollection as RuntimeParameters

    # Create patient SIF file
    call CreatePatientSIF {
        input:
            patient_id = patient_id,
            absolute_maf = absolute_maf,
            absolute_segtab = absolute_segtab,
            absolute_purity = absolute_purity,
            timepoints = timepoints,
            runtime_params = runtime_collection.create_patient_sif
    }

    # Run PhylogicNDT task
    call PhylogicNDTTask {
        input:
            patient_id = patient_id,
            sif_file = CreatePatientSIF.sif_file,
            run_with_BuildTree = run_with_BuildTree
            runtime_params = runtime_collection.phylogicndt_task
    }

    output {
        # Outputs from CreatePatientSIF
        File sif_file = CreatePatientSIF.sif_file
        
        # Outputs from PhylogicNDTTask
        Array[File] pie_plots = PhylogicNDTTask.pie_plots
        Array[File] one_d_mutation_plots = PhylogicNDTTask.one_d_mutation_plots
        Array[File] one_d_cluster_plots = PhylogicNDTTask.one_d_cluster_plots
        File cnvs = PhylogicNDTTask.cnvs
        File mut_ccfs = PhylogicNDTTask.mut_ccfs
        File unclustered = PhylogicNDTTask.unclustered
        File cluster_ccfs = PhylogicNDTTask.cluster_ccfs
        File? phylogic_report = PhylogicNDTTask.phylogic_report
        File? cell_population_abundances = PhylogicNDTTask.cell_population_abundances
        File? cell_population_mcmc_trace = PhylogicNDTTask.cell_population_mcmc_trace
        File? constrained_ccf = PhylogicNDTTask.constrained_ccf
        File? build_tree_posteriors = PhylogicNDTTask.build_tree_posteriors
    }
}


task CreatePatientSIF {
    input {
        String script = "https://github.com/phylyc/somatic_workflow/raw/master/python/create_patient_sif.py"

        String patient_id
        Array[File] absolute_maf
        Array[File]? absolute_segtab
        Array[Float] absolute_purity
        Array[Int]? timepoints

        Runtime runtime_params
    }

    command <<<
        set -e
        wget -O create_patient_sif.py ~{script}
        python create_patient_sif.py \
            --patient_id ~{patient_id} \
            --absolute_maf ~{sep=' ' absolute_maf} \
            ~{if defined(absolute_segtab) then "--absolute_segtab " + sep=' ' absolute_segtab else ""} \
            --absolute_purity ~{sep=' ' absolute_purity} \
            ~{if defined(timepoints) then "--timepoints " + sep=' ' timepoints else ""} \
            --outfile ~{patient_id}.sif
    >>>

    output {
        File sif_file = "~{patient_id}.sif"
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


task PhylogicNDTTask {
    input {
        String patient_id
        File sif_file
        Boolean run_with_BuildTree = true

        Runtime runtime_params
    }

    command <<<
        python PhylogicNDT.py Cluster \
            -i '~{patient_id}' \
            -sif '~{sif_file}' \
            ~{if run_with_BuildTree then "--run_with_BuildTree" else ""}
    >>>

    output {
        # CLUSTER outputs
        Array[File] pie_plots = glob("~{patient_id}_pie_plots/*.pieplot.svg")
        Array[File] one_d_mutation_plots = glob("~{patient_id}_1d_mutation_plots/*.mutations_ccfs.svg")
        Array[File] one_d_cluster_plots = glob("~{patient_id}_1d_cluster_plots/*.cluster_ccfs.svg")
        File cnvs = '~{patient_id}.cnvs.txt'
        File mut_ccfs = '~{patient_id}.mut_ccfs.txt'
        File unclustered = '~{patient_id}.unclustered.txt'
        File cluster_ccfs = '~{patient_id}.cluster_ccfs.txt'

        # BUILD_TREE outputs
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