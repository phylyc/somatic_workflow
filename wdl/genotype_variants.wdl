version development

import "runtime_collection.wdl" as rtc
import "runtimes.wdl" as rt


workflow GenotypeVariants {
    input {
        String script = "https://github.com/phylyc/genomics_workflows/raw/master/python/genotype.py"

        String individual_id
        Array[String] sample_names
        Array[File] pileups
        Array[File]? contamination_tables
        Array[File]? segmentation_tables
        File? common_germline_alleles  # SNP panel used to get pileups (and estimate contamination)
        File? common_germline_alleles_idx

        Int min_read_depth = 10
        Float min_genotype_likelihood = 0.90
        String format = "GT"
        Boolean select_hets = false
        Boolean save_sample_genotype_likelihoods = false
        Boolean verbose = true

        Boolean compress_output = false

        RuntimeCollection runtime_collection = GetRTC.rtc

        String genotype_docker = "civisanalytics/datascience-python:latest"
        Int preemptible = 1
        Int max_retries = 1
        Int disk_sizeGB = 1
        Int cpu = 1

        Int mem_genotype_variants = 8192
        Int time_startup = 10
        Int time_genotype_variants = 30
    }

    call rtc.DefineRuntimeCollection as GetRTC {
        input:
            genotype_docker = genotype_docker,
            preemptible = preemptible,
            max_retries = max_retries,
            disk_sizeGB = disk_sizeGB,
            cpu = cpu,
            mem_genotype_variants = mem_genotype_variants,
            time_startup = time_startup,
            time_genotype_variants = time_genotype_variants
    }

    call GenotypeVariantsTask {
        input:
            script = script,
            individual_id = individual_id,
            common_germline_alleles = common_germline_alleles,
            common_germline_alleles_idx = common_germline_alleles_idx,
            sample_names = sample_names,
            pileups = pileups,
            contamination_tables = contamination_tables,
            segmentation_tables = segmentation_tables,
            min_read_depth = min_read_depth,
            min_genotype_likelihood = min_genotype_likelihood,
            format = format,
            select_hets = select_hets,
            save_sample_genotype_likelihoods = save_sample_genotype_likelihoods,
            verbose = verbose,
            compress_output = compress_output,
            runtime_params = runtime_collection.genotype_variants
    }

    # todo: resort sample genotype likelihoods

    output {
        File vcf = GenotypeVariantsTask.vcf
        File vcf_idx = GenotypeVariantsTask.vcf_idx
        File ref_counts = GenotypeVariantsTask.ref_counts
        File alt_counts = GenotypeVariantsTask.alt_counts
        File other_alt_counts = GenotypeVariantsTask.other_alt_counts
        File sample_correlation = GenotypeVariantsTask.sample_correlation
        Array[File]? sample_genotype_likelihoods = GenotypeVariantsTask.sample_genotype_likelihoods
    }
}

task GenotypeVariantsTask {
    input {
        String script = "https://github.com/phylyc/genomics_workflows/raw/master/python/genotype.py"

        String individual_id
        Array[String] sample_names
        Array[File] pileups
        Array[File]? contamination_tables
        Array[File]? segmentation_tables
        File? common_germline_alleles  # SNP panel used to get pileups (and estimate contamination)
        File? common_germline_alleles_idx

        Int min_read_depth = 10
        Float min_genotype_likelihood = 0.90
        String model = "betabinom"
        String format = "GT"
        Boolean select_hets = false
        Boolean save_sample_genotype_likelihoods = false
        Boolean compress_output = false
        Boolean verbose = true

        Runtime runtime_params
    }

    String output_dir = "."

    String output_ref_counts = output_dir + "/" + individual_id + ".hets.ref_count.tsv" + (if compress_output then ".gz" else "")
    String output_alt_counts = output_dir + "/" + individual_id + ".hets.alt_count.tsv" + (if compress_output then ".gz" else "")
    String output_other_alt_counts = output_dir + "/" + individual_id + ".hets.other_alt_count.tsv" + (if compress_output then ".gz" else "")
    String output_sample_correlation = output_dir + "/" + individual_id + ".sample_correlation.tsv" + (if compress_output then ".gz" else "")
    String output_vcf = output_dir + "/" + individual_id + ".hets.vcf" + (if compress_output then ".gz" else "")
    String output_vcf_idx = output_vcf + (if compress_output then ".tbi" else ".idx")

    String sample_names_arg = sep("' ", prefix("--sample '", sample_names)) + "'"
    String pileups_arg = sep("' ", prefix("-P '", pileups)) + "'"
    String contamination_arg = if defined(contamination_tables) then sep("' ", prefix("-C '", select_first([contamination_tables]))) + "'" else ""
    String segmentation_arg = if defined(segmentation_tables) then sep("' ", prefix("-S '", select_first([segmentation_tables]))) + "'" else ""

    # Once "suffix" is implemented, we can use this over glob:
#    Array[String] possible_sample_outputs = suffix(".likelihoods.pileup" + (if compress_output then ".gz" else ""), sample_names)
#    Array[File]? output_sample_genotype_likelihoods = if save_sample_genotype_likelihoods then prefix(output_dir + "/", possible_sample_outputs) else None

    command <<<
        set -e
        wget -O genotype.py ~{script}
        python genotype.py \
            --output_dir '~{output_dir}' \
            ~{"--variant '" + common_germline_alleles + "'"} \
            --patient '~{individual_id}' \
            ~{sample_names_arg} \
            ~{pileups_arg} \
            ~{contamination_arg} \
            ~{segmentation_arg} \
            --min_read_depth ~{min_read_depth} \
            --min_genotype_likelihood ~{min_genotype_likelihood} \
            --model ~{model} \
            --format ~{format} \
            ~{if select_hets then "--select_hets" else ""} \
            ~{if save_sample_genotype_likelihoods then "--save_sample_genotype_likelihoods" else ""} \
            ~{if compress_output then "--compress_output" else ""} \
            ~{if verbose then "--verbose" else ""}

        # tabix not in docker
        touch '~{output_vcf_idx}'
    >>>

    output {
        File vcf = output_vcf
        File vcf_idx = output_vcf_idx
        File ref_counts = output_ref_counts
        File alt_counts = output_alt_counts
        File other_alt_counts = output_other_alt_counts
        File sample_correlation = output_sample_correlation
        Array[File]? sample_genotype_likelihoods = glob(output_dir + "/*.likelihoods.pileup" + (if compress_output then ".gz" else ""))
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