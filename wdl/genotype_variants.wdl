version development

import "runtimes.wdl"


workflow GenotypeVariants {
    input {
        String script = "https://github.com/phylyc/genomics_workflows/raw/master/python/genotype.py"

        String individual_id
        Array[String] sample_names
        Array[File] pileups
        Array[File]? contamination_tables
        Array[File]? segmentation_tables
        File? common_germline_alleles  # SNP panel used to get pileups and estimate contamination
        File? common_germline_alleles_idx

        Int min_read_depth = 10
        Float min_genotype_likelihood = 0.90
        String format = "GT"
        Boolean select_hets = false
        Boolean save_sample_genotype_likelihoods = false
        Boolean verbose = true

        Boolean compress_output = false

        Runtime genotype_variants_runtime = Runtimes.genotype_variants_runtime

        String genotype_docker = "civisanalytics/datascience-python:latest"
        Int preemptible = 1
        Int max_retries = 1
        Int disk_sizeGB = 1
        Int cpu = 1

        Int mem_genotype_variants = 8192
        Int time_startup = 10
        Int time_genotype_variants = 30
    }

    call runtimes.DefineRuntimes as Runtimes {
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

    call GenotypeVariants {
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
            runtime_params = genotype_variants_runtime
    }

    output {
        File vcf = GenotypeVariants.vcf
        File vcf_idx = GenotypeVariants.vcf_idx
        File ref_counts = GenotypeVariants.ref_counts
        File alt_counts = GenotypeVariants.alt_counts
        File other_alt_counts = GenotypeVariants.other_alt_counts
        File sample_correlation = GenotypeVariants.sample_correlation
        Array[File]? sample_genotype_likelihoods = GenotypeVariants.sample_genotype_likelihoods
    }
}

task GenotypeVariants {
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

    Array[String] possible_sample_outputs = suffix(".likelihoods.pileup" + (if compress_output then ".gz" else ""), sample_names)
    Array[File]? output_sample_genotype_likelihoods = if save_sample_genotype_likelihoods then prefix(output_dir + "/", possible_sample_outputs) else None

    String dollar = "$"

    command <<<
        set -e
        wget -O genotype.py ~{script}
        python genotype.py \
            --output_dir '~{output_dir}' \
            ~{"--variant '" + common_germline_alleles + "'"} \
            --patient '~{individual_id}' \
            ~{sep="' " prefix("--sample '", sample_names)}' \
            ~{sep="' " prefix("-P '", pileups)}' \
            ~{true="-S '" false="" defined(segmentation_tables)}~{default="" sep="' -S '" segmentation_tables}~{true="'" false="" defined(segmentation_tables)} \
            ~{true="-C '" false="" defined(contamination_tables)}~{default="" sep="' -C '" contamination_tables}~{true="'" false="" defined(contamination_tables)} \
            --min_read_depth ~{min_read_depth} \
            --min_genotype_likelihood ~{min_genotype_likelihood} \
            --model ~{model} \
            --format ~{format} \
            ~{true="--select_hets " false="" select_hets} \
            ~{true="--save_sample_genotype_likelihoods " false="" save_sample_genotype_likelihoods} \
            ~{true="--compress_output " false="" compress_output} \
            ~{true="--verbose " false="" verbose}

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
        Array[File]? sample_genotype_likelihoods = output_sample_genotype_likelihoods
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