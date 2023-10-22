version development

## Collection of Tasks

import "runtimes.wdl"


task GetSampleName {
    input {
        File bam

        Runtime runtime_params
    }

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.jar_override}
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            GetSampleName \
            -I '~{bam}' \
            -O bam_name.txt
    >>>

    output {
        String sample_name = read_string("bam_name.txt")
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

    parameter_meta {
        bam: {localization_optional: true}
    }
}

task PreprocessIntervals {
    input {
        File? interval_list
        Array[File]? interval_lists
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        Int bin_length = 0
        Int padding = 0
        String? preprocess_intervals_extra_args

        Runtime runtime_params
    }

    String preprocessed_intervals = "preprocessed.interval_list"

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.jar_override}
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            PreprocessIntervals \
            -R '~{ref_fasta}' \
            ~{"-L '" + interval_list + "'"} \
            ~{true="-I '" false="" defined(interval_lists)}~{default="" sep="' -I '" interval_lists}~{true="'" false="" defined(interval_lists)} \
            --bin-length ~{bin_length} \
            --padding ~{padding} \
            --interval-merging-rule OVERLAPPING_ONLY \
            -O '~{preprocessed_intervals}' \
            ~{preprocess_intervals_extra_args}
    >>>

    output {
        File preprocessed_interval_list = preprocessed_intervals
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

    parameter_meta {
        interval_list: {localization_optional: true}
        interval_lists: {localization_optional: true}
        ref_fasta: {localization_optional: true}
        ref_fasta_index: {localization_optional: true}
    }
}

task SplitIntervals {
    input {
        File? interval_list
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        Int scatter_count
        String? split_intervals_extra_args

        Runtime runtime_params
    }

    String extra_args = (
        select_first([split_intervals_extra_args, ""])
        # to avoid splitting intervals:
        # + " --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW"
        # Applied after inital scatter, so leads to more scattered intervals.
        # + " --dont-mix-contigs"
    )

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.jar_override}
        mkdir interval-files
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            SplitIntervals \
            -R '~{ref_fasta}' \
            ~{"-L '" + interval_list + "'"} \
            -scatter ~{scatter_count} \
            -O interval-files \
            ~{extra_args}
    >>>

    output {
        Array[File] interval_files = glob("interval-files/*.interval_list")
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

    parameter_meta {
        interval_list: {localization_optional: true}
        ref_fasta: {localization_optional: true}
        ref_fasta_index: {localization_optional: true}
        ref_dict: {localization_optional: true}
    }
}

task SelectVariants {
    input {
        String? interval_list
        File? ref_fasta
        File? ref_fasta_index
        File? ref_dict
        File vcf
        File vcf_idx
        Boolean select_passing = false
        Boolean keep_germline = false
        Boolean compress_output = false
        String? tumor_sample_name
        String? normal_sample_name
        String? select_variants_extra_args

        Runtime runtime_params
    }

    parameter_meta {
        ref_fasta: {localization_optional: true}
        ref_fasta_index: {localization_optional: true}
        ref_dict: {localization_optional: true}
        vcf: {localization_optional: true}
        vcf_idx: {localization_optional: true}
    }

    String uncompressed_input_vcf = basename(vcf, ".gz")
    String base_name = if defined(tumor_sample_name) then sub(select_first([tumor_sample_name, ""]), " ", "+") else basename(uncompressed_input_vcf, ".vcf")
    String output_base_name = base_name + ".selected"
    String uncompressed_output_vcf = output_base_name + ".tmp.vcf"
    String uncompressed_selected_vcf = output_base_name + ".vcf"
    String output_vcf = uncompressed_output_vcf + if compress_output then ".gz" else ""
    String output_vcf_idx = output_vcf + if compress_output then ".tbi" else ".idx"

    String dollar = "$"
    # Remove variants from multi-sample calling that are not present in the selected
    # tumor sample. We require that the total read depth is greater than the allelic
    # depth of the reference allele. This ensures that there is at least one alternate
    # allele present.
    String select_true_variants_arg = (
        if defined(tumor_sample_name)
        then "-select 'vc.getGenotype(\"" + tumor_sample_name + "\").getAD().0 < vc.getGenotype(\"" + tumor_sample_name + "\").getDP()'"
        else ""
    )

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.jar_override}
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            SelectVariants \
            ~{"-R '" + ref_fasta + "'"} \
            ~{"-L '" + interval_list + "'"} \
            -V '~{vcf}' \
            --output '~{uncompressed_output_vcf}' \
            --exclude-filtered false \
            ~{"--sample-name '" + tumor_sample_name + "'"} \
            ~{"--sample-name '" + normal_sample_name + "'"} \
            ~{select_true_variants_arg} \
            ~{select_variants_extra_args}

        # =======================================
        # We do the selection step using grep to also select germline variants.

        grep "^#" '~{uncompressed_output_vcf}' > '~{uncompressed_selected_vcf}'
        if ~{select_passing} ; then
            echo ">> Selecting PASSing variants ... "
            grep "PASS" '~{uncompressed_output_vcf}' >> '~{uncompressed_selected_vcf}'
            num_vars=~{dollar}(grep -v "^#" '~{uncompressed_output_vcf}' | wc -l)
            num_selected_vars=~{dollar}(grep -v "^#" '~{uncompressed_selected_vcf}' | wc -l)
            echo ">> Selected ~{dollar}num_selected_vars PASSing out of ~{dollar}num_vars variants."
        fi
        if ~{keep_germline} ; then
            echo ">> Selecting germline variants ... "
            grep -P "\tgermline\t" '~{uncompressed_output_vcf}' >> '~{uncompressed_selected_vcf}'
            num_vars=~{dollar}(grep -v "^#" '~{uncompressed_output_vcf}' | wc -l)
            num_selected_vars=~{dollar}(grep -P "\tgermline\t" '~{uncompressed_selected_vcf}' | wc -l)
            echo ">> Selected ~{dollar}num_selected_vars germline out of ~{dollar}num_vars variants."
        fi

        if ~{!select_passing} || ~{!keep_germline} ; then
            mv '~{uncompressed_output_vcf}' '~{uncompressed_selected_vcf}'
        fi

        # =======================================
        # Hack to correct a SelectVariants output bug. When selecting for samples, this
        # task only retains the first sample annotation in the header. Those annotations
        # are important for Funcotator to fill the t_alt_count and t_ref_count coverage
        # columns. This hack assumes that only one tumor sample and/or only one normal
        # sample have been selected.

        if ~{defined(tumor_sample_name)} ; then
            echo ">> Fixing tumor sample name in vcf header ... "
            input_header=~{dollar}(grep "##tumor_sample=" '~{uncompressed_selected_vcf}')
            corrected_header="##tumor_sample=~{tumor_sample_name}"
            sed -i "s/~{dollar}input_header/~{dollar}corrected_header/g" '~{uncompressed_selected_vcf}'
        fi
        if ~{defined(normal_sample_name)} ; then
            echo ">> Fixing normal sample name in vcf header ... "
            input_header=~{dollar}(grep "##normal_sample=" '~{uncompressed_selected_vcf}')
            corrected_header="##normal_sample=~{normal_sample_name}"
            sed -i "s/~{dollar}input_header/~{dollar}corrected_header/g" '~{uncompressed_selected_vcf}'
        fi

        grep -v "^#" '~{uncompressed_selected_vcf}' | wc -l > num_selected_vars.txt

        if ~{compress_output} ; then
            echo ">> Compressing selected vcf."
            bgzip -c '~{uncompressed_selected_vcf}' > '~{output_vcf}'
            gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
                IndexFeatureFile \
                --input '~{output_vcf}' \
                --output '~{output_vcf_idx}'
            rm '~{uncompressed_selected_vcf}'
            rm '~{uncompressed_output_vcf}.idx'
        fi
    >>>

    output {
        File selected_vcf = output_vcf
        File selected_vcf_idx = output_vcf_idx
        Int num_selected_variants = read_int("num_selected_vars.txt")
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

task MergeVCFs {
    # Consider replacing MergeVcfs with GatherVcfsCloud once the latter is out of beta.

	input {
        Array[File] vcfs
        Array[File] vcfs_idx
        String output_name
        Boolean compress_output = false

        Runtime runtime_params
    }

    # Optional localization leads to cromwell error.
    # parameter_meta {
    #     vcfs: {localization_optional: true}
    #     vcfs_idx: {localization_optional: true}
    # }

    String output_vcf = output_name + ".vcf" + if compress_output then ".gz" else ""
    String output_vcf_idx = output_vcf + if compress_output then ".tbi" else ".idx"

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.jar_override}
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            MergeVcfs \
            ~{sep="' " prefix("-I '", vcfs)}' \
            -O '~{output_vcf}'
    >>>

    output {
    	File merged_vcf = output_vcf
        File merged_vcf_idx = output_vcf_idx
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

task MergeBams {
    input {
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        Array[File]+ bams
        Array[File]+ bais
        String merged_bam_name

        Runtime runtime_params
    }

    parameter_meta {
        ref_fasta: {localization_optional: true}
        ref_fasta_index: {localization_optional: true}
        ref_dict: {localization_optional: true}
        # bams: {localization_optional: true}  # samtools requires localization
    }

    Int disk_spaceGB = 4 * ceil(size(bams, "GB")) + runtime_params.disk
    String output_bam_name = merged_bam_name + ".bam"
    String output_bai_name = merged_bam_name + ".bai"

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.jar_override}
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            GatherBamFiles \
            ~{sep="' " prefix("-I '", bams)}' \
            -O unsorted.out.bam \
            -R '~{ref_fasta}'

        # We must sort because adjacent scatters may have overlapping (padded) assembly
        # regions, hence overlapping bamouts

        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            SortSam \
            -I unsorted.out.bam \
            -O '~{output_bam_name}' \
            --SORT_ORDER coordinate \
            --VALIDATION_STRINGENCY LENIENT

        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            BuildBamIndex \
            -I '~{output_bam_name}' \
            --VALIDATION_STRINGENCY LENIENT
    >>>

    output {
        File merged_bam = output_bam_name
        File merged_bai = output_bai_name
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