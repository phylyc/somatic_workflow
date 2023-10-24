version development

## Call variants with Mutect2
##
## Known issues:
## force_call_alleles:
## genotype_germline_sites: Use with care! https://github.com/broadinstitute/gatk/issues/7391
##      @David Benjamin, 2021:
##      "The issue with the clustered events and haplotype filters when running in
##      -genotype-germline-sites is a real problem. These filters should only be
##      triggered by a cluster of technical artifacts or somatic variants, not by
##      germline variants. However, since the default mode of Mutect2 ignores most
##      germline sites, we overlooked this possibility. We need to fix those filters
##      so that they work as intended."
##      force-calling alleles of a large SNP panel influences the germline filtering
##      model in the same way.
## use_linked_de_bruijn_graph: This has trouble calling variants in complex regions.
##      Strongly recommended to use with recover_all_dangling_branches. This increases
##      compute cost though, and may still not guarantee that all variants are
##      being called. Ideally, run with and without and use the joint callset.

import "runtimes.wdl"
import "tasks.wdl"


workflow CallVariants {
    input {
        Array[File] scattered_interval_list
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        String individual_id
        Array[File]+ tumor_bams
        Array[File]+ tumor_bais
        Array[String]? tumor_bam_names
        Array[File]? normal_bams
        Array[File]? normal_bais
        Array[String]? normal_bam_names

        # resources
        File? force_call_alleles
        File? force_call_alleles_idx
        File? panel_of_normals
        File? panel_of_normals_idx
        File? germline_resource
        File? germline_resource_tbi

        Boolean run_orientation_bias_mixture_model = true

        Boolean compress_output = true
        Boolean make_bamout = false

        Boolean native_pair_hmm_use_double_precision = true
        Boolean use_linked_de_bruijn_graph = true
        Boolean recover_all_dangling_branches = true
        Boolean pileup_detection = true
        Boolean genotype_germline_sites = false  # use with care! (see above)
        Int downsampling_stride = 1
        Int max_reads_per_alignment_start = 50
        String? mutect2_extra_args

        Runtime variant_call_runtime = Runtimes.variant_call_runtime
        Runtime merge_vcfs_runtime = Runtimes.merge_vcfs_runtime
        Runtime merge_mutect_stats_runtime = Runtimes.merge_mutect_stats_runtime
        Runtime merge_bams_runtime = Runtimes.merge_bams_runtime
        Runtime learn_read_orientation_model_runtime = Runtimes.learn_read_orientation_model_runtime

        Int mem_additional_per_sample = 256  # this actually can depend on bam size (WES vs WGS)
        Int mem_variant_call_base = 4096
        Int mem_merge_vcfs = 512
        Int mem_merge_mutect_stats = 512 # 64
        Int mem_merge_bams = 8192  # wants at least 6G
        Int mem_learn_read_orientation_model_base = 6144

        Int time_startup = 10
        Int time_variant_call_total = 10000  # 6 d / scatter_count
        Int time_merge_vcfs = 10
        Int time_merge_mutect_stats = 1
        Int time_merge_bams = 60
        Int time_learn_read_orientation_model = 180  # 3 h

        String gatk_docker = "broadinstitute/gatk"
        File? gatk_override
        Int preemptible = 1
        Int max_retries = 1
        Int emergency_extra_diskGB = 0

        Int cpu_variant_call = 1  # good for PairHMM: 2
    }

    Int scatter_count = length(scattered_interval_list)
    Int gatk_override_size = if defined(gatk_override) then ceil(size(gatk_override, "GB")) else 0
    Int disk_padGB = 1 + gatk_override_size + emergency_extra_diskGB

    Int num_bams = length(tumor_bams) + length(select_first([normal_bams, []]))

    Int tumor_size = ceil(size(tumor_bams, "GB") + size(tumor_bais, "GB"))
    Int m2_output_size = if make_bamout then ceil(tumor_size / scatter_count) else 0
    Int m2_per_scatter_size = 1 + m2_output_size + disk_padGB

    call runtimes.DefineRuntimes as Runtimes {
        input:
            num_bams = num_bams,
            scatter_count = scatter_count,
            gatk_docker = gatk_docker,
            gatk_override = gatk_override,
            preemptible = preemptible,
            max_retries = max_retries,
            cpu_variant_call = cpu_variant_call,
            disk_variant_call = m2_per_scatter_size,
            mem_additional_per_sample = mem_additional_per_sample,
            mem_learn_read_orientation_model_base = mem_learn_read_orientation_model_base,
            mem_merge_vcfs = mem_merge_vcfs,
            mem_merge_mutect_stats = mem_merge_mutect_stats,
            mem_merge_bams = mem_merge_bams,
            mem_variant_call_base = mem_variant_call_base,
            time_startup = time_startup,
            time_learn_read_orientation_model = time_learn_read_orientation_model,
            time_merge_vcfs = time_merge_vcfs,
            time_merge_mutect_stats = time_merge_mutect_stats,
            time_merge_bams = time_merge_bams,
            time_variant_call_total = time_variant_call_total,
    }

    scatter (interval_list in scattered_interval_list) {
    	call Mutect2 {
            input:
                interval_list = interval_list,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                ref_dict = ref_dict,
                individual_id = individual_id,
                tumor_bams = tumor_bams,
                tumor_bais = tumor_bais,
                normal_bams = normal_bams,
                normal_bais = normal_bais,
                normal_sample_names = normal_bam_names,
                force_call_alleles = force_call_alleles,
                force_call_alleles_idx = force_call_alleles_idx,
                panel_of_normals = panel_of_normals,
                panel_of_normals_idx = panel_of_normals_idx,
                germline_resource = germline_resource,
                germline_resource_tbi = germline_resource_tbi,
                make_bamout = make_bamout,
                get_orientation_bias_priors = run_orientation_bias_mixture_model,
                compress_output = compress_output,
                genotype_germline_sites = genotype_germline_sites,
                native_pair_hmm_use_double_precision = native_pair_hmm_use_double_precision,
                use_linked_de_bruijn_graph = use_linked_de_bruijn_graph,
                recover_all_dangling_branches = recover_all_dangling_branches,
                pileup_detection = pileup_detection,
                downsampling_stride = downsampling_stride,
                max_reads_per_alignment_start = max_reads_per_alignment_start,
                m2_extra_args = mutect2_extra_args,
                runtime_params = variant_call_runtime,
		}
	}

    call tasks.MergeVCFs {
    	input:
            vcfs = Mutect2.vcf,
            vcfs_idx = Mutect2.vcf_idx,
            output_name = individual_id + ".unfiltered.merged",
            compress_output = compress_output,
            runtime_params = merge_vcfs_runtime
    }

    if (make_bamout) {
        call tasks.MergeBams {
            input:
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                ref_dict = ref_dict,
                bams = select_all(Mutect2.bam),
                bais = select_all(Mutect2.bai),
                merged_bam_name = individual_id + ".Mutect2.out",
                runtime_params = merge_bams_runtime
        }
    }

    call MergeMutectStats {
        input:
            stats = Mutect2.vcf_stats,
            individual_id = individual_id,
            runtime_params = merge_mutect_stats_runtime
    }

    if (run_orientation_bias_mixture_model) {
        call LearnReadOrientationModel {
            input:
                individual_id = individual_id,
                f1r2_counts = select_all(Mutect2.m2_artifact_priors),
                runtime_params = learn_read_orientation_model_runtime
        }
    }

    output {
        File vcf = MergeVCFs.merged_vcf
        File vcf_idx = MergeVCFs.merged_vcf_idx
        File mutect_stats = MergeMutectStats.merged_stats
        File? bam = MergeBams.merged_bam
        File? bai = MergeBams.merged_bai
        File? orientation_bias = LearnReadOrientationModel.orientation_bias
    }
}


task Mutect2 {
    input {
        File? interval_list
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        String individual_id
        Array[File] tumor_bams
        Array[File] tumor_bais
        Array[File]? normal_bams
        Array[File]? normal_bais
        Array[String]? normal_sample_names

        File? force_call_alleles
        File? force_call_alleles_idx
        File? panel_of_normals
        File? panel_of_normals_idx
        File? germline_resource
        File? germline_resource_tbi

        Boolean genotype_germline_sites = false
        Boolean native_pair_hmm_use_double_precision = true
        Boolean use_linked_de_bruijn_graph = true
        Boolean recover_all_dangling_branches = true
        Boolean pileup_detection = true

        # The linked de-Bruijn graph implementation has trouble calling variants
        # in complex regions, even when recovering all dangling branches.
        # Reducing the downsampling by increasing the following parameters might
        # solve the issue. It increases compute cost though.
        Int downsampling_stride = 1
        Int max_reads_per_alignment_start = 50

        String? m2_extra_args

        Boolean get_orientation_bias_priors = true
        Boolean compress_output = false
        Boolean make_bamout = false

        Runtime runtime_params
    }

    parameter_meta {
        interval_list: {localization_optional: true}
        ref_fasta: {localization_optional: true}
        ref_fasta_index: {localization_optional: true}
        ref_dict: {localization_optional: true}
        tumor_bams: {localization_optional: true}
        tumor_bais: {localization_optional: true}
        normal_bams: {localization_optional: true}
        normal_bais: {localization_optional: true}
        panel_of_normals: {localization_optional: true}
        panel_of_normals_idx: {localization_optional: true}
        germline_resource: {localization_optional: true}
        germline_resource_tbi: {localization_optional: true}
    }

    Boolean normal_is_present = defined(normal_bams) && (length(select_first([normal_bams])) > 0)

    String output_vcf = individual_id + if compress_output then ".vcf.gz" else ".vcf"
    String output_vcf_idx = output_vcf + if compress_output then ".tbi" else ".idx"
    String output_stats = output_vcf + ".stats"

    String output_bam = individual_id + ".bamout.bam"
    String output_bai = individual_id + ".bamout.bai"
    String make_bamout_arg = if make_bamout then "--bam-output " + output_bam else ""
    String output_artifact_priors = individual_id + ".f1r2_counts.tar.gz"
    String run_ob_filter_arg = if get_orientation_bias_priors then "--f1r2-tar-gz " + output_artifact_priors else ""

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.jar_override}
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            Mutect2 \
            --reference '~{ref_fasta}' \
            ~{sep="' " prefix("-I '", tumor_bams)}' \
            ~{true="-I '" false="" normal_is_present}~{default="" sep="' -I '" normal_bams}~{true="'" false="" normal_is_present} \
            ~{true="-normal '" false="" normal_is_present}~{default="" sep="' -normal '" normal_sample_names}~{true="'" false="" normal_is_present} \
            --output '~{output_vcf}' \
            ~{"--intervals '" + interval_list + "'"} \
            ~{"--alleles '" + force_call_alleles + "'"} \
            ~{"-pon '" + panel_of_normals + "'"} \
            ~{make_bamout_arg} \
            ~{run_ob_filter_arg} \
            ~{"--germline-resource '" + germline_resource + "'"} \
            ~{true="--genotype-germline-sites true" false="" genotype_germline_sites} \
            ~{true="--linked-de-bruijn-graph true" false="" use_linked_de_bruijn_graph} \
            ~{true="--recover-all-dangling-branches true" false="" recover_all_dangling_branches} \
            ~{true="--pileup-detection true" false="" pileup_detection} \
            --smith-waterman FASTEST_AVAILABLE \
            --pair-hmm-implementation FASTEST_AVAILABLE \
            ~{true="--native-pair-hmm-use-double-precision true" false="" native_pair_hmm_use_double_precision} \
            ~{"--downsampling-stride " + downsampling_stride} \
            ~{"--max-reads-per-alignment-start " + max_reads_per_alignment_start} \
            --seconds-between-progress-updates 300 \
            ~{m2_extra_args}
    >>>

    output {
        File vcf = output_vcf
        File vcf_idx = output_vcf_idx
        File vcf_stats = output_stats
        File? bam = output_bam
        File? bai = output_bai
        File? m2_artifact_priors = output_artifact_priors
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

task MergeMutectStats {
    input {
        Array[File]+ stats
        String individual_id

        Runtime runtime_params
    }

    # Optional localization leads to cromwell error.
    # parameter_meta {
    #     stats: {localization_optional: true}
    # }

    String output_name = individual_id + ".merged.stats"

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.jar_override}
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            MergeMutectStats \
            ~{sep="' " prefix("-stats '", stats)}' \
            --output '~{output_name}'
    >>>

    output {
        File merged_stats = output_name
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

task LearnReadOrientationModel {
    # This tool uses the forward and reverse read counts collected in the variant call
    # workflow. The basic idea of this tool is that true variants will be supported by
    # reads on both strands of DNA; whereas, artifacts will be supported by reads heavily
    # biased to one of the two strands. CollectF1R2 counts the number of reads on each
    # strand of DNA at a putatively variant site and stores this information along with
    # the nucleotide context. LearnReadOrientationModel uses the information from CollectF1R2
    # to estimate a prior probability that a site with a given context suffers from an
    # artifact. The two main artifacts of concern are OXOG (a result of sequencing technology)
    # and FFPE.

    input {
        String individual_id
        Array[File] f1r2_counts

        Runtime runtime_params
    }

    # Optional localization leads to cromwell error.
    # parameter_meta {
    #     f1r2_counts: {localization_optional: true}
    # }

    String output_name = individual_id + ".artifact_priors.tar.gz"
    Boolean f1r2_counts_empty = (length(f1r2_counts) == 0)

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.jar_override}

        if ~{f1r2_counts_empty} ; then
            echo "ERROR: f1r2_counts_tar_gz must be supplied and non empty."
            false
        fi

        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            LearnReadOrientationModel \
            ~{sep="' " prefix("-I '", f1r2_counts)}' \
            --output '~{output_name}'
    >>>

    output {
        File orientation_bias = output_name
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