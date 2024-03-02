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

import "patient.wdl" as p
import "workflow_arguments.wdl" as wfargs
import "runtime_collection.wdl" as rtc
import "runtimes.wdl" as rt
import "tasks.wdl"


workflow CallVariants {
    input {
        Patient patient
        WorkflowArguments args
        RuntimeCollection runtime_collection
    }

    scatter (tumor_sample in patient.tumor_samples) {
        scatter (seq_run in tumor_sample.sequencing_runs) {
            File seq_tumor_bams = seq_run.bam
            File seq_tumor_bais = seq_run.bai
        }
    }
    Array[File] tumor_bams = flatten(seq_tumor_bams)
    Array[File] tumor_bais = flatten(seq_tumor_bais)

    if (patient.has_normal) {
        scatter (normal_sample in patient.normal_samples) {
            scatter (seq_run in normal_sample.sequencing_runs) {
                File seq_normal_bams = seq_run.bam
                File seq_normal_bais = seq_run.bai
            }
            String normal_sample_names = normal_sample.bam_name
        }
        Array[File]? normal_bams = flatten(seq_normal_bams)
        Array[File]? normal_bais = flatten(seq_normal_bais)
    }

    scatter (interval_list in args.scattered_interval_list) {
    	call Mutect2 {
            input:
                interval_list = interval_list,
                ref_fasta = args.ref_fasta,
                ref_fasta_index = args.ref_fasta_index,
                ref_dict = args.ref_dict,
                individual_id = patient.name,
                tumor_bams = tumor_bams,
                tumor_bais = tumor_bais,
                normal_bams = normal_bams,
                normal_bais = normal_bais,
                normal_sample_names = normal_sample_names,
                force_call_alleles = args.force_call_alleles,
                force_call_alleles_idx = args.force_call_alleles_idx,
                panel_of_normals = args.snv_panel_of_normals,
                panel_of_normals_idx = args.snv_panel_of_normals_idx,
                germline_resource = args.germline_resource,
                germline_resource_tbi = args.germline_resource_tbi,
                make_bamout = args.make_bamout,
                get_orientation_bias_priors = args.run_orientation_bias_mixture_model,
                compress_output = args.compress_output,
                genotype_germline_sites = args.mutect2_genotype_germline_sites,
                native_pair_hmm_use_double_precision = args.mutect2_native_pair_hmm_use_double_precision,
                use_linked_de_bruijn_graph = args.mutect2_use_linked_de_bruijn_graph,
                recover_all_dangling_branches = args.mutect2_recover_all_dangling_branches,
                pileup_detection = args.mutect2_pileup_detection,
                downsampling_stride = args.mutect2_downsampling_stride,
                max_reads_per_alignment_start = args.mutect2_max_reads_per_alignment_start,
                m2_extra_args = args.mutect2_extra_args,
                runtime_params = runtime_collection.mutect2,
		}

        # TODO: add strelka
	}

    call tasks.MergeVCFs {
    	input:
            vcfs = Mutect2.vcf,
            vcfs_idx = Mutect2.vcf_idx,
            output_name = patient.name + ".unfiltered.merged",
            compress_output = args.compress_output,
            runtime_params = runtime_collection.merge_vcfs
    }

    if (args.make_bamout) {
        call tasks.MergeBams {
            input:
                ref_fasta = args.ref_fasta,
                ref_fasta_index = args.ref_fasta_index,
                ref_dict = args.ref_dict,
                bams = select_all(Mutect2.bam),
                bais = select_all(Mutect2.bai),
                merged_bam_name = patient.name + ".Mutect2.out",
                runtime_params = runtime_collection.merge_bams
        }
    }

    call MergeMutectStats {
        input:
            stats = Mutect2.vcf_stats,
            individual_id = patient.name,
            runtime_params = runtime_collection.merge_mutect_stats
    }

    if (args.run_orientation_bias_mixture_model) {
        call LearnReadOrientationModel {
            input:
                individual_id = patient.name,
                f1r2_counts = select_all(Mutect2.m2_artifact_priors),
                runtime_params = runtime_collection.learn_read_orientation_model
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
    String output_artifact_priors = individual_id + ".f1r2_counts.tar.gz"

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{select_first([runtime_params.jar_override, "/root/gatk.jar"])}
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            Mutect2 \
            --reference '~{ref_fasta}' \
            ~{sep="' " prefix("-I '", tumor_bams)}' \
            ~{sep="' " prefix("-I '", select_first([normal_bams, []]))}~{if normal_is_present then "'" else ""} \
            ~{sep="' " prefix("-normal '", select_first([normal_sample_names, []]))}~{if normal_is_present then "'" else ""} \
            --output '~{output_vcf}' \
            ~{"--intervals '" + interval_list + "'"} \
            ~{"--alleles '" + force_call_alleles + "'"} \
            ~{"-pon '" + panel_of_normals + "'"} \
            ~{if make_bamout then "--bam-output " + output_bam else ""} \
            ~{if get_orientation_bias_priors then "--f1r2-tar-gz " + output_artifact_priors else ""} \
            ~{"--germline-resource '" + germline_resource + "'"} \
            ~{if genotype_germline_sites then "--genotype-germline-sites true" else ""} \
            ~{if use_linked_de_bruijn_graph then "--linked-de-bruijn-graph true" else ""} \
            ~{if recover_all_dangling_branches then "--recover-all-dangling-branches true" else ""} \
            ~{if pileup_detection then "--pileup-detection true" else ""} \
            --smith-waterman FASTEST_AVAILABLE \
            --pair-hmm-implementation FASTEST_AVAILABLE \
            ~{if native_pair_hmm_use_double_precision then "--native-pair-hmm-use-double-precision true" else ""} \
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

    String output_name = individual_id + ".merged.stats"

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{select_first([runtime_params.jar_override, "/root/gatk.jar"])}
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

    # Optional localization leads to cromwell error.
    # parameter_meta {
    #     stats: {localization_optional: true}
    # }
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

    String output_name = individual_id + ".artifact_priors.tar.gz"

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{select_first([runtime_params.jar_override, "/root/gatk.jar"])}
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

    # Optional localization leads to cromwell error.
    # parameter_meta {
    #     f1r2_counts: {localization_optional: true}
    # }
}