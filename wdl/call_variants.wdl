version development

## Call variants with Mutect2
##
## Known issues:
## use_linked_de_bruijn_graph: This has trouble calling variants in complex regions.
##      Strongly recommended to use with recover_all_dangling_branches. This increases
##      compute cost though, and may still not guarantee that all variants are
##      being called. Ideally, run with and without and use the joint callset.
## --dont-use-soft-clipped-bases: https://gatk.broadinstitute.org/hc/en-us/community/posts/7259437599771-Evil-default
##      Using soft-clipped bases is necessary for indel calling. However, it can
##      introduce false positives, especially when DNA is sequenced from shorter
##      than 2xread_length fragments (like cfDNA), resulting in frequent read-through
##      into the adapters, which has lots of soft-clipped bases, and thus totally
##      unexpected mutation calls. The solution for that is to hard-clip the adapters
##      before variant calling.
##
## SOLVED IN GATK v4.6.0.0:
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

    # TODO: add Strelka2 and pipe via force-calling alleles into Mutect2

    scatter (interval_list in args.scattered_interval_list) {
        # Mutect1
         scatter (sample in patient.samples){ 
             scatter {seq_run in sample.sequencing_runs} {
                 if (sample.name != patient.matched_normal_sample.name) { # Mutect1 doesn't like same bam being referenced multiple times
                     call Mutect1 {
                         input:
                             interval_list = interval_list,
                             ref_fasta = args.files.ref_fasta,
                             ref_fasta_index = args.files.ref_fasta_index,
                             ref_dict = args.files.ref_dict,
                             germline_resource = args.files.germline_resource,
                             germline_resource_idx = args.files.germline_resource_idx,
                             panel_of_normals = args.files.snv_panel_of_normals,
                             panel_of_normals_idx = args.files.snv_panel_of_normals_idx,
                             contamination_table = sample.contamination,
 
                             individual_id = patient.name,
                             tumor_sample_name = sample.name,                                       
                             tumor_bam = seq_run.bam,
                             tumor_bai = seq_run.bai,
                             normal_sample_name = patient.matched_normal_sample.name,                
                             normal_bam = patient.matched_normal_sample.sequencing_runs[0].bam,
                             normal_bai = patient.matched_normal_sample.sequencing_runs[0].bai,
 
                             runtime_params = runtime_collection.mutect1,
                     }
                 }
             }
         }
         # Array[Array[File?]] Mutect1.vcf
         Array[File] mutect1_vcfs = select_all(flatten(Mutect1.mutect1_vcf))
 
         # TODO: Add Mutect1 stats + force call alleles TASK
 
        # Mutect2
    	call Mutect2 {
            input:
                interval_list = interval_list,
                ref_fasta = args.files.ref_fasta,
                ref_fasta_index = args.files.ref_fasta_index,
                ref_dict = args.files.ref_dict,
                individual_id = patient.name,
                tumor_bams = tumor_bams,
                tumor_bais = tumor_bais,
                normal_bams = normal_bams,
                normal_bais = normal_bais,
                normal_sample_names = normal_sample_names,
                force_call_alleles = args.files.force_call_alleles,                                      # TODO: replace with Mutect1 + force_call_alleles input
                force_call_alleles_idx = args.files.force_call_alleles_idx,                              # TODO: replace with Mutect1 + force_call_alleles input
                panel_of_normals = args.files.snv_panel_of_normals,
                panel_of_normals_idx = args.files.snv_panel_of_normals_idx,
                germline_resource = args.files.germline_resource,
                germline_resource_idx = args.files.germline_resource_idx,
                make_bamout = args.make_bamout,
                get_orientation_bias_priors = args.run_orientation_bias_mixture_model,
                compress_output = args.compress_output,
                genotype_germline_sites = args.mutect2_genotype_germline_sites,
                native_pair_hmm_use_double_precision = args.mutect2_native_pair_hmm_use_double_precision,
                dont_use_soft_clipped_bases = args.mutect2_dont_use_soft_clipped_bases,
                use_linked_de_bruijn_graph = args.mutect2_use_linked_de_bruijn_graph,
                recover_all_dangling_branches = args.mutect2_recover_all_dangling_branches,
                pileup_detection = args.mutect2_pileup_detection,
                downsampling_stride = args.mutect2_downsampling_stride,
                pcr_snv_qual = args.mutect2_pcr_snv_qual,
                pcr_indel_qual = args.mutect2_pcr_indel_qual,
                max_reads_per_alignment_start = args.mutect2_max_reads_per_alignment_start,
                m2_extra_args = args.mutect2_extra_args,
                runtime_params = runtime_collection.mutect2,
		}
	}

    call tasks.MergeVCFs {
    	input:
            vcfs = Mutect2.vcf,
            vcfs_idx = Mutect2.vcf_idx,
            output_name = patient.name,
            compress_output = args.compress_output,
            runtime_params = runtime_collection.merge_vcfs
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
        Array[File?] bams = select_all(Mutect2.bam)
        Array[File?] bais = select_all(Mutect2.bai)
        File? orientation_bias = LearnReadOrientationModel.orientation_bias
    }
}

task Mutect1 {
     input {
         File interval_list
         File ref_fasta
         File ref_fasta_index
         File ref_dict
         File germline_resource
         File germline_resource_idx
         File? panel_of_normals
         File? panel_of_normals_idx
 
         String individual_id
         String tumor_sample_name
         File tumor_bam
         File tumor_bai
         String? normal_sample_name
         File? normal_bam
         File? normal_bai
 
         Int downsample_to_coverage = 99999
         Int max_alt_alleles_in_normal_count = 10000
         Int max_alt_alleles_in_normal_qscore_sum = 10000
         File? contamination_table                         
         Int tumor_f_pretest = 0
         Float initial_tumor_lod = 0.5
         Int tumor_lod = 0
 
         Runtime runtime_params
     }
 
     command <<<
         set -euxo pipefail
 
         # Dynamically calculate memory usage for Cromwell's retry_with_more_memory feature:
         # This only works in a cloud VM! For HPC, you need to check cgroups.
         machine_mb=$( free -m | awk '{ print $NF }' | head -n 2 | tail -1 )   # change to -g for GB output from free
 
         ten_percent_machine_mb=$(( $machine_mb * 10 / 100 ))
         runtime_overhead_mb=~{runtime_params.machine_mem - runtime_params.command_mem}
 
         # Leave 10% of memory or runtime_overhead for overhead, whichever is larger
         if [ $ten_percent_machine_mb -gt $runtime_overhead_mb ]; then
             overhead_mb=$ten_percent_machine_mb
         else
             overhead_mb=$runtime_overhead_mb
         fi
         command_mb=$(( $machine_mb  - $overhead_mb ))
 
         echo ""
         echo "Available memory: $machine_mb MB."
         echo "Overhead: $overhead_mb MB."
         echo "Using $command_mb MB of memory for Mutect1."
         echo ""
 
         # Parse contamination_table
         if [ -f ~{contamination_table} ]; then
             fraction_contamination=$(tail -n 1 ~{contamination_table} | awk '{print $2}')
         else
             fraction_contamination=0
         fi
 
         # Run MuTect1 (without force-calling, will be done by Mutect2)
         java "-Xmx${command_mb}m" -jar /usr/local/bin/mutect-1.1.7.jar \
             --analysis_type MuTect \
             --reference_sequence ~{ref_fasta} \
             --intervals ~{interval_list} \
             --tumor_sample_name ~{tumor_sample_name} \
             -I:tumor ~{tumor_bam} \
             ~{"--normal_sample_name '" + normal_sample_name + "'"} \
             ~{"-I:normal '" + normal_bam + "'"} \         
             --normal_panel ~{panel_of_normals} \                                                # TODO: this must be VCF v4.1, currently ours is v4.2
             --dbsnp ~{germline_resource} \                                                      # TODO: this must be VCF v4.1, currently ours is v4.2
             --downsample_to_coverage ~{downsample_to_coverage} \
             --max_alt_alleles_in_normal_count ~{max_alt_alleles_in_normal_count} \
             --max_alt_alleles_in_normal_qscore_sum ~{max_alt_alleles_in_normal_qscore_sum} \
             --fraction_contamination ~{fraction_contamination} \                              
             --tumor_f_pretest ~{tumor_f_pretest} \
             --initial_tumor_lod ~{initial_tumor_lod} \
             --tumor_lod ~{tumor_lod} \
             --out ~{tumor_sample_name}.MuTect1.call_stats.txt \
             --coverage_file ~{tumor_sample_name}.MuTect1.coverage.wig.txt \
             --power_file ~{tumor_sample_name}.MuTect1.power.wig.txt \
             --vcf ~{tumor_sample_name}.MuTect1.vcf                                         #(not sure if needed) VCF output of mutation candidates
     >>>
     
     output {
         File mutect1_call_stats="~{tumor_sample_name}.MuTect1.call_stats.txt"      
         File mutect1_coverage_wig="~{tumor_sample_name}.MuTect1.coverage.wig.txt"
         File mutect1_power_wig="~{tumor_sample_name}.MuTect1.power.wig.txt"
         File mutect1_vcf="~{tumor_sample_name}.MuTect1.vcf"
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
 
 # TASK: mutect1 output, drop sample columns and union with force_call_alleles (COSMIC), subset to intervals to reduce size (move to tasks.wdl later)
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
        File? germline_resource_idx

        Boolean genotype_germline_sites = false
        Boolean native_pair_hmm_use_double_precision = true
        Boolean dont_use_soft_clipped_bases = false
        Boolean use_linked_de_bruijn_graph = true
        Boolean recover_all_dangling_branches = true
        Boolean pileup_detection = true

        # The linked de-Bruijn graph implementation has trouble calling variants
        # in complex regions, even when recovering all dangling branches.
        # Reducing the downsampling by increasing the following parameters might
        # solve the issue. It increases compute cost though.
        Int downsampling_stride = 1
        Int max_reads_per_alignment_start = 50

        # Increase for high quality panel sequencing data
        # 40 = 1e-2 error rate (standard WES)
        # 80 = 1e-4 error rate (duplex sequencing)
        Int pcr_snv_qual = 40
        Int pcr_indel_qual = 40

        String? m2_extra_args

        Boolean get_orientation_bias_priors = true
        Boolean compress_output = false
        Boolean make_bamout = false

        Runtime runtime_params
    }

    Boolean normal_is_present = defined(normal_bams) && (length(select_first([normal_bams])) > 0)

    String output_vcf = individual_id + if compress_output then ".vcf.gz" else ".vcf"
    String output_vcf_idx = output_vcf + if compress_output then ".tbi" else ".idx"
    String output_stats = output_vcf + ".stats"

    String output_bam = individual_id + ".bamout.bam"
    String output_bai = individual_id + ".bamout.bai"
    String output_artifact_priors = individual_id + ".f1r2_counts.tar.gz"

    String dollar = "$"

    command <<<
        set -e

        export GATK_LOCAL_JAR=~{select_first([runtime_params.jar_override, "/root/gatk.jar"])}

        # Dynamically calculate memory usage for Cromwell's retry_with_more_memory feature:
        # This only works in a cloud VM! For HPC, you need to check cgroups.
        machine_mb=$( free -m | awk '{ print $NF }' | head -n 2 | tail -1 )   # change to -g for GB output from free

        ten_percent_machine_mb=$(( $machine_mb * 10 / 100 ))
        runtime_overhead_mb=~{runtime_params.machine_mem - runtime_params.command_mem}

        # Leave 10% of memory or runtime_overhead for overhead, whichever is larger
        if [ $ten_percent_machine_mb -gt $runtime_overhead_mb ]; then
            overhead_mb=$ten_percent_machine_mb
        else
            overhead_mb=$runtime_overhead_mb
        fi
        command_mb=$(( $machine_mb  - $overhead_mb ))

        n_threads=$(nproc)

        echo ""
        echo "Available memory: $machine_mb MB."
        echo "Overhead: $overhead_mb MB."
        echo "Using $command_mb MB of memory for GATK and $n_threads threads."
        echo ""

        echo ""
        # This warning is from enabling the linked de-Bruijn graph implementation:
        echo "Suppressing the following warning messages: 'Dangling End recovery killed because of a loop (findPath)'"

        # This warning appears if multiple sequencing runs from the same sample are supplied:
        echo "Suppressing the following warning messages: 'More than two reads with the same name found. Using two reads randomly to combine as a fragment.'"
        echo ""

        gatk --java-options "-Xmx~{dollar}{command_mb}m" \
            Mutect2 \
            --reference '~{ref_fasta}' \
            ~{sep="' " prefix("-I '", tumor_bams)}' \
            ~{sep="' " prefix("--read-index '", tumor_bais)}' \
            ~{if normal_is_present then "-I '" else ""}~{default="" sep="' -I '" normal_bams}~{if normal_is_present then "'" else ""} \
            ~{if normal_is_present then "--read-index '" else ""}~{default="" sep="' --read-index '" normal_bais}~{if normal_is_present then "'" else ""} \
            ~{if normal_is_present then "-normal '" else ""}~{default="" sep="' -normal '" normal_sample_names}~{if normal_is_present then "'" else ""} \
            --output '~{output_vcf}' \
            ~{"--intervals '" + interval_list + "'"} \
            ~{"--alleles '" + force_call_alleles + "'"} \
            ~{"-pon '" + panel_of_normals + "'"} \
            ~{if make_bamout then "--bam-output '" + output_bam + "'" else ""} \
            ~{if get_orientation_bias_priors then "--f1r2-tar-gz '" + output_artifact_priors + "'" else ""} \
            ~{"--germline-resource '" + germline_resource + "'"} \
            ~{if genotype_germline_sites then "--genotype-germline-sites true" else ""} \
            ~{if dont_use_soft_clipped_bases then "--dont-use-soft-clipped-bases true" else ""} \
            ~{if use_linked_de_bruijn_graph then "--linked-de-bruijn-graph true" else ""} \
            ~{if recover_all_dangling_branches then "--recover-all-dangling-branches true" else ""} \
            ~{if pileup_detection then "--pileup-detection true" else ""} \
            --pcr-snv-qual ~{pcr_snv_qual} \
            --pcr-indel-qual ~{pcr_indel_qual} \
            --smith-waterman FASTEST_AVAILABLE \
            --pair-hmm-implementation FASTEST_AVAILABLE \
            --native-pair-hmm-threads $n_threads \
            ~{if native_pair_hmm_use_double_precision then "--native-pair-hmm-use-double-precision true" else ""} \
            ~{"--downsampling-stride " + downsampling_stride} \
            ~{"--max-reads-per-alignment-start " + max_reads_per_alignment_start} \
            --seconds-between-progress-updates 60 \
            ~{m2_extra_args} \
            2> >( \
                grep -v 'Dangling End recovery killed because of a loop (findPath)' | \
                grep -v 'More than two reads with the same name found' \
                >&2 \
            )
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

    parameter_meta {
        interval_list: {localization_optional: true}
        ref_fasta: {localization_optional: true}
        ref_fasta_index: {localization_optional: true}
        ref_dict: {localization_optional: true}
        tumor_bams: {localization_optional: true}
        tumor_bais: {localization_optional: true}
        normal_bams: {localization_optional: true}
        normal_bais: {localization_optional: true}
        force_call_alleles: {localization_optional: true}
        force_call_alleles_idx: {localization_optional: true}
        panel_of_normals: {localization_optional: true}
        panel_of_normals_idx: {localization_optional: true}
        germline_resource: {localization_optional: true}
        germline_resource_idx: {localization_optional: true}
    }
}

task MergeMutectStats {
    input {
        Array[File]+ stats
        String individual_id

        Runtime runtime_params
    }

    String output_name = individual_id + ".stats"

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
    # CAUTION: The EM algorithm may not converge for some context, however, it will still
    # output orientation bias priors, which are used in FilterMutectCalls to determine
    # variant status. If the EM did not converge, the artifact-classification may
    # lead to false positives and false negatives.

    input {
        String individual_id
        Array[File] f1r2_counts

        Int max_depth = 200
        Int max_num_em_iterations = 100  # default 20

        Runtime runtime_params
    }

    String output_name = individual_id + ".artifact_priors.tar.gz"

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{select_first([runtime_params.jar_override, "/root/gatk.jar"])}
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            LearnReadOrientationModel \
            ~{sep="' " prefix("-I '", f1r2_counts)}' \
            --output '~{output_name}' \
            --max-depth ~{max_depth} \
            --num-em-iterations ~{max_num_em_iterations}
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