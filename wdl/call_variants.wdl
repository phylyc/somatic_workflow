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

import "shard.wdl" as sh
import "sequencing_run.wdl" as seqrun
import "patient.wdl" as p
import "patient.update_samples.wdl" as p_update_s
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

    if (patient.has_normal) {
        scatter (normal_sample in patient.normal_samples) {
            String normal_sample_names = normal_sample.bam_name
        }
    }

    # TODO: add Strelka2

    # TODO: Estimate shard sizes based on bais and intervals

    scatter (shard in patient.shards) {
        if (size(shard.raw_calls_mutect2_vcf) == 0 && !shard.skip) {
            if (args.run_variant_calling_mutect1) {
                # Mutect1
                # Old GATK does not stream input files, so we have to localize them.
                # Localizing the full bams is very expensive, so we first subset them
                # to the sharded intervals.
                scatter (sample in patient.samples) {
                    scatter (seq_run in sample.sequencing_runs) {
                        call tasks.PrintReads as SubsetToShard {
                            input:
                                interval_list = shard.intervals,
                                ref_fasta = args.files.ref_fasta,
                                ref_fasta_index = args.files.ref_fasta_index,
                                ref_dict = args.files.ref_dict,
                                prefix = seq_run.sample_name,
                                bams = [seq_run.bam],
                                bais = [seq_run.bai],
                                runtime_params = runtime_collection.subset_bam_to_shard
                        }
                        call seqrun.UpdateSequencingRun as SeqRunShard {
                            input:
                                sequencing_run = seq_run,
                                bam = SubsetToShard.output_bam,
                                bai = SubsetToShard.output_bai
                        }
                    }
                }
                call p_update_s.UpdateSamples as PatientShard {
                    input:
                        patient = patient,
                        sequencing_runs = SeqRunShard.updated_sequencing_run
                }

                scatter (sample in PatientShard.updated_patient.samples) {
                    scatter (seq_run in sample.sequencing_runs) {
                        # Only supply matched normal sample for tumor samples as
                        # Mutect1 does not like the same bam for tumor and normal.
                        if (sample.is_tumor && defined(PatientShard.updated_patient.matched_normal_sample)) {
                            Sample matched_normal_sample = select_first([PatientShard.updated_patient.matched_normal_sample])
                            String matched_normal_sample_name = matched_normal_sample.name
                            SequencingRun matched_normal_sample_seq_run = select_first(matched_normal_sample.sequencing_runs)
                            File matched_normal_bam = matched_normal_sample_seq_run.bam
                            File matched_normal_bai = matched_normal_sample_seq_run.bai
                        }
                        call Mutect1 {
                            input:
                                interval_list = shard.intervals,
                                ref_fasta = args.files.ref_fasta,
                                ref_fasta_index = args.files.ref_fasta_index,
                                ref_dict = args.files.ref_dict,
                                germline_resource = args.files.germline_resource_v4_1,
                                germline_resource_idx = args.files.germline_resource_v4_1_idx,
                                panel_of_normals = args.files.snv_panel_of_normals_v4_1,
                                panel_of_normals_idx = args.files.snv_panel_of_normals_v4_1_idx,
                                contamination_table = sample.contamination_table,

                                patient_id = patient.name,
                                tumor_sample_name = seq_run.sample_name,
                                tumor_bam = seq_run.bam,
                                tumor_bai = seq_run.bai,
                                normal_sample_name = matched_normal_sample_name,
                                normal_bam = matched_normal_bam,
                                normal_bai = matched_normal_bai,

                                initial_tumor_lod = args.mutect1_initial_tumor_lod,
                                tumor_lod = args.mutect1_tumor_lod_to_emit,

                                runtime_params = runtime_collection.mutect1,
                        }
                    }
                }
                Array[File] mutect1_vcfs = select_all(flatten(Mutect1.mutect1_vcf))
                Array[File] mutect1_vcfs_idx = select_all(flatten(Mutect1.mutect1_vcf_idx))

                call MergeMutect1ForceCallVCFs {
                    input:
                        interval_list = shard.intervals,
                        ref_fasta = args.files.ref_fasta,
                        ref_fasta_index = args.files.ref_fasta_index,
                        ref_dict = args.files.ref_dict,

                        mutect1_vcfs = mutect1_vcfs,
                        mutect1_vcfs_idx = mutect1_vcfs_idx,
                        force_call_alleles = args.files.force_call_alleles,
                        force_call_alleles_idx = args.files.force_call_alleles_idx,
                        runtime_params = runtime_collection.merge_mutect1_forcecall_vcfs
                }

                Float shard_bams_size = size(flatten(SubsetToShard.output_bam), "GB")
            }

            scatter (m2_sample in patient.samples) {
                scatter (m2_seq_run in m2_sample.sequencing_runs) {
                    File full_bam = m2_seq_run.bam
                    File full_bai = m2_seq_run.bai
                }
                String bam_names = m2_sample.bam_name
                String sample_names = m2_sample.name
            }
            Float full_bams_size = size(flatten(full_bam), "GB")

            Int m2_diskGB = (
                runtime_collection.mutect2.disk
                + if args.make_bamout then ceil(1.2 * select_first([shard_bams_size, full_bams_size / length(patient.shards)])) else 0
            )

            if (shard.is_high_mem) {
                call rt.UpdateRuntimeParameters as Mutect2Runtime {
                    input:
                        runtime_params = runtime_collection.mutect2,
                        machine_mem = ceil(args.mutect2_high_mem_factor * runtime_collection.mutect2.machine_mem),
                        command_mem = ceil(args.mutect2_high_mem_factor * runtime_collection.mutect2.command_mem),
                }
            }

            call Mutect2 {
                input:
                    interval_list = shard.intervals,
                    ref_fasta = args.files.ref_fasta,
                    ref_fasta_index = args.files.ref_fasta_index,
                    ref_dict = args.files.ref_dict,
                    patient_id = patient.name,
                    bams = flatten(full_bam),
                    bais = flatten(full_bai),
                    normal_sample_names = normal_sample_names,
                    bam_names = bam_names,
                    sample_names = sample_names,
                    force_call_alleles = if defined(MergeMutect1ForceCallVCFs.merged_vcf) then MergeMutect1ForceCallVCFs.merged_vcf else args.files.force_call_alleles,
                    force_call_alleles_idx = if defined(MergeMutect1ForceCallVCFs.merged_vcf_idx) then MergeMutect1ForceCallVCFs.merged_vcf_idx else args.files.force_call_alleles_idx,
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
                    initial_tumor_lod = args.mutect2_initial_tumor_lod,
                    tumor_lod_to_emit = args.mutect2_tumor_lod_to_emit,
                    max_reads_per_alignment_start = args.mutect2_max_reads_per_alignment_start,
                    m2_extra_args = args.mutect2_extra_args,
                    diskGB = m2_diskGB,
                    runtime_params = select_first([Mutect2Runtime.params, runtime_collection.mutect2]),
            }

            call sh.UpdateShard as AddMutect2Calls {
                input:
                    shard = shard,
                    raw_calls_mutect2_vcf = Mutect2.vcf,
                    raw_calls_mutect2_vcf_idx = Mutect2.vcf_idx,
                    raw_mutect2_stats = Mutect2.vcf_stats,
                    raw_mutect2_bam_out = Mutect2.bam,
                    raw_mutect2_bai_out = Mutect2.bai,
                    raw_mutect2_artifact_priors = Mutect2.m2_artifact_priors,
            }
        }
        Shard updated_shard = select_first([AddMutect2Calls.updated_shard, shard])

        File? raw_mutect2_vcf = updated_shard.raw_calls_mutect2_vcf
        File? raw_mutect2_vcf_idx = updated_shard.raw_calls_mutect2_vcf_idx
        File? raw_mutect2_stats = updated_shard.raw_mutect2_stats
        File? raw_mutect2_bam = updated_shard.raw_mutect2_bam_out
        File? raw_mutect2_bai = updated_shard.raw_mutect2_bai_out
        File? raw_mutect2_artifact_priors = updated_shard.raw_mutect2_artifact_priors
	}

    if (!defined(patient.raw_snv_calls_vcf)) {
        call tasks.GatherVCFs {
            input:
                vcfs = select_all(raw_mutect2_vcf),
                vcfs_idx = select_all(raw_mutect2_vcf_idx),
                output_name = patient.name,
                compress_output = args.compress_output,
                runtime_params = runtime_collection.gather_vcfs
        }
    }

    if (!defined(patient.mutect2_stats)) {
        call MergeMutectStats {
            input:
                stats = select_all(raw_mutect2_stats),
                patient_id = patient.name,
                runtime_params = runtime_collection.merge_mutect_stats
        }
    }

    if (args.run_orientation_bias_mixture_model && (!defined(patient.orientation_bias))) {
        call LearnReadOrientationModel {
            input:
                patient_id = patient.name,
                f1r2_counts = select_all(raw_mutect2_artifact_priors),
                runtime_params = runtime_collection.learn_read_orientation_model
        }
    }

    call p.UpdatePatient {
        input:
            patient = patient,
            shards = updated_shard,
            raw_snv_calls_vcf = GatherVCFs.merged_vcf,
            raw_snv_calls_vcf_idx = GatherVCFs.merged_vcf_idx,
            mutect2_stats = MergeMutectStats.merged_stats,
            orientation_bias = LearnReadOrientationModel.orientation_bias,
    }

    output {
        Patient updated_patient = UpdatePatient.updated_patient
    }
}

task Mutect1 {
    # Documentation at https://gist.github.com/sbamin/041a4bdb7c5b184321b1468345bd2aa8
    # Tool defaults:
    #    downsample_to_coverage=1000
    #    baqGapOpenPenalty=40.0
    #    quantize_quals=0
    #    preserve_qscores_less_than=6
    #    globalQScorePrior=-1.0
    #    initial_tumor_lod=4.0
    #    tumor_lod=6.3
    #    minimum_mutation_cell_fraction=0.0
    #    normal_lod=2.2
    #    normal_artifact_lod=1.0
    #    strand_artifact_lod=2.0
    #    strand_artifact_power_threshold=0.9
    #    dbsnp_normal_lod=5.5
    #    minimum_normal_allele_fraction=0.0
    #    tumor_f_pretest=0.005
    #    min_qscore=5
    #    ap_events_threshold=3
    #    heavily_clipped_read_fraction=0.3
    #    fraction_mapq0_threshold=0.5
    #    pir_median_threshold=10.0
    #    pir_mad_threshold=3.0
    #    required_maximum_alt_allele_mapping_quality_score=20
    #    max_alt_alleles_in_normal_count=2
    #    max_alt_alleles_in_normal_qscore_sum=20
    #    max_alt_allele_in_normal_fraction=0.03
    input {
        File interval_list
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        File? germline_resource
        File? germline_resource_idx
        File? panel_of_normals
        File? panel_of_normals_idx

        String patient_id
        String tumor_sample_name
        File tumor_bam
        File tumor_bai
        String? normal_sample_name
        File? normal_bam
        File? normal_bai

        Int downsample_to_coverage = 99999
        File? contamination_table
        Int tumor_f_pretest = 0
        Float? initial_tumor_lod
        Float? tumor_lod

        Boolean only_passing_calls = true

        Runtime runtime_params
    }

    String call_stats = tumor_sample_name + ".MuTect1.call_stats.txt"
    String coverage_wig = tumor_sample_name + ".MuTect1.coverage.wig.txt"
    String power_wig = tumor_sample_name + ".MuTect1.power.wig.txt"
    String mutect1_vcf = tumor_sample_name + ".MuTect1.vcf"
    String mutect1_vcf_idx = tumor_sample_name + ".MuTect1.vcf.idx"

    # COMPUTE DISK SIZE
    Int diskGB = runtime_params.disk + ceil(size(ref_fasta, "GB") + size(germline_resource, "GB") + size(panel_of_normals, "GB") + size(tumor_bam, "GB") + size(normal_bam, "GB"))

    String dollar = "$"

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

        n_threads=$(nproc)

        echo ""
        echo "Available memory: $machine_mb MB."
        echo "Overhead: $overhead_mb MB."
        echo "Using $command_mb MB of memory for Mutect1 and 1/$n_threads threads."
        echo ""

        # Parse contamination_table
        if [ "~{defined(contamination_table)}" == "true" ]; then
            fraction_contamination=$(tail -n 1 '~{contamination_table}' | awk '{print $2}')
        else
            fraction_contamination=0
        fi

        # Run MuTect1 (without force-calling, will be done by Mutect2)
        java "-Xmx~{dollar}{command_mb}m" -jar /muTect-1.1.6.jar \
            --analysis_type MuTect \
            --reference_sequence '~{ref_fasta}' \
            --intervals '~{interval_list}' \
            --tumor_sample_name '~{tumor_sample_name}' \
            -I:tumor '~{tumor_bam}' \
            ~{"--normal_sample_name '" + normal_sample_name + "'"} \
            ~{"-I:normal '" + normal_bam + "'"} \
            ~{"--normal_panel '" + panel_of_normals + "'"} \
            ~{"--dbsnp '" + germline_resource + "'"} \
            --downsample_to_coverage ~{downsample_to_coverage} \
            --fraction_contamination $fraction_contamination \
            ~{"--tumor_f_pretest " + tumor_f_pretest} \
            ~{"--initial_tumor_lod " + initial_tumor_lod} \
            ~{"--tumor_lod " + tumor_lod} \
            --out '~{call_stats}' \
            --coverage_file '~{coverage_wig}' \
            --power_file '~{power_wig}' \
            ~{if only_passing_calls then "--only_passing_calls" else ""} \
            --vcf '~{mutect1_vcf}' \
            --filter_reads_with_N_cigar \
            --filter_mismatching_base_and_quals \
            --filter_bases_not_stored
    >>>

    output {
        File mutect1_call_stats = call_stats
        File mutect1_coverage_wig = coverage_wig
        File mutect1_power_wig = power_wig
        File mutect1_vcf = mutect1_vcf
        File mutect1_vcf_idx = mutect1_vcf_idx
    }

    runtime {
        docker: runtime_params.docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: runtime_params.machine_mem + " MB"
        runtime_minutes: runtime_params.runtime_minutes
        disks: "local-disk " + diskGB + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }
}

task MergeMutect1ForceCallVCFs {
    input {
        File interval_list
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        Array[File] mutect1_vcfs # not gzipped
        Array[File] mutect1_vcfs_idx # not gzipped
        File? force_call_alleles # gzipped
        File? force_call_alleles_idx # gzipped

        Runtime runtime_params
    }

    String base = "merged_mutect1"
    String mutect1_vcf = base + ".vcf.gz"
    String mutect1_vcf_idx = mutect1_vcf + ".tbi"
    String mutect1_pass_no_genotypes_vcf = base + ".pass.no_genotypes.vcf.gz"
    String mutect1_pass_no_genotypes_vcf_idx = mutect1_pass_no_genotypes_vcf + ".tbi"
    String force_call_alleles_subset_vcf = "forcecall.subset.vcf.gz"
    String force_call_alleles_subset_vcf_idx = force_call_alleles_subset_vcf + ".tbi"
    String mutect1_pass_no_genotypes_forcecall_vcf = base + ".pass.no_genotypes.forcecall_alleles.vcf.gz"
    String mutect1_pass_no_genotypes_forcecall_vcf_idx = mutect1_pass_no_genotypes_forcecall_vcf + ".tbi"
    String mutect1_pass_no_genotypes_forcecall_dedup_vcf = base + ".pass.no_genotypes.forcecall_alleles.dedup.vcf.gz"
    String mutect1_pass_no_genotypes_forcecall_dedup_vcf_idx = mutect1_pass_no_genotypes_forcecall_dedup_vcf + ".tbi"

    # COMPUTE DISK SIZE
    Int diskGB = runtime_params.disk + ceil(2 * size(mutect1_vcfs, "GB"))

    String dollar = "$"

    command <<<
        set -euxo pipefail

        # Loop over each VCF, drop genotypes, add to empty array
        no_gt_vcfs=()
        for vcf in ~{sep="' " prefix(" '", mutect1_vcfs)}'; do
            # Construct a new output file name, inputs: ".vcf", outputs: ".noGT.vcf.gz"
            out_vcf="${vcf%.vcf}.noGT.vcf.gz"

            # Drop genotype information and compress the output
            bcftools view --drop-genotypes "$vcf" -Oz > "$out_vcf"
            bcftools index -t "$out_vcf"

            # Append the processed file to the array
            no_gt_vcfs+=("$out_vcf")
        done

        # Concat all samples VCFs from shardX and compress the merged output
        bcftools concat -a "~{dollar}{no_gt_vcfs[@]}" -Oz > '~{mutect1_vcf}'

        # Drop FORMAT and sample name columns (i.e. only keep #CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO columns)
        bcftools view -i 'FILTER=="PASS"' '~{mutect1_vcf}' -Oz > '~{mutect1_pass_no_genotypes_vcf}'
        rm -f '~{mutect1_vcf}' '~{mutect1_vcf_idx}'

        gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
            IndexFeatureFile \
            --input '~{mutect1_pass_no_genotypes_vcf}' \
            --output '~{mutect1_pass_no_genotypes_vcf_idx}'

        if [ "~{defined(force_call_alleles)}" == "true" ]; then
            gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
                SelectVariants \
                ~{"-R '" + ref_fasta + "'"} \
                ~{"-L '" + interval_list + "'"} \
                -V '~{force_call_alleles}' \
                --lenient \
                -O '~{force_call_alleles_subset_vcf}'

            # Union with force_call_alleles vcf
            bcftools concat -a '~{mutect1_pass_no_genotypes_vcf}' '~{force_call_alleles_subset_vcf}' -Oz > '~{mutect1_pass_no_genotypes_forcecall_vcf}'
            rm -f '~{mutect1_pass_no_genotypes_vcf}' '~{mutect1_pass_no_genotypes_vcf_idx}' '~{force_call_alleles_subset_vcf}' '~{force_call_alleles_subset_vcf_idx}'

            # Deduplicate based on #CHROM, POS, REF, ALT (sorting is needed for --rm-dup to work properly)
            bcftools sort '~{mutect1_pass_no_genotypes_forcecall_vcf}' | \
                bcftools norm --rm-dup both -Oz \
                > '~{mutect1_pass_no_genotypes_forcecall_dedup_vcf}'
            rm -f '~{mutect1_pass_no_genotypes_forcecall_vcf}'

            gatk --java-options "-Xmx~{runtime_params.command_mem}m" \
                IndexFeatureFile \
                --input '~{mutect1_pass_no_genotypes_forcecall_dedup_vcf}' \
                --output '~{mutect1_pass_no_genotypes_forcecall_dedup_vcf_idx}'
        else
            mv '~{mutect1_pass_no_genotypes_vcf}' '~{mutect1_pass_no_genotypes_forcecall_dedup_vcf}'
            mv '~{mutect1_pass_no_genotypes_vcf_idx}' '~{mutect1_pass_no_genotypes_forcecall_dedup_vcf_idx}'
        fi
    >>>

    output {
        File merged_vcf = mutect1_pass_no_genotypes_forcecall_dedup_vcf
        File merged_vcf_idx = mutect1_pass_no_genotypes_forcecall_dedup_vcf_idx
    }

    runtime {
        docker: runtime_params.docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: runtime_params.machine_mem + " MB"
        runtime_minutes: runtime_params.runtime_minutes
        disks: "local-disk " + diskGB + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

    parameter_meta {
        interval_list: {localization_optional: true}
        ref_fasta: {localization_optional: true}
        ref_fasta_index: {localization_optional: true}
        ref_dict: {localization_optional: true}
        force_call_alleles: {localization_optional: true}
        force_call_alleles_idx: {localization_optional: true}
    }
}

task Mutect2 {
    input {
        File? interval_list
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        String patient_id
        Array[File] bams
        Array[File] bais
        Array[String]? normal_sample_names

        Array[String] bam_names
        Array[String] sample_names

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
        Int? downsampling_stride
        Int? max_reads_per_alignment_start

        # Increase for high quality panel sequencing data
        # 40 = 1e-2 error rate (standard WES)
        # 80 = 1e-4 error rate (duplex sequencing)
        Int pcr_snv_qual = 40
        Int pcr_indel_qual = 40

        Float? initial_tumor_lod
        Float? tumor_lod_to_emit

        String? m2_extra_args

        Boolean get_orientation_bias_priors = true
        Boolean compress_output = false
        Boolean make_bamout = false

        Int diskGB = runtime_params.disk + if make_bamout then ceil(1.2 * size(bams, "GB")) else 0
        Runtime runtime_params
    }

    Boolean normal_is_present = defined(normal_sample_names) && (length(select_first([normal_sample_names])) > 0)

    String output_vcf = patient_id + if compress_output then ".vcf.gz" else ".vcf"
    String output_vcf_idx = output_vcf + if compress_output then ".tbi" else ".idx"
    String output_stats = output_vcf + ".stats"

    String output_bam = patient_id + ".bamout.bam"
    String output_bai = patient_id + ".bamout.bai"
    String output_artifact_priors = patient_id + ".f1r2_counts.tar.gz"

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
            ~{sep="' " prefix("-I '", bams)}' \
            ~{sep="' " prefix("--read-index '", bais)}' \
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
            ~{"--initial-tumor-lod " + initial_tumor_lod} \
            ~{"--tumor-lod-to-emit " + tumor_lod_to_emit} \
            --seconds-between-progress-updates 60 \
            ~{m2_extra_args} \
            2> >( \
                grep -v 'Dangling End recovery killed because of a loop (findPath)' | \
                grep -v 'More than two reads with the same name found' \
                >&2 \
            )

        #####
        # Add map of bam names to chosen samples names to header. This makes it
        # easier for downstream tools to map called variants to samples.

        # Convert comma-separated strings to arrays
        IFS=',' read -r -a bam_names <<< "~{sep="," bam_names}"
        IFS=',' read -r -a sample_names <<< "~{sep="," sample_names}"
        # Create string to inset into header
        header=""
        for i in "~{dollar}{!bam_names[@]}"; do
            header+="##BAM_NAME:~{dollar}{bam_names[i]}=SAMPLE_NAME:~{dollar}{sample_names[i]}\n"
        done
        header="~{dollar}{header%??}"  # Trim the trailing newline
        # Insert string
        if [ "~{compress_output}" == "true" ]; then
            zcat '~{output_vcf}' | \
                awk -v header="$header" '/^#CHROM/ {print header} {print}' | \
                bgzip > 'tmp.~{output_vcf}'
        else
            awk -v header="$header" '/^#CHROM/ {print header} {print}' '~{output_vcf}' > 'tmp.~{output_vcf}'
        fi
        mv 'tmp.~{output_vcf}' '~{output_vcf}'
        bcftools index -f -t -o '~{output_vcf_idx}' '~{output_vcf}'
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
        disks: "local-disk " + diskGB + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

    parameter_meta {
        interval_list: {localization_optional: true}
        ref_fasta: {localization_optional: true}
        ref_fasta_index: {localization_optional: true}
        ref_dict: {localization_optional: true}
        bams: {localization_optional: true}
        bais: {localization_optional: true}
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
        String patient_id

        Runtime runtime_params
    }

    String output_name = patient_id + ".stats"

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
        String patient_id
        Array[File] f1r2_counts

        Int max_depth = 200
        Int max_num_em_iterations = 100  # default 20

        Runtime runtime_params
    }

    String output_name = patient_id + ".artifact_priors.tar.gz"

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