version development

import "sequencing_run.wdl" as seqrun
import "sample.wdl" as s
import "patient.wdl" as p
import "patient.define.wdl" as p_def
import "patient.merge.wdl" as p_merge
import "patient.update_samples.wdl" as p_update_s
import "patient.out.wdl" as p_out
import "workflow_arguments.wdl" as wfargs
import "workflow_resources.wdl" as wfres
import "runtime_collection.wdl" as rtc

import "collect_callable_loci.wdl" as ccl
import "collect_read_counts.wdl" as crc
import "collect_allelic_counts.wdl" as cac
import "harmonize_samples.wdl" as hs
import "calculate_contamination.wdl" as cc
import "genotype_variants.wdl" as gv
import "model_segments.wdl" as ms
import "call_variants.wdl" as cv
import "filter_variants.wdl" as fv
import "annotate_variants.wdl" as av
import "tasks.wdl"
#import "calculate_tumor_mutation_burden.wdl" as tmb
#import "filter_segments.wdl" as fs
import "absolute.wdl" as abs
import "absolute_extract.wdl" as abs_extract
import "phylogicndt.wdl" as phylogicndt


workflow MultiSampleSomaticWorkflow {
    input {
        # This string is used to label the outputs of the workflow.
        String patient_id
        String? sex
        # If defined, all arrays must have the same length. Each entry corresponds to a
        # sequencing run, and the same index in each array corresponds to the same run.
        # Several sequencing runs from the same physical sample are allowed and will be
        # grouped based on the sample name. If the sample name is not provided, it will
        # be inferred from the bam.
        Array[String]? sample_names
        # Ordinal timepoints of sample collection for phylogenetic inference.
        Array[Int]? timepoints
        Array[File]+ bams
        Array[File]+ bais
        # For targeted sequencing, the (possibly padded and ideally blacklist-removed)
        # target intervals must be supplied. For whole genome sequencing, the intervals
        # are just the chromosomal intervals (ideally blacklist-removed).
        Array[File]+ target_intervals
        # The target_intervals annotated with gc content, mappability, and segmental duplications.
        Array[File]? annotated_target_intervals
        # If a panel of normals is not available for the sequencing platform of a sample,
        # its corresponding path must point to an empty file (of size 0B). The
        # annotated_target_intervals will instead be used for denoising.
        Array[File]? cnv_panel_of_normals
        # Setting this avoids double counting evidence from paired-end reads. This is
        # particularly important for cell-free DNA samples, where the majority of
        # templates is shorter than twice the read length.
        Array[Boolean]? is_paired_end
        # Whether to use the sample for total copy ratio (tCR) and/or allelic copy ratio
        # (aCR) estimation. If not provided, all samples will be used.
        Array[Boolean]? use_sample_for_tCR  # Boolean inputs, ensure correct format!
        Array[Boolean]? use_sample_for_aCR  # Boolean inputs, ensure correct format!

        # A list of normal sample names. If not provided, all samples will be treated as
        # tumor samples.
        Array[String]? normal_sample_names

        Patient? input_patient
        WorkflowArguments? input_args
        WorkflowResources? input_resources
        RuntimeCollection? input_runtime_collection
    }


###############################################################################
#                                                                             #
#                              PREPROCESSING                                  #
#                                                                             #
###############################################################################

    if (!defined(input_runtime_collection)) {
        call rtc.DefineRuntimeCollection as RuntimeParameters {
            input:
                num_bams = length(bams),
        }
    }
    RuntimeCollection runtime_collection = select_first([input_runtime_collection, RuntimeParameters.rtc])

    if (!defined(input_resources)) {
        call wfres.DefineWorkflowResources as Files
    }
    WorkflowResources resources = select_first([input_resources, Files.resources])

    if (!defined(input_args)) {
        call wfargs.DefineWorkflowArguments as Parameters {
            input:
                resources = resources,
                runtime_collection = runtime_collection,
        }
    }
    WorkflowArguments args = select_first([input_args, Parameters.arguments])

    if (!defined(input_patient)) {
        call p_def.DefinePatient as Cache {
            input:
                name = patient_id,
                sex = sex,
                sample_names = sample_names,
                timepoints = timepoints,
                bams = bams,
                bais = bais,
                target_intervals = target_intervals,
                annotated_target_intervals = annotated_target_intervals,
                cnv_panel_of_normals = cnv_panel_of_normals,
                is_paired_end = is_paired_end,
                use_for_tCR = use_sample_for_tCR,
                use_for_aCR = use_sample_for_aCR,
                normal_sample_names = normal_sample_names,
                scattered_intervals_for_variant_calling = args.files.scattered_intervals_for_variant_calling,
                runtime_collection = runtime_collection,
        }
    }

    # TODO: add parse_input task to check for validity, then add "after parse_input" to all calls

    Patient patient = select_first([input_patient, Cache.patient])


###############################################################################
#                                                                             #
#                             COVERAGE WORKFLOW                               #
#                                                                             #
###############################################################################


    Patient coverage_workflow_patient = patient

    scatter (sample in coverage_workflow_patient.samples) {
        scatter (sequencing_run in sample.sequencing_runs) {
            if (args.run_collect_callable_loci && (size(sequencing_run.callable_loci) == 0)) {
                call ccl.CollectCallableLoci {
                    input:
                        ref_fasta = args.files.ref_fasta,
                        ref_fasta_index = args.files.ref_fasta_index,
                        ref_dict = args.files.ref_dict,
                        sample_name = sample.name,
                        bam = sequencing_run.bam,
                        bai = sequencing_run.bai,
                        is_paired_end = sequencing_run.is_paired_end,
                        runtime_collection = runtime_collection,
                }
            }

            if (args.run_collect_total_read_counts && (size(sequencing_run.total_read_counts) == 0) && (size(sequencing_run.denoised_total_copy_ratios) == 0)) {
                call crc.CollectReadCounts {
                    input:
                        ref_fasta = args.files.ref_fasta,
                        ref_fasta_index = args.files.ref_fasta_index,
                        ref_dict = args.files.ref_dict,
                        sample_name = sample.name,
                        bam = sequencing_run.bam,
                        bai = sequencing_run.bai,
                        interval_list = sequencing_run.target_intervals,
                        annotated_interval_list = sequencing_run.annotated_target_intervals,
                        read_count_panel_of_normals = sequencing_run.cnv_panel_of_normals,
                        is_paired_end = sequencing_run.is_paired_end,
                        max_soft_clipped_bases = args.collect_read_counts_max_soft_clipped_bases,
                        runtime_collection = runtime_collection,
                }
            }

            if (args.run_collect_allelic_read_counts && (size(sequencing_run.snppanel_allelic_pileup_summaries) == 0)) {
                call cac.CollectAllelicCounts {
                    input:
                        ref_dict = args.files.ref_dict,
                        bam = sequencing_run.bam,
                        bai = sequencing_run.bai,
                        is_paired_end = sequencing_run.is_paired_end,
                        sample_name = sample.name + ".snppanel",
                        interval_list = sequencing_run.target_intervals,
                        scattered_interval_list = args.files.scattered_intervals_for_pileups,
                        variants = args.files.common_germline_alleles,
                        variants_idx = args.files.common_germline_alleles_idx,
                        getpileupsummaries_extra_args = args.getpileupsummaries_extra_args,
                        minimum_population_allele_frequency = args.min_snppanel_pop_af,
                        maximum_population_allele_frequency = args.max_snppanel_pop_af,
                        minimum_read_depth = args.min_snppanel_read_depth,
                        padding = args.het_to_interval_mapping_max_distance,
                        runtime_collection = runtime_collection,
                }
            }

            call seqrun.UpdateSequencingRun as SeqAddCoverage {
                input:
                    sequencing_run = sequencing_run,
                    callable_loci = CollectCallableLoci.bed,
                    total_read_counts = CollectReadCounts.read_counts,
                    denoised_total_copy_ratios = CollectReadCounts.denoised_copy_ratios,
                    snppanel_allelic_pileup_summaries = CollectAllelicCounts.pileup_summaries,
            }
        }
    }

    call p_update_s.UpdateSamples as PatientAddCoverage {
        input:
            patient = coverage_workflow_patient,
            sequencing_runs = SeqAddCoverage.updated_sequencing_run,
    }

    # todo: FilterIntervals

    call hs.HarmonizeSamples {
        input:
            ref_dict = args.files.ref_dict,
            harmonize_copy_ratios_script = args.script_harmonize_copy_ratios,
            merge_pileups_script = args.script_merge_pileups,
            samples = PatientAddCoverage.updated_patient.samples,
            harmonize_min_target_length = args.harmonize_min_target_length,
            pileups_min_read_depth = args.min_snppanel_read_depth,
            compress_output = false,
            runtime_collection = runtime_collection,
    }

    call p_update_s.UpdateSamples as ConsensusPatient {
        input:
            patient = PatientAddCoverage.updated_patient,
            harmonized_callable_loci = HarmonizeSamples.harmonized_callable_loci,
            harmonized_denoised_total_copy_ratios = HarmonizeSamples.harmonized_denoised_copy_ratios,
            harmonized_snppanel_allelic_pileup_summaries = HarmonizeSamples.merged_allelic_counts,
            allelic_pileup_summaries = HarmonizeSamples.merged_allelic_counts,  # Will be overwritten later
    }

    if (defined(HarmonizeSamples.merged_allelic_counts) && args.run_contamination_model) {
        scatter (sample in ConsensusPatient.updated_patient.samples) {
            if (size(sample.contamination_table) == 0) {
                if (sample.is_tumor && defined(ConsensusPatient.updated_patient.matched_normal_sample)) {
                    Sample cov_matched_normal_sample = select_first([ConsensusPatient.updated_patient.matched_normal_sample])
                    File? matched_normal_pileups = cov_matched_normal_sample.harmonized_snppanel_allelic_pileup_summaries
                }

                call cc.CalculateContamination {
                    input:
                        tumor_pileups = sample.harmonized_snppanel_allelic_pileup_summaries,
                        normal_pileups = matched_normal_pileups,
                        runtime_collection = runtime_collection,
                }
            }

            File contam_table = select_first([CalculateContamination.contamination_table, sample.contamination_table])
        }

        call p_update_s.UpdateSamples as AddContaminationToSamples {
            input:
                patient = ConsensusPatient.updated_patient,
                contamination_table = contam_table,
        }

        # Perform a first-pass single-sample segmentation to get prior allelic
        # copy ratio segmentations for genotyping.
        call ms.ModelSegments as FirstPassSegmentation {
            input:
                patient = AddContaminationToSamples.updated_patient,
                args = args,
                runtime_collection = runtime_collection,
                pre_select_hets = false,
                gvcf = args.files.common_germline_alleles,
                gvcf_idx = args.files.common_germline_alleles_idx,
        }
    }

    Patient coverage_workflow_updated_patient = select_first([FirstPassSegmentation.updated_patient, ConsensusPatient.updated_patient])


###############################################################################
#                                                                             #
#                               SNV WORKFLOW                                  #
#                                                                             #
###############################################################################


    Patient snv_patient = coverage_workflow_updated_patient

    if (args.run_variant_calling) {
        call cv.CallVariants {
            input:
                patient = patient,
                args = args,
                runtime_collection = runtime_collection,
        }

        # So that the coverage workflow can run in parallel to calling SNVs
        call p_merge.MergePatients as AddSNVCallsToPatient {
            input:
                patient = snv_patient,
                other = CallVariants.updated_patient
        }
    }

    if (args.run_variant_filter) {
        call fv.FilterVariants {
            input:
                patient = select_first([AddSNVCallsToPatient.updated_patient, snv_patient]),
                args = args,
                runtime_collection = runtime_collection,
        }

        if (args.keep_germline && defined(FilterVariants.updated_patient.germline_vcf)) {
            # Collect allelic pileups for all putative germline sites that were
            # not yet collected via the coverage workflow, then merge them.
            # This allows for more sensitive aCR segmentation.
            # Only collect SNPs since Indels or MNVs are too likely misclassified.
            call tasks.SelectVariants as SelectGermlineNotInResource {
                input:
                    ref_fasta = args.files.ref_fasta,
                    ref_fasta_index = args.files.ref_fasta_index,
                    ref_dict = args.files.ref_dict,
                    vcf = select_first([FilterVariants.updated_patient.germline_vcf]),
                    vcf_idx = select_first([FilterVariants.updated_patient.germline_vcf_idx]),
                    interval_blacklist = args.files.common_germline_alleles,
                    interval_blacklist_idx = args.files.common_germline_alleles_idx,
                    compress_output = args.compress_output,
                    select_variants_extra_args = "--select-type-to-include SNP",
                    runtime_params = runtime_collection.select_variants
            }

            if (SelectGermlineNotInResource.num_selected_variants > 0) {
                scatter (sample in FilterVariants.updated_patient.samples) {
                    String bam_names = sample.bam_name
                    String these_sample_names = sample.name
                    File? allelic_pileups = sample.allelic_pileup_summaries
                }
                call cac.VcfToPileupVariants as GermlineVariantsNotInResource {
                    input:
                        vcf = SelectGermlineNotInResource.selected_vcf,
                        vcf_idx = SelectGermlineNotInResource.selected_vcf_idx,
                        sample_names = bam_names,
                        compress_output = args.compress_output,
                        runtime_params = runtime_collection.vcf_to_pileup_variants,
                }

                call hs.MergeAllelicCounts as MergeGermlineAllelicCounts {
                    input:
                        ref_dict = args.files.ref_dict,
                        script = args.script_merge_pileups,
                        sample_names = flatten([these_sample_names, these_sample_names]),
                        allelic_counts = select_all(flatten([allelic_pileups, GermlineVariantsNotInResource.pileups])),
                        compress_output = args.compress_output,
                        runtime_params = runtime_collection.merge_allelic_counts,
                }

                # sort output to match order of sample_names since glob doesn't guarantee order
                scatter (sample in FilterVariants.updated_patient.samples) {
                    scatter (allelic_count in MergeGermlineAllelicCounts.merged_allelic_counts) {
                        String this_sample_name = basename(basename(allelic_count, ".gz"), ".pileup")
                        if (sample.name == this_sample_name) {
                            File this_allelic_counts = allelic_count
                        }
                    }
                    Array[File] this_sample_allelic_counts = select_all(this_allelic_counts)
                }
                Array[File] sorted_allelic_counts = flatten(this_sample_allelic_counts)

                call p_update_s.UpdateSamples as ExtendAllelicPileups {
                    input:
                        patient = FilterVariants.updated_patient,
                        allelic_pileup_summaries = sorted_allelic_counts,
                }

                call p.UpdatePatient as AddGermlineAlleles {
                    input:
                        patient = ExtendAllelicPileups.updated_patient,
                        rare_germline_alleles = GermlineVariantsNotInResource.variants,
                        rare_germline_alleles_idx = GermlineVariantsNotInResource.variants_idx
                }
            }
        }

        Patient filtered_snv_patient = select_first([AddGermlineAlleles.updated_patient, FilterVariants.updated_patient])

        if (args.run_variant_annotation) {
            # The sample scatter needs to be outside of the call to AnnotateVariants
            # since cromwell shits the bed for piping optional inputs into a nested scatter.
            scatter (sample in filtered_snv_patient.samples) {
                if (size(sample.annotated_somatic_variants) == 0) {
                    if (sample.is_tumor && defined(filtered_snv_patient.matched_normal_sample)) {
                        Sample cnv_matched_normal_sample = select_first([filtered_snv_patient.matched_normal_sample])
                        String? matched_normal_sample_name = cnv_matched_normal_sample.name
                        String? matched_normal_bam_name = cnv_matched_normal_sample.bam_name
                    }

                    call av.AnnotateVariants {
                        input:
                            vcf = select_first([filtered_snv_patient.somatic_vcf]),
                            vcf_idx = select_first([filtered_snv_patient.somatic_vcf_idx]),
                            num_variants = filtered_snv_patient.num_somatic_variants,
                            individual_id = filtered_snv_patient.name,
                            tumor_sample_name = sample.name,
                            tumor_bam_name = sample.bam_name,
                            normal_sample_name = matched_normal_sample_name,
                            normal_bam_name = matched_normal_bam_name,
                            args = args,
                            runtime_collection = runtime_collection,
                    }
                }
                File annot_som_var = select_first([AnnotateVariants.annotated_variants, sample.annotated_somatic_variants])
                File? annot_som_var_idx = if defined(AnnotateVariants.annotated_variants_idx) then AnnotateVariants.annotated_variants_idx else sample.annotated_somatic_variants_idx
            }
            if (length(select_all(annot_som_var)) > 0) {
                Array[File] annotated_variants_idx = select_all(annot_som_var_idx)
            }

            ## Postpone to below, so CNV workflow can run in parallel to Funcotator.
            #call p_update_s.UpdateSamples as AddAnnotatedVariantsToSamples {
            #    input:
            #        patient = filtered_snv_patient,
            #        annotated_somatic_variants = annot_som_var,
            #        annotated_somatic_variants_idx = annotated_variants_idx,
            #}
        }
    }

    Patient snv_upated_patient = select_first([
        # AddAnnotatedVariantsToSamples.updated_patient,
        filtered_snv_patient,
        AddSNVCallsToPatient.updated_patient,
        snv_patient
    ])


###############################################################################
#                                                                             #
#                               CNV WORKFLOW                                  #
#                                                                             #
###############################################################################


    Patient cnv_patient = snv_upated_patient

    # ModelSegments requires the allelic counts to be pulled down at the same
    # set of loci for all samples. GetPileupSummaries does not guarantee this,
    # however, GenotypeVariants enforces this. Estimating the contamination is
    # helpful for genotyping.

    scatter (sample in cnv_patient.samples) {
        String gt_sample_names = sample.name
        File? pileups = sample.allelic_pileup_summaries
        File? contaminations = sample.contamination_table
        File? af_segmentations = select_first([sample.called_copy_ratio_segmentation, sample.af_segmentation_table])
        File? af_model_params = sample.af_model_parameters
    }
    Array[File] gt_pileups = select_all(pileups)
    if (length(select_all(contaminations)) > 0) {
        Array[File] contamination_tables = select_all(contaminations)
    }
    if (length(select_all(af_segmentations)) > 0) {
        Array[File] segmentation_tables = select_all(af_segmentations)
    }
    if (length(select_all(af_model_params)) > 0) {
        Array[File] af_pre_model_parameters = select_all(af_model_params)
    }

    if (cnv_patient.has_normal) {
        scatter (normal_sample in cnv_patient.normal_samples) {
            String? n_sample_name = normal_sample.name
        }
        Array[String] gt_normal_sample_names = select_all(n_sample_name)
    }

    if (length(gt_pileups) > 0) {
        call gv.GenotypeVariants {
            input:
                script = args.script_genotype_variants,
                patient_id = cnv_patient.name,
                sex = cnv_patient.sex,
                sample_names = gt_sample_names,
                normal_sample_names = gt_normal_sample_names,
                pileups = gt_pileups,
                contamination_tables = contamination_tables,
                segmentation_tables = segmentation_tables,
                af_model_parameters = af_pre_model_parameters,
                common_germline_alleles = args.files.common_germline_alleles,
                common_germline_alleles_idx = args.files.common_germline_alleles_idx,
                rare_germline_alleles = cnv_patient.rare_germline_alleles,
                rare_germline_alleles_idx = cnv_patient.rare_germline_alleles_idx,
                compress_output = args.compress_output,
                min_read_depth = args.min_snppanel_read_depth,
                normal_to_tumor_weight = args.genotype_variants_normal_to_tumor_weight,
                min_genotype_likelihood = args.genotype_variants_min_genotype_likelihood,
                outlier_prior = args.genotype_variants_outlier_prior,
                overdispersion = args.genotype_variants_overdispersion,
                ref_bias = args.genotype_variants_ref_bias,
                select_hets = false,
                save_sample_genotype_likelihoods = true,
                runtime_collection = runtime_collection,
        }

        # todo: phase gvcf

        call p_update_s.UpdateSamples as AddPileupsToSamples {
            input:
                patient = cnv_patient,
                allelic_pileup_summaries = GenotypeVariants.sample_genotype_likelihoods,  # Careful: This is not technically in a pileup format!
        }

        call p.UpdatePatient as AddGVCFtoPatient {
            input:
                patient = AddPileupsToSamples.updated_patient,
                gvcf = GenotypeVariants.vcf,
                gvcf_idx = GenotypeVariants.vcf_idx,
                snp_ref_counts = GenotypeVariants.ref_counts,
                snp_alt_counts = GenotypeVariants.alt_counts,
                snp_other_alt_counts = GenotypeVariants.other_alt_counts,
                snp_sample_correlation = GenotypeVariants.sample_correlation,
                snp_sample_correlation_min = GenotypeVariants.sample_correlation_min
        }
    }

    if (args.run_model_segments) {
        call ms.ModelSegments {
            input:
                patient = select_first([AddGVCFtoPatient.updated_patient, cnv_patient]),
                args = args,
                runtime_collection = runtime_collection,
        }

#        if (args.run_filter_segments) {
#            call fs.FilterSegments {
#                input:
#                    patient = ModelSegments.updated_patient,
#                    args = args,
#                    runtime_collection = runtime_collection,
#            }
#        }
    }

    # todo: FuncotateSegments

    # Only update here so Funcotator can run in parallel to CNV workflow.
    call p_update_s.UpdateSamples as AddAnnotatedVariantsToSamples {
        input:
            patient = select_first([ModelSegments.updated_patient, AddGVCFtoPatient.updated_patient, cnv_patient]),
            annotated_somatic_variants = annot_som_var,
            annotated_somatic_variants_idx = annotated_variants_idx,
    }

    Patient cnv_updated_patient = AddAnnotatedVariantsToSamples.updated_patient


###############################################################################
#                                                                             #
#                             CLONAL WORKFLOW                                 #
#                                                                             #
###############################################################################


    Patient clonal_patient = cnv_updated_patient

    if (args.run_clonal_decomposition) {
        scatter (sample in clonal_patient.samples) {
            if (defined(sample.called_copy_ratio_segmentation) && defined(sample.af_model_parameters) && !defined(sample.absolute_acr_rdata) && !defined(sample.absolute_acr_plot)) {
                call abs.Absolute {
                    input:
                        acs_conversion_script = args.script_acs_conversion,
                        sample_name = sample.name,
                        copy_ratio_segmentation = select_first([sample.called_copy_ratio_segmentation]),
                        af_model_parameters = select_first([sample.af_model_parameters]),
                        annotated_variants = sample.annotated_somatic_variants,
                        purity = sample.purity,
                        ploidy = sample.ploidy,
                        sex = clonal_patient.sex,
                        min_hets = args.absolute_min_hets,
                        min_probes = args.absolute_min_probes,
                        maf90_threshold = args.absolute_maf90_threshold,
                        genome_build = args.absolute_genome_build,
                        runtime_collection = runtime_collection
                }
            }

            File acs_cr_segmentation = select_first([Absolute.acs_copy_ratio_segmentation, sample.acs_copy_ratio_segmentation])
            Float acs_cr_skew = select_first([Absolute.acs_copy_ratio_skew, sample.acs_copy_ratio_skew])
            File? snv_maf = if defined(Absolute.snv_maf) then Absolute.snv_maf else sample.absolute_snv_maf
            File? indel_maf = if defined(Absolute.indel_maf) then Absolute.indel_maf else sample.absolute_indel_maf
            File acr_rdata = select_first([Absolute.acr_rdata, sample.absolute_acr_rdata])
            File acr_plot = select_first([Absolute.acr_plot, sample.absolute_acr_plot])

            if (defined(sample.absolute_solution)) {
                call abs_extract.AbsoluteExtract {
                    input:
                        map_to_absolute_copy_number_script = args.script_map_to_absolute_copy_number,
                        sample_name = sample.name,
                        sex = clonal_patient.sex,
                        rdata = acr_rdata,
                        called_solution = select_first([sample.absolute_solution]),
                        analyst_id = args.analyst_id,
                        copy_ratio_type = "allelic",
                        acs_copy_ratio_segmentation = acs_cr_segmentation,
                        acs_copy_ratio_skew = acs_cr_skew,
                        snv_maf = snv_maf,
                        indel_maf = indel_maf,
                        gvcf = clonal_patient.gvcf,
                        genome_build = args.absolute_genome_build,
                        runtime_collection = runtime_collection
                }
            }
        }

        if (length(select_all(snv_maf)) > 0) {
            Array[File] abs_snv_maf = select_all(snv_maf)
        }
        if (length(select_all(indel_maf)) > 0) {
            Array[File] abs_indel_maf = select_all(indel_maf)
        }
        if (length(select_all(AbsoluteExtract.absolute_maf)) > 0) {
            Array[File] abs_maf = select_all(AbsoluteExtract.absolute_maf)
        }
        if (length(select_all(AbsoluteExtract.absolute_segtab)) > 0) {
            Array[File] abs_segtab = select_all(AbsoluteExtract.absolute_segtab)
        }
        if (length(select_all(AbsoluteExtract.absolute_table)) > 0) {
            Array[File] abs_table = select_all(AbsoluteExtract.absolute_table)
        }
        if (length(select_all(AbsoluteExtract.absolute_purity)) > 0) {
            Array[Float] abs_purity = select_all(AbsoluteExtract.absolute_purity)
        }
        if (length(select_all(AbsoluteExtract.absolute_ploidy)) > 0) {
            Array[Float] abs_ploidy = select_all(AbsoluteExtract.absolute_ploidy)
        }

        call p_update_s.UpdateSamples as AddAbsoluteResultsToSamples {
            input:
                patient = clonal_patient,
                acs_copy_ratio_segmentation = acs_cr_segmentation,
                acs_copy_ratio_skew = acs_cr_skew,
                absolute_acr_rdata = acr_rdata,
                absolute_acr_plot = acr_plot,
                absolute_snv_maf = abs_snv_maf,
                absolute_indel_maf = abs_indel_maf,
                absolute_maf = abs_maf,
                absolute_segtab = abs_segtab,
                absolute_table = abs_table,
                purity = abs_purity,
                ploidy = abs_ploidy,
        }
    }
    
    Patient clonal_updated_patient = select_first([AddAbsoluteResultsToSamples.updated_patient, clonal_patient])

    # todo: add phylogicNDT
    if (args.run_phylogicndt) {
        scatter (t_sample in clonal_updated_patient.tumor_samples){
            File? t_sample_absolute_maf = t_sample.absolute_maf
            File? t_sample_absolute_segtab = t_sample.absolute_segtab
            Float? t_sample_absolute_purity = t_sample.purity
            Int? t_sample_absolute_timepoint = t_sample.absolute_timepoint
        }

        Array[File] phylogic_absolute_maf = select_all(t_sample_absolute_maf)
        Array[File]? phylogic_absolute_segtab = select_all(t_sample_absolute_segtab)
        Array[Float] phylogic_absolute_purity = select_all(t_sample_absolute_purity)
        Array[Int]? phylogic_absolute_timepoint = select_all(t_sample_absolute_timepoint)

        call phylogicndt.PhylogicNDT {
            input:
                patient_id = clonal_updated_patient.name,
                absolute_maf = phylogic_absolute_maf,
                absolute_segtab = phylogic_absolute_segtab,
                absolute_purity = phylogic_absolute_purity,
                timepoints = phylogic_absolute_timepoint,
                runtime_collection = runtime_collection
        }

        # if (length(select_all(PhylogicNDT.phylogic_pie_plots)) > 0) {
        #     Array[File] phylogic_pie_plots = select_all(PhylogicNDT.phylogic_pie_plots)
        # }
        # if (length(select_all(PhylogicNDT.phylogic_mutation_plots)) > 0) {
        #     Array[File] phylogic_mutation_plots = select_all(PhylogicNDT.phylogic_mutation_plots)
        # }
        # if (length(select_all(PhylogicNDT.phylogic_cluster_plots)) > 0) {
        #     Array[File] phylogic_cluster_plots = select_all(PhylogicNDT.phylogic_cluster_plots)
        # }

        # File phylogic_cnvs = if defined(PhylogicNDT.phylogic_cnvs) then PhylogicNDT.phylogic_cnvs else clonal_updated_patient.phylogic_cnvs

    }

###############################################################################
#                                                                             #
#                                  OUTPUT                                     #
#                                                                             #
###############################################################################

    Patient out_patient = clonal_updated_patient

    call p_out.Output {
        input:
            patient = out_patient
    }

    output {
        # for each sequencing run:
        # CACHE (as returned by the workflow)
        Array[Array[File]]? callable_loci = Output.callable_loci
        Array[Array[File]]? total_read_counts = Output.total_read_counts
        Array[Array[File]]? denoised_total_copy_ratios = Output.denoised_total_copy_ratios
        Array[Array[File]]? snppanel_allelic_pileup_summaries = Output.snppanel_allelic_pileup_summaries

        # for each sample:
        # CACHE (as returned by the workflow)
        Array[String]? sample_names_ordered = Output.sample_names_ordered
        Array[File]? harmonized_callable_loci = Output.harmonized_callable_loci
        Array[File]? harmonized_denoised_total_copy_ratios = Output.harmonized_denoised_total_copy_ratios
        Array[File]? harmonized_snppanel_allelic_pileup_summaries = Output.harmonized_snppanel_allelic_pileup_summaries
        Array[File]? contamination_table = Output.contamination_table
        Array[File]? af_segmentation_table = Output.af_segmentation_table
        Array[File]? allelic_pileup_summaries = Output.allelic_pileup_summaries
        Array[File]? aggregated_allelic_read_counts = Output.aggregated_allelic_read_counts
        Array[Float]? genotype_error_probabilities = Output.genotype_error_probabilities
        Array[File]? af_model_parameters = Output.af_model_parameters
        Array[File]? cr_model_parameters = Output.cr_model_parameters
        Array[File]? called_copy_ratio_segmentation = Output.called_copy_ratio_segmentation
        Array[File]? cr_plot = Output.cr_plot
        Array[File]? acs_copy_ratio_segmentation = Output.acs_copy_ratio_segmentation
        Array[Float]? acs_copy_ratio_skew = Output.acs_copy_ratio_skew
        Array[File]? annotated_somatic_variants = Output.annotated_somatic_variants
        Array[File?]? annotated_somatic_variants_idx = Output.annotated_somatic_variants_idx
        Array[File]? absolute_acr_rdata = Output.absolute_acr_rdata
        Array[File]? absolute_acr_plot = Output.absolute_acr_plot
        Array[File]? absolute_snv_maf = Output.absolute_snv_maf
        Array[File]? absolute_indel_maf = Output.absolute_indel_maf
        Array[Int]? absolute_solution = Output.absolute_solution
        Array[Int]? absolute_timepoint = Output.absolute_timepoint
        Array[File]? absolute_maf = Output.absolute_maf
        Array[File]? absolute_segtab = Output.absolute_segtab
        Array[File]? absolute_table = Output.absolute_table
        Array[Float]? purity = Output.purity
        Array[Float]? ploidy = Output.ploidy

        # Array[File]? phylogic_pie_plots = PhylogicNDT.phylogic_pie_plots
        # Array[File]? phylogic_mutation_plots = PhylogicNDT.phylogic_mutation_plots
        # Array[File]? phylogic_cluster_plots = PhylogicNDT.phylogic_cluster_plots
        # File? phylogic_cnvs = PhylogicNDT.phylogic_cnvs
        # File? phylogic_mut_ccfs = PhylogicNDT.phylogic_mut_ccfs
        # File? phylogic_unclustered = PhylogicNDT.phylogic_unclustered
        # File? phylogic_cluster_ccfs = PhylogicNDT.phylogic_cluster_ccfs
        # File? phylogic_report = PhylogicNDT.phylogic_report
        # File? phylogic_cell_population_abundances = PhylogicNDT.phylogic_cell_population_abundances
        # File? phylogic_cell_population_mcmc_trace = PhylogicNDT.phylogic_cell_population_mcmc_trace
        # File? phylogic_constrained_ccf = PhylogicNDT.phylogic_constrained_ccf
        # File? phylogic_build_tree_posteriors = PhylogicNDT.phylogic_build_tree_posteriors

        Array[File]? first_pass_cr_segmentations = FirstPassSegmentation.called_copy_ratio_segmentations
        Array[File]? first_pass_cr_plots = FirstPassSegmentation.cr_plots
        Array[File]? first_pass_af_model_parameters = FirstPassSegmentation.af_model_final_parameters
        Array[File]? first_pass_cr_model_parameters = FirstPassSegmentation.cr_model_final_parameters

        # for each interval shard:
        # CACHE (as returned by the workflow)
        Array[File]? raw_calls_mutect2_vcf_scattered = Output.raw_calls_mutect2_vcf_scattered
        Array[File]? raw_calls_mutect2_vcf_idx_scattered = Output.raw_calls_mutect2_vcf_idx_scattered
        Array[File]? raw_mutect2_stats_scattered = Output.raw_mutect2_stats_scattered
        Array[File]? raw_mutect2_bam_out_scattered = Output.raw_mutect2_bam_out_scattered
        Array[File]? raw_mutect2_bai_out_scattered = Output.raw_mutect2_bai_out_scattered
        Array[File]? raw_mutect2_artifact_priors_scattered = Output.raw_mutect2_artifact_priors_scattered

        # for patient
        File? raw_snv_calls_vcf = out_patient.raw_snv_calls_vcf
        File? raw_snv_calls_vcf_idx = out_patient.raw_snv_calls_vcf_idx
        File? mutect2_stats = out_patient.mutect2_stats
        File? orientation_bias = out_patient.orientation_bias
        File? filtered_vcf = out_patient.filtered_vcf
        File? filtered_vcf_idx = out_patient.filtered_vcf_idx
        File? filtering_stats = out_patient.filtering_stats
        File? somatic_vcf = out_patient.somatic_vcf
        File? somatic_vcf_idx = out_patient.somatic_vcf_idx
        Int? num_somatic_variants = out_patient.num_somatic_variants
        File? germline_vcf = out_patient.germline_vcf
        File? germline_vcf_idx = out_patient.germline_vcf_idx
        File? rare_germline_alleles = out_patient.rare_germline_alleles
        File? rare_germline_alleles_idx = out_patient.rare_germline_alleles_idx
        File? somatic_calls_bam = out_patient.somatic_calls_bam
        File? somatic_calls_bai = out_patient.somatic_calls_bai
        File? gvcf = out_patient.gvcf
        File? gvcf_idx = out_patient.gvcf_idx
        File? snp_ref_counts = out_patient.snp_ref_counts
        File? snp_alt_counts = out_patient.snp_alt_counts
        File? snp_other_alt_counts = out_patient.snp_other_alt_counts
        File? snp_sample_correlation = out_patient.snp_sample_correlation
        Float? snp_sample_correlation_min = out_patient.snp_sample_correlation_min
        File? modeled_segments = out_patient.modeled_segments
    }
}
