version development

# Compute half of the pipeline (CNV, SNV, copy-ratio modeling, clonal analysis,
# output). It takes a fully assembled patient and the resolved argument/resource/
# runtime structs and contains no preprocessing or input construction of its own;
# MultiSampleSomaticWorkflow performs all of that and calls this workflow. Keeping
# this split lets the entry point present a clean input surface on Terra.

import "sequencing_run.wdl" as seqrun
import "sample.wdl" as s
import "patient.wdl" as p
import "patient.merge.wdl" as p_merge
import "patient.update_samples.wdl" as p_update_s
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
import "filter_segments.wdl" as fs
import "absolute.wdl" as abs
import "absolute_extract.wdl" as abs_extract
import "phylogicndt.wdl" as phylogicndt


workflow MultiSampleSomaticWorkflowRun {
    input {
        Patient patient
        WorkflowArguments args
        WorkflowResources resources
        RuntimeCollection runtime_collection
    }


    scatter (sample in patient.samples) {
        File? cached_called_copy_ratio_segmentation = sample.called_copy_ratio_segmentation
        File? cached_af_model_parameters = sample.af_model_parameters
        File? cached_absolute_rdata = sample.absolute_acr_rdata
    }
    Boolean skip_to_clonal_decomposition = (
        (
            (length(select_all(cached_called_copy_ratio_segmentation)) > 0)
            && (length(select_all(cached_af_model_parameters)) > 0)
        ) || (length(select_all(cached_absolute_rdata)) > 0)
    )

    if (!skip_to_clonal_decomposition) {


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
                        sex_genotype = coverage_workflow_patient.sex,
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

    # todo: infer sex genotype

    call p_update_s.UpdateSamples as PatientAddCoverage {
        input:
            patient = coverage_workflow_patient,
            sequencing_runs = SeqAddCoverage.updated_sequencing_run,
    }

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

        # Perform a first-pass single-sample segmentation to detect segments that
        # need to be filtered.
        call ms.ModelSegments as PreFirstPassSegmentation {
            input:
                patient = AddContaminationToSamples.updated_patient,
                args = args,
                runtime_collection = runtime_collection,
                pre_select_hets = false,
                gvcf = args.files.common_germline_alleles,
                gvcf_idx = args.files.common_germline_alleles_idx,
        }

        if (args.filter_segments_min_probes > 1) {
            # Remove copy ratio observations that are outliers which induce artificial
            # segmentation boundaries.
            call fs.FilterSegments {
                input:
                    patient = PreFirstPassSegmentation.updated_patient,
                    args = args,
                    runtime_collection = runtime_collection,
            }

            # Perform a second-pass single-sample segmentation to get prior allelic
            # copy ratio segmentations for genotyping.
            call ms.ModelSegments as FirstPassSegmentation {
                input:
                    patient = FilterSegments.updated_patient,
                    args = args,
                    runtime_collection = runtime_collection,
                    pre_select_hets = false,
                    gvcf = args.files.common_germline_alleles,
                    gvcf_idx = args.files.common_germline_alleles_idx,
            }
        }
    }

    Patient coverage_workflow_updated_patient = select_first([FirstPassSegmentation.updated_patient, PreFirstPassSegmentation.updated_patient, ConsensusPatient.updated_patient])


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
                if (defined(filtered_snv_patient.somatic_vcf) && (size(sample.annotated_somatic_variants) == 0)) {
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
        if (length(select_all([sample.called_copy_ratio_segmentation, sample.af_segmentation_table])) > 1) {
            File? af_segmentations = select_first([sample.called_copy_ratio_segmentation, sample.af_segmentation_table])
        }
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
                phasing_log_ratio_cap = args.genotype_variants_phasing_log_ratio_cap,
                phasing_sample_llr_threshold = args.genotype_variants_phasing_sample_llr_threshold,
                phasing_consensus_fdr = args.genotype_variants_phasing_consensus_fdr,
                phasing_max_num_contig_segs = args.genotype_variants_phasing_max_num_contig_segs,
                select_hets = false,
                save_sample_genotype_likelihoods = true,
                runtime_collection = runtime_collection,
        }

        if ((args.genome_build == "hg38" || args.genome_build == "hg19")) {
            call tasks.RunPeddy {
                input:
                    patient_id = cnv_patient.name,
                    gvcf = GenotypeVariants.vcf,
                    gvcf_idx = GenotypeVariants.vcf_idx,
                    genome_build = args.genome_build,
                    runtime_params = runtime_collection.ancestry
            }
        }

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
                snp_sample_correlation_min = GenotypeVariants.sample_correlation_min,
                ancestry_pred = RunPeddy.ancestry_pred,
                ancestry_prob = RunPeddy.ancestry_prob,
                ancestry_background_pca_table = RunPeddy.ancestry_background_pca_table,
                ancestry_pca_plot = RunPeddy.ancestry_pca_plot
        }
    }

    if (args.run_model_segments) {
        call ms.ModelSegments {
            input:
                patient = select_first([AddGVCFtoPatient.updated_patient, cnv_patient]),
                args = args,
                runtime_collection = runtime_collection,
        }
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


    } # skip_to_clonal_decomposition

    Patient clonal_patient = select_first([cnv_updated_patient, patient])

    if (args.run_clonal_decomposition) {
        scatter (sample in clonal_patient.samples) {
            # Force ABSOLUTE's purity (alpha) and ploidy (tau) only when the user supplied
            # both; a lone value is ignored so ABSOLUTE infers both (alpha and tau are
            # forced together or not at all).
            if (defined(sample.user_purity) && defined(sample.user_ploidy)) {
                Float forced_purity = select_first([sample.user_purity])
                Float forced_ploidy = select_first([sample.user_ploidy])
            }
            if (defined(sample.called_copy_ratio_segmentation) && defined(sample.af_model_parameters) && !defined(sample.absolute_acr_rdata) && !defined(sample.absolute_acr_plot)) {
                call abs.Absolute {
                    input:
                        acs_conversion_script = args.script_acs_conversion,
                        sample_name = sample.name,
                        copy_ratio_segmentation = select_first([sample.called_copy_ratio_segmentation]),
                        af_model_parameters = sample.af_model_parameters,
                        annotated_variants = sample.annotated_somatic_variants,
                        purity = forced_purity,
                        ploidy = forced_ploidy,
                        sex = clonal_patient.sex,
                        min_hets = args.absolute_min_hets,
                        min_probes = args.absolute_min_probes,
                        maf90_threshold = args.absolute_maf90_threshold,
                        genome_build = args.genome_build,
                        runtime_collection = runtime_collection
                }
            }

            File acs_cr_segmentation = select_first([Absolute.acs_copy_ratio_segmentation, sample.acs_copy_ratio_segmentation])
            Float acs_cr_skew = select_first([Absolute.acs_copy_ratio_skew, sample.acs_copy_ratio_skew])
            File? snv_maf = if defined(Absolute.snv_maf) then Absolute.snv_maf else sample.absolute_snv_maf
            File? indel_maf = if defined(Absolute.indel_maf) then Absolute.indel_maf else sample.absolute_indel_maf
            File acr_rdata = select_first([Absolute.acr_rdata, sample.absolute_acr_rdata])
            File acr_plot = select_first([Absolute.acr_plot, sample.absolute_acr_plot])
            File tcr_rdata = select_first([Absolute.tcr_rdata, sample.absolute_tcr_rdata])
            File tcr_plot = select_first([Absolute.tcr_plot, sample.absolute_tcr_plot])

            if (defined(sample.absolute_acr_solution)) {
                call abs_extract.AbsoluteExtract as AbsoluteACRExtract {
                    input:
                        map_to_absolute_copy_number_script = args.script_map_to_absolute_copy_number,
                        calculate_cancer_cell_fraction_script = args.script_calculate_cancer_cell_fraction,
                        sample_name = sample.name,
                        sex = clonal_patient.sex,
                        rdata = acr_rdata,
                        called_solution = select_first([sample.absolute_acr_solution]),
                        analyst_id = args.analyst_id,
                        copy_ratio_type = "allelic",
                        acs_copy_ratio_segmentation = acs_cr_segmentation,
                        acs_copy_ratio_skew = acs_cr_skew,
                        snv_maf = snv_maf,
                        indel_maf = indel_maf,
                        gvcf = clonal_patient.gvcf,
                        genome_build = args.genome_build,
                        runtime_collection = runtime_collection
                }
            }

            if (defined(sample.absolute_tcr_solution)) {
                call abs_extract.AbsoluteExtract as AbsoluteTCRExtract {
                    input:
                        map_to_absolute_copy_number_script = args.script_map_to_absolute_copy_number,
                        calculate_cancer_cell_fraction_script = args.script_calculate_cancer_cell_fraction,
                        sample_name = sample.name,
                        sex = clonal_patient.sex,
                        rdata = tcr_rdata,
                        called_solution = select_first([sample.absolute_tcr_solution]),
                        analyst_id = args.analyst_id,
                        copy_ratio_type = "total",
                        acs_copy_ratio_segmentation = acs_cr_segmentation,
                        acs_copy_ratio_skew = acs_cr_skew,
                        snv_maf = snv_maf,
                        indel_maf = indel_maf,
                        gvcf = clonal_patient.gvcf,
                        genome_build = args.genome_build,
                        runtime_collection = runtime_collection
                }
            }
        }

        # ABSOLUTE input
        if (length(select_all(snv_maf)) > 0) {
            Array[File] abs_snv_maf = select_all(snv_maf)
        }
        if (length(select_all(indel_maf)) > 0) {
            Array[File] abs_indel_maf = select_all(indel_maf)
        }

        # ACR
        if (length(select_all(AbsoluteACRExtract.absolute_maf)) > 0) {
            Array[File] abs_acr_maf = select_all(AbsoluteACRExtract.absolute_maf)
        }
        if (length(select_all(AbsoluteACRExtract.absolute_segtab)) > 0) {
            Array[File] abs_acr_segtab = select_all(AbsoluteACRExtract.absolute_segtab)
        }
        if (length(select_all(AbsoluteACRExtract.absolute_segtab_igv)) > 0) {
            Array[File] abs_acr_segtab_igv = select_all(AbsoluteACRExtract.absolute_segtab_igv)
        }
        if (length(select_all(AbsoluteACRExtract.absolute_table)) > 0) {
            Array[File] abs_acr_table = select_all(AbsoluteACRExtract.absolute_table)
        }
        if (length(select_all(AbsoluteACRExtract.absolute_purity)) > 0) {
            Array[Float] abs_acr_purity = select_all(AbsoluteACRExtract.absolute_purity)
        }
        if (length(select_all(AbsoluteACRExtract.absolute_ploidy)) > 0) {
            Array[Float] abs_acr_ploidy = select_all(AbsoluteACRExtract.absolute_ploidy)
        }

        # TCR
        if (length(select_all(AbsoluteTCRExtract.absolute_maf)) > 0) {
            Array[File] abs_tcr_maf = select_all(AbsoluteTCRExtract.absolute_maf)
        }
        if (length(select_all(AbsoluteTCRExtract.absolute_segtab)) > 0) {
            Array[File] abs_tcr_segtab = select_all(AbsoluteTCRExtract.absolute_segtab)
        }
        if (length(select_all(AbsoluteTCRExtract.absolute_segtab_igv)) > 0) {
            Array[File] abs_tcr_segtab_igv = select_all(AbsoluteTCRExtract.absolute_segtab_igv)
        }
        if (length(select_all(AbsoluteTCRExtract.absolute_table)) > 0) {
            Array[File] abs_tcr_table = select_all(AbsoluteTCRExtract.absolute_table)
        }
        if (length(select_all(AbsoluteTCRExtract.absolute_purity)) > 0) {
            Array[Float] abs_tcr_purity = select_all(AbsoluteTCRExtract.absolute_purity)
        }
        if (length(select_all(AbsoluteTCRExtract.absolute_ploidy)) > 0) {
            Array[Float] abs_tcr_ploidy = select_all(AbsoluteTCRExtract.absolute_ploidy)
        }

        call p_update_s.UpdateSamples as AddAbsoluteResultsToSamples {
            input:
                patient = clonal_patient,
                acs_copy_ratio_segmentation = acs_cr_segmentation,
                acs_copy_ratio_skew = acs_cr_skew,
                absolute_snv_maf = abs_snv_maf,
                absolute_indel_maf = abs_indel_maf,

                absolute_acr_rdata = acr_rdata,
                absolute_acr_plot = acr_plot,
                absolute_acr_maf = abs_acr_maf,
                absolute_acr_segtab = abs_acr_segtab,
                absolute_acr_segtab_igv = abs_acr_segtab_igv,
                absolute_acr_table = abs_acr_table,
                absolute_acr_purity = abs_acr_purity,
                absolute_acr_ploidy = abs_acr_ploidy,

                absolute_tcr_rdata = tcr_rdata,
                absolute_tcr_plot = tcr_plot,
                absolute_tcr_maf = abs_tcr_maf,
                absolute_tcr_segtab = abs_tcr_segtab,
                absolute_tcr_segtab_igv = abs_tcr_segtab_igv,
                absolute_tcr_table = abs_tcr_table,
                absolute_tcr_purity = abs_tcr_purity,
                absolute_tcr_ploidy = abs_tcr_ploidy,
        }

        # Only run PhylogicNDT if there are MAFs with ccf annotation
        if (length(select_first([abs_acr_maf, []])) > 0) {
            scatter (sample in AddAbsoluteResultsToSamples.updated_patient.samples) {
                if (defined(sample.absolute_acr_maf) && defined(sample.absolute_acr_purity) && (sample.absolute_acr_purity > 0)) {
                    String? phylogic_sample_name = sample.name
                    File? sample_absolute_maf = sample.absolute_acr_maf
                    File? sample_absolute_segtab = sample.absolute_acr_segtab
                    Float? sample_purity = sample.absolute_acr_purity
                    Int? sample_timepoint = sample.timepoint
                }
            }
            if (args.phylogic_use_segtab && length(select_all(sample_absolute_segtab)) > 0) {
                Array[File]? phylogic_absolute_segtabs = select_all(sample_absolute_segtab)
            }
            if (length(select_all(sample_timepoint)) > 0) {
                Array[Int]? phylogic_timepoints = select_all(sample_timepoint)
            }

            if (length(select_all(phylogic_sample_name)) > 0) {
                call phylogicndt.PhylogicNDT {
                    input:
                        patient_id = AddAbsoluteResultsToSamples.updated_patient.name,
                        sample_names = select_all(phylogic_sample_name),
                        absolute_mafs = select_all(sample_absolute_maf),
                        absolute_segtabs = phylogic_absolute_segtabs,
                        absolute_purities = select_all(sample_purity),
                        timepoints = phylogic_timepoints,
                        use_indels = args.phylogic_use_indels,
                        impute_missing_snvs = args.phylogic_impute_missing_snvs,
                        min_coverage = args.phylogic_min_coverage,
                        driver_genes_file = args.files.phylogic_driver_genes_file,
                        focal_cnv_intervals = args.files.phylogic_focal_cnv_intervals,
                        genome_build = args.genome_build,
                        runtime_collection = runtime_collection
                }

                call p.UpdatePatient as AddPhylogicToPatient {
                    input:
                        patient = AddAbsoluteResultsToSamples.updated_patient,
                        phylogic_sif_file = PhylogicNDT.sif_file,
                        phylogic_report = PhylogicNDT.report,
                        phylogic_ccfs_cnvs = PhylogicNDT.ccfs_cnvs,
                        phylogic_ccfs_snvs = PhylogicNDT.ccfs_snvs,
                        phylogic_constrained_ccf = PhylogicNDT.constrained_ccf,
                        phylogic_cluster_ccfs = PhylogicNDT.cluster_ccfs,
                        phylogic_build_tree_posteriors = PhylogicNDT.build_tree_posteriors,
                        phylogic_growth_rates = PhylogicNDT.growth_rates,
                        phylogic_growth_rate_plot = PhylogicNDT.growth_rate_plot,
                        phylogic_timing_report = PhylogicNDT.timing_report,
                        phylogic_timing_wgd_supporting_events = PhylogicNDT.timing_wgd_supporting_events,
                        phylogic_timing_graph = PhylogicNDT.timing_graph,
                        phylogic_timing_comparison = PhylogicNDT.timing_comparison,
                        phylogic_timing_table = PhylogicNDT.timing_table,
                }
            }
        }
    }

    # TODO: add mutational signature decomposition
    # TODO: calculate TMB
    
    Patient clonal_updated_patient = select_first([AddPhylogicToPatient.updated_patient, AddAbsoluteResultsToSamples.updated_patient, clonal_patient])


###############################################################################
#                                                                             #
#                                  OUTPUT                                     #
#                                                                             #
###############################################################################


    Patient out_patient = clonal_updated_patient

    output {
        Patient output_patient = out_patient
        WorkflowArguments output_args = args
        WorkflowResources output_resources = resources
        RuntimeCollection output_runtime_collection = runtime_collection
    }
}
