"""Single source of truth for the SequencingRun / Sample / Shard / Patient struct fields.

Edit the field lists here, then run `generate.py` to (re)emit the WDL struct
definitions and the bulk single-pass constructors. The point is that each field is
declared exactly ONCE; the long `object {...}` blocks that are tedious and drift-prone
to maintain by hand (the struct definitions, their per-element Update* constructors,
and patient.out.wdl) become generated output, so a rename can never again leave one
consumer pointing at a field that no longer exists. patient.define.wdl stays hand-
written — it remaps input names and exposes only a curated subset — so adding a field
here that DefinePatient should accept still requires a matching edit there.

Field semantics
---------------
- wdl_type : the field's WDL type, e.g. "File?", "Float?", "Array[SequencingRun]".
- overlaid : True  -> exposed as an optional overlay input `Array[<elem>]?` in the bulk
                      constructor and set per element via
                      `if length(<f>_arr) > 0 then <f>_arr[i] else existing.<f>`.
                      Covers both derived "cache" fields AND structural fields that a
                      module legitimately replaces (e.g. Sample.sequencing_runs, which
                      call_variants swaps for shard-subset runs).
             False -> a fixed identity field, set once at construction and carried over
                      unchanged by the bulk constructor (no overlay input).
- origin   : optional trailing provenance comment.

`overlaid` is the only per-field policy the generator needs.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class F:
    name: str
    wdl_type: str
    overlaid: bool
    origin: str = ""


def fixed(name, wdl_type, origin=""):
    """Identity field: set at construction, carried over by the bulk constructor."""
    return F(name, wdl_type, overlaid=False, origin=origin)


def over(name, wdl_type, origin=""):
    """Overlaid field: exposed as an overlay input and set per element."""
    return F(name, wdl_type, overlaid=True, origin=origin)


SEQUENCING_RUN = [
    fixed("name", "String"),
    over("sample_name", "String"),
    over("timepoint", "Int?", "manual"),
    fixed("bam", "File"),
    fixed("bai", "File"),
    fixed("target_intervals", "File"),
    over("annotated_target_intervals", "File?"),
    over("cnv_panel_of_normals", "File?"),
    over("is_paired_end", "Boolean?"),
    over("use_for_tCR", "Boolean"),
    over("use_for_aCR", "Boolean"),
    over("callable_loci", "File?", "GATK CollectCallableLoci"),
    over("total_read_counts", "File?", "GATK CollectReadCounts"),
    over("denoised_total_copy_ratios", "File?", "GATK DenoiseReadCounts"),
    over("snppanel_allelic_pileup_summaries", "File?", "GATK GetPileupSummaries"),
]

SAMPLE = [
    fixed("name", "String"),
    fixed("bam_name", "String"),
    over("sequencing_runs", "Array[SequencingRun]"),
    fixed("is_tumor", "Boolean"),
    over("harmonized_callable_loci", "File?", "GATK CallableLoci | harmonization"),
    over("harmonized_denoised_total_copy_ratios", "File?", "GATK DenoiseReadCounts | harmonization"),
    over("harmonized_snppanel_allelic_pileup_summaries", "File?", "GATK GetPileupSummaries | harmonization"),
    over("contamination_table", "File?", "GATK CalculateContamination"),
    over("af_segmentation_table", "File?", "GATK CalculateContamination / ModelSegments"),
    over("allelic_pileup_summaries", "File?", "snppanel + rare germline pileups"),
    over("aggregated_allelic_read_counts", "File?", "PileupToAllelicCounts"),
    over("genotype_error_probabilities", "Float?", "PileupToAllelicCounts"),
    over("af_model_parameters", "File?", "GATK ModelSegments"),
    over("cr_model_parameters", "File?", "GATK ModelSegments"),
    over("called_copy_ratio_segmentation", "File?", "GATK CallCopyRatioSegments + merge"),
    over("cr_plot", "File?", "GATK PlotModeledSegments"),
    over("acs_copy_ratio_segmentation", "File?", "ModelSegmentsToACSConversion"),
    over("acs_copy_ratio_skew", "Float?", "ModelSegmentsToACSConversion"),
    over("annotated_somatic_variants", "File?", "GATK Funcotator"),
    over("annotated_somatic_variants_idx", "File?", "GATK Funcotator"),
    over("absolute_snv_maf", "File?", "ABSOLUTE input"),
    over("absolute_indel_maf", "File?", "ABSOLUTE input"),
    over("absolute_acr_rdata", "File?", "ABSOLUTE"),
    over("absolute_acr_plot", "File?", "ABSOLUTE"),
    over("absolute_acr_solution", "Int?", "manual"),
    over("absolute_acr_maf", "File?", "ABSOLUTE + Postprocess"),
    over("absolute_acr_segtab", "File?", "ABSOLUTE + Postprocess"),
    over("absolute_acr_segtab_igv", "File?", "ABSOLUTE + Postprocess"),
    over("absolute_acr_table", "File?", "ABSOLUTE"),
    over("absolute_acr_purity", "Float?", "ABSOLUTE"),
    over("absolute_acr_ploidy", "Float?", "ABSOLUTE"),
    over("absolute_tcr_rdata", "File?", "ABSOLUTE"),
    over("absolute_tcr_plot", "File?", "ABSOLUTE"),
    over("absolute_tcr_solution", "Int?", "manual"),
    over("absolute_tcr_maf", "File?", "ABSOLUTE + Postprocess"),
    over("absolute_tcr_segtab", "File?", "ABSOLUTE + Postprocess"),
    over("absolute_tcr_segtab_igv", "File?", "ABSOLUTE + Postprocess"),
    over("absolute_tcr_table", "File?", "ABSOLUTE"),
    over("absolute_tcr_purity", "Float?", "ABSOLUTE"),
    over("absolute_tcr_ploidy", "Float?", "ABSOLUTE"),
    over("timepoint", "Int?", "manual"),
    over("user_purity", "Float?", "manual: forces ABSOLUTE alpha/tau (both or neither)"),
    over("user_ploidy", "Float?", "manual: forces ABSOLUTE alpha/tau (both or neither)"),
]

# Patient's sample-set fields are computed inside UpdateSamples (not overlay inputs and
# not plain carry-overs), so the generator special-cases them by name; every other field
# is carried over from the input patient. The overlaid/fixed flags here would drive a
# generated bulk UpdatePatient / DefinePatient, which follow the same recipe.
PATIENT_COMPUTED = {"samples", "tumor_samples", "normal_samples", "matched_normal_sample"}

PATIENT = [
    fixed("name", "String"),
    fixed("sex", "String?"),
    fixed("samples", "Array[Sample]"),
    fixed("tumor_samples", "Array[Sample]"),
    fixed("normal_samples", "Array[Sample]"),
    fixed("has_tumor", "Boolean"),
    fixed("has_normal", "Boolean"),
    fixed("matched_normal_sample", "Sample?"),
    over("shards", "Array[Shard]"),
    over("raw_snv_calls_vcf", "File?"),
    over("raw_snv_calls_vcf_idx", "File?"),
    over("mutect2_stats", "File?"),
    over("orientation_bias", "File?"),
    over("filtered_vcf", "File?"),
    over("filtered_vcf_idx", "File?"),
    over("filtering_stats", "File?"),
    over("mask_vcf", "File?"),
    over("mask_vcf_idx", "File?"),
    over("mask_name", "String?"),
    over("somatic_vcf", "File?"),
    over("somatic_vcf_idx", "File?"),
    over("num_somatic_variants", "Int?"),
    over("germline_vcf", "File?"),
    over("germline_vcf_idx", "File?"),
    over("num_germline_variants", "Int?"),
    over("somatic_calls_bam", "File?"),
    over("somatic_calls_bai", "File?"),
    over("rare_germline_alleles", "File?"),
    over("rare_germline_alleles_idx", "File?"),
    over("gvcf", "File?"),
    over("gvcf_idx", "File?"),
    over("snp_ref_counts", "File?"),
    over("snp_alt_counts", "File?"),
    over("snp_other_alt_counts", "File?"),
    over("snp_sample_correlation", "File?"),
    over("snp_sample_correlation_min", "Float?"),
    over("modeled_segments", "File?"),
    over("phylogic_sif_file", "File?"),
    over("phylogic_report", "File?"),
    over("phylogic_ccfs_cnvs", "File?"),
    over("phylogic_ccfs_snvs", "File?"),
    over("phylogic_constrained_ccf", "File?"),
    over("phylogic_cluster_ccfs", "File?"),
    over("phylogic_build_tree_posteriors", "File?"),
    over("phylogic_growth_rates", "File?"),
    over("phylogic_growth_rate_plot", "File?"),
    over("phylogic_timing_report", "File?"),
    over("phylogic_timing_wgd_supporting_events", "File?"),
    over("phylogic_timing_graph", "File?"),
    over("phylogic_timing_comparison", "File?"),
    over("phylogic_timing_table", "File?"),
]

SHARD = [
    fixed("id", "Int"),
    fixed("intervals", "File"),
    fixed("skip", "Boolean"),
    fixed("is_high_mem", "Boolean"),
    over("raw_calls_mutect2_vcf", "File?"),
    over("raw_calls_mutect2_vcf_idx", "File?"),
    over("raw_mutect2_stats", "File?"),
    over("raw_mutect2_bam_out", "File?"),
    over("raw_mutect2_bai_out", "File?"),
    over("raw_mutect2_artifact_priors", "File?"),
]

STRUCTS = {"SequencingRun": SEQUENCING_RUN, "Sample": SAMPLE, "Shard": SHARD, "Patient": PATIENT}

