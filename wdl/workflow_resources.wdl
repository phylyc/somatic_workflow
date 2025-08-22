version development


struct WorkflowResources {
    # Intervals for short variant calling. The interval_list and interval_lists will
    # be combined and interval_blacklist will be subtracted from the result.
    File? interval_list
    File? interval_blacklist
    Array[File]? interval_lists
    File? preprocessed_intervals
    Array[File]? scattered_intervals_for_variant_calling
    Array[File]? scattered_intervals_for_pileups

    # Reference fasta, index, and dictionary files.
    File ref_fasta
    File ref_fasta_index
    File ref_dict

    # VCF file of variants to force mutation calling at.
    File? force_call_alleles
    File? force_call_alleles_idx
    # VCF file of common sequencing artifacts to filter out.
    File? snv_panel_of_normals
    File? snv_panel_of_normals_idx
    # VCF file of germline alleles with population allele frequencies.
    File? germline_resource
    File? germline_resource_idx
    # VCF file of common germline alleles (population allele frequency > 5%) to collect
    # allelic counts at for allelic copy ratio (aCR) and contamination estimation.
    File? common_germline_alleles
    File? common_germline_alleles_idx
    # BWA index image file for realignment task (to hg38).
    File? realignment_bwa_mem_index_image
    # List of transcript names to use for the Funcotator annotation.
    File? funcotator_transcript_list
    # Tarball of data sources for Funcotator. If not provided, the tarball will automatically
    # be downloaded from the GATK resource bundle, which is much slower.
    File? funcotator_data_sources_tar_gz
    # File containing HUGO symbols of driver genes, one per line, to be annotated in the tree.
    File? phylogic_driver_genes_file

    # Mutect1 specific resources
    File? germline_resource_v4_1
    File? germline_resource_v4_1_idx
    File? snv_panel_of_normals_v4_1
    File? snv_panel_of_normals_v4_1_idx
}


workflow DefineWorkflowResources {
    input {
        File? interval_list
        File? interval_blacklist
        Array[File]? interval_lists
        File? preprocessed_intervals
        Array[File]? scattered_intervals_for_variant_calling
        Array[File]? scattered_intervals_for_pileups
        File? ref_fasta
        File? ref_fasta_index
        File? ref_dict
        File? force_call_alleles
        File? force_call_alleles_idx
        File? snv_panel_of_normals
        File? snv_panel_of_normals_idx
        File? germline_resource
        File? germline_resource_idx
        File? common_germline_alleles
        File? common_germline_alleles_idx
        File? realignment_bwa_mem_index_image
        File? funcotator_transcript_list
        File? funcotator_data_sources_tar_gz
        File? phylogic_driver_genes_file
        
        # Mutect1 specific resources
        File? germline_resource_v4_1
        File? germline_resource_v4_1_idx
        File? snv_panel_of_normals_v4_1
        File? snv_panel_of_normals_v4_1_idx
    }

    WorkflowResources files = object {
        interval_list: interval_list,
        interval_blacklist: interval_blacklist,
        interval_lists: interval_lists,
        preprocessed_intervals: preprocessed_intervals,
        scattered_intervals_for_variant_calling: scattered_intervals_for_variant_calling,
        scattered_intervals_for_pileups: scattered_intervals_for_pileups,
        ref_fasta: select_first([ref_fasta]),
        ref_fasta_index: select_first([ref_fasta_index]),
        ref_dict: select_first([ref_dict]),
        force_call_alleles: force_call_alleles,
        force_call_alleles_idx: force_call_alleles_idx,
        snv_panel_of_normals: snv_panel_of_normals,
        snv_panel_of_normals_idx: snv_panel_of_normals_idx,
        germline_resource: germline_resource,
        germline_resource_idx: germline_resource_idx,
        common_germline_alleles: common_germline_alleles,
        common_germline_alleles_idx: common_germline_alleles_idx,
        realignment_bwa_mem_index_image: realignment_bwa_mem_index_image,
        funcotator_transcript_list: funcotator_transcript_list,
        funcotator_data_sources_tar_gz: funcotator_data_sources_tar_gz,
        phylogic_driver_genes_file: phylogic_driver_genes_file,
        germline_resource_v4_1: germline_resource_v4_1,
        germline_resource_v4_1_idx: germline_resource_v4_1_idx,
        snv_panel_of_normals_v4_1: snv_panel_of_normals_v4_1,
        snv_panel_of_normals_v4_1_idx: snv_panel_of_normals_v4_1_idx,
    }

    output {
        WorkflowResources resources = files
    }
}