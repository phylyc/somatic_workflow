version development


import "workflow_resources.wdl" as wfres


workflow UpdateWorkflowResources {
    input {
        WorkflowResources resources
        File? interval_list
        File? interval_blacklist
        Array[File]? interval_lists
        File? preprocessed_intervals
        Array[File]? scattered_intervals_for_variant_calling_m1
        Array[File]? scattered_intervals_for_variant_calling_m2
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
    }

    WorkflowResources files = object {
        interval_list: if (defined(interval_list)) then interval_list else resources.interval_list,
        interval_blacklist: if (defined(interval_blacklist)) then interval_blacklist else resources.interval_blacklist,
        interval_lists: if (defined(interval_lists)) then interval_lists else resources.interval_lists,
        preprocessed_intervals: if (defined(preprocessed_intervals)) then preprocessed_intervals else resources.preprocessed_intervals,
        scattered_intervals_for_variant_calling_m1: if (defined(scattered_intervals_for_variant_calling_m1)) then scattered_intervals_for_variant_calling_m1 else resources.scattered_intervals_for_variant_calling_m1,
        scattered_intervals_for_variant_calling_m2: if (defined(scattered_intervals_for_variant_calling_m2)) then scattered_intervals_for_variant_calling_m2 else resources.scattered_intervals_for_variant_calling_m2,
        scattered_intervals_for_pileups: if (defined(scattered_intervals_for_pileups)) then scattered_intervals_for_pileups else resources.scattered_intervals_for_pileups,
        ref_fasta: if (defined(ref_fasta)) then select_first([ref_fasta]) else resources.ref_fasta,
        ref_fasta_index: if (defined(ref_fasta_index)) then select_first([ref_fasta_index]) else resources.ref_fasta_index,
        ref_dict: if (defined(ref_dict)) then select_first([ref_dict]) else resources.ref_dict,
        force_call_alleles: if (defined(force_call_alleles)) then force_call_alleles else resources.force_call_alleles,
        force_call_alleles_idx: if (defined(force_call_alleles_idx)) then force_call_alleles_idx else resources.force_call_alleles_idx,
        snv_panel_of_normals: if (defined(snv_panel_of_normals)) then snv_panel_of_normals else resources.snv_panel_of_normals,
        snv_panel_of_normals_idx: if (defined(snv_panel_of_normals_idx)) then snv_panel_of_normals_idx else resources.snv_panel_of_normals_idx,
        germline_resource: if (defined(germline_resource)) then germline_resource else resources.germline_resource,
        germline_resource_idx: if (defined(germline_resource_idx)) then germline_resource_idx else resources.germline_resource_idx,
        common_germline_alleles: if (defined(common_germline_alleles)) then common_germline_alleles else resources.common_germline_alleles,
        common_germline_alleles_idx: if (defined(common_germline_alleles_idx)) then common_germline_alleles_idx else resources.common_germline_alleles_idx,
        realignment_bwa_mem_index_image: if (defined(realignment_bwa_mem_index_image)) then realignment_bwa_mem_index_image else resources.realignment_bwa_mem_index_image,
        funcotator_transcript_list: if (defined(funcotator_transcript_list)) then funcotator_transcript_list else resources.funcotator_transcript_list,
        funcotator_data_sources_tar_gz: if (defined(funcotator_data_sources_tar_gz)) then funcotator_data_sources_tar_gz else resources.funcotator_data_sources_tar_gz
    }

    output {
        WorkflowResources updated_resources = files
    }
}