version development


struct WorkflowResources {
    File? interval_list
    File? interval_blacklist
    Array[File]? interval_lists

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    File? force_call_alleles
    File? force_call_alleles_idx
    File? snv_panel_of_normals
    File? snv_panel_of_normals_idx
    File? germline_resource
    File? germline_resource_tbi
    File? common_germline_alleles
    File? common_germline_alleles_idx
    File? realignment_bwa_mem_index_image
    File? funcotator_transcript_list
    File? funcotator_data_sources_tar_gz
}


workflow DefineWorkflowResources {
    input {
        File? interval_list
        File? interval_blacklist
        Array[File]? interval_lists
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        File? force_call_alleles
        File? force_call_alleles_idx
        File? snv_panel_of_normals
        File? snv_panel_of_normals_idx
        File? germline_resource
        File? germline_resource_tbi
        File? common_germline_alleles
        File? common_germline_alleles_idx
        File? realignment_bwa_mem_index_image
        File? funcotator_transcript_list
        File? funcotator_data_sources_tar_gz
    }

    WorkflowResources files = object {
        interval_list: interval_list,
        interval_blacklist: interval_blacklist,
        interval_lists: interval_lists,
        ref_fasta: ref_fasta,
        ref_fasta_index: ref_fasta_index,
        ref_dict: ref_dict,
        force_call_alleles: force_call_alleles,
        force_call_alleles_idx: force_call_alleles_idx,
        snv_panel_of_normals: snv_panel_of_normals,
        snv_panel_of_normals_idx: snv_panel_of_normals_idx,
        germline_resource: germline_resource,
        germline_resource_tbi: germline_resource_tbi,
        common_germline_alleles: common_germline_alleles,
        common_germline_alleles_idx: common_germline_alleles_idx,
        realignment_bwa_mem_index_image: realignment_bwa_mem_index_image,
        funcotator_transcript_list: funcotator_transcript_list,
        funcotator_data_sources_tar_gz: funcotator_data_sources_tar_gz
    }

    output {
        WorkflowResources resources = files
    }
}