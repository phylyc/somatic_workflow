version development


struct Shard {
    Int id
    File intervals
    Boolean skip
    Boolean is_high_mem
    File? raw_calls_mutect2_vcf
    File? raw_calls_mutect2_vcf_idx
    File? raw_mutect2_stats
    File? raw_mutect2_bam_out
    File? raw_mutect2_bai_out
    File? raw_mutect2_artifact_priors
}


workflow UpdateShard {
    input {
        Shard shard
        Int? id
        Boolean? skip
        Boolean? is_high_mem
        File? intervals
        File? raw_calls_mutect2_vcf
        File? raw_calls_mutect2_vcf_idx
        File? raw_mutect2_stats
        File? raw_mutect2_bam_out
        File? raw_mutect2_bai_out
        File? raw_mutect2_artifact_priors
    }

    Shard s = object {
        id: if (defined(id)) then id else shard.id,
        intervals: if (defined(intervals)) then intervals else shard.intervals,
        skip: if (defined(skip)) then skip else shard.skip,
        is_high_mem: if (defined(is_high_mem)) then is_high_mem else shard.is_high_mem,
        raw_calls_mutect2_vcf: if (defined(raw_calls_mutect2_vcf)) then raw_calls_mutect2_vcf else shard.raw_calls_mutect2_vcf,
        raw_calls_mutect2_vcf_idx: if (defined(raw_calls_mutect2_vcf_idx)) then raw_calls_mutect2_vcf_idx else shard.raw_calls_mutect2_vcf_idx,
        raw_mutect2_stats: if (defined(raw_mutect2_stats)) then raw_mutect2_stats else shard.raw_mutect2_stats,
        raw_mutect2_bam_out: if (defined(raw_mutect2_bam_out)) then raw_mutect2_bam_out else shard.raw_mutect2_bam_out,
        raw_mutect2_bai_out: if (defined(raw_mutect2_bai_out)) then raw_mutect2_bai_out else shard.raw_mutect2_bai_out,
        raw_mutect2_artifact_priors: if (defined(raw_mutect2_artifact_priors)) then raw_mutect2_artifact_priors else shard.raw_mutect2_artifact_priors
    }

    output {
        Shard updated_shard = s
    }
}