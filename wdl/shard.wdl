version development


struct Shard {
    File intervals
    Boolean skip
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
        Boolean? skip
        File? intervals
        File? raw_calls_mutect2_vcf
        File? raw_calls_mutect2_vcf_idx
        File? raw_mutect2_stats
        File? raw_mutect2_bam_out
        File? raw_mutect2_bai_out
        File? raw_mutect2_artifact_priors
    }

    Shard s = object {
        intervals: if (defined(intervals)) then intervals else shard.intervals,
        skip: if (defined(skip)) then skip else shard.skip,
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