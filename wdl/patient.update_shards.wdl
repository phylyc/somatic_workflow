version development

import "shard.wdl" as sh
import "patient.wdl" as p


workflow UpdateShards {
    input {
        Patient patient
        Array[File]? raw_calls_mutect2_vcf_scattered
        Array[File]? raw_calls_mutect2_vcf_idx_scattered
        Array[File]? raw_mutect2_stats_scattered
        Array[File]? raw_mutect2_bam_out_scattered
        Array[File]? raw_mutect2_bai_out_scattered
        Array[File]? raw_mutect2_artifact_priors_scattered
    }

    Array[File] rcmvs = select_first([raw_calls_mutect2_vcf_scattered, []])
    if (length(rcmvs) > 0) {
        scatter (pair in zip(patient.shards, rcmvs)) {
            call sh.UpdateShard as UpdateRawMutect2vcf {
                input:
                    shard = pair.left,
                    raw_calls_mutect2_vcf = pair.right,
            }
        }
    }
    Array[Shard] shards_m2_vcf = select_first([UpdateRawMutect2vcf.updated_shard, patient.shards])

    Array[File] rcmvixs = select_first([raw_calls_mutect2_vcf_idx_scattered, []])
    if (length(rcmvixs) > 0) {
        scatter (pair in zip(shards_m2_vcf, rcmvixs)) {
            call sh.UpdateShard as UpdateRawMutect2vcfIdx {
                input:
                    shard = pair.left,
                    raw_calls_mutect2_vcf_idx = pair.right,
            }
        }
    }
    Array[Shard] shards_m2_vcf_idx = select_first([UpdateRawMutect2vcfIdx.updated_shard, shards_m2_vcf])

    Array[File] rmss = select_first([raw_mutect2_stats_scattered, []])
    if (length(rmss) > 0) {
        scatter (pair in zip(shards_m2_vcf_idx, rmss)) {
            call sh.UpdateShard as UpdateRawMutect2Stats {
                input:
                    shard = pair.left,
                    raw_mutect2_stats = pair.right,
            }
        }
    }
    Array[Shard] shards_m2_stats = select_first([UpdateRawMutect2Stats.updated_shard, shards_m2_vcf_idx])

    Array[File] rmbams = select_first([raw_mutect2_bam_out_scattered, []])
    if (length(rmbams) > 0) {
        scatter (pair in zip(shards_m2_stats, rmbams)) {
            call sh.UpdateShard as UpdateRawMutect2BamOut {
                input:
                    shard = pair.left,
                    raw_mutect2_bam_out = pair.right,
            }
        }
    }
    Array[Shard] shards_m2_bam_out = select_first([UpdateRawMutect2BamOut.updated_shard, shards_m2_stats])

    Array[File] rmbais = select_first([raw_mutect2_bai_out_scattered, []])
    if (length(rmbais) > 0) {
        scatter (pair in zip(shards_m2_bam_out, rmbais)) {
            call sh.UpdateShard as UpdateRawMutect2BaiOut {
                input:
                    shard = pair.left,
                    raw_mutect2_bai_out = pair.right,
            }
        }
    }
    Array[Shard] shards_m2_bai_out = select_first([UpdateRawMutect2BaiOut.updated_shard, shards_m2_bam_out])

    Array[File] rmaps = select_first([raw_mutect2_artifact_priors_scattered, []])
    if (length(rmaps) > 0) {
        scatter (pair in zip(shards_m2_bai_out, rmaps)) {
            call sh.UpdateShard as UpdateRawMutect2ArtifactPriors {
                input:
                    shard = pair.left,
                    raw_mutect2_artifact_priors = pair.right,
            }
        }
    }
    Array[Shard] shards_m2_artifact_priors = select_first([UpdateRawMutect2ArtifactPriors.updated_shard, shards_m2_bai_out])

    Array[Shard] shards = shards_m2_artifact_priors

    call p.UpdatePatient {
        input:
            patient = patient,
            shards = shards,
    }

    output {
        Patient updated_patient = UpdatePatient.updated_patient
    }
}