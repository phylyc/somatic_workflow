version development


struct Sample {
    File bam
    File bai
    File target_intervals
    Boolean? paired_end
    String bam_sample_name
    String assigned_sample_name
    Boolean is_tumor
}


workflow Patient {
    input {
        String individual_id
        Array[File]+ tumor_bams
        Array[File]+ tumor_bais
        Array[File]+ tumor_target_intervals
        Array[String]+ tumor_bam_names
        Array[String]+ tumor_sample_names
        Array[File]? normal_bams
        Array[File]? normal_bais
        Array[File]? normal_target_intervals
        Array[String]? normal_bam_names
        Array[String]? normal_sample_names
    }

    Array[File] non_optional_normal_bams = select_first([normal_bams, []])
    Array[File] non_optional_normal_bais = select_first([normal_bais, []])
    Array[File] non_optional_normal_target_intervals = select_first([normal_target_intervals, []])
    Array[String] non_optional_normal_bam_names = select_first([normal_bam_names, []])
    Array[String] non_optional_normal_sample_names = select_first([normal_sample_names, []])

    # The object syntax is likely going to be deprecated in newer WDL versions!
    # Definition via {"bam": sample.left.left, ..., "is_tumor": true} leads to
    # Error(s): No coercion defined from wom value(s) '"true"' of type 'String' to 'Boolean'.
    # :(
    scatter (sample in transpose(
        [
            tumor_bams,
            tumor_bais,
            tumor_target_intervals,
            tumor_bam_names,
            tumor_sample_names
        ]
    )) {
        Sample tumors = object {
            bam: sample[0],
            bai: sample[1],
            target_intervals: sample[2],
            bam_sample_name: sample[3],
            assigned_sample_name: sample[4],
            is_tumor: true
        }
    }

    scatter (sample in transpose(
        [
            non_optional_normal_bams,
            non_optional_normal_bais,
            non_optional_normal_target_intervals,
            non_optional_normal_bam_names,
            non_optional_normal_sample_names
        ]
    )) {
        Sample normals = object {
            bam: sample[0],
            bai: sample[1],
            target_intervals: sample[2],
            bam_sample_name: sample[3],
            assigned_sample_name: sample[4],
            is_tumor: false
        }
    }

    output {
        String name = individual_id
        Array[Sample] samples = flatten([tumors, normals])
        Array[Sample] tumor_samples = select_first([tumors, []])
        Array[Sample] normal_samples = select_first([normals, []])
        Boolean has_tumor = defined(tumor_bams) && (length(tumor_bams) > 0)
        Boolean has_normal = defined(normal_bams) && (length(non_optional_normal_bams) > 0)

        Array[File] bams = flatten([tumor_bams, non_optional_normal_bams])
        Array[File] bais = flatten([tumor_bais, non_optional_normal_bais])
        Array[File] target_intervals = flatten([tumor_target_intervals, non_optional_normal_target_intervals])
        Array[String] bam_names = flatten([tumor_bam_names, non_optional_normal_bam_names])
        Array[String] tumor_bam_names = tumor_bam_names
        Array[String] normal_bam_names = non_optional_normal_bam_names
        Array[String] sample_names = flatten([tumor_sample_names, non_optional_normal_sample_names])
    }
}