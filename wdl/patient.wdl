version development

import "sample.wdl" as s


struct Patient {
    String name
    Array[Sample] samples
    Array[Sample] tumor_samples
    Array[Sample] normal_samples

    Boolean has_tumor
    Boolean has_normal
    Sample? matched_normal_sample
    File? rare_germline_alleles
    File? rare_germline_alleles_idx
}


workflow UpdatePatient {
    input {
        Patient patient
        String? name
        Array[Sample]? samples
        Array[Sample]? tumor_samples
        Array[Sample]? normal_samples
        Boolean? has_tumor
        Boolean? has_normal
        Sample? matched_normal_sample
        File? rare_germline_alleles
        File? rare_germline_alleles_idx
    }

    Patient pat = object {
        name: select_first([name, patient.name]),
        samples: select_first([samples, patient.samples]),
        tumor_samples: select_first([tumor_samples, patient.tumor_samples]),
        normal_samples: select_first([normal_samples, patient.normal_samples]),
        has_tumor: select_first([has_tumor, patient.has_tumor]),
        has_normal: select_first([has_normal, patient.has_normal]),
        # cannot use select_first for optional fields:
        matched_normal_sample: if defined(matched_normal_sample) then matched_normal_sample else patient.matched_normal_sample,
        rare_germline_alleles: if defined(rare_germline_alleles) then rare_germline_alleles else patient.rare_germline_alleles,
        rare_germline_alleles_idx: if defined(rare_germline_alleles_idx) then rare_germline_alleles_idx else patient.rare_germline_alleles_idx
    }

    output {
        Patient updated_patient = pat
    }
}