version development

import "runtime_collection.wdl" as rtc
import "runtimes.wdl" as rt
import "patient.wdl"


workflow CallAncestry {
    input {
        Patient patient
        String genome_build
        RuntimeCollection runtime_collection
    }

    if (defined(patient.gvcf) && defined(patient.gvcf_idx) && (genome_build == "hg38" || genome_build == "hg19")) {
        call RunPeddy {
            input:
                patient_id = patient.name,
                gvcf = select_first([patient.gvcf]),
                gvcf_index = select_first([patient.gvcf_idx]),
                genome_build = genome_build,
                runtime_params = runtime_collection.ancestry
        }
    }

    output {
        String? ancestry_pred = RunPeddy.ancestry_pred
        Float? ancestry_prob = RunPeddy.ancestry_prob
        File? background_pca_table = RunPeddy.background_pca_table
        File? pca_plot = RunPeddy.pca_plot
    }
}

task RunPeddy {
    input {
        String patient_id
        File gvcf
        File gvcf_index
        String genome_build
        Runtime runtime_params
    }
    command <<<
        set -e

        # Filter the VCF to only retain autosomes
        if [[ "~{genome_build}" == "hg38" ]]; then
            regions=$(echo chr{1..22} | tr ' ' ',')
        else
            regions=$(echo {1..22} | tr ' ' ',')
        fi
        bcftools view -r "$regions" "~{gvcf}" -Oz -o "~{patient_id}.filt.vcf.gz"
        bcftools index --tbi "~{patient_id}.filt.vcf.gz"

        # Add synthetic AD and DP fields required by peddy
        zcat "~{patient_id}.filt.vcf.gz" | awk -F'\t' 'BEGIN{OFS="\t"}
          /^##/{print;next}
          /^#CHROM/{
            print "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"synthetic allele depths\">";
            print "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"synthetic depth\">";
            print; next}
          {
            $9=$9":AD:DP";
            for(i=10;i<=NF;i++){
              split($i,a,":"); g=a[1]; gsub(/\|/,"/",g);
              if(g=="0/0") ad="30,0"; else if(g=="1/1") ad="0,30";
              else if(g=="0/1"||g=="1/0") ad="15,15"; else ad="0,0";
              $i=$i":"ad":"(ad=="0,0"?0:30);
            }
            print}' | bgzip > "~{patient_id}.peddy_input.vcf.gz"
        bcftools index --tbi "~{patient_id}.peddy_input.vcf.gz"

        # Create dummy ped file
        printf "~{patient_id}\t~{patient_id}\t0\t0\t0\t-9\n" > "~{patient_id}.ped"

        # Run peddy
        if [ "~{genome_build}" == "hg38" ]; then
            peddy --plot -p 4 --sites ~{genome_build} --prefix ~{patient_id} "~{patient_id}.peddy_input.vcf.gz" "~{patient_id}.ped"
        else
            peddy --plot -p 4 --prefix ~{patient_id} "~{patient_id}.peddy_input.vcf.gz" "~{patient_id}.ped"
        fi

        # Getting the prediction and associated probability
        awk -F',' 'NR==2 {print $12}' "~{patient_id}.het_check.csv" > pred.txt
        awk -F',' 'NR==2 {print $13}' "~{patient_id}.het_check.csv" > prob.txt
    >>>
    output {
        String ancestry_pred = read_string("pred.txt")
        Float ancestry_prob = read_float("prob.txt")
        File background_pca_table = "~{patient_id}.background_pca.json"
        File pca_plot = "~{patient_id}.pca_check.png"
    }
    runtime {
        docker: runtime_params.docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: runtime_params.machine_mem + " MB"
        runtime_minutes: runtime_params.runtime_minutes
        disks: "local-disk " + runtime_params.disk + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }
}
