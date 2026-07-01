# Reference resources

This page provides links to some publicly available reference files for the workflow


| Variable | b37 | hg38 |
|---------------------------------|------------------|---------------------|
| ref_fasta | gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta | gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta |
| ref_fasta_index | gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai | gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai |
| ref_fasta_dict | gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.dict | gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict |
| funcotator_data_sources_tar_gz | gs://gcp-public-data--broad-references/funcotator/funcotator_dataSources.v1.8.hg19.20230908s.tar.gz | gs://gcp-public-data--broad-references/funcotator/funcotator_dataSources.v1.8.hg38.20230908s.tar.gz |
| funcotator_transcript_list | gs://broad-public-datasets/funcotator/transcriptList.exact_uniprot_matches.AKT1_CRLF2_FGFR1.txt |  |


### Force-call alleles

A force-calling resource may be useful for known drivers or targeted mutation lists.

Example source family:

- COSMIC targeted screens mutant VCF

`https://cancer.sanger.ac.uk/cosmic/download/cosmic/v103/completetargetedscreensmutantvcf`

## Outbound network access

Outbound internet access is only needed for a few situations:

- container retrieval
- helper scripts provided by URL instead of local files
- automatic Funcotator bundle download if not supplied locally
- HTML formatting steps for PhylogicNDT output

Where reproducibility matters, prefer local or pinned resources over implicit downloads.
