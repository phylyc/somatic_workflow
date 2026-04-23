# Reference resources

This page summarizes the external reference data and auxiliary resources used by the workflow.

## Absolutely mandatory

Only the reference genome triplet is truly mandatory:

- reference FASTA
- FASTA index
- sequence dictionary

Example public bucket root:

- `gs://gcp-public-data--broad-references`

Example hg19 FASTA path:

- `gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta`

## Strongly recommended resources

For a sensible run, users should usually also provide:

- germline resource VCF with population allele frequencies
- common germline alleles VCF
- SNV panel of normals
- CNV panel(s) of normals
- force-call allele resource, where appropriate
- Funcotator data sources
- realignment BWA-MEM index image if the realignment artifact filter is used

## Suggested public sources

### Reference genomes

Use exact public reference files from:

- `gs://gcp-public-data--broad-references/hg19/...`
- `gs://gcp-public-data--broad-references/hg38/...`

`<TODO_EXACT_HG19_RESOURCE_TABLE>`
`<TODO_EXACT_HG38_RESOURCE_TABLE>`

### Germline resource and common germline alleles

Recommended source family:

- gnomAD

`<TODO_GNOMAD_RESOURCE_EXAMPLES>`

### Funcotator resources

Recommended source family:

- GATK Funcotator bundle

`<TODO_FUNCOTATOR_BUNDLE_LINK>`

### Force-call alleles

A force-calling resource may be useful for known drivers or targeted mutation lists.

Example source family:

- COSMIC targeted screens mutant VCF

`https://cancer.sanger.ac.uk/cosmic/download/cosmic/v103/completetargetedscreensmutantvcf`

## Panels of normals

### SNV PoN

A public SNV PoN may exist, but users should not assume that an old generic PoN is optimal for a modern assay. The preferred practice is to provide an assay-appropriate PoN when possible.

The repository contains workflows for generating PoNs:

- `wdl/resources/create_snv_pon.wdl`
- `wdl/resources/create_cnv_pon.wdl`

### CNV PoN

Preferred practice:

- platform-specific
- sex-specific when chrX/chrY behavior matters

If no suitable CNV PoN is available for a run, an empty file can be supplied as a fallback.

## Outbound network access

Outbound internet access is only needed for a few situations:

- container retrieval
- helper scripts provided by URL instead of local files
- automatic Funcotator bundle download if not supplied locally
- HTML formatting steps for PhylogicNDT output

Where reproducibility matters, prefer local or pinned resources over implicit downloads.
