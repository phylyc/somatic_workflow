# Supported use cases

## Assay-specific notes

### Whole-exome sequencing

This is the main target use case. Default runtime settings are optimized for WES.

### Whole-genome sequencing

WGS is supported, but is less extensively tested than WES. Runtime settings may need to be adjusted for scale.

### Targeted panels

Targeted panel data below roughly a few hundred genes is best treated as SNV-only input. CNV inference is generally too sparse to be reliable for small targeted panels.

### ULP WGS

ULP WGS is somewhat supported for CNV calling. Purity/ploidy estimation is not considered supported in this mode.

## Matched normals versus tumor-only

Matched normals improve:

- artifact suppression
- germline discrimination
- confidence in somatic interpretation

However, matched normals are not mandatory. If no normal sample is provided, the workflow treats all samples as tumor samples. Rare germline SNVs that were misclassified as somatic will likely cluster on the trunk of the clonal phylogenetic tree.

## Mixed sequencing platforms

Mixed platforms are technically allowed, but should be used with care.

For CNV analysis, target harmonization reduces the usable signal to the intersection of target intervals across runs/samples. This can seriously degrade CNV quality if platforms target very different regions.

Examples:

- WES + WES from similar capture designs: usually acceptable
- WES + targeted panel: often problematic for CNV analysis
- ULP WGS + WES for the same biospecimen: useful in some contexts, but often one modality should be excluded from specific copy-ratio calculations

For SNV calling, mixed platforms can still be useful, but error-model differences matter. In such cases, consider adjusting `Parameters.mutect2_pcr_snv_qual` and related settings appropriately.
