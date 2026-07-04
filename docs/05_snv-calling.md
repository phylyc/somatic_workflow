# SNV calling

The SNV workflow is organized around four ideas:

1. call likely short variants in each sample
2. use those calls to seed multi-sample calling
3. apply both statistical and hard filters
4. annotate the final somatic calls

## Goal

Starting from aligned reads and reference resources, this stage produces filtered and annotated somatic SNV/indel calls and related germline-aware outputs.

## Inputs

Main inputs include:

- BAM/BAI files
- reference FASTA, FASTA index, and sequence dictionary
- target intervals
- germline resource
- common germline alleles
- SNV panel of normals
- force-call alleles
- matched normal sample designation
- realignment index image

## Outputs

Key outputs from this part of the workflow include:

- somatic VCF
- germline VCF
- gVCF
- filtered/annotated mutation calls
- rare germline allele set used downstream
- mutation-support BAM + BAI when enabled

## Panel of normals
Recommended: Use **CreateSNVPanelOfNormals** workflow to create `snv_panel_of_normals` and `snv_panel_of_normals_idx`. This is needed for the SNV workflow

    -   Dockstore: <https://dockstore.org/workflows/github.com/phylyc/somatic_workflow/CreateSNVPanelOfNormals:v2.0.0?tab=info>
    -   WDL: <https://github.com/phylyc/somatic_workflow/blob/v2.0.0/wdl/resources/create_snv_pon.wdl>

## Force calling and filtering

A Mutect1-based callset is used as an additional source of candidate alleles for force-calling with Mutect2. This is intentional: Mutect1 and Mutect2 have different behavior, and the force-call strategy helps recover variants that are obvious in the reads but may be missed by Mutect2 graph assembly in some regions.

An additional force-calling resource may be supplied, for example a known mutation set or a cancer hotspot resource (see docs/10_resources.md)

Downstream filtering uses information from:

- read orientation bias
- contamination estimates
- copy-ratio segmentation
- germline resources
- sequencing-artifact resources
- hard thresholds on read depth, alternate counts, and read-quality metrics
- optional realignment-artifact filtering

## Matched normal versus tumor-only

If a matched normal is available, it improves filtering and interpretation. If no matched normal exists, the workflow still runs in tumor-only mode.

## Filtering parameters users are most likely to adjust

These are the main SNV-related settings that users may wish to change:

- `z3_Parameters.hard_filter_min_total_depth`
- `z3_Parameters.hard_filter_min_total_alt_count`
- `z3_Parameters.min_read_depth`
- `z3_Parameters.mutect2_pcr_indel_qual`
- `z3_Parameters.mutect2_pcr_snv_qual`

Most other settings should be treated as advanced.

## Important caveats for mixed platforms

When combining data from different platforms, SNV calling can still be informative, but platform-specific error behavior matters. In those settings, `mutect2_pcr_snv_qual` may need to be adapted to the effective mixed-platform error profile.

## References

SNV calling, filtering and annotation are implemented with GATK. The underlying
methods and reference resources:

- Mutect — Cibulskis K, Lawrence M, Carter S, et al. *Nat Biotechnol* 31, 213–219 (2013). <https://doi.org/10.1038/nbt.2514>
- Mutect2 — Benjamin D, Sato T, Cibulskis K, Getz G, Stewart C, Lichtenstein L. *bioRxiv* 861054 (2019). <https://doi.org/10.1101/861054>
- GATK (incl. Funcotator annotation) — <https://github.com/broadinstitute/gatk>
- BWA-MEM (realignment-artifact filter) — Li H. *arXiv*:1303.3997 (2013). <https://arxiv.org/abs/1303.3997>
- gnomAD (germline resource) — Karczewski KJ, et al. *Nature* 581, 434–443 (2020). <https://doi.org/10.1038/s41586-020-2308-7>
- COSMIC (example force-call resource) — Tate JG, et al. *Nucleic Acids Res* 47(D1):D941–D947 (2019). <https://doi.org/10.1093/nar/gky1015>

See [References](11_references.md) for the full bibliography.

