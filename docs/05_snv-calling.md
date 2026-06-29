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

    -   Dockstore: <https://dockstore.org/workflows/github.com/phylyc/somatic_workflow/CreateSNVPanelOfNormals:master?tab=info>
    -   WDL: <https://github.com/phylyc/somatic_workflow/blob/master/wdl/resources/create_snv_pon.wdl>

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

- `Parameters.hard_filter_min_total_depth`
- `Parameters.hard_filter_min_total_alt_count`
- `Parameters.min_read_depth`
- `Parameters.mutect2_pcr_indel_qual`
- `Parameters.mutect2_pcr_snv_qual`

Most other settings should be treated as advanced.

## Important caveats for mixed platforms

When combining data from different platforms, SNV calling can still be informative, but platform-specific error behavior matters. In those settings, `mutect2_pcr_snv_qual` may need to be adapted to the effective mixed-platform error profile.
