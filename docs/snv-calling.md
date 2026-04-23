# SNV calling

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
- optional force-call alleles
- optional matched normal sample designation
- optional realignment index image

## Outputs

Key outputs from this part of the workflow include:

- somatic VCF
- germline VCF
- gVCF
- filtered/annotated mutation calls
- rare germline allele set used downstream
- mutation-support BAM/BAM index when enabled

## User-facing method summary

The SNV workflow is organized around four ideas:

1. call likely short variants in each sample
2. use those calls to seed multi-sample calling
3. apply both statistical and hard filters
4. annotate the final somatic calls

A Mutect1-based callset is used as an additional source of candidate alleles for force-calling with Mutect2. This is intentional: Mutect1 and Mutect2 have different behavior, and the force-call strategy helps recover variants that are obvious in the reads but may be missed by Mutect2 graph assembly in some regions.

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

## Force-calling resource

A force-calling resource may be supplied, for example a known mutation set or a cancer hotspot resource.

`<TODO_FORCE_CALL_RESOURCE_LINK>`

## Parameters users are most likely to adjust

These are the main SNV-related settings that ordinary users may need to touch:

- `Parameters.hard_filter_min_total_depth`
- `Parameters.hard_filter_min_total_alt_count`
- `Parameters.min_read_depth`
- `Parameters.mutect2_pcr_indel_qual`
- `Parameters.mutect2_pcr_snv_qual`
- `Parameters.run_*` toggles for enabling/disabling specific workflow blocks

Most other settings should be treated as advanced.

## Notes for mixed platforms

When combining data from different platforms, SNV calling can still be informative, but platform-specific error behavior matters. In those settings, `mutect2_pcr_snv_qual` may need to be adapted to the effective mixed-platform error profile.
