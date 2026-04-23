# ABSOLUTE review and rerun

## Why this step is separate

ABSOLUTE is intentionally handled as a two-pass process:

1. run the workflow to generate candidate ABSOLUTE solutions and review plots
2. manually choose one solution per sample
3. rerun the workflow with cached outputs plus `Cache.absolute_solution`

This manual checkpoint is expected normal behavior, not an exceptional workaround.

## What to review

The main file to inspect is:

- `absolute_acr_plot`

Review the candidate solutions for each sample and decide which one is biologically and technically most plausible.

## Practical review criteria

Use the plots to assess:

1. alignment of copy-number levels with visible segmental plateaus
2. alignment of the purity/2 line with the VAF mode of somatic SNVs
3. consistency across all samples, especially with respect to whole-genome doubling support

This page is not intended to be a full ABSOLUTE tutorial; it is a workflow-operational guide.

## How to encode the selected solution

Provide:

- `Cache.absolute_solution`

This is an array of integers in the same order as the PDF plot pages.

Interpretation:

- positive integer: selected solution number on the plot page
- `0`: sample has effectively zero or too little tumor content for purity estimation
- `-1`: sample failed CNV calling for practical purposes, for example because segmentation is too fragmented despite evident aneuploidy

These values are expected and should be documented as standard workflow behavior.

## Minimal rerun JSON snippet

```json
{
  "MultiSampleSomaticWorkflow.patient_id": "Patient",
  "Cache.absolute_solution": [2, 0],
  "Cache.sample_names_ordered": [
    "<TODO_SAMPLE_ORDER_1>",
    "<TODO_SAMPLE_ORDER_2>"
  ],
  "Cache.acs_copy_ratio_segmentation": [
    "<TODO_PREVIOUS_ACS_SEGMENTATION_1>",
    "<TODO_PREVIOUS_ACS_SEGMENTATION_2>"
  ],
  "Cache.absolute_acr_rdata": [
    "<TODO_PREVIOUS_ABSOLUTE_RDATA_1>",
    "<TODO_PREVIOUS_ABSOLUTE_RDATA_2>"
  ],
  "Cache.absolute_acr_plot": [
    "<TODO_PREVIOUS_ABSOLUTE_PLOT_1>",
    "<TODO_PREVIOUS_ABSOLUTE_PLOT_2>"
  ]
}
```

Adjust the cached inputs to match the actual outputs of the first run.

## What happens after rerun

After the rerun, the workflow performs ABSOLUTE extraction and rescues:

- dropped copy-number segments
- dropped somatic SNVs/indels

This can be viewed as a completion layer around the ABSOLUTE output, using the selected purity/ploidy solution.
