# ABSOLUTE review and rerun

[ABSOLUTE](https://www.nature.com/articles/nbt.2203) is an algorithm that infers tumor purity (fraction of tumor cells in the sample) and ploidy (average genome-wide copy number) from segmented copy-ratio data and point mutation variant allele fractions. Given a purity/ploidy solution, it converts relative copy ratios into integer (absolute) copy numbers and estimates the cancer cell fraction of each somatic mutation.

ABSOLUTE initially outputs multiple purity/ploidy combinations, that can explain the same data — it produces a ranked list of candidate solutions that must be manually reviewed to select the most plausible one.

## Two-pass workflow

1. Run the workflow to generate candidate ABSOLUTE solutions and review plots
2. Manually choose one solution per sample
3. Rerun the workflow with cached outputs plus `Cache.absolute_solution`. You can find a minimal rerun input json below. Note that this only contains the additional Cache inputs. Ensure that other inputs from the first run remain filled out

## What to review

The main file to inspect is:

- `absolute_acr_plot`

Review the candidate solutions for each sample and decide which one is biologically and technically most plausible.

## Practical review suggestions

- Whole genome Allelic copy number plot:
  1. The goal is to align the copy number segments to the most plausible integer copy number grid.
  2. To achieve this, look for modal peaks on the "Seg clonality" part of the plot (right side). The goal is to align these peaks with the integer copy number gridlines.

- Point mutation VAF Density plot
  1. The goal is to align the purity/2 line with the VAF peak of somatic SNVs

- ABSOLUTE solutions are ranked by a composite likelihood score. In many cases, the correct solution is among the first few. However, treat this only as a heuristic for prioritizing review order, not as a guarantee of correctness - especially for low-purity or highly aneuploid samples.

- If running multiple samples per patient, it can be helpful to look at acr plots for different samples side by side. Look for shared events between samples and try to assign them to the same integer copy number

- When in doubt, pick the solution that is biologically most probable

This page is not intended to be a full ABSOLUTE tutorial; it is a workflow-operational guide. For more details, refer to the ABSOLUTE paper linked above

## How to encode the selected solution

Provide `Cache.absolute_solution`. This is an array of integers in the same order as `Cache.sample_names_ordered`

- positive integer: selected solution number on the plot page
- `0`: sample has effectively zero or too little tumor content for purity estimation
- `-1`: sample failed CNV calling for practical purposes, for example because segmentation is too fragmented despite evident aneuploidy

## Minimal example rerun JSON

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

Adjust the other cached inputs to match the actual outputs of the first run.

## What happens after rerun

After the rerun, the workflow performs ABSOLUTE extraction and rescues:

- dropped copy-number segments
- dropped somatic SNVs/indels

This can be viewed as a completion layer around the ABSOLUTE output, using the selected purity/ploidy solution.
