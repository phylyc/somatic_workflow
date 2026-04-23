# PhylogicNDT

## Goal

This stage performs multi-sample clonal clustering, phylogenetic reconstruction, and molecular timing using somatic variants and ABSOLUTE-informed copy-number context.

## When it runs

PhylogicNDT is most useful after a successful ABSOLUTE selection/extraction step. In practice, it is part of the second-pass workflow after the ABSOLUTE manual review.

## User-facing outputs

The most user-facing deliverable is:

- `phylogic_report`

Additional outputs may be important for expert downstream analysis, including cluster- and event-level posterior summaries.

## User-facing method summary

PhylogicNDT clusters mutations across samples from the same patient using cancer cell fraction profiles, reconstructs likely phylogenetic trees, and estimates relative molecular timing.

In this workflow, the implementation additionally supports focal CNAs as timing events. Events with strong local support can be timed directly; events with insufficient support may still be retained as prior-only or otherwise non-estimable timing events rather than being forced into an overconfident posterior.

## Interpretation guidance

For most users, the HTML report is the right entry point. Expert users may then inspect detailed posterior files and cluster-level outputs for downstream analysis.

## Caveat

This page intentionally keeps the mathematical details light. For a scientific description, see [Methodology](methodology.md).
