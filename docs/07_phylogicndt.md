# PhylogicNDT

[PhylogicNDT](https://github.com/broadinstitute/PhylogicNDT) is a Bayesian framework for analyzing the clonal architecture and evolutionary history of tumors. Given cancer cell fractions of somatic mutations (from ABSOLUTE), it clusters co-occurring mutations into clonal populations, infers phylogenetic tree topologies relating those clones, and estimates the relative timing of somatic events (both point mutations and focal copy-number alterations)

## Goal

This stage performs multi-sample clonal clustering, phylogenetic reconstruction, and molecular timing using somatic variants and ABSOLUTE-informed copy-number context. In this workflow, the implementation additionally supports focal CNAs as timing events. Events with strong local support can be timed directly; events with insufficient support may still be retained as prior-only or otherwise non-estimable timing events rather than being forced into an overconfident posterior.

## When it runs

PhylogicNDT is most useful after a successful ABSOLUTE selection/extraction step. In practice, it is part of the second-pass workflow after the ABSOLUTE manual review.

## Outputs

PhylogicNDT is run `phylogic_n_replicates` times (independent clustering and MCMC replicates, patient id suffixed with `.REP<i>`). Every `phylogic_*` field on the `Patient` is therefore an `Array[File]` holding one entry per replicate.

Key outputs include:

- `phylogic_report` — HTML report with tree visualizations and cluster summaries
- `phylogic_sif_file` — sample information file used as PhylogicNDT input
- `phylogic_growth_rates` — estimated clonal growth rates (if timepoints provided)
- `phylogic_growth_rate_plot` — growth rate visualization
- `phylogic_timing_report` — molecular timing report
- `phylogic_timing_table` — timing estimates per event
- `phylogic_timing_graph` — timing DAG visualization
- `phylogic_timing_comparison` — cross-sample timing comparison
- `phylogic_timing_wgd_supporting_events` — events supporting whole-genome doubling timing

## References

- PhylogicNDT — Leshchiner I, Livitz D, Gainor JF, et al. Comprehensive analysis of tumour initiation, spatial and temporal progression under multiple lines of treatment. *bioRxiv* 508127 (2019). <https://doi.org/10.1101/508127>
  - Implementation used: <https://github.com/phylyc/PhylogicNDT> (a fork of the original with bugfixes and runtime improvements).

See [References](11_references.md) for the full bibliography.
