# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.0] - 2026-07-04

First supported release of the MultiSampleSomaticWorkflow (MSSW). This is a
complete rewrite and should be treated as a new baseline; it is **not**
compatible with the v1.0.0 (2024) prototype. See the release notes for the full
overview.

### Highlights
- **Multi-sample somatic pipeline** for WES, WGS, and panel DNA data, organized
  around a patient → samples → sequencing-runs model, with matched normal
  optional and tumor-only supported.
- **SNV calling**: Mutect1 single-sample force-calling into multi-sample
  Mutect2, followed by statistical, hard, and realignment filtering, plus
  functional annotation.
- **CNV calling**: multi-sample allelic + total copy-ratio segmentation, SNP
  genotyping across samples with HET phasing, per-sample copy-ratio inference,
  event calling, and plotting.
- **Clonal analysis**: ABSOLUTE (with a documented manual solution-review and
  re-run step) and PhylogicNDT phylogenetic, growth-kinetics, and timing
  analysis.
- **Genetic ancestry calling** via PCA against reference populations.
- **Contamination estimation**, coverage/callable-loci aggregation, and target
  harmonization across mixed panels.
- **Failure recovery** controls for retrying or skipping individual Mutect2
  shards.
- Full browsable documentation at <https://phylyc.github.io/somatic_workflow/>,
  a `CITATION.cff`, and Dockstore publication of the workflow and its resource
  sub-workflows.

### Requirements
- GATK 4.6.2.0
- Cromwell (tested on v86)
- Docker (tags listed under `docker/`)

[2.0.0]: https://github.com/phylyc/somatic_workflow/compare/v1.0.0...v2.0.0
