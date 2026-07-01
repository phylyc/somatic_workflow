# References

This page collects the primary methods and tools the workflow relies on, grouped
by stage. Please cite the relevant works when publishing results produced with
this pipeline.

## Core framework — GATK / WDL / Cromwell

- Van der Auwera GA, O'Connor BD. *Genomics in the Cloud: Using Docker, GATK, and
  WDL in Terra* (1st Edition). O'Reilly Media, 2020.
- GATK (v4.6.x), Broad Institute. GitHub repository:
  <https://github.com/broadinstitute/gatk>

Most preprocessing, SNV calling, contamination estimation, copy-ratio collection,
segmentation, and functional annotation (Funcotator) steps are implemented with
GATK tools.

## SNV calling

- **Mutect** — Cibulskis K, Lawrence M, Carter S, et al. Sensitive
  detection of somatic point mutations in impure and heterogeneous cancer
  samples. *Nat Biotechnol* 31, 213–219 (2013).
  <https://doi.org/10.1038/nbt.2514>
- **Mutect2** — Benjamin D, Sato T, Cibulskis K, Getz G, Stewart C,
  Lichtenstein L. Calling Somatic SNVs and Indels with Mutect2. *bioRxiv* 861054
  (2019). <https://doi.org/10.1101/861054>

## CNV calling

- **Tangent normalization** (denoising of total read counts) — Gao GF, Oh C,
  Saksena G, et al. Tangent normalization for somatic copy-number inference in
  cancer genome analysis. *Bioinformatics* 38(20):4677–4686 (2022).
  <https://doi.org/10.1093/bioinformatics/btac586>
- **Kernel change-point segmentation** (ModelSegments) — Celisse A, Marot G,
  Pierre-Jean M, Rigaill G. New efficient algorithms for multiple change-point
  detection with kernels. 2016. hal-01413230v2.
  <https://hal.science/hal-01413230v2>
- **GATK** (CollectReadCounts / DenoiseReadCounts / ModelSegments) — GATK (v4.6.x),
  Broad Institute. <https://github.com/broadinstitute/gatk>

## Purity, ploidy and absolute copy number — ABSOLUTE

- Carter SL, Cibulskis K, Helman E, et al. Absolute quantification of somatic DNA
  alterations in human cancer. *Nat Biotechnol* 30, 413–421 (2012).
  <https://doi.org/10.1038/nbt.2203>
  - Implementation used: <https://github.com/phylyc/absolute> (a fork of the
    original with bugfixes and runtime improvements).

## Clonal analysis — PhylogicNDT

- Leshchiner I, Livitz D, Gainor JF, et al. Comprehensive analysis of tumour
  initiation, spatial and temporal progression under multiple lines of treatment.
  *bioRxiv* 508127 (2019). <https://doi.org/10.1101/508127>
  - Implementation used: <https://github.com/phylyc/PhylogicNDT> (a fork of the
    original with bugfixes and runtime improvements).

## Ancestry inference — Peddy

- Pedersen BS, Quinlan AR. Who's Who? Detecting and Resolving Sample Anomalies in
  Human DNA Sequencing Studies with Peddy. *Am J Hum Genet* 100:406–413 (2017).
  <https://doi.org/10.1016/j.ajhg.2017.01.017>

## Supporting tools

- **SAMtools / BCFtools** (VCF manipulation, `samtools quickcheck` input
  validation) — Danecek P, Bonfield JK, Liddle J, et al. Twelve years of SAMtools
  and BCFtools. *GigaScience* 10(2):giab008 (2021).
  <https://doi.org/10.1093/gigascience/giab008>
- **BWA-MEM** (realignment-artifact filter) — Li H. Aligning sequence reads, clone
  sequences and assembly contigs with BWA-MEM. *arXiv*:1303.3997 (2013).
  <https://arxiv.org/abs/1303.3997>

## Reference data resources

- **gnomAD** (`germline_resource`, `common_germline_alleles`) — Karczewski KJ,
  Francioli LC, Tiao G, et al. The mutational constraint spectrum quantified from
  variation in 141,456 humans. *Nature* 581, 434–443 (2020).
  <https://doi.org/10.1038/s41586-020-2308-7>
- **COSMIC** (example `force_call_alleles` resource) — Tate JG, Bamford S, Jubb HC,
  et al. COSMIC: the Catalogue Of Somatic Mutations In Cancer. *Nucleic Acids Res*
  47(D1):D941–D947 (2019). <https://doi.org/10.1093/nar/gky1015>
- **1000 Genomes** (Peddy reference populations) — 1000 Genomes Project Consortium.
  A global reference for human genetic variation. *Nature* 526, 68–74 (2015).
  <https://doi.org/10.1038/nature15393>

