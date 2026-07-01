# Ancestry Calling

## Goal

The ancestry calling step uses [Peddy](https://peddy.readthedocs.io/) to:

1. **Filter to autosomal SNPs** - Uses BCFtools to extract autosomal regions from the gVCF, removing sex chromosomes to avoid sex-specific biases
2. **Ancestry inference** - Compares samples against a 1000 Genomes reference populations - peddy first projects the query sample onto the Principal Components of a 1000 Genomes reference population using ~25,000 varaint sites. Next, it fits an SVM to predict ancestry
3. **Generate reports** - Produces QC plots and ancestry assignments

## Inputs

The ancestry calling stage requires:

- A germline VCF file (produced during genotyping)
- A reference genome build designation (e.g., "hg38" or "hg19")

## Outputs

Key outputs from ancestry calling include:

- `ancestry_pred` - Predicted ancestry super population label (e.g., "AFR", "EUR", "EAS", "SAS", "AMR")
- `ancestry_prob`- Confidence probability for the ancestry prediction
- `ancestry_background_pca_table` - JSON file with PCA projections for the reference population (2504 samples) from the 1000 genomes project
- `ancestry_pca_plot` - PNG image showing PCA plots (PC1 vs PC2 and PC1 vs PC3) of samples against reference

## When to use ancestry calling

Ancestry calling is useful for:

- **Sample verification**: Confirm samples are correctly labeled and unrelated
- **Population stratification**: Account for ancestry in downstream analyses
- **QC dashboards**: Monitor ancestry distribution across a cohort of samples

## When ancestry calling runs

Ancestry calling runs automatically as part of the genotyping step whenever a
germline gVCF is produced and the `genome_build` is `hg19` or `hg38` (peddy's
reference panels are only defined for those builds). 

## Interpreting ancestry results

- `ancestry_pred`: Three-letter ancestry code (EUR = European; AFR = African American; AMR = Ad-mixed American; EAS = East Asian; SAS = South Asian; UNKNOWN = `ancestry_prob` is below 0.65)
- `ancestry_prob`: Confidence in the prediction (0.0 to 1.0)
  - Typically > 0.9 is high confidence
  - < 0.7 may indicate admixed ancestry or non-reference populations

## References

- Peddy — Pedersen BS, Quinlan AR. Who's Who? Detecting and Resolving Sample Anomalies in Human DNA Sequencing Studies with Peddy. *Am J Hum Genet* 100:406–413 (2017). <https://doi.org/10.1016/j.ajhg.2017.01.017>
- 1000 Genomes (reference populations) — 1000 Genomes Project Consortium. A global reference for human genetic variation. *Nature* 526, 68–74 (2015). <https://doi.org/10.1038/nature15393>
- SAMtools / BCFtools (autosome filtering / VCF prep) — Danecek P, et al. *GigaScience* 10(2):giab008 (2021). <https://doi.org/10.1093/gigascience/giab008>

See [References](11_references.md) for the full bibliography.

