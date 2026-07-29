# Outputs

### Per **sequencing run**

| Name                                | Type                 | Meaning                          |
| ----------------------------------- |----------------------| -------------------------------- |
| `callable_loci`                     | Array\[Array\[File]] | Regions passing coverage filters |
| `total_read_counts`                 | Array\[Array\[File]] | Counts in targets (pre-denoise)  |
| `denoised_total_copy_ratios`        | Array\[Array\[File]] | tCR after tangent-normalization  |
| `snppanel_allelic_pileup_summaries` | Array\[Array\[File]] | Allelic counts summaries         |

### Per **sample**

| Name                                           | Type               | Meaning                                                                          |
| ---------------------------------------------- |--------------------|----------------------------------------------------------------------------------|
| `sample_names_ordered`                         | Array\[String]     | Sample order used downstream                                                     |
| `harmonized_callable_loci`                     | Array\[File]       | Union of callable loci per seq run                                               |
| `harmonized_denoised_total_copy_ratios`        | Array\[File]       | Harmonized tCR per sample                                                        |
| `harmonized_snppanel_allelic_pileup_summaries` | Array\[File]       | Harmonized allelic summaries                                                     |
| `contamination_table`                          | Array\[File]       | Per-sample contamination estimates                                               |
| `af_segmentation_table`                        | Array\[File]       | Allelic copy-ratio segmentation (intervals annotated with minor allele fraction) |
| `allelic_pileup_summaries`                     | Array\[File]       | Per-sample allelic pileups                                                       |
| `aggregated_allelic_read_counts`               | Array\[File]       | Per target interval aggregated allelic counts                                    |
| `genotype_error_probabilities`                 | Array\[Float]      | Genotyping error metrics                                                         |
| `af_model_parameters` / `cr_model_parameters`  | Array\[File]       | Model params (AF & copy-ratio)                                                   |
| `called_copy_ratio_segmentation`               | Array\[File]       | Segmentation calls (amp/del)                                                     |
| `cr_plot`                                      | Array\[File]       | Segmentation plot                                                                |
| `acs_copy_ratio_segmentation`                  | Array\[File]       | Allelic copy-ratio segmentation                                                  |
| `acs_copy_ratio_skew`                          | Array\[Float]      | Allelic skew per sample                                                          |
| `annotated_somatic_variants` / `_idx`          | Array\[File]       | Annotated somatic VCF(+ index)                                                   |
| `absolute_*` (incl. `purity`, `ploidy`)        | Array\[File/Float] | ABSOLUTE outputs/metrics                                                         |

### Per **interval shard**

| Name                                                   | Type         | Meaning                                |
| ------------------------------------------------------ | ------------ |----------------------------------------|
| `raw_calls_mutect2_vcf_scattered` / `_idx_scattered`   | Array\[File] | Shard VCF + index                      |
| `raw_mutect2_stats_scattered`                          | Array\[File] | Mutect2 stats per shard                |
| `raw_mutect2_bam_out_scattered` / `_bai_out_scattered` | Array\[File] | BAMs with evidence of called mutations |
| `raw_mutect2_artifact_priors_scattered`                | Array\[File] | Artifact priors                        |

### **Patient** level

| Name                                                         | Type | Meaning                                        |
| ------------------------------------------------------------ | ---- |------------------------------------------------|
| `raw_snv_calls_vcf` / `_idx`                                 | File | Merged raw SNV calls                           |
| `mutect2_stats`                                              | File | Combined stats                                 |
| `orientation_bias`                                           | File | Read-orientation artifacts                     |
| `filtered_vcf` / `_idx`                                      | File | After statistical/hard/realignment filters     |
| `filtering_stats`                                            | File | Filtering summaries                            |
| `somatic_vcf` / `_idx`                                       | File | Final selected somatic calls                   |
| `germline_vcf` / `_idx`                                      | File | Germline calls retained                        |
| `rare_germline_alleles` / `_idx`                             | File | Rare germline sites used for CNV               |
| `somatic_calls_bam` / `_bai`                                 | File | Evidence BAM for selected somatic calls        |
| `gvcf` / `_idx`                                              | File | Genotyping VCF                                 |
| `snp_ref_counts` / `snp_alt_counts` / `snp_other_alt_counts` | File | Allele-specific count tables                   |
| `snp_sample_correlation`                                     | File | Cross-sample variant SNP genotype correlations |
| `modeled_segments`                                           | File | CNV modeled segments                           |
| `phylogic_*` (report/timing/graph/growth)                    | Array[File] | PhylogicNDT outputs, one entry per replicate (`phylogic_n_replicates`) |
| `ancestry_pred`                                              | String | Ancestry classification                      |
| `ancestry_prob`                                              | Float | Ancestry prediction confidence                |
| `ancestry_background_pca_table`                              | File | Peddy reference-population PCA projections (JSON) |
| `ancestry_pca_plot`                                          | File | Peddy PCA plot (PNG)                           |


All outputs can be provided as `Cache` input in which case the corresponding task to generate that output is skipped if possible.

If `Array[File] z2_Cache.absolute_acr_rdata` is supplied as input, the workflow skips straight to the clonal decomposition step. Use that together with `Array[Int] z2_Cache.absolute_solution` to extract the chosen purity/ploidy mode for each sample and run phylogenetic inference. (Other required cached inputs for this to run successfully are `Array[File] z2_Cache.acs_copy_ratio_segmentation` and `Array[Float] z2_Cache.acs_copy_ratio_skew`). See docs/06_absolute-review-and-rerun.md for more details

## Output review guide

The main deliverables are:

- `absolute_acr/tcr_segtab`
- `absolute_acr/tcr_maf`
- `contamination_table`
- `gvcf`
- `snp_sample_correlation`
- PhylogicNDT outputs, especially `phylogic_report`

Depending on the project, somatic and germline VCF outputs may also be primary deliverables.

### Important QC checks

#### 1. Sample identity / consistency

Inspect:

- `snp_sample_correlation`
- `snp_sample_correlation_min`

Suggested thresholds for further investigation:

- WES: investigate below `0.9`
- WGS: investigate below `0.99`

Values below these levels should prompt a review of whether all samples truly come from the same individual.

#### 2. Contamination

Inspect `contamination_table`

A contamination estimate above roughly `0.04` (4% contamination) should prompt closer review.

#### 3. Segmentation quality

Inspect `cr_plot` and related segmentation outputs.

Potential warning signs:

- extremely fragmented segmentation
- more than roughly `1000` segments
- noisy profiles with weak plateau structure

#### 4. ABSOLUTE plausibility

Inspect `absolute_acr_plot`

Potential warning signs:

- poor agreement between proposed copy states and visible segmental plateaus
- poor agreement between purity expectations and SNV VAF modes
- inconsistent whole-genome doubling interpretation across samples

#### 5. Normal sample behavior

If a designated normal sample shows visible aneuploidy in the copy-ratio plots, that does not automatically invalidate the run. In this workflow, a "normal" sample is only weakly privileged as germline normal; if the data support substantial tumor content, it can effectively behave like a tumor sample.

### Other intermediate or cache-style outputs

Many other workflow outputs are primarily intermediate data products or cacheable artifacts. They remain useful for debugging and reruns, but are not the first files most users should inspect.
