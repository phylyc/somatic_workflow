# Inputs

The workflow uses a hierarchical biological and technical model:

- **patient** → one or more samples
- **sample** → one or more sequencing runs

## Inputs table

The input arguments are grouped into 
```wdl
MultiSampleSomaticWorkflow # (core patient-related inputs)
z1_Files # (reference data resources)
z2_Cache # (workflow output)
z3_Parameters # (workflow parameters)
z4_RuntimeParameters # (per-task runtime settings)
```
> **Note** - Please note that the "Required" tags are specified for a full, sensible run. Some of these might not be required in certain use cases - for example, if users choose to only run certain parts of the workflow or to pass in intermediate outputs using Cache input fields.

> **Note on the `z1_`–`z4_` prefixes.** These groups are sub-workflow calls, and
> the numeric prefixes are only there to force a sensible top-to-bottom ordering in
> the Terra inputs UI. Use the prefixed names in input keys, e.g.
> `MultiSampleSomaticWorkflow.z1_Files.ref_fasta` and
> `MultiSampleSomaticWorkflow.z2_Cache.absolute_solution`.

### A) Core input (`MultiSampleSomaticWorkflow`)

| Variable Name | Required/Recommended/Optional | Description | Type |
|----------------|-------------------------|----------------|----------------|
| `patient_id` | Required | Patient identifier | String |
| `sample_names` | Required | Sample identifiers for each BAM, ideally with no spaces. Runs with the same sample name are grouped into one physical sample. Whitespace is acceptable. In practice, avoid unusual shell-sensitive characters. | Array[String] |
| `bams` | Required | Sample's BAM file | Array[File] |
| `bais` | Required | Sample's BAI file | Array[File] |
| `target_intervals` | Required | Per-run capture kit intervals used for CNV read-count collection and allelic pileups. For WES, these are the assay target regions (ideally padded and blacklist-removed). For WGS, use binned chromosomal intervals. Different runs may use different interval files if sequenced on different platforms. | Array[File] |
| `annotated_target_intervals` | Optional | Target intervals annotated with GC content, mappability, and segmental duplications | Array[File] |
| `cnv_panel_of_normals` | Recommended | **(Highly recommended)** Denoised total read counts by filtering out sequencing platform artifacts. If not available, needs to point to an empty file (0B) and `annotated_target_intervals` will be used to denoise. Refer to docs/04_cnv-calling.md for more details on creation. | Array[File] |
| `is_paired_end` | Recommended | Whether paired-end sequencing was used to process the sample | Array[Boolean] |
| `use_sample_for_tCR` | Recommended | Whether to use the sample for denoised copy ratio (dCR) estimation. Example for when to change: If one biospecimen has both ULP WGS and WES, exclude the ULP WGS run from total-copy-ratio estimation when the WES signal is much more informative. Or, if WES and targeted-panel data are analyzed together, exclude panel data from total-copy-ratio estimation to avoid collapsing the signal to the narrow panel interval set | Array[Boolean] |
| `use_sample_for_aCR` | Recommended | Whether to use the sample for allelic copy ratio (aCR) estimation | Array[Boolean] |
| `normal_sample_names` | Recommended | A list of sample names that should be treated as normal samples. The first listed normal is treated as the matched normal when a single matched normal is needed. If several normal samples exist, the preferred choice is typically the one with highest sequencing depth, although similarly sequenced normals often behave comparably. For tumor-only analysis, leave `normal_sample_names` empty or omit it. | Array[String] |
| `sex` | Required | Patient sex genotype (female, male or XX, XY) | String |
| `timepoints` | Recommended | Relative collection times for samples. The workflow is invariant to overall translation and rescaling of these values, so they may be represented in different units as long as the ordering and spacing make sense. Used for phylogenetic analysis. Especially useful in cfDNA settings. | Array[String] |
| `total_mean_read_depth` | Recommended | Total mean read depth across all bams. Used to dynamically split the genome into shards for variant calling with Mutect2 to keep overall runtime reasonable. | Int |

### B) Reference data (`z1_Files`)

- Refer to [resources](10_resources.md) for links to available public resources

| Variable Name | Required/Recommended/Optional | Description | Type |
|----------------|-------------------------|----------------|----------------|
| `ref_fasta` | Required | Reference FASTA (hg38 or b37) | File |
| `ref_fasta_index` | Required | Reference FASTA index (hg38 or b37) | File |
| `ref_fasta_dict` | Required | Reference FASTA dictionary (hg38 or b37) | File |
| `interval_list` / `interval_lists` | Recommended | Genomic intervals for SNV calling. Preprocessed into scatter-gather shards for Mutect1/Mutect2. Supplying either `interval_list` or `preprocessed_intervals` is recommended for WES / panel sequencing. Otherwise, SNV calling is performed across the full genome. | File |
| `interval_blacklist` | Optional | Regions to subtract from `interval_list` before variant calling (e.g. centromeres, low-complexity) | File |
| `preprocessed_intervals` | Recommended | Pre-computed output of `PreprocessIntervals` using `interval_list`, `interval_lists`, and `interval_blacklist` (interval_list already binned/padded/blacklist-removed). If supplied, skips preprocessing. Supplying either `interval_list` or `preprocessed_intervals` is recommended for WES / panel sequencing. | File |
| `force_call_alleles` | Recommended | VCF file of sites to force mutation calling at (e.g. COSMIC or known driver mutations) | File |
| `force_call_alleles_idx` | Recommended | VCF index file for force_call_alleles | File |
| `funcotator_data_sources_tar_gz` | Recommended | Tarball of data sources for Funcotator. If not provided, the tarball will be automatically downloaded from the GATK resource bundle (slower) | File |
| `funcotator_transcript_list` | Recommended | List of transcript names to use for Funcotator annotation | File |
| `germline_resource` | Recommended | VCF file with AF field for annotating/filtering germline alleles (gnomAD) | File |
| `germline_resource_idx` | Recommended | VCF index file for germline_resource | File |
| `germline_resource_v4_1` | Optional | Same as germline_resource in VCF v4.1 format (for Mutect1) | File |
| `germline_resource_v4_1_idx` | Optional | Same as germline_resource_idx in VCF v4.1 format (for Mutect1) | File |
| `common_germline_alleles` | Recommended | VCF file of common biallelic germline alleles (e.g. population allele frequency \> 5%) to collect allelic counts for aCR and contamination estimation | File |
| `common_germline_alleles_idx` | Recommended | VCF index file for common_germline_alleles | File |
| `snv_panel_of_normals` | Recommended | VCF file of common SNV sequencing artifacts to filter out; ideally matches the sequencing platform(s) of the samples | File |
| `snv_panel_of_normals_idx` | Recommended | VCF index file for snv_panel_of_normals | File |
| `snv_panel_of_normals_v4_1` | Optional | Same as snv_panel_of_normals in VCF v4.1 format (for Mutect1) | File |
| `snv_panel_of_normals_v4_1_idx` | Optional | Same as snv_panel_of_normals_idx in VCF v4.1 format (for Mutect1) | File |
| `realignment_bwa_mem_index_image` | Optional | Only needed if realigning to a different genome (e.g. hg38) | File |

### C) Workflow Parameters (`z3_Parameters`)

For defaults, please check out: <https://github.com/phylyc/somatic_workflow/blob/master/wdl/workflow_arguments.wdl>

| Argument | Required/Recommended/Optional | Description / Example | Type |
|----------------|-----------------------|------------------|----------------|
| `genome_build` | Required | "hg19" or "hg38" | String |
| `high_mem_shards` | Optional | Indices of shards that require additional memory allocation, e.g. `[0, 5, 99, 102]` | Array[Int] |
| `mutect2_high_mem_factor` | Optional | Multiplier applied to the default memory allocation for high_mem_shards, e.g. `2.0` will use 2x memory | Float |
| `huge_mem_shards` | Optional | Indices of shards that require the highest memory tier. This tier takes precedence if a shard is also listed in `high_mem_shards` | Array[Int] |
| `mutect2_huge_mem_factor` | Optional | Multiplier applied to the default memory allocation for `huge_mem_shards`; defaults to `4.0` | Float |
| `scatter_count_base_for_variant_calling` | Optional | Number of shards to split the BAMs into for parallelized variant calling | Int |
| `run_reorder_bam_contigs` | Optional | true or false | Boolean |
| `run_collect_callable_loci` | Optional | true or false | Boolean |
| `run_collect_allelic_read_counts` | Optional | true or false | Boolean |
| `run_collect_total_read_counts` | Optional | true or false | Boolean |
| `run_contamination_model` | Optional | true or false | Boolean |
| `run_model_segments` | Optional | true or false | Boolean |
| `run_orientation_bias_mixture_model` | Optional | true or false | Boolean |
| `run_variant_annotation` | Optional | true or false | Boolean |
| `run_variant_annotation_scattered` | Optional | true or false | Boolean |
| `run_variant_calling` | Optional | true or false | Boolean |
| `run_variant_calling_mutect1` | Optional | true or false | Boolean |
| `run_variant_filter` | Optional | true or false | Boolean |
| `run_realignment_filter` | Optional | true or false | Boolean |
| `run_clonal_decomposition` | Optional | true or false | Boolean |

### D) Runtime parameters (`z4_RuntimeParameters`)

The runtime input group intentionally exposes only settings that have demonstrated case-specific operational value. These include the Mutect1 and Mutect2 memory/CPU/retry controls, GetPileupSummaries, HarmonizeCopyRatios, genotype-variants, Funcotator, PhylogicNDT, panel-of-normals construction, dynamic scatter dimensions, and a small set of shared retry/disk/preemptibility controls.

`JUST_RUN_IM_WILLING_TO_PAY` is an explicit user policy switch for prioritizing successful completion over compute cost; `mem_BIG_MACHINE` controls the corresponding high-memory target. Individual ABSOLUTE runtime settings, Docker images, startup and machine-overhead constants, and other stable task settings are fixed in `wdl/runtime_collection.wdl`. All `time_*` values are source-controlled rather than Terra inputs because they are used for HPC scheduler limits; an HPC deployment should adjust them in its checked-out workflow source.

---

## Minimal example JSON

- This example shows a matched tumor-normal WES-style structure. Replace the placeholder file paths with real locations. 
- Note that this is a minimal input template. Review other docs and add inputs as necessary
- Note that the below json files have hard-coded filenames. If you are running the workflow from a Terra data table, replace the relevant file paths with Terra column names, for example, you might replace `["file1.bam", "file2.bam"]` with `this.bams` where the `bams` column in your Terra data table contains an `Array[File]` of your bam files.

```json
{
  "MultiSampleSomaticWorkflow.patient_id": "Patient",
  "MultiSampleSomaticWorkflow.sex": "XX",
  "MultiSampleSomaticWorkflow.bams": [
    "file1.bam",
    "file2.bam"
  ],
  "MultiSampleSomaticWorkflow.bais": [
    "file1.bai",
    "file2.bai"
  ],
  "MultiSampleSomaticWorkflow.sample_names": [
    "Sample1",
    "Sample2"
  ],
  "MultiSampleSomaticWorkflow.normal_sample_names": [
    "Sample2"
  ],
  "MultiSampleSomaticWorkflow.target_intervals": [
    "WES_target_intervals.interval_list",
    "WES_target_intervals.interval_list"
  ],
  "MultiSampleSomaticWorkflow.cnv_panel_of_normals": [
    "WES_pon.hdf5",
    "WES_pon.hdf5"
  ],
  "MultiSampleSomaticWorkflow.z1_Files.ref_dict": "gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.dict",
  "MultiSampleSomaticWorkflow.z1_Files.ref_fasta": "gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta",
  "MultiSampleSomaticWorkflow.z1_Files.ref_fasta_index": "gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai",
  "MultiSampleSomaticWorkflow.z1_Files.germline_resource": "gnomad_sites.af_only.vcf.gz",
  "MultiSampleSomaticWorkflow.z1_Files.germline_resource_idx": "gnomad_sites.af_only.vcf.gz.tbi",
  "MultiSampleSomaticWorkflow.z1_Files.common_germline_alleles": "gnomad_sites.AFgt0.05.vcf.gz",
  "MultiSampleSomaticWorkflow.z1_Files.common_germline_alleles_idx": "gnomad_sites.AFgt0.05.vcf.gz.tbi",
  "MultiSampleSomaticWorkflow.z1_Files.snv_panel_of_normals": "snv_pon.vcf.gz",
  "MultiSampleSomaticWorkflow.z1_Files.snv_panel_of_normals_idx": "snv_pon.vcf.gz.tbi",
  "MultiSampleSomaticWorkflow.z1_Files.realignment_bwa_mem_index_image": "reference.fasta.img"
}
```
