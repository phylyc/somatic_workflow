# Multi-sample somatic variant discovery (GATK4+)
The MultiSampleSomaticWorkflow (MSSW) is a pipeline for somatic analysis of tissue or liquid biopsy samples from cancer patients. The modality of the data is assumed to be DNA, supporting WES, WGS, and panel sequencing data.
This repository contains a WDL (Workflow Description Language) workflow using GATK4. The workflow covers preprocessing, single nucleotide variant (SNV) calling, copy number variation (CNV) calling, and clonal analysis.


## 📖 Documentation
Full documentation is published as a browsable site (with sidebar navigation and search):

**➡️ https://phylyc.github.io/somatic_workflow/**

You can also read the same pages directly in the [`docs/`](docs/) folder:

| Page | What it covers |
| --- | --- |
| [Getting started](docs/00_getting-started.md) | First-run walkthrough and orientation |
| [Supported use cases](docs/01_supported-use-cases.md) | What data and scenarios are supported |
| [Inputs](docs/02_inputs.md) | Full input reference and minimal run template |
| [Outputs](docs/03_outputs.md) | Reviewing and understanding outputs |
| [CNV calling](docs/04_cnv-calling.md) | Copy number variation calling details |
| [SNV calling](docs/05_snv-calling.md) | Single nucleotide variant calling details |
| [ABSOLUTE review and rerun](docs/06_absolute-review-and-rerun.md) | Manual solution selection and re-run |
| [PhylogicNDT](docs/07_phylogicndt.md) | Phylogenetic and clonal dynamics analysis |
| [Ancestry calling](docs/08_ancestry-calling.md) | Genetic ancestry inference |
| [Failure recovery](docs/09_failure-recovery.md) | Recovering from failed shards/tasks |
| [Reference resources](docs/10_resources.md) | Reference and resource bundles |
| [References](docs/11_references.md) | Citations for underlying tools |


## Why multi-sample analysis?
Multi-sample analysis improves sensitivity by borrowing evidence across samples from the same patient. This is especially useful in low-purity settings and for subclonal variants that may be weak in one sample but clearer in another.


## Quickstart
1) **Import the workflow**:
   - The workflow is available on [Dockstore](https://dockstore.org/workflows/github.com/phylyc/somatic_workflow/MultiSampleSomaticWorkflow:master) to import into e.g. a [Terra](https://app.terra.bio/) workspace.
   - Alternatively, clone this repository to run the workflow locally using [Cromwell](https://cromwell.readthedocs.io/en/stable/)
2) **Prepare inputs**: Prepare an input json file. See `docs/02_inputs.md` for a full description of inputs and for a minimal run template. For Terra, upload the json file under the INPUTS tab on the Terra workflow page to fill out the inputs (or fill out manually).
3) **Run 1**: Run the full workflow - CNV calling, SNV calling and Clonal Analysis with ABSOLUTE. If running locally using Cromwell:
   ```bash
      java -jar cromwell.jar run multi-sample_somatic_workflow.wdl --inputs multi-sample_somatic_workflow.inputs.json
   ```
4) **Manual review**: ABSOLUTE outputs multiple possible solutions (purity and ploidy estimates) per sample. These solutions can be found in the `absolute_acr_plot` output pdf. Review all solutions manually and select one per sample to pass into `z2_Cache.absolute_acr_solution`. See `docs/06_absolute-review-and-rerun.md` for a full guide.
5) **Populate cache inputs**: Prepare an updated version of the input json file with Cache inputs added. For Terra, upload this json under the INPUTS tab on workflow page to fill out the inputs.
6) **Run 2**: Run the workflow again with `z2_Cache.absolute_acr_solution` + cached outputs to complete ABSOLUTE extraction and downstream clonal analysis.
7) **Understanding outputs**: Refer to `docs/03_outputs.md` for a guide on reviewing and understanding outputs.

**For a comprehensive understanding of the workflow:**
1) First, please take a look at the Workflow Structure section below.
2) Please refer to the `docs` folder, starting with `docs/00_getting-started.md`


## Basic concepts
**Patient → samples → runs.** A patient has ≥1 sample(s), and a sample may have multiple sequencing runs. Runs are grouped by `sample_names` which are inferred from bam read group names if not provided. Matched normal optional.

## Common troubleshooting strategies and suggestions

### A few Mutect2 shards failed while most succeeded
- Retry problematic shards with more memory: set `z2_Cache.high_mem_shards` to the list of shard IDs and re-run (call-caching will reuse prior successes if enabled).
- If needed, increase `z3_Parameters.mutect2_high_mem_factor` (e.g., 3).
- If only some high-memory shards fail again, list those IDs in `z2_Cache.huge_mem_shards`; they use `z3_Parameters.mutect2_huge_mem_factor` (4x the default memory allocation by default).
- As a last resort, skip shards with `z2_Cache.skip_shards` (acknowledging lost calls in those regions).

Also see `docs/09_failure-recovery.md`

### Mutect2 consideration for complex regions
- Mutect2 `use_linked_de_bruijn_graph`, while increasing sensitivity, has trouble calling variants in complex regions. It is strongly recommended (necessary) to use it together with `recover_all_dangling_branches` (both turned on by default). Pre-calling with Mutect1 also helps.

### Other issues
- For assistance running workflows on GCP or locally, refer to the [GATK tutorial](https://gatk.broadinstitute.org/hc/en-us/articles/360035530952).
- Access necessary reference and resources bundles via the [GATK Resource Bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle). Also see `docs/10_resources.md`.



## Workflow Structure

The workflow is organized into the following main tasks:

### 1. Preprocessing
- **1.1 Preprocess and Split Calling Intervals** Prepare the intervals for scatter-gather tasks in shards for SNV calling.
- **1.2 Define Patient/Sample/SequencingRun sets**: Define a patient as a set of samples, and each sample as a set of sequencing runs. Sequencing runs are grouped by sample name.

### 2. Coverage Aggregation
- **2.1 Callable loci per sample**: Collect genomic regions per sample with sufficient coverage for SNV detection
- **2.2 Total read counts**: Collect total read counts in target intervals and perform denoising of total copy ratios via tangent-normalization.
- **2.3 Allelic read count**: Collect allelic read counts at common germline sites (SNP panel: `common_germline_alleles`).
- **2.4 Harmonization of target intervals and allelic counts**: Harmonize target intervals across samples, subset to intersection; merge read count data from multiple sequencing runs per sample.
- **2.5 Contamination estimation**: Estimate *out-of-patient* contamination in each sample.
- **2.6 First-pass copy ratio segmentation**: Get a prior single-sample allelic copy ratio segmentation for filtering SNVs (3.3) and genotyping of germline sites (4.1).
  - **2.6.1 Filter total read counts**: Filters total read counts based on first-pass segmentation (turned off by default).

### 3. SNV Calling
- **3.1 Mutect1 single-sample calling**: Use Mutect1 for single-sample mutation calling in tumor-normal mode if a matched normal sample is available, otherwise in tumor-only mode.
- **3.2 Mutect2 multi-sample calling**: Use Mutect2 for multi-sample mutation calling. Force-call alleles that were called via Mutect1.
- **3.3 Filter Variant Calls**: Annotate and select somatic vs germline vs artifactual variants based on various filters.
  - **3.3a Filter**: Apply statistical filters for sequencing artifacts, germline variants, read orientation bias, and contamination, among others.
  - **3.3b Hard Filter**: Apply hard filters based on base quality, mapping quality, fragment length, read depth, read orientation quality, position on the read, and population allele frequency.
  - **3.3c Realignment Filter** Filter based on realignment success (to hg38 or whichever reference is given by `realignment_bwa_mem_index_image`).
- **3.4 Annotate SNVs**: Annotate SNVs with functional information.
- **3.5 Tumor mutational burden (TMB) estimation**: (coming soon)

### 4. CNV Calling
- **4.1 SNP genotyping**: Genotype allelic count data at common (from 2.3) and rare (from 3.3) germline sites using evidence across all samples; harmonize loci across samples; phase HETs using allelic imbalance.
- **4.2 Multi-sample segmentation**: Segment denoised total copy ratios and allelic copy ratio across multiple samples
- **4.3 Per-sample total and allelic copy ratio inference**: Infer total and allelic copy ratios for each sample.
- **4.4 Per-sample event calling**: Call amplifications and deletions for each sample. 
- **4.5 Per-sample segmentation plotting**: Plot the segmented denoised copy ratios and allelic copy ratios for each sample.

### 5. Clonal Analysis
- **5.1 ABSOLUTE** Perform per-sample clonal analysis of the identified somatic variants and estimate tumor purity, ploidy, and absolute copy number.
- **5.2 ABSOLUTE extraction**: Extract results for one chosen solution (needs manual input).
  - **5.2a**: Rescue dropped segments and SNVs.
- **5.3 PhylogicNDT**: Build phylogenetic trees from all samples, perform growth kinetics, and timing analysis.

### 6. Ancestry Calling
- Infer the patient's genetic ancestry from the genotyped gVCF. Runs PCA against reference populations and reports the most likely ancestry assignment with a confidence probability. Supports hg38 and hg19.

Please remember to always review the intermediate results to ensure that the final results are as expected. Inappropriate filtering or parameter settings can lead to misleading output.

## Runtime
Runtime parameters are optimized for implementations on Google Cloud Platform (GCP) and HPC cluster with SLURM scheduler and for applications to whole exome sequencing data. Provide RAM/disk for Mutect2/CNV steps.


## FAQ
**Q: Can I run tumor-only (no matched normal)?**
- Yes, Mutect1/2 and CNV tumor-only calling are supported. Expect more rare germline calls in the somatic SNV call set and consider stronger filtering (e.g. GERMQ >= 40) for SNV burden tests. 

**Q: No CNV panel of normals (PoN) for my sequencing platform?**
- Provide a **0-byte** file for that entry. The workflow will fall back to `annotated_target_intervals` (using gc-content) for denoising. This is still provides decent results for WGS samples, but WES samples tend to have more batch structure. 

**Q: My targeted panels differ across samples. Is this OK?**
- Harmonization intersects targets across samples and mixed panels reduce effective target space. Prefer consistent panels for multi-sample CNV detection. Using highly dissimilar panels across samples will still work, but might lead to unreliable CNV calls and loss of SNV calls in regions that are not common between the panels.


## Software version requirements :
- **GATK**: Version 4.6.2.0. 
- **Cromwell**: Tested successfully on version 86.
- **Docker**: Docker tags listed in `docker/`.

## Citation
This workflow builds on GATK (Mutect/Mutect2, copy-ratio tools, Funcotator),
ABSOLUTE, PhylogicNDT and Peddy, among other tools. If you use it in a
publication, please cite the relevant methods listed in
[`docs/11_references.md`](docs/11_references.md).
