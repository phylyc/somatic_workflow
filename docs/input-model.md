# Input model

The workflow uses a hierarchical biological and technical model:

- **patient** → one or more samples
- **sample** → one or more sequencing runs

Each entry in the run-level input arrays corresponds to one sequencing run. Arrays that describe sequencing runs must therefore have matching length and matching order.

## Core identifiers

### `patient_id`

A user-defined identifier for the case. This label is used in workflow outputs.

### `sample_names`

A user-defined name per sequencing run. Runs with the same sample name are grouped into one physical sample.

Whitespace is acceptable. In practice, avoid unusual shell-sensitive characters.

### `normal_sample_names`

A list of sample names that should be treated as normal samples. The first listed normal is treated as the matched normal when a single matched normal is needed.

If several normal samples exist, the preferred choice is typically the one with highest sequencing depth, although similarly sequenced normals often behave comparably.

## Timepoints

`timepoints` represent relative collection times for samples. The workflow is invariant to overall translation and rescaling of these values, so they may be represented in different units as long as the ordering and spacing make sense for the analysis.

Timepoints are mainly relevant for growth-kinetics interpretations and are especially useful in some cfDNA settings. They are generally less meaningful when comparing geographically distinct lesions or separate primaries.

## Target intervals

Each sequencing run must have a corresponding `target_intervals` file.

- For WES or targeted sequencing, these are assay target intervals, ideally padded and cleaned of problematic regions.
- For WGS, these are binned whole-genome interval definitions, again ideally curated for problematic regions.

## CNV panel of normals

A `cnv_panel_of_normals` entry may be supplied for each sequencing run.

Preferred practice:

- platform-specific PoN
- optionally sex-specific PoN if accurate chrX/chrY behavior matters

If no suitable CNV PoN is available for a run, an empty placeholder file may be used and the workflow will fall back to target-annotation-based denoising behavior using GC content.

## `is_paired_end`

Set this to `true` for paired-end platforms. This is especially important for cfDNA-like libraries where overlapping mates can otherwise bias read counting.

## `use_sample_for_tCR` and `use_sample_for_aCR`

These per-run/sample flags allow excluding samples from specific copy-ratio calculations.

Typical examples:

- if one biospecimen has both ULP WGS and WES, exclude the ULP WGS run from total-copy-ratio estimation when the WES signal is much more informative
- if WES and targeted-panel data are analyzed together, exclude panel data from total-copy-ratio estimation to avoid collapsing the signal to the narrow panel interval set

## Tumor-only input

Tumor-only analysis is represented simply by omitting `normal_sample_names` or leaving it empty.
