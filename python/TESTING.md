# Python test suite

Assertion-based pytest suite for the eight scripts under `python/`. 

## Running

One-time environment (a venv with the runtime deps + pytest):

```bash
python -m venv --system-site-packages .venv-test
.venv-test/bin/pip install -r python/requirements.txt
```

Run from the `python/` directory:

```bash
cd python
../.venv-test/bin/pytest                       # everything
../.venv-test/bin/pytest -m unit               # fast pure-function tier only
../.venv-test/bin/pytest -m "not regression"   # skip the end-to-end script runs
../.venv-test/bin/pytest tests/test_genotype.py # one script
```

CI runs the same thing (`.github/workflows/python-tests.yml`, Python 3.11).

## Layout

All tests for a script live in **one** file, `tests/test_<script>.py`:

```
tests/
  conftest.py        # sys.path shim, the `synth`/`synth_chr` cohort fixtures, --update-golden
  helpers.py         # df builders, assert_frame_close, run_script, synthetic-chain runners
  test_<script>.py   # unit / io / stat tests + a "black-box" (synthetic cohort) section
  test_cli_smoke.py  # every script: `--help` + import (argparse/wiring guard)
  test_synthetic_generator.py  # self-tests for the cohort generator
  fixtures/          # tiny committed inputs (e.g. the map_to_absolute characterization seg)
  snapshots/         # committed map_to_absolute characterization outputs
  synthetic/         # the cohort generator (code); its data/ output is gitignored
```

Test tiers are pytest markers: `unit` (pure functions, no IO), `io` (file
parsing/writing), `stat` (numerical correctness of the statistical cores),
`regression` (end-to-end script runs against the synthetic cohort), `deep`
(needs external binaries; currently unused). Floats are compared with tolerances
(`assert_frame_close`, `pytest.approx`); ints / labels / flags are exact.
End-to-end runs use `--threads 1` for determinism.

## What is tested

Each script has the pure-logic tiers (`unit`/`io`/`stat`) plus a black-box
regression section that drives the real script on the synthetic cohort and checks
the output against the known ground truth.

| Script | unit/io/stat highlights | black-box (truth-based) highlights |
|---|---|---|
| `merge_pileups` | contig ordering, locus sort, header round-trip, count-sum & max-AF, dedup, depth filter | exact run-sum; population AF carried not recomputed; merged-depth boundary; empty & chr-prefixed inputs |
| `pileup_to_allelic_counts` | gvcf read, interval remap/aggregation, error-prob clip, het selection | alleles from gvcf; hom + het site handling; min-depth filter; determinism; sex-chromosome hets like autosomal |
| `acs_conversion` | sex→ploidy `tau`, `f` snapping, σ error-propagation, CNLOH labels, skew, min_hets/probes | output schema; f recovered from truth CN; skew from MEAN_BIAS; sex-chromosome presence; σ NA-iff-f-NA; probe/het/length pass-through |
| `genotype` | ploidies, error-prob, likelihoods (hom/het/argmax), phasing helpers, IO classes, joint likelihood, pop-AF prior | genotype invariant to sample subset & to `--normal`; phasing matches truth haplotype; count tables vs pileup; biallelic-SNV-only; empty-pileup sample degenerate not fatal |
| `harmonize_copy_ratios` | interval algebra, agg funcs, INTERSECTION/UNION, germline-CNV filter | preserves value on matching grid; intersection grid recovery; subset-contig dropped; empty/single/disjoint/partial cases; sorted-disjoint-unique output |
| `calculate_cancer_cell_fraction` | sex/chrom normalize, beta-binom vs scipy, CCF posterior/CI, power calc | ccf_hat vs truth (clonal/subclonal/absent); MT/unmapped/zero-cov dropped; per-sample isolation; CI calibration & width-vs-coverage; cross-sample clonal consistency; multiplicity |
| `map_to_absolute_copy_number` | extracted pure helpers (`chromosomal_ploidy`, `wmode`, `map_cn_to_cluster`, `split_alleles_for_segment`, `allele_ccf`, `segment_ccf`) + characterization snapshot | total/allelic/LOH/HZ/amp recovery; WGD=1/2 baselines; male/XXY/female sex-chromosome CN; allelic-sum identity; every ABSOLUTE segment present |
| `validate_inputs` | `parse_bool`, `sq_contigs`, ordered-subsequence, full rule-engine matrix | generated cohort manifest passes; bai-count mismatch rejected |

## The synthetic golden cohort

The regression tier is driven by one deterministic generator
(`tests/synthetic/make_synthetic_testing_data.py`). From a single ground-truth
tumor model it forward-simulates every upstream tool output (pileups, gVCFs,
denoised copy ratios, ModelSegments, ABSOLUTE segtabs, MAFs) **and** emits the
ground-truth answer key, so outputs are validated against known truth rather than
a blessed snapshot.

The cohort covers: a female (XX), a male (XY), and an XXY patient; WES-like,
ULP-WGS-like (shallow, no allelics), and deep-panel-like (no segmentation) assays;
samples with multiple sequencing runs; WGD counts 0/1/2; clonal & subclonal
events; and degenerate edge cases (empty / single / duplicate / no-normal inputs,
disjoint & partial copy-ratio grids, chr-prefixed contigs).

**The cohort is not committed** — the `synth` / `synth_chr` fixtures in
`conftest.py` regenerate it (seed 0) into a temp directory for every test session.
To produce an inspectable copy locally:

```bash
cd python
# plain contig names
python tests/synthetic/make_synthetic_testing_data.py --outdir tests/synthetic/data
# chr-prefixed variant
python tests/synthetic/make_synthetic_testing_data.py --outdir tests/synthetic/data_chr --chr
```

`--seed N` changes the RNG seed (tests use `0`). The output is plain TSV (gVCFs
plain text) so it diffs cleanly; `tests/synthetic/data/` is gitignored.
