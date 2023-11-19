#!/bin/bash

individual_id="test"
samples=( \
  "test_N" \
  "test_T" \
)
#samples=( \
#  "test_empty"
#)

output_dir="test_data/output"
pileup_dir="test_data"
segments_dir="test_data"
contamination_dir="test_data"

mkdir -p $output_dir


sample_args=""
pileup_args=""
segments_args=""
contamination_args=""
for sample in "${samples[@]}" ; do
  sample_args="$sample_args--sample $sample "
  pileup_args="$pileup_args--pileup $pileup_dir/$sample.pileup "
  segments_args="$segments_args--segments $segments_dir/$sample.segments "
  contamination_args="$contamination_args--contamination $contamination_dir/$sample.contamination "
done

python -u genotype.py \
  --output_dir $output_dir \
  --patient $individual_id \
  --variant "test_data/test_reference.vcf" \
  $sample_args \
  $pileup_args \
  $segments_args \
  $contamination_args \
  --min_read_depth 10 \
  --min_genotype_likelihood 0.95 \
  --model "betabinom" \
  --format "GT:AD:DP:PL" \
  --threads 4 \
  --save_sample_genotype_likelihoods \
  --verbose \
  > $output_dir/$individual_id.log \
  2>&1
