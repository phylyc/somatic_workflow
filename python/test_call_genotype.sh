#!/bin/bash

individual_id=$1
samples_string=$2
names_string=$3
normal_sample=$4

# Convert delimited strings into arrays
IFS=',' read -r -a names <<< "$names_string"
IFS=',' read -r -a samples <<< "$samples_string"

output_dir="test_data/output/genotype/$individual_id"
pileup_dir="test_data/input"
segments_dir="test_data/input"
contamination_dir="test_data/input"

mkdir -p $output_dir


sample_args=""
for name in "${names[@]}" ; do
  sample_args="$sample_args--sample $name "
done
pileup_args=""
segments_args=""
contamination_args=""
for sample in "${samples[@]}" ; do
  pileup_args="$pileup_args--pileup $pileup_dir/$sample.pileup "
  segments_args="$segments_args--segments $segments_dir/$sample.segments "
  contamination_args="$contamination_args--contamination $contamination_dir/$sample.contamination "
done
normal_sample_arg=""
if [ ! "$normal_sample" == "" ]; then
  normal_sample_arg="--normal_sample $normal_sample"
fi

python -u genotype.py \
  --output_dir $output_dir \
  --patient $individual_id \
  --variant "test_data/input/test_reference.vcf" \
  $sample_args \
  $pileup_args \
  $segments_args \
  $contamination_args \
  $normal_sample_arg \
  --min_read_depth 10 \
  --min_genotype_likelihood 0.95 \
  --model "betabinom" \
  --format "GT:AD:DP:PL" \
  --threads 4 \
  --save_sample_genotype_likelihoods \
  --verbose \
  > $output_dir/$individual_id.log \
  2>&1
