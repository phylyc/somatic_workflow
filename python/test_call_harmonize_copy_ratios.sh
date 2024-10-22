#!/bin/bash

individual_id=$1
samples_string=$2
names_string=$3

# Convert delimited strings into arrays
IFS=',' read -r -a names <<< "$names_string"
IFS=',' read -r -a samples <<< "$samples_string"

output_dir="test_data/output/harmonize_copy_ratios/$individual_id"
cr_dir="test_data/input"

mkdir -p $output_dir


sample_args=""
for name in "${names[@]}" ; do
  sample_args="$sample_args--sample $name "
done
cr_args=""
for sample in "${samples[@]}" ; do
  cr_args="$cr_args--copy_ratio $cr_dir/$sample.denoised_CR.tsv "
done

python -u harmonize_copy_ratios.py \
  --output_dir $output_dir \
  $sample_args \
  $cr_args \
  --suffix ".harmonized.CR.tsv" \
  --verbose \
  --threads 4 \
  --min_target_length 2 \
  > $output_dir/$individual_id.log \
  2>&1
