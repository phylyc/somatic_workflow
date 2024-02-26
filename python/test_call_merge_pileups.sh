#!/bin/bash

individual_id=$1
samples_string=$2
names_string=$3

# Convert delimited strings into arrays
IFS=',' read -r -a names <<< "$names_string"
IFS=',' read -r -a samples <<< "$samples_string"

output_dir="test_data/output/merge_pileups/$individual_id"
pileup_dir="test_data/input"

mkdir -p $output_dir


sample_args=""
for name in "${names[@]}" ; do
  sample_args="$sample_args--sample $name "
done
pileup_args=""
for sample in "${samples[@]}" ; do
  pileup_args="$pileup_args--pileup $pileup_dir/$sample.pileup "
done

python -u merge_pileups.py \
  --output_dir $output_dir \
  $sample_args \
  $pileup_args \
  --verbose \
  > $output_dir/$individual_id.log \
  2>&1
