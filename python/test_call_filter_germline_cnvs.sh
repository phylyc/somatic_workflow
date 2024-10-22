#!/bin/bash

individual_id=$1
samples_string=$2
names_string=$3
normal_sample=$4

# Convert delimited strings into arrays
IFS=',' read -r -a names <<< "$names_string"
IFS=',' read -r -a samples <<< "$samples_string"

output_dir="test_data/output/filter_germline_cnvs/$individual_id"
cr_dir="test_data/input"

mkdir -p $output_dir

sample_args=""
for name in "${names[@]}" ; do
  sample_args="$sample_args--sample $name "
done
cr_args=""
for sample in "${samples[@]}" ; do
  cr_args="$cr_args--copy_ratio $cr_dir/$sample.modelFinal.called.seg "
done
normal_sample_arg=""
if [ ! "$normal_sample" == "" ]; then
  normal_sample_arg="--normal_sample $normal_sample"
fi

python -u harmonize_copy_ratios.py \
  --output_dir $output_dir \
  $sample_args \
  $cr_args \
  $normal_sample_arg \
  --suffix ".modelFinal.called.germline_filtered.seg" \
  --column_names \
    CONTIG START END \
    NUM_POINTS_COPY_RATIO NUM_POINTS_ALLELE_FRACTION \
    LOG2_COPY_RATIO_POSTERIOR_10 LOG2_COPY_RATIO_POSTERIOR_50 LOG2_COPY_RATIO_POSTERIOR_90 \
    MINOR_ALLELE_FRACTION_POSTERIOR_10 MINOR_ALLELE_FRACTION_POSTERIOR_50 MINOR_ALLELE_FRACTION_POSTERIOR_90 \
    CALL \
  --column_types \
    str int int \
    int int \
    float float float \
    float float float \
    str \
  --agg_col NUM_POINTS_COPY_RATIO \
  --agg_func mean \
  --agg_col NUM_POINTS_ALLELE_FRACTION \
  --agg_func mean \
  --agg_col LOG2_COPY_RATIO_POSTERIOR_10 \
  --agg_func mean \
  --agg_col LOG2_COPY_RATIO_POSTERIOR_50 \
  --agg_func mean \
  --agg_col LOG2_COPY_RATIO_POSTERIOR_90 \
  --agg_func mean \
  --agg_col MINOR_ALLELE_FRACTION_POSTERIOR_10 \
  --agg_func mean \
  --agg_col MINOR_ALLELE_FRACTION_POSTERIOR_50 \
  --agg_func mean \
  --agg_col MINOR_ALLELE_FRACTION_POSTERIOR_90 \
  --agg_func mean \
  --agg_col CALL \
  --agg_func first \
  --filter_germline_calls \
  --threads 1 \
  --min_target_length 10 \
  --verbose \
  > $output_dir/$individual_id.log \
  2>&1
