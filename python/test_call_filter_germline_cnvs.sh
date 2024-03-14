#!/bin/bash

tumor_sample=$1
normal_sample=$2

output_dir="test_data/output/filter_germline_cnvs"
cr_dir="test_data/input"

mkdir -p $output_dir

python -u filter_germline_cnvs.py \
  --output $output_dir/$tumor_sample.germline_filtered.seg \
  --tumor_called_copy_ratio_segmentation $cr_dir/$tumor_sample.modelFinal.called.seg \
  --normal_called_copy_ratio_segmentation $cr_dir/$normal_sample.modelFinal.called.seg \
  --verbose \
  --threads 1 \
  --min_segment_length 10 \
  > $output_dir/$tumor_sample.log \
  2>&1
