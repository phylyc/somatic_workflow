#!/bin/bash

# conda activate wdl

validate_wdl () {
  local wdl=$1
  echo "> miniwdl check:"
  miniwdl check --no-quant-check $wdl
  echo "> womtool validate:"
  java -jar ../cromwell/womtool-85.jar validate $wdl
  echo ""
}

for wdl in wdl/*.wdl; do
  echo $wdl
  validate_wdl $wdl
done

# conda deactivate