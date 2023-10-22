#!/bin/bash

validate_wdl () {
  local wdl=$1
  java -jar ../cromwell/womtool-85.jar validate $wdl
}

for wdl in wdl/*.wdl; do
  echo $wdl
  validate_wdl $wdl
done
