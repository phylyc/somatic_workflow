#!/bin/bash

#                                SAMPLES (files)         NAMES
bash call_genotype.sh test       "test_N,test_T"         "test_N,test_T"
bash call_genotype.sh test_empty "test_empty"            "test_empty"
bash call_genotype.sh test_dup   "test_N,test_T,test_T2" "test_N,test_T,test_T"
