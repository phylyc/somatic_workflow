#!/bin/bash

#                                        INDIVIDUAL SAMPLES (files)         NAMES
bash test_call_genotype.sh               test       "test_N,test_T"         "test_N,test_T"         "test_N"
bash test_call_genotype.sh               test_one   "test_N"                "test_N"                "test_N"
bash test_call_genotype.sh               test_empty "test_empty"            "test_empty"            ""
bash test_call_genotype.sh               test_emptyN "test_empty,test_N"     "test_empty,test_N"     "test_N"
bash test_call_genotype.sh               test_dup   "test_N,test_T,test_T2" "test_N,test_T,test_T"  "test_N"

bash test_call_merge_pileups.sh          test       "test_N,test_T"         "test_N,test_T"
bash test_call_merge_pileups.sh          test_one   "test_N"                "test_N"
bash test_call_merge_pileups.sh          test_empty "test_empty"            "test_empty"
bash test_call_merge_pileups.sh          test_dup   "test_N,test_T,test_T2" "test_N,test_T,test_T"
bash test_call_merge_pileups.sh          test_chr   "test_T.chr"            "test_T.chr"

bash test_call_harmonize_copy_ratios.sh  test       "test_N,test_T"         "test_N,test_T"
bash test_call_harmonize_copy_ratios.sh  test_one   "test_N"                "test_N"
bash test_call_harmonize_copy_ratios.sh  test_empty "test_empty"            "test_empty"
bash test_call_harmonize_copy_ratios.sh  test_dup   "test_N,test_T,test_T2" "test_N,test_T,test_T"
bash test_call_harmonize_copy_ratios.sh  test_none  "test_T,test_T3"        "test_T,test_T3"        # empty intersection of intervals

bash test_call_filter_germline_cnvs.sh   test       "test_N,test_T"         "test_N,test_T"         "test_N"
bash test_call_filter_germline_cnvs.sh   test_no_N  "test_N,test_T"         "test_N,test_T"         ""
bash test_call_filter_germline_cnvs.sh   test_one   "test_N"                "test_N"                "test_N"
bash test_call_filter_germline_cnvs.sh   test_empty "test_empty"            "test_empty"            ""
bash test_call_filter_germline_cnvs.sh   test_dup   "test_N,test_T,test_T2" "test_N,test_T,test_T"  "test_N"
bash test_call_filter_germline_cnvs.sh   test_none  "test_T,test_T3"        "test_T,test_T3"        ""        # empty intersection of intervals
