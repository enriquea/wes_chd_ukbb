#!/usr/bin/env bash

nohup python  ../sample_qc/relatedness_inference.py \
                 --mt_input_path 'file:///home/ubuntu/data/hail_data/mts/chd_ukbb_split_v2_09092020.mt' \
                 --ht_output_path 'hdfs://spark-master:9820/dir/chd_ukbb_samples_relatedness_kinship_24092020.ht' \
                 --ld_pruning \
                 --write_to_file \
                 > nohup.log 2>&1 &



