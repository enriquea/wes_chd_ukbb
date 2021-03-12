#!/usr/bin/env bash

nohup python  ../sample_qc/relatedness_inference.py \
                 --mt_input_path 'file:///home/ubuntu/data/hail_data/mts/chd_ukbb_split_v2_09092020.mt' \
                 --ht_output_path 'hdfs://spark-master:9820/dir/hail_data/chd_ukbb.relatedness_kinship.ht' \
                 --skip_filter_data \
                 --skip_prune_ld \
                 --write_to_file \
                 --overwrite \
                 > relatedness_inference_maf1.log 2>&1 &



