#!/usr/bin/env bash

nohup python -u sample_qc.sex_imputation.py \
                 --mt_input_path 'file:///home/ubuntu/data/hail_data/mts/chd_ukbb_split_v2_09092020.mt' \
                 --ht_output_path 'hdfs://spark-master:9820/dir/chd_ukbb_samples_sex_annotation_24092020.ht' \
                 --write_to_file \
                 > nohup.log 2>&1 &

