#!/bin/bash

nohup python -u data/scripts/sample_qc.py \
                 --mt_input_path 'file:///home/ubuntu/data/hail_data/mts/chd_ukbb_split_v2_09092020.mt' \
                 --ht_output_path 'hdfs://spark-master:9820/dir/chd_ukbb_split_v2_09092020_raw' \
                 --write_to_file \
                 > nohup.log 2>&1 &

