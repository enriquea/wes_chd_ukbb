#!/usr/bin/env bash

nohup python -m sample_qc.relatedness_inference \
                 --mt_input_path 'file:///home/ubuntu/data/hail_data/mts/chd_ukbb_split_v2_09092020.mt' \
                 --skip_compute_pc_relate \
                 --write_to_file \
                 --overwrite \
                 > relatedness_inference_29042021-v1.log 2>&1 &



