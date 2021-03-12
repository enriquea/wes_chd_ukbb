#!/usr/bin/env bash

nohup python -m sample_qc.platform_pca --mt_input_path 'file:///home/ubuntu/data/hail_data/mts/chd_ukbb_split_v2_09092020.mt' \
                                       --ht_output_path 'file:///home/ubuntu/data/hail_data/sample_qc/chd_ukbb.platform_pca.ht'  \
                                       --ht_intervals 'file:///home/ubuntu/data/resources/intervals/ssv3_ssv4_ssv5_idt.intervals_list.union.ht' \
                                       --write_to_file \
                                       --overwrite \
                                       > platform_pca.log 2>&1 &
