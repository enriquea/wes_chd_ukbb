#!/bin/bash

nohup python -u data/scripts/ancestry_inference_pca.py \
                 --mt_input_path 'file:///home/ubuntu/data/hail_data/mts/chd_ukbb_split_v2_09092020.mt' \
                 --mt_input_1kg 'file:///home/ubuntu/data/hail_data/mts/phase3_1kg_snp_biallelic_hg38.mt' \
                 --ht_interval 'file:///home/ubuntu/data/resources/intervals/purcell_5k_intervals/purcell5k.ht' \
                 --ht_output_path 'hdfs://spark-master:9820/dir/chd_ukb_purcell5k_1kg_merged_pca_maf0.05_14092020.ht' \
                 --write_to_file \
                 > nohup.log 2>&1 &

