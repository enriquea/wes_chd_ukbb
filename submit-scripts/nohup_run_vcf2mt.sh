#!/bin/bash

nohup python -u data/scripts/vcf2mt.py \
      --vcf_path 'file:///home/ubuntu/data/chd_ukbb_vcfs/all_chromosomes/*.vcf.gz' \
      --output_path 'file:///home/ubuntu/data/hail_data/mts/chd_ukbb_08092020.mt' \
      --force_bgz \
      --overwrite \
      > nohup.log 2>&1 &

