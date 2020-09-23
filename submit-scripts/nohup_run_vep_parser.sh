#!/bin/bash

nohup python data/scripts/vep_parser.py \
             --vcf_vep_path 'file:///home/ubuntu/data/chd_ukbb_vep/recal_snp_recal_indelgenome.sorted.vcf.gz' \
             --tb_output_path 'file:///home/ubuntu/data/hail_data/hts/chd_ukbb_split_vep_09092020.ht' \
             --force_bgz \
             --write_to_file \
             --keep_info_fields VQSLOD POSITIVE_TRAIN_SITE NEGATIVE_TRAIN_SITE \
             > nohup.log 2>&1 &

