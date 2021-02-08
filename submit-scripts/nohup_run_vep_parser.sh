#!/bin/bash

nohup python -m vep.vep_parser \
             --vcf_vep_path file:///home/ubuntu/data/chd_ukbb_vep/recal_snp_recal_indelgenome.sorted.test.vcf.bgz \
             --tb_output_path file:///home/ubuntu/data/hail_data/hts/chd_ukbb.vep.split.ht \
             --force_bgz \
             --write_to_file \
             --overwrite \
             --split_multi_allelic \
             > vep_parser.log 2>&1 &

