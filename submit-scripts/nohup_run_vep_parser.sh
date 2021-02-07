#!/bin/bash

nohup python -m vep.vep_parser \
             --vcf_vep_path ./testdata/vcf/recal_snp_recal_indelgenome.sorted.test.vcf.bgz \
             --tb_output_path ./testdata/hail_hts/vep.ht \
             --force_bgz \
             --write_to_file \
             --overwrite \
             --split_multi_allelic \
             > vep_parser.log 2>&1 &

