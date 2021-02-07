#!/usr/bin/env bash

nohup python -m sample_qc.platform_pca --mt_input_path ./testdata/hail_mts/chd_ukbb_split_chr20XY_testdata.mt \
                                       --ht_output_path ./testdata/hail_hts/platfrom_pca.tets.ht \
                                       --ht_intervals ./testdata/intervals/ssv5_idt_intersect.intervals.GRCh38.ht \
                                       --write_to_file \
                                       --overwrite \
                                       > platform_pca.log 2>&1 &
