#!/usr/bin/env bash


nohup python -m sample_qc.sex_chrom_coverage \
                 --exome_cohort 'chd_ukbb' \
                 --interval_to_include 'path/to/ht/interval' \
                 --overwrite \
                 --write_to_file \
                 > sex_chrom_coverage.log 2>&1 &