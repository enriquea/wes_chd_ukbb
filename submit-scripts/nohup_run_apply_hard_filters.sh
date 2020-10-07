#!/usr/bin/env bash

nohup python -m sample_qc.apply_hard_filters \
                 --exome_cohort 'chd_ukbb' \
                 --skip_filter_step \
                 --overwrite \
                 --write_to_file \
                 > apply_hard_filters.log 2>&1 &