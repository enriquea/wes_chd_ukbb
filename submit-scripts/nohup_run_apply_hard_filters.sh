#!/usr/bin/env bash

nohup python -m sample_qc.apply_hard_filters \
                 --exome_cohort 'chd_ukbb' \
                 --overwrite \
                 --write_to_file \
                 > apply_hard_filters.log 2>&1 &