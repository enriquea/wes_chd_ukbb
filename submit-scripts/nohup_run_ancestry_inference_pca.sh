#!/bin/bash

nohup python -m sample_qc.ancestry_inference \
                 --exome_cohort 'chd_ukbb' \
                 --n_pcs 20 \
                 --overwrite \
                 --write_to_file \
                 --default_reference 'GRCh38' \
                 > pca_ancestry_inference.log 2>&1 &

