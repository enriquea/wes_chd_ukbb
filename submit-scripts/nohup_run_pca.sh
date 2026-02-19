#!/bin/bash

nohup python -m pipelines.pca \
                 -i 'chd_ukbb' \
                 -o 'file:///home/ubuntu/data/hail_data/sample_qc/pca_european_samples_03112021.ht' \
                 --n_pcs 10 \
                 --maf_threshold 0.01 \
                 --overwrite \
                 --write_to_file \
                 --default_reference 'GRCh38' \
                 > pca_pipeline.log 2>&1 &

