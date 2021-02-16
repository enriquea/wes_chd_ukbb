import argparse
import hail as hl
import os

"""
Simple script to generate a fully bi-allelic MatrixTable
from a MT contating multi-allelic variants. 

All Call fields are also split/downcoded.

See https://hail.is/docs/0.2/methods/genetics.html#hail.methods.split_multi_hts
for details.
"""

# init Hail
hl.init(default_reference='GRCh38')

# path to MatrixTable
mt_path = 'file:///home/ubuntu/data/hail_data/mts/chd_ukbb_08092020.mt'

# read matrix
mt = hl.read_matrix_table(mt_path)

# split multi-allelic variant into bi-allelic. It will also downcode the genotypes.
mt = hl.split_multi_hts(mt)

# write MT to disk
mt.write("file:///home/ubuntu/data/hail_data/mts/chd_ukbb_split_09092020.mt")

# stop Hail
hl.stop()
