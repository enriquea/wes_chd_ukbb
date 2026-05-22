"""
Simple script to generate a fully bi-allelic MatrixTable
from a MT containing multi-allelic variants.

All Call fields are also split/downcoded.

See https://hail.is/docs/0.2/methods/genetics.html#hail.methods.split_multi_hts
for details.
"""

import hail as hl

from utils.config import NFS_DIR


def main() -> None:
    """Split multi-allelic variants in the CHD-UKBB MatrixTable into bi-allelic sites."""
    # init Hail
    hl.init(default_reference='GRCh38')

    # path to MatrixTable
    mt_path = f'{NFS_DIR}/hail_data/mts/chd_ukbb_08092020.mt'

    # read matrix
    mt = hl.read_matrix_table(mt_path)

    # split multi-allelic variant into bi-allelic. It will also downcode the genotypes.
    mt = hl.split_multi_hts(mt)

    # write MT to disk
    mt.write(f"{NFS_DIR}/hail_data/mts/chd_ukbb_split_09092020.mt")

    # stop Hail
    hl.stop()


if __name__ == '__main__':
    main()
