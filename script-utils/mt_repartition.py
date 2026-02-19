# Decrease the number of original partitions in the MatrixTable

import hail as hl

from utils.config import NFS_DIR

# Init Hail on cluster
hl.init(default_reference='GRCh38')

# read MT
mt = hl.read_matrix_table(f'{NFS_DIR}/hail_data/mts/chd_ukbb_split_09092020.mt')

# downsample number of partitions
mt = mt.repartition(n_partitions=2000,
                    shuffle=False)
# write new MT
mt.write(f'{NFS_DIR}/hail_data/mts/chd_ukbb_split_v2_09092020.mt')

# stop Hail
hl.stop()