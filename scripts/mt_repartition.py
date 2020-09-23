# Decrease the number of original partitions in the MatrixTable

import hail as hl

# Init Hail on cluster
hl.init(default_reference='GRCh38')

# read MT
mt = hl.read_matrix_table('file:///home/ubuntu/data/hail_data/mts/chd_ukbb_split_09092020.mt')

# downsample number of partitions
mt = mt.repartition(n_partitions=2000,
                    shuffle=False)
# write new MT
mt.write('file:///home/ubuntu/data/hail_data/mts/chd_ukbb_split_v2_09092020.mt')

# stop Hail
hl.stop()