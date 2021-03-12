"""

Generate intervals HT from different exome capture products

"""

import hail as hl

from utils.liftover import liftover_intervals
from utils.intervals import (import_intervals_from_bed,
                             generate_interval_list_ht,
                             write_intervals_ht)


hl.init()

nfs_dir = 'file:///home/ubuntu/data'

#### SSV2 ####

# SureSelect Clinical Exomes V2 (hg37)
ssv2_hg37 = import_intervals_from_bed(
    bed_path=f'{nfs_dir}/resources/intervals/SSV2/S30409818_hs_hg19/S30409818_Padded.bed',
    platform_label='ssv2',
    genome_ref='GRCh37')
write_intervals_ht(ssv2_hg37, overwrite=True)

# SureSelect Clinical Exomes V2 (hg38)
ssv2_hg38 = import_intervals_from_bed(
    bed_path=f'{nfs_dir}/resources/intervals/SSV2/S30409818_hs_hg38/S30409818_Padded.bed',
    platform_label='ssv2',
    genome_ref='GRCh38')
write_intervals_ht(ssv2_hg38, overwrite=True)

#### SSV3 ####

# SureSelect Exomes V3 (hg37)
ssv3_hg37 = import_intervals_from_bed(
    bed_path=f'{nfs_dir}/resources/intervals/SSV3/ddd_exome_v3_probes_fixed_plus100bp_nonredundant.bed',
    platform_label='ssv3',
    genome_ref='GRCh37')
write_intervals_ht(ssv3_hg37, overwrite=True)

# SureSelect Exomes V3 (lifted over hg38)
ssv3_lifted = liftover_intervals(t=ssv3_hg37)
write_intervals_ht(ssv3_lifted, overwrite=True)


#### SSV4 ####

# SureSelect Exomes V4 (hg37)
ssv4_hg37 = import_intervals_from_bed(
    bed_path=f'{nfs_dir}/resources/intervals/SSV4/S03723314/S03723314_Padded.bed',
    platform_label='ssv4',
    genome_ref='GRCh37')
write_intervals_ht(ssv4_hg37, overwrite=True)

# SureSelect Exomes V4 (lifted over hg38)
ssv4_lifted = liftover_intervals(t=ssv4_hg37)
write_intervals_ht(ssv4_lifted, overwrite=True)


#### SSV5 ####

# SureSelect Exomes V5 (hg37)
ssv5_hg37 = import_intervals_from_bed(
    bed_path=f'{nfs_dir}/resources/intervals/SSV5/S04380110_hs_hg19/S04380110_Padded.bed',
    platform_label='ssv5',
    genome_ref='GRCh37')
write_intervals_ht(ssv5_hg37, overwrite=True)

# SureSelect Exomes V5 (hg38)
ssv5_hg38 = import_intervals_from_bed(
    bed_path=f'{nfs_dir}/resources/intervals/SSV5/S04380110_hs_hg38/S04380110_Padded.bed',
    platform_label='ssv5',
    genome_ref='GRCh38')
write_intervals_ht(ssv5_hg38, overwrite=True)


#### IDT_xGen ####

# IDT xGen (UKBB) (hg38)
idt_xgen = import_intervals_from_bed(
    bed_path=f'{nfs_dir}/resources/intervals/ukbb_capture_kit/xgen_plus_spikein.b38.bed',
    platform_label='idt_xgen',
    genome_ref='GRCh38')
write_intervals_ht(idt_xgen, overwrite=True)

#### Intervals union list ####

ht_intervals = generate_interval_list_ht()
ht_intervals.write(
    f'{nfs_dir}/resources/intervals/list.intervals.GRCh38.ht',
    overwrite=True)


