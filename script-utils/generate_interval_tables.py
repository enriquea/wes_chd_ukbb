"""Generate intervals HT from different exome capture products."""

import hail as hl

from utils.liftover import liftover_intervals
from utils.intervals import (import_intervals_from_bed,
                              generate_interval_list_ht,
                              write_intervals_ht)

from utils.config import NFS_DIR


def process_ssv2(nfs_dir: str) -> None:
    """Import and write SSV2 interval tables for hg37 and hg38."""
    ssv2_hg37 = import_intervals_from_bed(
        bed_path=f'{nfs_dir}/resources/intervals/SSV2/S30409818_hs_hg19/S30409818_Padded.bed',
        platform_label='ssv2',
        genome_ref='GRCh37')
    write_intervals_ht(ssv2_hg37, overwrite=True)

    ssv2_hg38 = import_intervals_from_bed(
        bed_path=f'{nfs_dir}/resources/intervals/SSV2/S30409818_hs_hg38/S30409818_Padded.bed',
        platform_label='ssv2',
        genome_ref='GRCh38')
    write_intervals_ht(ssv2_hg38, overwrite=True)


def process_ssv3(nfs_dir: str) -> None:
    """Import SSV3 hg37 intervals, write them, then liftover to hg38 and write."""
    ssv3_hg37 = import_intervals_from_bed(
        bed_path=f'{nfs_dir}/resources/intervals/SSV3/ddd_exome_v3_probes_fixed_plus100bp_nonredundant.bed',
        platform_label='ssv3',
        genome_ref='GRCh37')
    write_intervals_ht(ssv3_hg37, overwrite=True)

    ssv3_lifted = liftover_intervals(t=ssv3_hg37)
    write_intervals_ht(ssv3_lifted, overwrite=True)


def process_ssv4(nfs_dir: str) -> None:
    """Import SSV4 hg37 intervals, write them, then liftover to hg38 and write."""
    ssv4_hg37 = import_intervals_from_bed(
        bed_path=f'{nfs_dir}/resources/intervals/SSV4/S03723314/S03723314_Padded.bed',
        platform_label='ssv4',
        genome_ref='GRCh37')
    write_intervals_ht(ssv4_hg37, overwrite=True)

    ssv4_lifted = liftover_intervals(t=ssv4_hg37)
    write_intervals_ht(ssv4_lifted, overwrite=True)


def process_ssv5(nfs_dir: str) -> None:
    """Import and write SSV5 interval tables for hg37 and hg38."""
    ssv5_hg37 = import_intervals_from_bed(
        bed_path=f'{nfs_dir}/resources/intervals/SSV5/S04380110_hs_hg19/S04380110_Padded.bed',
        platform_label='ssv5',
        genome_ref='GRCh37')
    write_intervals_ht(ssv5_hg37, overwrite=True)

    ssv5_hg38 = import_intervals_from_bed(
        bed_path=f'{nfs_dir}/resources/intervals/SSV5/S04380110_hs_hg38/S04380110_Padded.bed',
        platform_label='ssv5',
        genome_ref='GRCh38')
    write_intervals_ht(ssv5_hg38, overwrite=True)


def process_idt_xgen(nfs_dir: str) -> None:
    """Import and write IDT xGen (UKBB) hg38 interval table."""
    idt_xgen = import_intervals_from_bed(
        bed_path=f'{nfs_dir}/resources/intervals/ukbb_capture_kit/xgen_plus_spikein.b38.bed',
        platform_label='idt_xgen',
        genome_ref='GRCh38')
    write_intervals_ht(idt_xgen, overwrite=True)


def build_interval_list(nfs_dir: str) -> hl.Table:
    """Generate the union interval list HT and write it to disk."""
    ht_intervals = generate_interval_list_ht()
    ht_intervals.write(
        f'{nfs_dir}/resources/intervals/list.intervals.GRCh38.ht',
        overwrite=True)
    return ht_intervals


def main() -> None:
    """Generate interval tables for all exome capture platforms."""
    hl.init()

    nfs_dir = NFS_DIR

    #### SSV2 ####
    process_ssv2(nfs_dir)

    #### SSV3 ####
    process_ssv3(nfs_dir)

    #### SSV4 ####
    process_ssv4(nfs_dir)

    #### SSV5 ####
    process_ssv5(nfs_dir)

    #### IDT_xGen ####
    process_idt_xgen(nfs_dir)

    #### Intervals union list ####
    build_interval_list(nfs_dir)


if __name__ == '__main__':
    main()
