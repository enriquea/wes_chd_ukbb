# eam
# 2021-03-08

"""

A set of utils function to import HailTable-like capture intervals from
different sources.

Example table intervals:

----------------------------------------
Global fields:
    'date': str
    'source': str
----------------------------------------
Row fields:
    'interval': interval<locus<GRCh38>>
    'target': str
----------------------------------------
Key: ['interval']
----------------------------------------

"""

import hail as hl

nfs_dir = 'file:///home/ubuntu/data'


def get_ssv3_intervals_ht(genome_ref: str = 'GRCh38') -> hl.Table:
    """
    Get exome capture intervals SureSelect V3
    Source: ddd_exome_v3_probes_fixed_plus100bp_nonredundant.bed

    :param genome_ref: One of GRCh37 or GRCh38 (lifted over)
    :return: A Table of intervals
    """
    return hl.read_table(
        f'{nfs_dir}/resources/intervals/ssv3.intervals.{genome_ref}.ht'
    )


def get_ssv4_intervals_ht(genome_ref: str = 'GRCh38') -> hl.Table:
    """
    Get exome capture intervals SureSelect V4
    Source: ss_all_exon_covered_v4_plus100bp.bed

    :param genome_ref: One of GRCh37 or GRCh38 (lifted over)
    :return: A Table of intervals
    """
    return hl.read_table(
        f'{nfs_dir}/resources/intervals/ssv4.intervals.{genome_ref}.ht'
    )


def get_ssv5_intervals_ht(genome_ref: str = 'GRCh38') -> hl.Table:
    """
    Get exome capture intervals SureSelect V5
    Source: S04380110/S04380110_Padded.bed

    :param genome_ref: GRCh38
    :return: A Table of intervals
    """
    return hl.read_table(
        f'{nfs_dir}/resources/intervals/ssv5.intervals.{genome_ref}.ht'
    )


def get_idt_xgen_intervals_ht(genome_ref: str = 'GRCh38') -> hl.Table:
    """
    Get exome capture intervals IDT xGen (UKBB)
    Source: xgen_plus_spikein.b38.bed

    :param genome_ref: GRCh38
    :return: A Table of intervals
    """
    return hl.read_table(
        f'{nfs_dir}/resources/intervals/idt_xgen.intervals.{genome_ref}.ht'
    )


def generate_interval_list_ht() -> hl.Table:
    """
    Generate a list of intervals (union) from ssv3, ssv4, ssv5 and IDT xGen.

    :return: A joint table (union) of intervals
    """
    intervals = [get_ssv3_intervals_ht(),
                 get_ssv4_intervals_ht(),
                 get_ssv5_intervals_ht(),
                 get_idt_xgen_intervals_ht()]

    ht_interval = (hl.Table.union(*intervals)
                   .select_globals()
                   )
    assert ht_interval.key.interval.dtype == hl.dtype('interval<locus<GRCh38>>')

    return ht_interval

