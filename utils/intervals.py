# eam
# 2021-03-08

"""

A set of utils function to parse/import HailTable-like capture intervals from
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
from utils.generic import current_date


nfs_dir = 'file:///home/ubuntu/data'


# required global fields to be annotated in the interval HT
GLOBAL_ANNOTATION_FIELDS = ('date',  # date table was created
                            'source',  # bed file generating the interval HT
                            'reference_genome',  # genome build
                            'platform_label'  # unique label
                            )


def _get_interval_ht_path(platform_label: str,
                          genome_ref: str) -> str:
    """
    Generate path to write HT intervals.
    Note: `platform_label` and `genome_ref` are required.

    :param platform_label: Unique capture interval identifier (e.g. 'ssv3')
    :param genome_ref: genome reference
    :return: Path to HT intervals
    """
    return f"{nfs_dir}/resources/intervals/{platform_label}.intervals.{genome_ref}.ht"


def write_intervals_ht(ht_interval: hl.Table,
                       overwrite: bool = False) -> None:
    """
    Handle writing HT interval to disk. Check if required globals fields exists.

    :param ht_interval: HT intervals
    :param overwrite: Overwrite file if already exists.
    :return: None
    """

    assert all([f in ht_interval.index_globals() for f in GLOBAL_ANNOTATION_FIELDS])

    platform_label = ht_interval.platform_label.collect()[0]
    genome_ref = ht_interval.reference_genome.collect()[0]

    return ht_interval.write(_get_interval_ht_path(platform_label=platform_label,
                                                   genome_ref=genome_ref),
                             overwrite=overwrite)


def import_intervals_from_bed(bed_path: str,
                              platform_label: str,
                              genome_ref: str) -> hl.Table:
    """
    Handle importing BED files as intervals. Recode contig if necessary and
    annotate global meta-info.
    Note: `platform_label` and `genome_ref` are required, since these info
           will be used as global annotations.

    :param bed_path: Path to capture interval BED file
    :param platform_label: Unique capture interval identifier (e.g. 'ssv3')
    :param genome_ref: Either 'GRCh37' or 'GRCh38

    :return: HailTable keyed by interval
    """
    
    # genome references
    rg37 = hl.get_reference('GRCh37')
    rg38 = hl.get_reference('GRCh38')

    # dict contig recode from rg38 -> rg37.
    # only autosomal and sex chromosomes
    CONTIG_RECODING_HG38_TO_HG37 = {contig: contig.replace('chr', '')
                                    for contig in rg38.contigs[:24]}

    # dict contig recode from rg37 -> rg38.
    # only autosomal and sex chromosomes
    CONTIG_RECODING_HG37_TO_HG38 = {CONTIG_RECODING_HG38_TO_HG37.get(k): k
                                    for k in CONTIG_RECODING_HG38_TO_HG37.keys()}


    # Recode contig if chromosome field in BED file miss-match with genome reference.
    if genome_ref == 'GRCh37':
        contig_recoding = CONTIG_RECODING_HG38_TO_HG37
    elif genome_ref == 'GRCh38':
        contig_recoding = CONTIG_RECODING_HG37_TO_HG38
    else:
        contig_recoding = None

    ht_intervals = hl.import_bed(bed_path,
                                 reference_genome=genome_ref,
                                 contig_recoding=contig_recoding)

    global_ann_expr = dict(zip(GLOBAL_ANNOTATION_FIELDS,
                               (current_date(), bed_path, genome_ref, platform_label)))

    ht_intervals = (ht_intervals
                    .annotate_globals(**global_ann_expr)
                    .key_by('interval')
                    .repartition(100)
                    )
    return ht_intervals


def get_ssv2_intervals_ht(genome_ref: str = 'GRCh38') -> hl.Table:
    """
    Get exome capture intervals SureSelect V2 Clinical Exomes
    Source: S30409818_hs_<hg19 or hg38>/S30409818_Padded.bed

    :param genome_ref: One of GRCh37 or GRCh38
    :return: A Table of intervals
    """
    return hl.read_table(
        f'{nfs_dir}/resources/intervals/ssv2.intervals.{genome_ref}.ht'
    )


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
    Source: S03723314_Padded.bed

    :param genome_ref: One of GRCh37 or GRCh38 (lifted over)
    :return: A Table of intervals
    """
    return hl.read_table(
        f'{nfs_dir}/resources/intervals/ssv4.intervals.{genome_ref}.ht'
    )


def get_ssv5_intervals_ht(genome_ref: str = 'GRCh38') -> hl.Table:
    """
    Get exome capture intervals SureSelect V5
    Source: S04380110_hs_<hg19 or hg38>/S04380110_Padded.bed

    :param genome_ref: One of GRCh37 or GRCh38
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


def get_ssv5_idt_intervals_intersect_ht(genome_ref: str = 'GRCh38') -> hl.Table:
    """
    Get exome capture intervals IDT xGen (UKBB) intersected with SSV5
    with bad regions removed
    Source: Generated by HGI

    :param genome_ref: GRCh38
    :return: A Table of intervals
    """
    return hl.read_table(
        f'{nfs_dir}/resources/intervals/ssv5_idt_intersect.intervals.{genome_ref}.ht'
    )


def generate_interval_list_ht(genome_ref: str = 'GRCh38') -> hl.Table:
    """
    Generate a list of intervals (union)

    :return: A joint table (union) of intervals
    """
    
    intervals = [get_ssv2_intervals_ht(),
                 get_ssv3_intervals_ht(),
                 get_ssv4_intervals_ht(),
                 get_ssv5_intervals_ht(),
                 get_idt_xgen_intervals_ht()
                ]
    
    # get global annotation(s) from input tables
    sources = [t.source.collect()[0] for t in intervals]
    platform_labels = [t.platform_label.collect()[0] for t in intervals]
    
    global_ann_expr = dict(zip(GLOBAL_ANNOTATION_FIELDS,
                               (current_date(), 
                                sources, 
                                genome_ref, 
                                platform_labels)))
    
    # keep only the interval <key> field for all tables
    intervals = [ht.key_by('interval').select()
                 for ht in intervals]

    ht_interval = (hl.Table.union(*intervals)
                   .select_globals()
                   )
    
    ht_interval = ht_interval.annotate_globals(**global_ann_expr)
    
    assert ht_interval.key.interval.dtype == hl.dtype(f'interval<locus<{genome_ref}>>')

    return ht_interval
