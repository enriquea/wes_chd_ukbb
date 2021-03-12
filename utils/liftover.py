# eam
# 2021-03-12

from gnomad.utils.liftover import *

import logging
import hail as hl

from utils.generic import current_date

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

hl.init()

project_dir = 'file:///home/ubuntu/data/projects/wes_chd_ukbb'


def liftover_intervals(t: hl.Table,
                       keep_missing_interval: bool = False) -> hl.Table:
    """
    Liftover locus in intervals from one coordinate system (hg37) to another (hg38)

    # Example input table description
    #
    # ----------------------------------------
    # Global fields:
    #     None
    # ----------------------------------------
    # Row fields:
    #     'interval': interval<locus<GRCh37>>
    # ----------------------------------------
    # Key: ['interval']
    # ----------------------------------------


    :param t: Table of intervals on GRCh37
    :param keep_missing_interval: If True, keep missing (non-lifted) intervals in the output Table.
    :return: Table with intervals lifted over GRCh38 added.
    """

    rg37 = hl.get_reference("GRCh37")
    rg38 = hl.get_reference("GRCh38")

    if not rg37.has_liftover("GRCh38"):
        rg37.add_liftover(
            f'{project_dir}/data/resources/grch37_to_grch38.over.chain.gz', rg38
        )

    t = t.annotate(
        start=hl.liftover(t.interval.start, "GRCh38"),
        end=hl.liftover(t.interval.end, "GRCh38"),
    )

    t = t.filter(
        (t.start.contig == "chr" + t.interval.start.contig)
        & (t.end.contig == "chr" + t.interval.end.contig)
    )

    t = t.key_by()

    t = (t
         .select(interval=hl.locus_interval(t.start.contig,
                                            t.start.position,
                                            t.end.position,
                                            reference_genome=rg38,
                                            invalid_missing=True),
                 interval_hg37=t.interval
                 )
         )

    # bad intervals
    missing = t.aggregate(hl.agg.counter(~hl.is_defined(t.interval)))
    logger.info(
        f"Number of missing intervals: {missing[True]} out of {t.count()}..."
    )

    # update globals annotations
    global_ann_expr = {'date': current_date(),
                       'reference_genome': 'GRCh38',
                       'was_lifted': True
                       }
    t = t.annotate_globals(**global_ann_expr)

    if not keep_missing_interval:
        logger.info(
            f"Filtering out {missing[True]} missing intervals..."
        )
        t = t.filter(hl.is_defined(t.interval), keep=True)

    return t.key_by("interval")

