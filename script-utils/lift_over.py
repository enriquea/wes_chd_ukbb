# eam
# 2021-02-16

from gnomad.utils.liftover import *

import logging
import hail as hl


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

hl.init()

project_dir = 'file:///home/ubuntu/data/projects/wes_chd_ukbb'

# prepare genome references
rg37 = hl.get_reference('GRCh37')
rg38 = hl.get_reference('GRCh38')
rg37.add_liftover(f'{project_dir}/data/resources/grch37_to_grch38.over.chain.gz', rg38)


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

    if not keep_missing_interval:
        logger.info(
            f"Filtering out {missing[True]} missing intervals..."
        )
        t = t.filter(hl.is_defined(t.interval), keep=True)

    return t.key_by("interval")


# lift over german AFs (hg37 -> hg38)
path_ht_ger_af_hg37 = f'{project_dir}/data/resources/german_pop_af.ht'
output_path = f'{project_dir}/data/resources/german_af_hg38'

lift_data(t=hl.read_table(path_ht_ger_af_hg37),
          gnomad=False,
          data_type='exomes',
          rg=rg38,
          path=output_path,
          overwrite=False)


# lift over RUMC exomes AFs (hg37 -> hg38)
path_ht_rumc_af_hg37 = f'{project_dir}/data/resources/rumc_af_18032020.ht'
output_path = f'{project_dir}/data/resources/rumc_af_hg38'

lift_data(t=hl.read_table(path_ht_rumc_af_hg37),
          gnomad=False,
          data_type='exomes',
          rg=rg38,
          path=output_path,
          overwrite=False)

# lift over SureSelect V3 intervals
ssv3 = hl.import_bed(f"{project_dir}/data/resources/intervals/ddd_exome_v3_probes_fixed_plus100bp_nonredundant.bed")
ssv3_lifted = liftover_intervals(t=ssv3)
ssv3_lifted.show()
(ssv3_lifted
 .repartition(100)
 .write(f"{project_dir}/data/resources/intervals/ssv3.intervals.hg38.ht")
 )

# lift over SureSelect V4 intervals
ssv4 = hl.import_bed(f"{project_dir}/data/resources/intervals/ss_all_exon_covered_v4_plus100bp.bed")
ssv4_lifted = liftover_intervals(t=ssv4)
ssv4_lifted.show()
(ssv4_lifted
 .repartition(100)
 .write(f"{project_dir}/data/resources/intervals/ssv4.intervals.hg38.ht")
 )

hl.stop()
