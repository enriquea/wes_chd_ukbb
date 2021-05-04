# eam
# 2021-04-29

"""
Finalise sample QC.
 Pipeline actions:
         1. Flag samples with hard filters
         2. Flag samples with population filters
         3. Flag related samples
         4. Flag pop/platform-specific outliers
         5. Release sample QCed MT with adjusted genotypes


usage: finalise_sample_qc.py [-h] [--exome_cohort EXOME_COHORT]
                             [--write_to_file] [--overwrite]
                             [--default_ref_genome DEFAULT_REF_GENOME]

optional arguments:
  -h, --help            show this help message and exit
  --exome_cohort EXOME_COHORT
                        One of <chd_ukbb> or <chd_ddd>
  --write_to_file       Write output to BGZ-compressed file
  --overwrite           Overwrite pre-existing data
  --default_ref_genome DEFAULT_REF_GENOME
                        Default reference genome to start Hail


"""
import argparse
import logging

import hail as hl
from gnomad.utils.annotations import annotate_adj

from utils.data_utils import (get_mt_data,
                              get_sample_qc_ht_path,
                              get_qc_mt_path)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("Finalise Sample QC")
logger.setLevel(logging.INFO)

# hdfs_dir = 'hdfs://spark-master:9820'
nfs_dir = 'file:///home/ubuntu/data'


def get_related_samples_to_drop() -> hl.Table:
    return hl.read_table(
        f'{nfs_dir}/hail_data/sample_qc/chd_ukbb.related_samples_to_remove.ht'
    )


def main(args):
    # Start Hail
    hl.init(default_reference=args.default_ref_genome)

    # Import unfiltered split MT
    mt = (get_mt_data(dataset=args.exome_cohort, part='unfiltered')
          .select_cols()
          )

    # Annotate samples filters
    sample_qc_filters = {}

    # 1. Add sample hard filters annotation expr
    sample_qc_hard_filters_ht = hl.read_table(
        get_sample_qc_ht_path(dataset=args.exome_cohort,
                              part='hard_filters')
    )

    sample_qc_filters.update(
        {'hard_filters': sample_qc_hard_filters_ht[mt.s]['hard_filters']}
    )

    # 2. Add population qc filters annotation expr
    sample_qc_pop_ht = hl.read_table(
        get_sample_qc_ht_path(dataset=args.exome_cohort,
                              part='population_qc')
    )

    sample_qc_filters.update(
        {'predicted_pop': sample_qc_pop_ht[mt.s]['predicted_pop']}
    )

    # 3. Add relatedness filters annotation expr
    related_samples_to_drop = get_related_samples_to_drop()
    related_samples = hl.set(related_samples_to_drop
                             .aggregate(hl.agg.collect_as_set(related_samples_to_drop.node.id))
                             )

    sample_qc_filters.update(
        {'is_related': related_samples.contains(mt.s)}
    )

    # 4. Add stratified sample qc (population/platform) annotation expr
    sample_qc_pop_platform_filters_ht = hl.read_table(
        get_sample_qc_ht_path(dataset=args.exome_cohort,
                              part='stratified_metrics_filter')
    )

    sample_qc_filters.update(
        {'pop_platform_filters': sample_qc_pop_platform_filters_ht[mt.s]['pop_platform_filters']}
    )

    mt = (mt
          .annotate_cols(**sample_qc_filters)
          )

    # Final sample qc filter joint expression
    final_sample_qc_ann_expr = {'pass_filters': hl.cond(
        (hl.len(mt.hard_filters) == 0) &
        (hl.len(mt.pop_platform_filters) == 0) &
        (mt.predicted_pop == 'EUR') &
        ~mt.is_related,
        True, False)}
    mt = (mt
          .annotate_cols(**final_sample_qc_ann_expr)
          )

    # Export final sample QC annotations to file
    if args.write_to_file:
        (mt.cols().export(
            f'{nfs_dir}/hail_data/sample_qc/chd_ukbb.sample_qc.final_joint_filters.tsv.bgz')
         )

    # Release final sample QCed MT with adjusted genotypes filtered
    mt = annotate_adj(mt)
    mt = mt.filter_entries(
        mt.adj
    ).select_entries('GT', 'DP', 'GQ', 'adj')

    # write MT
    mt.write(
        get_qc_mt_path(dataset=args.exome_cohort,
                       part='sample_qc_adj_genotypes',
                       split=True),
        overwrite=args.overwrite
    )

    # Stop Hail
    hl.stop()

    print("Finished!")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exome_cohort', help="One of <chd_ukbb> or <chd_ddd>",
                        type=str, default=None)
    parser.add_argument('--write_to_file', help='Write output to BGZ-compressed file',
                        action='store_true')
    parser.add_argument('--overwrite', help='Overwrite pre-existing data',
                        action='store_true')
    parser.add_argument('--default_ref_genome', help='Default reference genome to start Hail',
                        type=str, default='GRCh38')

    args = parser.parse_args()

    main(args)
