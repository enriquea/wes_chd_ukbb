# eam
# 16.05.22

"""

Annotate cohort-specific allelic frequencies (internal)
stratified by major continental ancestries.

"""

import argparse
import logging

import hail as hl

from utils.data_utils import (get_mt_data,
                              get_sample_qc_ht_path,
                              get_qc_mt_path,
                              get_variant_af_pops_ht_path)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("Annotate AF")
logger.setLevel(logging.INFO)

nfs_dir = 'file:///home/ubuntu/data'
nfs_tmp = 'file:///home/ubuntu/data/tmp'
hdfs_dir = 'hdfs://spark-master:9820/dir'


def annotate_release_samples(ht: hl.Table):
    ht=ht.annotate(release_sample=hl.cond(
        (hl.len(ht.hard_filters) == 0) &
        (hl.len(ht.pop_platform_filters) == 0) &
        ~ht.is_related,
        True, False))
    return ht


def main(args):
    # Start Hail
    hl.init(default_reference=args.default_ref_genome)

    ds = args.exome_cohort

    if args.run_test_mode:
        logger.info('Running pipeline on test data...')
        mt = (get_mt_data(part='raw_chr20')
              .sample_rows(0.1)
              )
    else:
        logger.info('Running pipeline on MatrixTable wih adjusted genotypes...')
        mt = hl.read_matrix_table(get_qc_mt_path(dataset=ds,
                                                 part='unphase_adj_genotypes',
                                                 split=True))

    # 1. Annotate sample qc filters
    sample_qc_path_ht = get_sample_qc_ht_path(dataset=args.exome_cohort,
                                              part='final_qc')
    sample_qc_ht = hl.read_table(sample_qc_path_ht)
    sample_qc_ht = (annotate_release_samples(sample_qc_ht)
                    .select('predicted_pop', 'release_sample'))

    mt = (mt
          .annotate_cols(**sample_qc_ht[mt.s])
          )

    # 2. keep release samples
    mt = (mt
          .filter_cols(mt.release_sample, keep=True)
          )

    # 3. compute AFs stratified by population
    pops = mt.aggregate_cols(hl.agg.collect_as_set(mt.predicted_pop))

    # compute call stats by pop groups
    ht_cs = (mt
             .annotate_rows(**{f'call_stats_{p}': hl.agg.filter(mt.predicted_pop == p,
                                                                hl.agg.call_stats(mt.GT, mt.alleles))
                               for p in pops}
                            )
             .rows()
             )

    # annotate (unnest) AFs, ACs and Homozygote counts per population groups
    field_call_stats = ['AF', 'AC', 'homozygote_count']
    ht_cs = (ht_cs
             .annotate(**{f'{f}_{p}': ht_cs[f'call_stats_{p}'].get(f)[1]
                          for p in pops
                          for f in field_call_stats})
             )

    # export table
    ht_cs = (ht_cs
             .checkpoint(get_variant_af_pops_ht_path(dataset=ds),
                         overwrite=args.overwrite)
             )

    # Stop Hail
    hl.stop()

    print("Finished...")


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
    parser.add_argument('--run_test_mode', help='Run pipeline on smaller chunk of data (chr20) for testing propose',
                        action='store_true')

    args = parser.parse_args()

    main(args)
