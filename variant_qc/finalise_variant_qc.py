# eam
# 2021-05-13

"""
Finalize variant QC

Actions:
 - Apply hard filters
 - Apply VQSR filter
 - Apply RF filter


usage: finalise_variant_qc.py [-h] [--exome_cohort EXOME_COHORT]
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

from utils.data_utils import (get_qc_mt_path,
                              get_variant_qc_ht_path,
                              get_vep_annotation_ht)

# from utils.constants import *

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# hdfs_dir = 'hdfs://spark-master:9820'
nfs_dir = 'file:///home/ubuntu/data'

INBREEDING_COEFFICIENT_CUTOFF = -0.3

RF_PROBABILITY_SNV_CUTOFF = 0.2  # this cutoff could be lower for SNVs (0.1?)

RF_PROBABILITY_INDEL_CUTOFF = 0.2


def main(args):
    # Start Hail
    hl.init(default_reference=args.default_ref_genome)

    # Import adj genotype MT and remove
    mt = hl.read_matrix_table(get_qc_mt_path(dataset=args.exome_cohort,
                                             part='sample_qc_adj_genotypes',
                                             split=True))

    # keep samples passing QC filtering
    mt = (mt
          .filter_cols(mt.pass_filters)
          .select_cols()
          .select_rows()
          )

    # import variant info fields (vcf info)
    variant_info_ht = (get_vep_annotation_ht()
                       .drop('vep')
                       )

    # Add useful annotation for variant hard filter
    ht = (mt
          .annotate_rows(inbreeding_coeff=variant_info_ht[mt.row_key].info.InbreedingCoeff,
                         vqsr_filter=variant_info_ht[mt.row_key].filters,
                         VQSLOD=variant_info_ht[mt.row_key].info.VQSLOD,
                         gt_counts=hl.agg.count_where(hl.is_defined(mt.GT))  # expected MT filtered to high-quality GT
                         )
          .rows()
          )

    # 1. Apply variant hard filters
    # hard filter expression
    variant_hard_filter_expr = {'fail_inbreeding_coeff': ht.inbreeding_coeff < INBREEDING_COEFFICIENT_CUTOFF,
                                'AC0': ht.gt_counts == 0}

    ht = (ht
          .annotate(**variant_hard_filter_expr)
          )

    # 2. Apply VQSR filter
    ht = (ht
          .annotate(fail_vqsr=hl.len(ht.vqsr_filter) != 0)
          )

    # 3. Apply RF filter

    # import/parse rf final HT
    ht_rf = hl.read_table(
        get_variant_qc_ht_path(part='rf_result')
    )

    ht_rf = (ht_rf
             .select(rf_probability_tp=ht_rf.rf_probability['TP'],
                     variant_type=ht_rf.variant_type)
             )

    ht = (ht
          .annotate(**ht_rf[ht.key])
          )

    ht = (ht
          .annotate(fail_rf=hl.case()
                    .when((ht.rf_probability_tp < RF_PROBABILITY_SNV_CUTOFF) & (ht.variant_type == 'snv'), True)
                    .when((ht.rf_probability_tp < RF_PROBABILITY_INDEL_CUTOFF) & (ht.variant_type == 'indel'), True)
                    .default(False)
                    )
          )

    # 4. Summary final variant QC

    # final variant qc filter joint expression
    final_variant_qc_ann_expr = {
        'pass_variant_qc_filters': hl.cond(
            ~ht.fail_inbreeding_coeff &
            ~ht.AC0 &
            ~ht.fail_vqsr &
            ~ht.fail_rf,
            True, False)}
    ht = (ht
          .annotate(**final_variant_qc_ann_expr)
          )

    # Counts the number of variants (snv and indels) affected by every filter and add as global field
    filter_flags = ['fail_inbreeding_coeff',
                    'AC0',
                    'fail_vqsr',
                    'fail_rf',
                    'pass_variant_qc_filters']

    summary_filter_expr = {v: hl.struct(**{f: hl.agg.filter(ht.variant_type == v, hl.agg.counter(ht[f]))
                                           for f in filter_flags})
                           for v in ['snv', 'indel']
                           }

    ht = ht.annotate_globals(summary_filter=ht.aggregate(summary_filter_expr, _localize=False))

    # write HT variant QC final table
    output_path = get_variant_qc_ht_path(dataset=args.exome_cohort,
                                         part='final_qc')
    ht = ht.checkpoint(
        output_path,
        overwrite=args.overwrite
    )

    # print filter summary
    logger.info(f'Variant QC filter summary: {ht.summary_filter.collect()}')

    # export HT to file
    if args.write_to_file:
        ht.export(
            f'{output_path}.tsv.bgz'
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
