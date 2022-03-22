# eam
# 22.03.22


"""

Retrieving variant/genotype/sample information
after extensive filtering from a Hail MatrixTable.

usage: mt_query.py [-h] [-i EXOME_COHORT] [-o OUTPUT_FILE] [-s GENESET_FILE]
                   [--apply_sample_qc_filtering]
                   [--apply_variant_qc_filtering] [--apply_af_filtering]
                   [--af_max_threshold AF_MAX_THRESHOLD] [--filter_biallelic]
                   [--default_ref_genome DEFAULT_REF_GENOME]

optional arguments:
  -h, --help            show this help message and exit
  -i EXOME_COHORT, --exome_cohort EXOME_COHORT
                        One of <chd_ukbb> or <chd_ddd>
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        Path to output file
  -s GENESET_FILE, --geneset_file GENESET_FILE
                        Path to gene-set file. Expected one-column TSV file
                        without header
  --apply_sample_qc_filtering
                        Apply sample QC-filtering
  --apply_variant_qc_filtering
                        Apply variant QC-filtering
  --apply_af_filtering  Apply AF-based variant filtering
  --af_max_threshold AF_MAX_THRESHOLD
                        Allelic frequency cutoff to filter variants (max)
  --filter_biallelic    Run burden test on bi-allelic variants only
  --default_ref_genome DEFAULT_REF_GENOME
                        Default reference genome to start Hail


"""

import argparse
import functools
import logging
import operator

import hail as hl

from utils.qc import (apply_sample_qc_filtering,
                      apply_variant_qc_filtering)
from utils.annotation import annotate_variant_id
from utils.data_utils import (get_af_annotation_ht,
                              get_sample_meta_data,
                              get_qc_mt_path,
                              get_vep_scores_ht,
                              get_vep_annotation_ht)
from utils.expressions import (af_filter_expr,
                               bi_allelic_expr)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("Burden testing pipeline")
logger.setLevel(logging.INFO)

nfs_dir = 'file:///home/ubuntu/data'
nfs_tmp = 'file:///home/ubuntu/data/tmp'
hdfs_dir = 'hdfs://spark-master:9820/dir'
project_dir = f'{nfs_dir}/projects/wes_chd_ukbb'


def parse_geneset(geneset_file: str) -> hl.expr.SetExpression:
    # parse geneset file
    geneset = hl.import_table(paths=geneset_file,
                              no_header=True,
                              delimiter="\t",
                              min_partitions=50,
                              impute=False
                              )
    geneset = geneset.aggregate(hl.agg.collect_as_set(geneset.f0))
    return geneset


def main(args):
    ## Init Hail
    hl.init(default_reference=args.default_ref_genome)

    ## Import unfiltered MT with adjusted genotypes
    ds = args.exome_cohort
    mt = hl.read_matrix_table(get_qc_mt_path(dataset=ds,
                                             part='unphase_adj_genotypes',
                                             split=True))

    ## Sample-QC filtering
    if args.apply_sample_qc_filtering:
        logger.info('Applying per sample QC filtering...')

        mt = apply_sample_qc_filtering(mt)

        logger.info('Writing sample qc-filtered mt with rare variants (internal maf 0.01) to disk...')
        mt = (mt
              .write(f'{hdfs_dir}/chd_ukbb.sample_qc_filtered.mt',
                     overwrite=True)
              )

    ## Variant-QC filtering
    if args.apply_variant_qc_filtering:
        logger.info('Applying per variant QC filtering...')

        mt = apply_variant_qc_filtering(mt)

        # write hard filtered MT to disk
        logger.info('Writing variant qc-filtered mt with rare variants (internal maf 0.01) to disk...')
        mt = (mt
              .write(f'{hdfs_dir}/chd_ukbb.variant_qc_filtered.mt',
                     overwrite=True)
              )

    ## Filtering by AFs

    # allelic frequency cut-off
    maf_cutoff = args.af_max_threshold

    if args.apply_af_filtering:
        # Annotate allelic frequencies from external source,
        # and compute internal AF on samples passing QC
        af_ht = get_af_annotation_ht()

        mt = (mt
              .annotate_rows(**af_ht[mt.row_key])
              )

        filter_expressions = [af_filter_expr(mt, 'internal_af', af_cutoff=maf_cutoff),
                              af_filter_expr(mt, 'gnomad_genomes_af', af_cutoff=maf_cutoff),
                              af_filter_expr(mt, 'gnomAD_AF', af_cutoff=maf_cutoff),
                              af_filter_expr(mt, 'ger_af', af_cutoff=maf_cutoff),
                              af_filter_expr(mt, 'rumc_af', af_cutoff=maf_cutoff),
                              af_filter_expr(mt, 'bonn_af', af_cutoff=maf_cutoff)
                              ]

        mt = (mt
              .filter_rows(functools.reduce(operator.iand, filter_expressions), keep=True)
              )

        logger.info('Writing AF-filtered MT to disk...')
        mt = (mt
              .write(f'{hdfs_dir}/chd_ukbb.qc_final.rare.mt',
                     overwrite=True)
              )

    ## Filter to bi-allelic variants
    if args.filter_biallelic:
        logger.info('Running burden test on biallelic variants...')
        mt = mt.filter_rows(bi_allelic_expr(mt))

    ## Add VEP-annotated fields
    vep_ht = get_vep_annotation_ht()

    mt = (mt
          .annotate_rows(LoF=vep_ht[mt.row_key].vep.LoF,
                         Consequence=vep_ht[mt.row_key].vep.Consequence,
                         DOMAINS=vep_ht[mt.row_key].vep.DOMAINS,
                         SYMBOL=vep_ht[mt.row_key].vep.SYMBOL)
          )

    ## Parse geneset
    geneset = parse_geneset(args.geneset_file)

    ## Filter to geneset
    mt = (mt
          .filter_rows(hl.set(geneset).contains(mt.SYMBOL))
          )

    ## Generate blind sample IDs
    mt = mt.add_col_index()

    mt = (mt
          .annotate_cols(BIID=hl.str('BLIND_ID_') + hl.str(mt.col_idx))
          )

    ## Add cases/controls sample annotations
    tb_sample = get_sample_meta_data()
    mt = (mt
          .annotate_cols(**tb_sample[mt.s])
          )

    mt = (mt
          .filter_cols(mt['phe.is_case'] | mt['phe.is_control'])
          )

    ## Annotate pathogenic scores
    ht_scores = get_vep_scores_ht()
    mt = mt.annotate_rows(**ht_scores[mt.row_key])

    ## Annotate variants ID
    mt = annotate_variant_id(mt)

    # annotate samples
    ann_expr = {'n_het_cases': hl.agg.filter(mt.GT.is_het() & mt['phe.is_case'], hl.agg.count()),
                'n_hom_cases': hl.agg.filter(mt.GT.is_hom_var() & mt['phe.is_case'], hl.agg.count()),
                'n_het_syndromic': hl.agg.filter(mt.GT.is_het() & mt['phe.is_syndromic'], hl.agg.count()),
                'n_hom_syndromic': hl.agg.filter(mt.GT.is_hom_var() & mt['phe.is_syndromic'], hl.agg.count()),
                'n_het_nonsyndromic': hl.agg.filter(mt.GT.is_het() & mt['phe.is_nonsyndromic'], hl.agg.count()),
                'n_hom_nonsyndromic': hl.agg.filter(mt.GT.is_hom_var() & mt['phe.is_nonsyndromic'], hl.agg.count()),
                'n_het_controls': hl.agg.filter(mt.GT.is_het() & ~mt['phe.is_case'], hl.agg.count()),
                'n_hom_controls': hl.agg.filter(mt.GT.is_hom_var() & ~mt['phe.is_case'], hl.agg.count()),
                'het_case_ids': hl.agg.filter(mt.GT.is_het() & mt['phe.is_case'],
                                              hl.delimit(hl.agg.collect_as_set(mt.BIID), '|')),
                'hom_case_ids': hl.agg.filter(mt.GT.is_hom_var() & mt['phe.is_case'],
                                              hl.delimit(hl.agg.collect_as_set(mt.BIID), '|')),
                'het_control_ids': hl.agg.filter(mt.GT.is_het() & ~mt['phe.is_case'],
                                                 hl.delimit(hl.agg.collect_as_set(mt.BIID), '|')),
                'hom_control_ids': hl.agg.filter(mt.GT.is_hom_var() & ~mt['phe.is_case'],
                                                 hl.delimit(hl.agg.collect_as_set(mt.BIID), '|'))
                }
    ht = (mt
          .annotate_rows(**ann_expr)
          .rows()
          .key_by()
          .select(*list(['vid', 'Consequence', 'SYMBOL', 'internal_af', 'gnomAD_AF', 'vep.MVP_score', 'vep.REVEL_score',
                         'vep.MPC_score', 'vep.CADD_PHRED']) + list(ann_expr.keys()))
          )

    # export results
    (ht
     .export(args.output_file)
     )


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--exome_cohort', help="One of <chd_ukbb> or <chd_ddd>",
                        type=str, default=None)

    parser.add_argument('-o', '--output_file', help='Path to output file',
                        type=str, default=None)

    parser.add_argument('-s',
                        '--geneset_file',
                        help="Path to gene-set file. Expected one-column TSV file without header",
                        default=None)

    parser.add_argument('--apply_sample_qc_filtering', help='Apply sample QC-filtering',
                        action='store_true')

    parser.add_argument('--apply_variant_qc_filtering', help='Apply variant QC-filtering',
                        action='store_true')

    parser.add_argument('--apply_af_filtering', help='Apply AF-based variant filtering',
                        action='store_true')

    parser.add_argument('--af_max_threshold', help='Allelic frequency cutoff to filter variants (max)',
                        type=float, default=0.001)

    parser.add_argument('--filter_biallelic', help='Run burden test on bi-allelic variants only',
                        action='store_true')

    parser.add_argument('--default_ref_genome', help='Default reference genome to start Hail',
                        type=str, default='GRCh38')

    args = parser.parse_args()

    main(args)
