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

from utils.config import NFS_DIR, HDFS_DIR

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("Burden testing pipeline")
logger.setLevel(logging.INFO)

nfs_dir = NFS_DIR
nfs_tmp = f'{NFS_DIR}/tmp'
hdfs_dir = f'{HDFS_DIR}/dir'
project_dir = f'{nfs_dir}/projects/wes_chd_ukbb'


def parse_geneset(geneset_file: str) -> hl.expr.SetExpression:
    """Parse a one-column TSV gene-set file and return its contents as a Hail set."""
    geneset = hl.import_table(paths=geneset_file,
                              no_header=True,
                              delimiter="\t",
                              min_partitions=50,
                              impute=False
                              )
    geneset = geneset.aggregate(hl.agg.collect_as_set(geneset.f0))
    return geneset


def load_mt_with_vep(dataset: str) -> hl.MatrixTable:
    """Load the unfiltered adjusted-genotype MT and annotate VEP fields from the VEP HailTable."""
    mt = hl.read_matrix_table(get_qc_mt_path(dataset=dataset,
                                             part='unphase_adj_genotypes',
                                             split=True))

    vep_ht = get_vep_annotation_ht()

    mt = (mt
          .annotate_rows(LoF=vep_ht[mt.row_key].vep.LoF,
                         Consequence=vep_ht[mt.row_key].vep.Consequence,
                         DOMAINS=vep_ht[mt.row_key].vep.DOMAINS,
                         SYMBOL=vep_ht[mt.row_key].vep.SYMBOL)
          )
    return mt


def filter_to_geneset(mt: hl.MatrixTable,
                      geneset: hl.expr.SetExpression,
                      nfs_tmp: str) -> hl.MatrixTable:
    """Filter MT rows to variants in the gene set and checkpoint to a temp path."""
    mt = (mt
          .filter_rows(hl.set(geneset).contains(mt.SYMBOL))
          .checkpoint(f'{nfs_tmp}/tmp.mt',
                      overwrite=True)
          )
    return mt


def apply_sample_qc(mt: hl.MatrixTable, hdfs_dir: str) -> hl.MatrixTable:
    """Apply per-sample QC filtering and checkpoint the result."""
    logger.info('Applying per sample QC filtering...')

    mt = apply_sample_qc_filtering(mt)

    logger.info('Writing sample qc-filtered MT to disk...')
    mt = (mt
          .checkpoint(f'{hdfs_dir}/chd_ukbb.sample_qc_filtered.mt',
                      overwrite=True)
          )
    return mt


def apply_variant_qc(mt: hl.MatrixTable, hdfs_dir: str) -> hl.MatrixTable:
    """Apply per-variant QC filtering and checkpoint the result."""
    logger.info('Applying per variant QC filtering...')

    mt = apply_variant_qc_filtering(mt)

    # write hard filtered MT to disk
    logger.info('Writing variant qc-filtered mt with rare variants (internal maf 0.01) to disk...')
    mt = (mt
          .checkpoint(f'{hdfs_dir}/chd_ukbb.variant_qc_filtered.mt',
                      overwrite=True)
          )
    return mt


def apply_af_filter(mt: hl.MatrixTable,
                    maf_cutoff: float,
                    hdfs_dir: str) -> hl.MatrixTable:
    """Annotate AF fields from external source and filter rows by allele-frequency cutoff."""
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
          .checkpoint(f'{hdfs_dir}/chd_ukbb.qc_final.rare.mt',
                      overwrite=True)
          )
    return mt


def filter_biallelic(mt: hl.MatrixTable) -> hl.MatrixTable:
    """Filter MT to bi-allelic variants only."""
    logger.info('Running burden test on biallelic variants...')
    mt = mt.filter_rows(bi_allelic_expr(mt))
    return mt


def add_blind_ids_and_sample_meta(mt: hl.MatrixTable) -> hl.MatrixTable:
    """Add blind sample IDs, annotate sample metadata, and filter to cases/controls."""
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
    return mt


def annotate_scores_and_variant_id(mt: hl.MatrixTable) -> hl.MatrixTable:
    """Annotate pathogenic scores from VEP scores HailTable and add variant ID."""
    ht_scores = get_vep_scores_ht()
    mt = mt.annotate_rows(**ht_scores[mt.row_key])

    ## Annotate variants ID
    mt = annotate_variant_id(mt)
    return mt


def aggregate_and_export(mt: hl.MatrixTable, output_file: str) -> None:
    """Aggregate per-genotype counts into row annotations, select columns, and export as TSV."""
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
     .export(output_file)
     )


def main(args):
    ## Init Hail
    hl.init(default_reference=args.default_ref_genome)

    ## Import unfiltered MT with adjusted genotypes and add VEP-annotated fields
    mt = load_mt_with_vep(args.exome_cohort)

    ## Parse geneset and filter to geneset
    geneset = parse_geneset(args.geneset_file)
    mt = filter_to_geneset(mt, geneset, nfs_tmp)

    ## Sample-QC filtering
    if args.apply_sample_qc_filtering:
        mt = apply_sample_qc(mt, hdfs_dir)

    ## Variant-QC filtering
    if args.apply_variant_qc_filtering:
        mt = apply_variant_qc(mt, hdfs_dir)

    ## Filtering by AFs
    # allelic frequency cut-off
    maf_cutoff = args.af_max_threshold

    if args.apply_af_filtering:
        mt = apply_af_filter(mt, maf_cutoff, hdfs_dir)

    ## Filter to bi-allelic variants
    if args.filter_biallelic:
        mt = filter_biallelic(mt)

    ## Generate blind sample IDs and annotate sample metadata
    mt = add_blind_ids_and_sample_meta(mt)

    ## Annotate pathogenic scores and variant ID
    mt = annotate_scores_and_variant_id(mt)

    ## Aggregate genotype counts and export
    aggregate_and_export(mt, args.output_file)


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
