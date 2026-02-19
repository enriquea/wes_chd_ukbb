# eam
# 2021-09-11


"""

Run case-control burden test (logistic regression) on gene-sets, stratified by variant functional category and proband
syndromic status.

usage: geneset_burden_logreg.py [-h] [-i EXOME_COHORT] [-s SET_FILE]
                                [-o OUTPUT_DIR] [--skip_sample_qc_filtering]
                                [--skip_variant_qc_filtering]
                                [--skip_af_filtering]
                                [--af_max_threshold AF_MAX_THRESHOLD]
                                [--filter_biallelic] [--filter_protein_domain]
                                [--hets] [--homs] [--chets] [--homs_chets]
                                [--write_to_file] [--overwrite]
                                [--default_ref_genome DEFAULT_REF_GENOME]
                                [--run_test_mode]

optional arguments:
  -h, --help            show this help message and exit
  -i EXOME_COHORT, --exome_cohort EXOME_COHORT
                        One of <chd_ukbb> or <chd_ddd>
  -s SET_FILE, --set_file SET_FILE
                        Path to gene-set file. Expected two-column TSV file
                        without header (first column denotes cluster name/id,
                        second column denotes gene names)
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Path to output directory
  --skip_sample_qc_filtering
                        Skip the sample QC filtering step
  --skip_variant_qc_filtering
                        Skip the variant QC filtering step
  --skip_af_filtering   Skip the allelic frequency filtering step
  --af_max_threshold AF_MAX_THRESHOLD
                        Allelic frequency cutoff to filter variants (max)
  --filter_biallelic    Run burden test on bi-allelic variants only
  --filter_protein_domain
                        Run burden test on variants within protein domain(s)
                        only
  --hets                Aggregate hets genotypes to run the burden test
  --homs                Aggregate homs genotypes to run the burden test
  --chets               Aggregate compound hets genotypes to run the burden
                        test
  --homs_chets          Aggregate compound hets and/or homs genotypes to run
                        the burden test
  --write_to_file       Write output to BGZ-compressed file
  --overwrite           Overwrite pre-existing data
  --default_ref_genome DEFAULT_REF_GENOME
                        Default reference genome to start Hail
  --run_test_mode       Run pipeline on smaller chunk of data (chr20) for
                        testing propose


"""

import argparse
import functools
import logging
import operator
import sys
import uuid

import hail as hl

from utils.data_utils import (get_af_annotation_ht,
                              get_sample_meta_data,
                              get_qc_mt_path,
                              get_vep_scores_ht,
                              get_mt_data, get_vep_annotation_ht)
from utils.expressions import (af_filter_expr,
                               bi_allelic_expr)
from utils.generic import current_date
from utils.qc import (apply_sample_qc_filtering,
                      apply_variant_qc_filtering)
from utils.stats import logistic_regression
from utils.vep import vep_protein_domain_filter_expr

from utils.config import NFS_DIR, HDFS_DIR

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("Burden testing pipeline")
logger.setLevel(logging.INFO)

nfs_dir = NFS_DIR
nfs_tmp = f'{NFS_DIR}/tmp'
hdfs_dir = f'{HDFS_DIR}/dir'

# Deleterious scores cutoff
MVP_THRESHOLD = 0.8
REVEL_THRESHOLD = 0.5
CADD_THRESHOLD = 25
MPC_THRESHOLD = 2


def generate_clusters_map(ht: hl.Table) -> hl.Table:
    """
    Generate a table mapping gene/features -> cluster_ids.
    Expected as input a two-field HailTable: f0: cluster id, f1: gene/feature symbol.

    :param ht: hl.HailTable
    :return:  hl.HailTable
    """

    # rename fields
    ht = (ht
          .rename({'f0': 'cluster_id', 'f1': 'gene'})
          )
    clusters = (ht
                .group_by('gene')
                .aggregate(cluster_id=hl.agg.collect_as_set(ht.cluster_id))
                .repartition(50)
                .key_by('gene')
                )

    return clusters


def load_mt(run_test_mode, exome_cohort):
    """Load raw (test mode: chr20 sample) or adjusted-genotype MatrixTable."""
    if run_test_mode:
        logger.info('Running pipeline on test data...')
        mt = (get_mt_data(part='raw_chr20')
              .sample_rows(0.1)
              )
    else:
        logger.info('Running pipeline on MatrixTable wih adjusted genotypes...')
        mt = hl.read_matrix_table(get_qc_mt_path(dataset=exome_cohort,
                                                 part='unphase_adj_genotypes',
                                                 split=True))
    return mt


def apply_qc_filters(mt, skip_sample_qc, skip_variant_qc, skip_af_filtering, maf_cutoff):
    """Apply sample QC, variant QC, and AF filters, checkpointing between steps."""
    if not skip_sample_qc:
        logger.info('Applying per sample QC filtering...')
        mt = apply_sample_qc_filtering(mt)

        logger.info('Writing sample qc-filtered mt with rare variants (internal maf 0.01) to disk...')
        mt = (mt
              .write(f'{hdfs_dir}/chd_ukbb.sample_qc_filtered.mt',
                     overwrite=True)
              )

    if not skip_variant_qc:
        logger.info('Applying per variant QC filtering...')

        if hl.hadoop_is_file(f'{hdfs_dir}/chd_ukbb.sample_qc_filtered.mt/_SUCCESS'):
            logger.info('Reading pre-existing sample qc-filtered MT...')
            mt = hl.read_matrix_table(
                f'{hdfs_dir}/chd_ukbb.sample_qc_filtered.mt'
            )
        mt = apply_variant_qc_filtering(mt)

        logger.info('Writing variant qc-filtered mt with rare variants (internal maf 0.01) to disk...')
        mt = (mt
              .write(f'{hdfs_dir}/chd_ukbb.variant_qc_filtered.mt',
                     overwrite=True)
              )

    if not skip_af_filtering:
        if hl.hadoop_is_file(f'{hdfs_dir}/chd_ukbb.variant_qc_filtered.mt/_SUCCESS'):
            logger.info('Reading pre-existing sample/variant qc-filtered MT...')
            mt = hl.read_matrix_table(
                f'{hdfs_dir}/chd_ukbb.variant_qc_filtered.mt'
            )

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

        logger.info(f'Writing sample/variant QCed MT with rare variants at maf: {maf_cutoff}.')
        mt = (mt
              .write(f'{hdfs_dir}/chd_ukbb.qc_final.rare.mt',
                     overwrite=True)
              )

    return mt


def annotate_vep_fields(mt):
    """Load the VEP annotation table and annotate LoF, Consequence, DOMAINS, SYMBOL onto rows."""
    if hl.hadoop_is_file(f'{hdfs_dir}/chd_ukbb.qc_final.rare.mt/_SUCCESS'):
        logger.info('Reading pre-existing sample/variant qc-filtered MT with rare variants...')
        mt = hl.read_matrix_table(
            f'{hdfs_dir}/chd_ukbb.qc_final.rare.mt'
        )

    vep_ht = get_vep_annotation_ht()
    mt = (mt
          .annotate_rows(LoF=vep_ht[mt.row_key].vep.LoF,
                         Consequence=vep_ht[mt.row_key].vep.Consequence,
                         DOMAINS=vep_ht[mt.row_key].vep.DOMAINS,
                         SYMBOL=vep_ht[mt.row_key].vep.SYMBOL)
          )
    return mt


def apply_variant_filters(mt, filter_biallelic, filter_protein_domain):
    """Optionally filter to bi-allelic variants or protein-domain variants."""
    if filter_biallelic:
        logger.info('Running burden test on biallelic variants...')
        mt = mt.filter_rows(bi_allelic_expr(mt))

    if filter_protein_domain:
        logger.info('Running burden test on variants within protein domain(s)...')
        mt = mt.filter_rows(
            vep_protein_domain_filter_expr(mt.DOMAINS),
            keep=True
        )

    return mt


def annotate_pathogenic_scores(mt):
    """Annotate rows with pathogenic/deleteriousness scores from the VEP scores table."""
    ht_scores = get_vep_scores_ht()
    mt = mt.annotate_rows(**ht_scores[mt.row_key])
    return mt


def build_consequence_groups(mt):
    """Build fixed consequence-group boolean expressions and explode rows by group."""
    score_expr_ann = {'hcLOF': mt.LoF == 'HC',
                      'syn': mt.Consequence == 'synonymous_variant',
                      'miss': mt.Consequence == 'missense_variant'
                      }

    score_expr_ann.update(
        {'missC': (hl.sum([(mt['vep.MVP_score'] >= MVP_THRESHOLD),
                           (mt['vep.REVEL_score'] >= REVEL_THRESHOLD),
                           (mt['vep.CADD_PHRED'] >= CADD_THRESHOLD)]) >= 2) & score_expr_ann.get('miss')}
    )

    score_expr_ann.update(
        {'hcLOF_missC': score_expr_ann.get('hcLOF') | score_expr_ann.get('missC')}
    )

    mt = (mt
          .annotate_rows(csq_group=score_expr_ann)
          )

    mt = (mt
          .transmute_rows(csq_group=hl.set(hl.filter(lambda x:
                                                     mt.csq_group.get(x),
                                                     mt.csq_group.keys())))
          )

    mt = (mt
          .filter_rows(hl.len(mt.csq_group) > 0)
          )

    mt = (mt
          .explode_rows(mt.csq_group)
          )

    return mt


def aggregate_gene_level(mt):
    """First-stage aggregation: group by gene/csq_group and aggregate hets, homs, chets."""
    mt_grouped = (mt
                  .group_rows_by(mt['SYMBOL'], mt['csq_group'])
                  .aggregate(hets=hl.agg.any(mt.GT.is_het()),
                             homs=hl.agg.any(mt.GT.is_hom_var()),
                             chets=hl.agg.count_where(mt.GT.is_het()) >= 2)
                  .repartition(100)
                  .persist()
                  )
    return mt_grouped


def annotate_gene_sets(mt_grouped, set_file):
    """Import gene-set clusters, annotate rows, explode by cluster_id, and filter defined entries."""
    clusters = hl.import_table(set_file,
                               no_header=True,
                               delimiter="\t",
                               min_partitions=50,
                               impute=False
                               )
    clusters = generate_clusters_map(clusters)

    mt_grouped = (mt_grouped
                  .annotate_rows(**clusters[mt_grouped.SYMBOL])
                  )

    mt_grouped = (mt_grouped
                  .explode_rows(mt_grouped.cluster_id)
                  )

    mt_grouped = (mt_grouped
                  .filter_rows(hl.is_defined(mt_grouped.csq_group) &
                               hl.is_defined(mt_grouped.cluster_id))
                  )

    return mt_grouped


def aggregate_by_genotype_mode(mt_grouped, hets, homs, chets, homs_chets):
    """Second-stage aggregation: group by gene-set/csq_group and union selected genotype modes."""
    mts = []

    if homs:
        mt_homs = (mt_grouped
                   .group_rows_by(mt_grouped.csq_group, mt_grouped.cluster_id)
                   .aggregate(mac=hl.int(hl.agg.sum(mt_grouped.homs)))
                   .repartition(100)
                   .persist()
                   .annotate_rows(agg_genotype='homs')
                   )
        mts.append(mt_homs)

    if chets:
        mt_chets = (mt_grouped
                    .group_rows_by(mt_grouped.csq_group, mt_grouped.cluster_id)
                    .aggregate(mac=hl.int(hl.agg.sum(mt_grouped.chets)))
                    .repartition(100)
                    .persist()
                    .annotate_rows(agg_genotype='chets')
                    )
        mts.append(mt_chets)

    if homs_chets:
        mt_homs_chets = (mt_grouped
                         .group_rows_by(mt_grouped.csq_group, mt_grouped.cluster_id)
                         .aggregate(mac=hl.int(hl.agg.count_where(mt_grouped.chets | mt_grouped.homs)))
                         .repartition(100)
                         .persist()
                         .annotate_rows(agg_genotype='homs_chets')
                         )
        mts.append(mt_homs_chets)

    if hets:
        mt_hets = (mt_grouped
                   .group_rows_by(mt_grouped.csq_group, mt_grouped.cluster_id)
                   .aggregate(mac=hl.int(hl.agg.sum(mt_grouped.hets)))
                   .repartition(100)
                   .persist()
                   .annotate_rows(agg_genotype='hets')
                   )
        mts.append(mt_hets)

    mt_joint = hl.MatrixTable.union_rows(*mts)
    return mt_joint


def annotate_sample_covariates(mt_joint):
    """Annotate sample columns with covariates and case/control metadata, filter to cases/controls."""
    covariates = hl.read_table(
        f'{nfs_dir}/hail_data/sample_qc/chd_ukbb.sample_covariates.ht'
    )
    mt_joint = (mt_joint
                .annotate_cols(**covariates[mt_joint.s])
                )

    tb_sample = get_sample_meta_data()
    mt_joint = (mt_joint
                .annotate_cols(**tb_sample[mt_joint.s])
                )

    mt_joint = (mt_joint
                .filter_cols(mt_joint['phe.is_case'] | mt_joint['phe.is_control'])
                )

    return mt_joint


def run_logreg_stratified(mt_joint, maf_cutoff):
    """Run logistic regression stratified by proband type (all_cases, syndromic, nonsyndromic)."""
    analysis = ['all_cases', 'syndromic', 'nonsyndromic']
    tbs = []
    covs = ['sex', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5']

    for proband in analysis:
        logger.info(f'Running burden test for {proband}...')

        mt_tmp = hl.MatrixTable

        if proband == 'all_cases':
            mt_tmp = mt_joint
        if proband == 'syndromic':
            mt_tmp = mt_joint.filter_cols(~mt_joint['phe.is_nonsyndromic'])
        if proband == 'nonsyndromic':
            mt_tmp = mt_joint.filter_cols(~mt_joint['phe.is_syndromic'])

        tb_logreg = logistic_regression(mt=mt_tmp,
                                        x_expr='mac',
                                        response='phe.is_case',
                                        covs=covs,
                                        pass_through=['agg_genotype'],
                                        extra_fields={'analysis': proband,
                                                      'maf': maf_cutoff,
                                                      'covs': '|'.join(covs)})

        tbs.append(tb_logreg)

    tb_final = hl.Table.union(*tbs)
    return tb_final


def export_results(tb, output_dir, exome_cohort, write_to_file):
    """Checkpoint results table and optionally export as TSV."""
    date = current_date()
    run_hash = str(uuid.uuid4())[:6]
    output_path = f'{output_dir}/{date}/{exome_cohort}.logreg_burden.{run_hash}.ht'

    tb = (tb
          .checkpoint(output=output_path)
          )

    if write_to_file:
        (tb
         .export(f'{output_path}.tsv')
         )

    return tb


def main(args):
    hl.init(default_reference=args.default_ref_genome)

    maf_cutoff = args.af_max_threshold

    mt = load_mt(args.run_test_mode, args.exome_cohort)

    mt = apply_qc_filters(mt,
                          args.skip_sample_qc_filtering,
                          args.skip_variant_qc_filtering,
                          args.skip_af_filtering,
                          maf_cutoff)

    mt = annotate_vep_fields(mt)

    mt = apply_variant_filters(mt, args.filter_biallelic, args.filter_protein_domain)

    mt = annotate_pathogenic_scores(mt)

    mt = build_consequence_groups(mt)

    mt_grouped = aggregate_gene_level(mt)

    mt_grouped = annotate_gene_sets(mt_grouped, args.set_file)

    mt_joint = aggregate_by_genotype_mode(mt_grouped, args.hets, args.homs, args.chets, args.homs_chets)

    mt_joint = annotate_sample_covariates(mt_joint)

    tb_final = run_logreg_stratified(mt_joint, maf_cutoff)

    export_results(tb_final, args.output_dir, args.exome_cohort, args.write_to_file)

    hl.stop()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--exome_cohort', help="One of <chd_ukbb> or <chd_ddd>",
                        type=str, default=None)

    parser.add_argument('-s',
                        '--set_file',
                        help="Path to gene-set file. Expected two-column TSV file without header "
                             "(first column denotes cluster name/id, second column denotes gene names)",
                        default=None)

    parser.add_argument('-o', '--output_dir', help='Path to output directory',
                        type=str, default=f'{nfs_dir}/hail_data/burden')

    parser.add_argument('--skip_sample_qc_filtering', help='Skip the sample QC filtering step',
                        action='store_true')

    parser.add_argument('--skip_variant_qc_filtering', help='Skip the variant QC filtering step',
                        action='store_true')

    parser.add_argument('--skip_af_filtering', help='Skip the allelic frequency filtering step',
                        action='store_true')

    parser.add_argument('--af_max_threshold', help='Allelic frequency cutoff to filter variants (max)',
                        type=float, default=0.001)

    parser.add_argument('--filter_biallelic', help='Run burden test on bi-allelic variants only',
                        action='store_true')

    parser.add_argument('--filter_protein_domain', help='Run burden test on variants within protein domain(s) only',
                        action='store_true')

    parser.add_argument('--hets', help='Aggregate hets genotypes to run the burden test',
                        action='store_true')

    parser.add_argument('--homs', help='Aggregate homs genotypes to run the burden test',
                        action='store_true')

    parser.add_argument('--chets', help='Aggregate compound hets genotypes to run the burden test',
                        action='store_true')

    parser.add_argument('--homs_chets', help='Aggregate compound hets and/or homs genotypes to run the burden test',
                        action='store_true')

    parser.add_argument('--write_to_file', help='Write output to BGZ-compressed file',
                        action='store_true')

    parser.add_argument('--overwrite', help='Overwrite pre-existing data',
                        action='store_true')

    parser.add_argument('--default_ref_genome', help='Default reference genome to start Hail',
                        type=str, default='GRCh38')

    parser.add_argument('--run_test_mode', help='Run pipeline on smaller chunk of data (chr20) for testing propose',
                        action='store_true')

    args = parser.parse_args()

    if args.hets + args.homs + args.chets + args.homs_chets == 0:
        sys.exit('Specifies at least one of hets, homs, chets or homs_chets options...')

    main(args)
