# eam
# 2021-09-08

"""

Run case-control burden test (Fisher Exact) stratified by variant functional category and proband
syndromic status.

usage: gene_burden_fet.py [-h] [-i EXOME_COHORT] [-o OUTPUT_DIR]
                          [--write_to_file] [--overwrite]
                          [--default_ref_genome DEFAULT_REF_GENOME]
                          [--run_test_mode] [--skip_sample_qc_filtering]
                          [--skip_variant_qc_filtering] [--skip_af_filtering]
                          [--af_max_threshold AF_MAX_THRESHOLD]
                          [--filter_biallelic] [--filter_protein_domain]
                          [--hets] [--homs] [--chets] [--homs_chets]
                          [--test_cadd_miss] [--test_revel_miss]
                          [--test_mvp_miss] [--test_pav_hc]
                          [--cadd_threshold CADD_THRESHOLD]
                          [--revel_threshold REVEL_THRESHOLD]
                          [--mvp_threshold MVP_THRESHOLD]

optional arguments:
  -h, --help            show this help message and exit
  -i EXOME_COHORT, --exome_cohort EXOME_COHORT
                        One of <chd_ukbb> or <chd_ddd>
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Path to output directory
  --write_to_file       Write output to BGZ-compressed file
  --overwrite           Overwrite pre-existing data
  --default_ref_genome DEFAULT_REF_GENOME
                        Default reference genome to start Hail
  --run_test_mode       Run pipeline on smaller chunk of data (chr20) for
                        testing propose

Options for QC and filtering:
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

Options to aggregate genotypes for burden testing:
  --hets                Aggregate hets genotypes to run the burden test
  --homs                Aggregate homs genotypes to run the burden test
  --chets               Aggregate compound hets genotypes to run the burden
                        test
  --homs_chets          Aggregate compound hets and/or homs genotypes to run
                        the burden test

Type of variant to be included in burden testing:
  --test_cadd_miss      Test missense variants passing a predefined CADD score
                        threshold
  --test_revel_miss     Test missense variants passing a predefined REVEL
                        score threshold
  --test_mvp_miss       Test missense variants passing a predefined MVP score
                        threshold
  --test_pav_hc         Test high confidence protein-altering variants

Pathogenic score threshold options:
  --cadd_threshold CADD_THRESHOLD
                        CADD score threshold
  --revel_threshold REVEL_THRESHOLD
                        REVEL score threshold
  --mvp_threshold MVP_THRESHOLD
                        MVP score threshold


"""

import argparse
import functools
import logging
import operator
import uuid
import sys

import hail as hl

from utils.data_utils import (get_af_annotation_ht,
                              get_sample_meta_data,
                              get_qc_mt_path,
                              get_vep_scores_ht,
                              get_mt_data, get_vep_annotation_ht)
from utils.expressions import (af_filter_expr,
                               bi_allelic_expr)
from utils.filter import filter_ccr
from utils.generic import current_date
from utils.qc import (apply_sample_qc_filtering,
                      apply_variant_qc_filtering)
from utils.stats import compute_fisher_exact
from utils.vep import vep_protein_domain_filter_expr

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("Burden testing pipeline")
logger.setLevel(logging.INFO)

nfs_dir = 'file:///home/ubuntu/data'
nfs_tmp = 'file:///home/ubuntu/data/tmp'
hdfs_dir = 'hdfs://spark-master:9820/dir'


def main(args):
    hl.init(default_reference=args.default_ref_genome)

    # Deleterious scores cutoff
    mvp_threshold = args.mvp_threshold
    revel_threshold = args.revel_threshold
    cadd_threshold = args.cadd_threshold
    # MPC_THRESHOLD = 2

    if args.run_test_mode:
        logger.info('Running pipeline on test data...')
        mt = (get_mt_data(part='raw_chr20')
              .sample_rows(0.1)
              )
    else:
        logger.info('Running pipeline on MatrixTable wih adjusted genotypes...')
        ds = args.exome_cohort
        mt = hl.read_matrix_table(get_qc_mt_path(dataset=ds,
                                                 part='unphase_adj_genotypes',
                                                 split=True))

    # 1. Sample-QC filtering
    if not args.skip_sample_qc_filtering:
        logger.info('Applying per sample QC filtering...')

        mt = apply_sample_qc_filtering(mt)

        logger.info('Writing sample qc-filtered MT to disk...')
        mt = (mt
              .write(f'{hdfs_dir}/chd_ukbb.sample_qc_filtered.mt',
                     overwrite=True)
              )

    # 2. Variant-QC filtering
    if not args.skip_variant_qc_filtering:

        logger.info('Applying per variant QC filtering...')

        if hl.hadoop_is_file(f'{hdfs_dir}/chd_ukbb.sample_qc_filtered.mt/_SUCCESS'):
            logger.info('Reading pre-existing sample qc-filtered MT...')
            mt = hl.read_matrix_table(
                f'{hdfs_dir}/chd_ukbb.sample_qc_filtered.mt'
            )
        mt = apply_variant_qc_filtering(mt)

        # write hard filtered MT to disk
        logger.info('Writing variant qc-filtered MT to disk...')
        mt = (mt
              .write(f'{hdfs_dir}/chd_ukbb.variant_qc_filtered.mt',
                     overwrite=True)
              )

    # 3. Annotate AFs

    # allelic frequency cut-off
    maf_cutoff = args.af_max_threshold

    if not args.skip_af_filtering:

        if hl.hadoop_is_file(f'{hdfs_dir}/chd_ukbb.variant_qc_filtered.mt/_SUCCESS'):
            logger.info('Reading pre-existing sample/variant qc-filtered MT...')
            mt = hl.read_matrix_table(
                f'{hdfs_dir}/chd_ukbb.variant_qc_filtered.mt'
            )

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

        logger.info('Writing qc-filtered MT filtered to external maf with to disk...')
        mt = (mt
              .write(f'{hdfs_dir}/chd_ukbb.qc_final.rare.mt',
                     overwrite=True)
              )

    # 4. ##### Burden Test ######

    logger.info('Running burden test...')

    if hl.hadoop_is_file(f'{hdfs_dir}/chd_ukbb.qc_final.rare.mt/_SUCCESS'):
        logger.info('Reading pre-existing sample/variant qc-filtered MT with rare variants...')
        mt = hl.read_matrix_table(
            f'{hdfs_dir}/chd_ukbb.qc_final.rare.mt'
        )

    ## Add VEP-annotated fields
    vep_ht = get_vep_annotation_ht()

    mt = (mt
          .annotate_rows(LoF=vep_ht[mt.row_key].vep.LoF,
                         Consequence=vep_ht[mt.row_key].vep.Consequence,
                         DOMAINS=vep_ht[mt.row_key].vep.DOMAINS,
                         SYMBOL=vep_ht[mt.row_key].vep.SYMBOL)
          )

    ## Filter to bi-allelic variants
    if args.filter_biallelic:
        logger.info('Running burden test on biallelic variants...')
        mt = mt.filter_rows(bi_allelic_expr(mt))

    ## Filter to variants within protein domain(s)
    if args.filter_protein_domain:
        logger.info('Running burden test on variants within protein domain(s)...')
        mt = mt.filter_rows(
            vep_protein_domain_filter_expr(mt.DOMAINS),
            keep=True
        )

    ## Filter to variants within CCRs
    if args.filter_ccr:
        logger.info('Running burden test on variants within CCRs...')
        mt = filter_ccr(mt, ccr_pct=args.ccr_pct_cutoff)

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

    ## Classify variant into (major) consequence groups
    score_expr_ann = {'lof': hl.set(['LC', 'HC']).contains(mt.LoF),
                      'syn': mt.Consequence == 'synonymous_variant',
                      'miss': mt.Consequence == 'missense_variant',
                      'hcLOF': mt.LoF == 'HC'
                      }

    # Update dict expr annotations with combinations of variant consequences categories
    score_expr_ann.update(
        {'missC': (hl.sum([(mt['vep.MVP_score'] >= mvp_threshold),
                           (mt['vep.REVEL_score'] >= revel_threshold),
                           (mt['vep.CADD_PHRED'] >= cadd_threshold)]) >= 2) & score_expr_ann.get('miss')}
    )

    if args.test_cadd_miss:
        score_expr_ann.update(
            {'cadd_miss': (mt['vep.CADD_PHRED'] >= cadd_threshold) & score_expr_ann.get('miss')}
        )

    if args.test_mvp_miss:
        score_expr_ann.update(
            {'mvp_miss': (mt['vep.MVP_score'] >= mvp_threshold) & score_expr_ann.get('miss')}
        )

    if args.test_revel_miss:
        score_expr_ann.update(
            {'revel_miss': (mt['vep.REVEL_score'] >= revel_threshold) & score_expr_ann.get('miss')}
        )

    if args.test_pav_hc:
        score_expr_ann.update(
            {'hcPAV': score_expr_ann.get('hcLOF') | score_expr_ann.get('missC')}
        )

    mt = (mt
          .annotate_rows(csq_group=score_expr_ann)
          )

    # Transmute csq_group and convert dict to set where the group is defined
    # (easier to explode and grouping later)
    mt = (mt
          .transmute_rows(csq_group=hl.set(hl.filter(lambda x:
                                                     mt.csq_group.get(x),
                                                     mt.csq_group.keys())))
          )

    mt = (mt
          .filter_rows(hl.len(mt.csq_group) > 0)
          )

    # Explode nested csq_group before grouping
    mt = (mt
          .explode_rows(mt.csq_group)
          )

    # print('Number of samples/variants: ')
    # print(mt.count())

    # Group mt by gene/csq_group.
    mt_grouped = (mt
                  .group_rows_by(mt['SYMBOL'], mt['csq_group'])
                  .aggregate(hets=hl.agg.any(mt.GT.is_het()),
                             homs=hl.agg.any(mt.GT.is_hom_var()),
                             chets=hl.agg.count_where(mt.GT.is_het()) >= 2,
                             homs_chets=(hl.agg.count_where(mt.GT.is_het()) >= 2) | (hl.agg.any(mt.GT.is_hom_var()))
                             )
                  .repartition(100)
                  .persist()
                  )
    mts = []

    if args.homs:
        # select homs genotypes.

        mt_homs = (mt_grouped
                   .select_entries(mac=mt_grouped.homs)
                   .annotate_rows(agg_genotype='homs')
                   )

        mts.append(mt_homs)

    if args.chets:
        # select compound hets (chets) genotypes.
        mt_chets = (mt_grouped
                    .select_entries(mac=mt_grouped.chets)
                    .annotate_rows(agg_genotype='chets')
                    )

        mts.append(mt_chets)

    if args.homs_chets:
        # select chets and/or homs genotypes.
        mt_homs_chets = (mt_grouped
                         .select_entries(mac=mt_grouped.homs_chets)
                         .annotate_rows(agg_genotype='homs_chets')
                         )

        mts.append(mt_homs_chets)

    if args.hets:
        # select hets genotypes
        mt_hets = (mt_grouped
                   .select_entries(mac=mt_grouped.hets)
                   .annotate_rows(agg_genotype='hets')
                   )

        mts.append(mt_hets)

    ## Joint MatrixTables
    mt_grouped = hl.MatrixTable.union_rows(*mts)

    # Generate table of counts
    tb_gene = (
        mt_grouped
        .annotate_rows(
            n_cases=hl.agg.filter(mt_grouped['phe.is_case'], hl.agg.sum(mt_grouped.mac)),
            n_syndromic=hl.agg.filter(mt_grouped['phe.is_syndromic'], hl.agg.sum(mt_grouped.mac)),
            n_nonsyndromic=hl.agg.filter(mt_grouped['phe.is_nonsyndromic'], hl.agg.sum(mt_grouped.mac)),
            n_controls=hl.agg.filter(mt_grouped['phe.is_control'], hl.agg.sum(mt_grouped.mac)),
            n_total_cases=hl.agg.filter(mt_grouped['phe.is_case'], hl.agg.count()),
            n_total_syndromic=hl.agg.filter(mt_grouped['phe.is_syndromic'], hl.agg.count()),
            n_total_nonsyndromic=hl.agg.filter(mt_grouped['phe.is_nonsyndromic'], hl.agg.count()),
            n_total_controls=hl.agg.filter(mt_grouped['phe.is_control'], hl.agg.count()))
        .rows()
    )

    # run fet stratified by proband type
    analysis = ['all_cases', 'syndromic', 'nonsyndromic']

    tbs = []
    for proband in analysis:
        logger.info(f'Running test for {proband}...')
        colCases = None
        colTotalCases = None
        colControls = 'n_controls'
        colTotalControls = 'n_total_controls'
        if proband == 'all_cases':
            colCases = 'n_cases'
            colTotalCases = 'n_total_cases'
        if proband == 'syndromic':
            colCases = 'n_syndromic'
            colTotalCases = 'n_total_syndromic'
        if proband == 'nonsyndromic':
            colCases = 'n_nonsyndromic'
            colTotalCases = 'n_total_nonsyndromic'

        tb_fet = compute_fisher_exact(tb=tb_gene,
                                      n_cases_col=colCases,
                                      n_control_col=colControls,
                                      total_cases_col=colTotalCases,
                                      total_controls_col=colTotalControls,
                                      correct_total_counts=True,
                                      root_col_name='fet',
                                      extra_fields={'analysis': proband,
                                                    'maf': maf_cutoff})

        # filter out zero-count genes
        tb_fet = (tb_fet
                  .filter(hl.sum([tb_fet[colCases], tb_fet[colControls]]) > 0,
                          keep=True)
                  )

        tbs.append(tb_fet)

    tb_final = hl.Table.union(*tbs)

    tb_final.describe()

    # export results
    date = current_date()
    run_hash = str(uuid.uuid4())[:6]
    output_path = f'{args.output_dir}/{date}/{args.exome_cohort}.fet_burden.{run_hash}.ht'

    tb_final = (tb_final
                .checkpoint(output=output_path)
                )

    if args.write_to_file:
        # write table to disk as TSV file
        (tb_final
         .export(f'{output_path}.tsv')
         )

    hl.stop()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--exome_cohort', help="One of <chd_ukbb> or <chd_ddd>",
                        type=str, default=None)

    parser.add_argument('-o', '--output_dir', help='Path to output directory',
                        type=str, default=f'{nfs_dir}/hail_data/burden')

    parser.add_argument('--write_to_file', help='Write output to BGZ-compressed file',
                        action='store_true')

    parser.add_argument('--overwrite', help='Overwrite pre-existing data',
                        action='store_true')

    parser.add_argument('--default_ref_genome', help='Default reference genome to start Hail',
                        type=str, default='GRCh38')

    parser.add_argument('--run_test_mode', help='Run pipeline on smaller chunk of data (chr20) for testing propose',
                        action='store_true')

    # QC filtering option
    qc_filtering_group = parser.add_argument_group('Options for QC and filtering')

    qc_filtering_group.add_argument('--skip_sample_qc_filtering', help='Skip the sample QC filtering step',
                                    action='store_true')

    qc_filtering_group.add_argument('--skip_variant_qc_filtering', help='Skip the variant QC filtering step',
                                    action='store_true')

    qc_filtering_group.add_argument('--skip_af_filtering', help='Skip the allelic frequency filtering step',
                                    action='store_true')

    qc_filtering_group.add_argument('--af_max_threshold', help='Allelic frequency cutoff to filter variants (max)',
                                    type=float, default=0.001)

    qc_filtering_group.add_argument('--filter_biallelic', help='Run burden test on bi-allelic variants only',
                                    action='store_true')

    qc_filtering_group.add_argument('--filter_ccr',
                                    help='Run burden test on variants within CCRs only',
                                    action='store_true')

    qc_filtering_group.add_argument('--ccr_pct_cutoff', help='CCRs percentile cutoff',
                                    type=float, default=95.0)

    qc_filtering_group.add_argument('--filter_protein_domain',
                                    help='Run burden test on variants within protein domain(s) only',
                                    action='store_true')

    # Aggregate genotype options
    agg_genotype_group = parser.add_argument_group('Options to aggregate genotypes for burden testing')

    agg_genotype_group.add_argument('--hets', help='Aggregate hets genotypes to run the burden test',
                                    action='store_true')

    agg_genotype_group.add_argument('--homs', help='Aggregate homs genotypes to run the burden test',
                                    action='store_true')

    agg_genotype_group.add_argument('--chets', help='Aggregate compound hets genotypes to run the burden test',
                                    action='store_true')

    agg_genotype_group.add_argument('--homs_chets',
                                    help='Aggregate compound hets and/or homs genotypes to run the burden test',
                                    action='store_true')

    # Type of variants for burden testing
    variant_type_group = parser.add_argument_group('Type of variant to be included in burden testing')

    variant_type_group.add_argument('--test_cadd_miss',
                                    help='Test missense variants passing a predefined CADD score threshold',
                                    action='store_true')

    variant_type_group.add_argument('--test_revel_miss',
                                    help='Test missense variants passing a predefined REVEL score threshold',
                                    action='store_true')

    variant_type_group.add_argument('--test_mvp_miss',
                                    help='Test missense variants passing a predefined MVP score threshold',
                                    action='store_true')

    variant_type_group.add_argument('--test_pav_hc',
                                    help='Test high confidence protein-altering variants',
                                    action='store_true')

    # Pathogenic score cutoff
    score_cutoff_group = parser.add_argument_group('Pathogenic score threshold options')

    score_cutoff_group.add_argument('--cadd_threshold', help='CADD score threshold', type=float, default=20.0)

    score_cutoff_group.add_argument('--revel_threshold', help='REVEL score threshold', type=float, default=0.5)

    score_cutoff_group.add_argument('--mvp_threshold', help='MVP score threshold', type=float, default=0.8)

    args = parser.parse_args()

    if args.hets + args.homs + args.chets + args.homs_chets == 0:
        sys.exit('Specifies at least one of hets, homs, chets or homs_chets options...')

    main(args)
