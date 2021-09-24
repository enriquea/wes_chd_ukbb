# eam
# 2021-09-08

"""

Run case-control burden test (Fisher Exact) stratified by variant functional category and proband
syndromic status.

usage: gene_burden_fet.py [-h] [-i EXOME_COHORT] [-o OUTPUT_DIR]
                          [--skip_sample_qc_filtering]
                          [--skip_variant_qc_filtering] [--skip_af_filtering]
                          [--af_max_threshold AF_MAX_THRESHOLD]
                          [--filter_biallelic] [--filter_protein_domain]
                          [--write_to_file] [--overwrite]
                          [--default_ref_genome DEFAULT_REF_GENOME]
                          [--run_test_mode]

optional arguments:
  -h, --help            show this help message and exit
  -i EXOME_COHORT, --exome_cohort EXOME_COHORT
                        One of <chd_ukbb> or <chd_ddd>
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
  --write_to_file       Write output to BGZ-compressed file
  --overwrite           Overwrite pre-existing data
  --default_ref_genome DEFAULT_REF_GENOME
                        Default reference genome to start Hail
  --run_test_mode       Run pipeline on smaller chunk of data (chr20) for
                        testing propose


"""

import functools
import operator
import argparse
import logging
import uuid

import hail as hl

from utils.data_utils import (get_af_annotation_ht,
                              get_sample_meta_data,
                              get_qc_mt_path,
                              get_vep_scores_ht,
                              get_variant_qc_ht_path,
                              get_sample_qc_ht_path,
                              get_mt_data, get_vep_annotation_ht)

from utils.filter import (filter_low_conf_regions,
                          filter_to_autosomes,
                          filter_to_cds_regions,
                          filter_capture_intervals,
                          remove_telomeres_centromes
                          )

from utils.generic import current_date
from utils.stats import compute_fisher_exact
from utils.expressions import (af_filter_expr,
                               bi_allelic_expr)

from utils.vep import vep_protein_domain_filter_expr


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("Burden testing pipeline")
logger.setLevel(logging.INFO)

nfs_dir = 'file:///home/ubuntu/data'
nfs_tmp = 'file:///home/ubuntu/data/tmp'
hdfs_dir = 'hdfs://spark-master:9820/dir'

# Deleterious scores cutoff
MVP_THRESHOLD = 0.8
REVEL_THRESHOLD = 0.5
CADD_THRESHOLD = 25
MPC_THRESHOLD = 2


def apply_sample_qc_filtering(mt: hl.MatrixTable,
                              keep_rare_variants: bool = True,
                              maf_threshold: float = 0.01) -> hl.MatrixTable:
    """
    Apply sample QC filtering, compute internal allelic frequencies on samples passing qc and
    adjusted phenotypes. Optionally, return MT filtered to rare variants.

    :param mt: hl.MatrixTable
    :param keep_rare_variants: Filter MT to rare variants
    :param maf_threshold: allelic frequency cutoff
    :return: hl.MatrixTable
    """
    # import variant qc final table
    sample_qc_ht = hl.read_table(
        get_sample_qc_ht_path(part='final_qc')
    )
    sample_qc_ht = (sample_qc_ht
                    .filter(sample_qc_ht.pass_filters)
                    )
    mt = (mt
          .filter_cols(hl.is_defined(sample_qc_ht[mt.col_key]))
          )
    # compute cohort-specific (internal) allelic frequencies on samples passing qc
    mt = (mt
          .annotate_rows(gt_stats=hl.agg.call_stats(mt.GT, mt.alleles))
          )
    mt = (mt
          .annotate_rows(internal_af=mt.gt_stats.AF[1],
                         internal_ac=mt.gt_stats.AC[1])
          )
    # filter out common variants base don internal af
    if keep_rare_variants:
        mt = (mt
              .filter_rows(af_filter_expr(mt, 'internal_af', maf_threshold))
              )

    return mt


def apply_variant_qc_filtering(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Apply variant QC filtering

    :param mt: hl.MatrixTable
    :return: hl.MatrixTable
    """
    # import variant qc final table
    variant_qc_ht = hl.read_table(
        get_variant_qc_ht_path(part='final_qc')
    )
    mt = (mt
          .annotate_rows(**variant_qc_ht[mt.row_key])
          )
    mt = (mt
          .filter_rows(~mt.fail_inbreeding_coeff &
                       ~mt.AC0 &
                       ~mt.fail_vqsr &
                       ~mt.fail_rf &
                       mt.is_coveraged_gnomad_genomes &
                       mt.is_defined_capture_intervals)
          )
    # filter low conf regions
    mt = filter_low_conf_regions(
        mt,
        filter_lcr=True,  # TODO: include also decoy and low coverage exome regions
        filter_segdup=True
    )
    # filter telomeres/centromes
    mt = remove_telomeres_centromes(mt)

    return mt


def main(args):
    hl.init(default_reference=args.default_ref_genome)

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

        logger.info('Writing sample qc-filtered mt with rare variants (internal maf 0.01) to disk...')
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
        logger.info('Writing variant qc-filtered mt with rare variants (internal maf 0.01) to disk...')
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
    score_expr_ann = {'hcLOF': mt.LoF == 'HC',
                      'syn': mt.Consequence == 'synonymous_variant',
                      'miss': mt.Consequence == 'missense_variant'
                      }

    # Update dict expr annotations with combinations of variant consequences categories
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
                  .group_rows_by(mt.csq_group, mt.SYMBOL)
                  .aggregate(n_het=hl.agg.count_where(mt.GT.is_het()))
                  .repartition(500)
                  .persist()
                  )

    # Generate table of counts
    tb_gene = (
        mt_grouped
            .annotate_rows(
            n_het_cases=hl.agg.filter(mt_grouped['phe.is_case'], hl.agg.count_where(mt_grouped.n_het > 0)),
            n_het_syndromic=hl.agg.filter(mt_grouped['phe.is_syndromic'], hl.agg.count_where(mt_grouped.n_het > 0)),
            n_het_nonsyndromic=hl.agg.filter(mt_grouped['phe.is_nonsyndromic'],
                                             hl.agg.count_where(mt_grouped.n_het > 0)),
            n_het_controls=hl.agg.filter(mt_grouped['phe.is_control'], hl.agg.count_where(mt_grouped.n_het > 0)),
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
        colControls = 'n_het_controls'
        colTotalControls = 'n_total_controls'
        if proband == 'all_cases':
            colCases = 'n_het_cases'
            colTotalCases = 'n_total_cases'
        if proband == 'syndromic':
            colCases = 'n_het_syndromic'
            colTotalCases = 'n_total_syndromic'
        if proband == 'nonsyndromic':
            colCases = 'n_het_nonsyndromic'
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
