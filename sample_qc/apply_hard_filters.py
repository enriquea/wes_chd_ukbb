import argparse
import logging

import hail as hl
from typing import Union

# from resources.sample_qc import get_sample_metadata
# from utils.expressions import bi_allelic_expr
from resources.data_utils import get_mt_data, get_qc_mt_path, get_sample_qc_ht_path
from utils.expressions import bi_allelic_expr
from utils.generic import unphase_mt

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("unified_sample_qc_a")
logger.setLevel(logging.INFO)

# dir to write temp files
hdfs_tmp_dir = 'hdfs://spark-master:9820/tmp'


def annotate_sex(mt: hl.MatrixTable,
                 male_threshold: float = 0.8,
                 female_threshold: float = 0.5) -> hl.MatrixTable:
    """
    Imputes sex, exports data, and annotates mt with this data
    NOTE: Evaluated in R (plots) and decided on cutoff of F<0.5 for females and F>0.6 for males (default) for exomes
    :param female_threshold:
    :param male_threshold:
    :param MatrixTable mt: MT containing samples to be ascertained for sex
    :rtype: MatrixTable
    """
    mt = hl.filter_intervals(mt, [hl.parse_locus_interval('chrX')])
    sex_ht = hl.impute_sex(mt.GT,
                           aaf_threshold=0.01,
                           female_threshold=female_threshold,
                           male_threshold=male_threshold)
    sex_colnames = ['f_stat', 'is_female']
    sex_ht = sex_ht.select(*sex_colnames)
    mt = mt.annotate_cols(**sex_ht[mt.col_key])
    return mt


def make_hard_filters_expr(ht: hl.Table) -> hl.expr.SetExpression:
    """
    NOTE: additional metadata imported from 'sex_chrom_coverage' output
    :param: Table ht: input MT
    :return: SetExpression
    """
    hard_filters = {
        # 'contamination': ht.freemix > 0.05,
        # 'chimera': ht.pct_chimeras > 0.05,
        'callrate': ht.callrate < 0.85,
        'coverage': ht.chr20_mean_dp == 0,
        'ambiguous_sex': ht.ambiguous_sex,
        'sex_aneuploidy': ht.sex_aneuploidy,
    }

    return hl.set(hl.filter(lambda x: hl.is_defined(x),
                            [hl.or_missing(filter_expr, name) for name, filter_expr in hard_filters.items()]))


def main(args):

    hl.init(args.default_reference)

    logger.info("Importing data...")

    # read MT
    if not args.skip_filter_step:

        # import unfiltered MT
        mt = get_mt_data(dataset=args.exome_cohort, part='unfiltered')

        logger.info("Filtering to bi-allelic, high-callrate, common SNPs for sample QC...")

        mt = mt.filter_rows(bi_allelic_expr(mt) &
                            hl.is_snp(mt.alleles[0], mt.alleles[1]) &
                            (hl.agg.mean(mt.GT.n_alt_alleles()) / 2 > 0.001) &
                            (hl.agg.fraction(hl.is_defined(mt.GT)) > 0.99))

        mt = (mt
              .annotate_cols(callrate=hl.agg.fraction(hl.is_defined(mt.GT)))
              .naive_coalesce(1000)
              )

        if not args.skip_write_qc_mt:
            (mt
             .write(get_qc_mt_path(part='high_callrate_common_snp_biallelic', split=True),
                    overwrite=args.overwrite)
             )

    qc_mt = hl.read_matrix_table(get_qc_mt_path(part='high_callrate_common_snp_biallelic', split=True))

    # unphase MT. required to impute_sex...
    qc_mt = unphase_mt(qc_mt)

    logger.info("Importing metadata...")
    meta_ht = (hl.read_table(get_sample_qc_ht_path(part='sex_chromosome_coverage'))
               .key_by('s')
               )

    meta_ht = meta_ht.annotate(normalized_y_coverage=meta_ht.chrY_ploidy)

    qc_mt = qc_mt.annotate_cols(**meta_ht[qc_mt.s])

    logger.info("Inferring sex...")
    qc_ht = annotate_sex(qc_mt,
                         female_threshold=0.5,
                         male_threshold=0.6).cols()

    # Flag individuals with ambiguous sex and aneuploidies
    qc_ht = qc_ht.annotate(ambiguous_sex=((qc_ht.f_stat >= 0.5) & (hl.is_defined(qc_ht.normalized_y_coverage) &
                                          (qc_ht.normalized_y_coverage <= 0.1))) |
                                         (hl.is_missing(qc_ht.f_stat)) |
                                         ((qc_ht.f_stat >= 0.4) & (qc_ht.f_stat <= 0.6) &
                                          (hl.is_defined(qc_ht.normalized_y_coverage) & (
                                                 qc_ht.normalized_y_coverage > 0.1))),
                           sex_aneuploidy=(qc_ht.f_stat < 0.4) & hl.is_defined(qc_ht.normalized_y_coverage) & (
                                   qc_ht.normalized_y_coverage > 0.1))

    logger.info("Annotating samples failing hard filters...")

    sex_expr = (hl.case()
                .when(qc_ht.ambiguous_sex, "ambiguous_sex")
                .when(qc_ht.sex_aneuploidy, "sex_aneuploidy")
                .when(qc_ht.is_female, "female")
                .default("male"))

    qc_ht = qc_ht.annotate(hard_filters=make_hard_filters_expr(qc_ht),
                           sex=sex_expr).key_by('s')

    qc_ht.describe()

    qc_ht.write(get_sample_qc_ht_path(part='hard_filters'),
                overwrite=args.overwrite)

    if args.write_to_file:
        qc_ht.export(get_sample_qc_ht_path(part='hard_filters') + '.tsv.bgz')

    # Export annotations to make rank list for relatedness (in final sample QC)

    # Check numbers:
    qc_ht = hl.read_table(get_sample_qc_ht_path(part='hard_filters'))
    sample_count = qc_ht.count()
    checkpoint1a = qc_ht.aggregate(hl.agg.count_where(hl.len(qc_ht['hard_filters']) == 0))
    logger.info('{} samples found before filtering'.format(sample_count))
    logger.info('{} samples found after checkpoint 1a (hard filters)'.format(checkpoint1a))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--exome_cohort', help="One of <chd_ukbb> or <chd_ddd>",
                        type=str, default=None)
    parser.add_argument('--skip_filter_step',
                        help='Skip filtering the input MT...will load already existing filtered MT with high-quality variants',
                        action='store_true')
    parser.add_argument('--overwrite', help='Overwrite pre-existing data', action='store_true')
    parser.add_argument('--default_reference', help='One of GRCh37 and GRCh38', type=str, default='GRCh38')
    parser.add_argument('--skip_write_qc_mt',
                        help='Skip writing pre-calculated MatrixTable containing high-quality variants',
                        action='store_true')
    parser.add_argument('--write_to_file', help='Optionally, write output HT to TSV file.', action='store_true')

    args = parser.parse_args()

    main(args)
