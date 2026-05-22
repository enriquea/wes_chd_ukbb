"""
Apply hard filters to samples for sample QC.

Filters the input MatrixTable to high-callrate, common, bi-allelic SNPs,
imputes sex, flags ambiguous sex / aneuploidies, annotates each sample with
the set of hard filters it fails, and writes the result to a Hail Table.
"""
import argparse
import logging

import hail as hl

from utils.data_utils import (get_mt_data,
                              get_qc_mt_path,
                              get_sample_qc_ht_path)
from utils.expressions import bi_allelic_expr
from utils.generic import unphase_mt

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("unified_sample_qc_a")
logger.setLevel(logging.INFO)


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
        'coverage': ht.chr20_mean_coverage == 0,
        'ambiguous_sex': ht.ambiguous_sex,
        'sex_aneuploidy': ht.sex_aneuploidy,
    }

    return hl.set(hl.filter(lambda x: hl.is_defined(x),
                            [hl.or_missing(filter_expr, name) for name, filter_expr in hard_filters.items()]))


def filter_and_write_qc_mt(exome_cohort: str, overwrite: bool) -> None:
    """Filter input MT to bi-allelic, high-callrate, common SNPs and write to disk."""
    # import unfiltered MT
    mt = get_mt_data(dataset=exome_cohort, part='raw')

    logger.info("Filtering to bi-allelic, high-callrate, common SNPs for sample QC...")

    mt = mt.filter_rows(bi_allelic_expr(mt) &
                        hl.is_snp(mt.alleles[0], mt.alleles[1]) &
                        (hl.agg.mean(mt.GT.n_alt_alleles()) / 2 > 0.001) &
                        (hl.agg.fraction(hl.is_defined(mt.GT)) > 0.99))

    logger.info(f"Writing filtered MT with {mt.count_rows()} variants...")
    (mt
     .annotate_cols(callrate=hl.agg.fraction(hl.is_defined(mt.GT)))
     .naive_coalesce(1000)
     .write(get_qc_mt_path(dataset=exome_cohort, part='high_callrate_common_snp_biallelic', split=True),
            overwrite=overwrite)
     )


def load_qc_mt_with_metadata(exome_cohort: str) -> hl.MatrixTable:
    """Load QC MatrixTable, unphase it, and annotate columns with sex-chrom coverage metadata."""
    qc_mt = hl.read_matrix_table(get_qc_mt_path(dataset=exome_cohort, part='high_callrate_common_snp_biallelic', split=True))

    # unphase MT. required for imputing sex...
    qc_mt = unphase_mt(qc_mt)

    logger.info("Importing sample metadata...")
    meta_ht = (hl.read_table(get_sample_qc_ht_path(part='sex_chrom_coverage'))
               .key_by('s')
               )

    qc_mt = qc_mt.annotate_cols(**meta_ht[qc_mt.s])
    return qc_mt


def compute_sex_annotations(qc_mt: hl.MatrixTable) -> hl.Table:
    """Infer sex and annotate ambiguous_sex and sex_aneuploidy flags; return cols Table."""
    logger.info("Inferring sex...")
    qc_ht = annotate_sex(qc_mt,
                         female_threshold=0.5,
                         male_threshold=0.6).cols()

    # Flag individuals with ambiguous sex and aneuploidies
    qc_ht = qc_ht.annotate(ambiguous_sex=((qc_ht.f_stat >= 0.5) & (hl.is_defined(qc_ht.chrY_normalized_coverage) &
                                                                   (qc_ht.chrY_normalized_coverage <= 0.1))) |
                                         (hl.is_missing(qc_ht.f_stat)) |
                                         ((qc_ht.f_stat >= 0.4) & (qc_ht.f_stat <= 0.6) &
                                          (hl.is_defined(qc_ht.chrY_normalized_coverage) & (
                                                  qc_ht.chrY_normalized_coverage > 0.1))),
                           sex_aneuploidy=(qc_ht.f_stat < 0.4) & hl.is_defined(qc_ht.chrY_normalized_coverage) & (
                                   qc_ht.chrY_normalized_coverage > 0.1))
    return qc_ht


def annotate_hard_filters(qc_ht: hl.Table) -> hl.Table:
    """Annotate table with hard_filters set and inferred sex label, keyed by sample ID."""
    logger.info("Annotating samples failing hard filters...")

    sex_expr = (hl.case()
                .when(qc_ht.ambiguous_sex, "ambiguous_sex")
                .when(qc_ht.sex_aneuploidy, "sex_aneuploidy")
                .when(qc_ht.is_female, "female")
                .default("male"))

    qc_ht = qc_ht.annotate(hard_filters=make_hard_filters_expr(qc_ht),
                           sex=sex_expr).key_by('s')
    return qc_ht


def write_qc_ht(qc_ht: hl.Table, overwrite: bool, write_to_file: bool) -> None:
    """Describe, write hard-filters Table, and optionally export as BGZ-compressed TSV."""
    qc_ht.describe()

    qc_ht.write(get_sample_qc_ht_path(part='hard_filters'),
                overwrite=overwrite)

    if write_to_file:
        qc_ht.export(get_sample_qc_ht_path(part='hard_filters') + '.tsv.bgz')


def log_filter_counts() -> None:
    """Read back the hard-filters Table and log pre/post-filter sample counts."""
    # Export annotations to make rank list for relatedness (in final sample QC)

    # Check numbers:
    qc_ht = hl.read_table(get_sample_qc_ht_path(part='hard_filters'))
    sample_count = qc_ht.count()
    checkpoint1a = qc_ht.aggregate(hl.agg.count_where(hl.len(qc_ht['hard_filters']) == 0))
    logger.info('{} samples found before filtering'.format(sample_count))
    logger.info('{} samples found after checkpoint 1a (hard filters)'.format(checkpoint1a))


def main(args):
    hl.init(default_reference=args.default_reference)

    logger.info("Importing data...")
    # read MT
    if not args.skip_filter_step:
        filter_and_write_qc_mt(exome_cohort=args.exome_cohort, overwrite=args.overwrite)

    qc_mt = load_qc_mt_with_metadata(exome_cohort=args.exome_cohort)

    qc_ht = compute_sex_annotations(qc_mt)

    qc_ht = annotate_hard_filters(qc_ht)

    write_qc_ht(qc_ht, overwrite=args.overwrite, write_to_file=args.write_to_file)

    log_filter_counts()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--exome_cohort', help="One of <chd_ukbb> or <chd_ddd>",
                        type=str, default=None)
    parser.add_argument('--skip_filter_step',
                        help='Skip filtering the input MT...will load already existing filtered MT with high-quality variants',
                        action='store_true')
    parser.add_argument('--overwrite', help='Overwrite pre-existing data', action='store_true')
    parser.add_argument('--default_reference', help='One of GRCh37 and GRCh38', type=str, default='GRCh38')
    parser.add_argument('--write_to_file', help='Optionally, write output HT to TSV file.', action='store_true')

    args = parser.parse_args()

    main(args)
