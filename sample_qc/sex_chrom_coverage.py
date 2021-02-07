"""

pipeline modified from gnomad.methods.impute_sex_ploidy.
It computes per sample normalized mean coverage on
chromosome X and Y using an autosomal contig (e.g. chr20)

"""

import argparse
import logging
from typing import Optional, Union

import hail as hl

from resources.data_utils import get_mt_data, get_sample_qc_ht_path, get_interval_ht
from utils.expressions import bi_allelic_expr
from utils.reference_genome import get_reference_genome

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def compute_mean_coverage(
    mt: hl.MatrixTable,
    excluded_calling_intervals: Optional[hl.Table] = None,
    included_calling_intervals: Optional[hl.Table] = None,
    normalization_contig: str = "chr20",
    chr_x: Optional[str] = None,
    chr_y: Optional[str] = None,
) -> hl.Table:
    """
    Computes sex chromosome mean coverage on dense Matrix Table and normalizes using the coverage of an autosomal
    chromosome (by default chr20).

    :param mt: Input sparse Matrix Table
    :param excluded_calling_intervals: Optional table of intervals to exclude from the computation. If defined,
           variant falling in the interval are excluded when computing chromosome mean coverage.
    :param included_calling_intervals: Optional table of intervals to use in the computation. If defined,
           variant falling in the interval are included when computing chromosome mean coverage.
    :param normalization_contig: Which chromosome to normalize by (default: chr20)
    :param chr_x: Optional X Chromosome contig name (by default uses the X contig in the reference)
    :param chr_y: Optional Y Chromosome contig name (by default uses the Y contig in the reference)
    :return: Table with mean coverage over chromosomes 20, X and Y and normalized coverage.
    """

    ref = get_reference_genome(mt.locus, add_sequence=False)

    if chr_x is None:
        if len(ref.x_contigs) != 1:
            raise NotImplementedError(
                "Found {0} X chromosome contigs ({1}) in Genome reference. sparse_impute_sex_ploidy currently only "
                "supports a single X chromosome contig. Please use the `chr_x` argument to  specify which X chromosome"
                " contig to use ".format(
                    len(ref.x_contigs), ",".join(ref.x_contigs)
                )
            )
        chr_x = ref.x_contigs[0]
    if chr_y is None:
        if len(ref.y_contigs) != 1:
            raise NotImplementedError(
                "Found {0} Y chromosome contigs ({1}) in Genome reference. sparse_impute_sex_ploidy currently only "
                "supports a single Y chromosome contig. Please use the `chr_y` argument to  specify which Y chromosome "
                "contig to use ".format(
                    len(ref.y_contigs), ",".join(ref.y_contigs)
                )
            )
        chr_y = ref.y_contigs[0]

    def get_chr_coverage_ann(chrom: str) -> hl.Table:

        logger.info(f"Working on {chrom}...")

        chr_mt = hl.filter_intervals(mt, [hl.parse_locus_interval(chrom)])

        # filtering if exclude/include intervals are defined
        if included_calling_intervals is not None:
            logger.info(f"Filtering variants defined in calling interval {args.interval_to_include}...")
            total_variants = chr_mt.count_rows()
            chr_mt = chr_mt.filter_rows(
                hl.is_defined(included_calling_intervals[chr_mt.locus])
            )
            variant_in_interval = chr_mt.count_rows()
            logger.info(f"Including {variant_in_interval} out of {total_variants} defined in intervals...")

        if excluded_calling_intervals is not None:
            logger.info(f"Excluding variants defined in interval {args.interval_to_exclude}...")
            total_variants = chr_mt.count_rows()
            chr_mt = chr_mt.filter_rows(
                hl.is_missing(excluded_calling_intervals[chr_mt.locus])
            )
            excluded_variants = total_variants - chr_mt.count_rows()
            logger.info(f"Excluding {excluded_variants} out of {total_variants} defined in intervals...")

        # exclude sex chromosome pseudo-autosomal-region (PAR) from the computation
        if chrom in ref.x_contigs:
            chr_mt = chr_mt.filter_rows(chr_mt.locus.in_x_nonpar())

        if chrom in ref.y_contigs:
            chr_mt = chr_mt.filter_rows(chr_mt.locus.in_y_nonpar())

        return chr_mt.select_cols(
            **{
               f"{chrom}_mean_coverage": hl.agg.mean(chr_mt.DP),
               f"{chrom}_callrate": hl.agg.fraction(hl.is_defined(chr_mt.GT))
            }
        ).cols()

    normalization_chrom_dp = get_chr_coverage_ann(normalization_contig)
    chrX_dp = get_chr_coverage_ann(chr_x)
    chrY_dp = get_chr_coverage_ann(chr_y)

    ht = normalization_chrom_dp.annotate(
        **chrX_dp[normalization_chrom_dp.key], **chrY_dp[normalization_chrom_dp.key],
    )

    return ht.annotate(
        **{
            f"{chr_x}_normalized_coverage": ht[f"{chr_x}_mean_coverage"]
            / (ht[f"{normalization_contig}_mean_coverage"]),
            f"{chr_y}_normalized_coverage": ht[f"{chr_y}_mean_coverage"]
            / (ht[f"{normalization_contig}_mean_coverage"]),
        }
    )


def main(args):

    # nfs_dir = 'file:///home/ubuntu/data'

    hl.init(default_reference=args.default_reference)

    logger.info("Importing data...")

    # import unfiltered MT
    mt = get_mt_data(dataset=args.exome_cohort, part='unfiltered')

    # keep bi-allelic variants
    mt = (mt
          .filter_rows(bi_allelic_expr(mt), keep=True)
          )

    # read intervals for filtering variants (used mainly for exomes)
    def _get_interval_table(interval: str) -> Union[None, hl.Table]:
        return get_interval_ht(name=interval,
                               reference=args.default_reference) if interval is not None else interval

    ht = compute_mean_coverage(mt=mt,
                               normalization_contig=args.normalization_contig,
                               included_calling_intervals=_get_interval_table(args.interval_to_include),
                               excluded_calling_intervals=_get_interval_table(args.interval_to_exclude),
                               chr_x=args.chr_x,
                               chr_y=args.chr_y)

    logger.info("Exporting data...")

    # write HT
    output_ht_path = get_sample_qc_ht_path(part='sex_chrom_coverage')
    ht.write(output=output_ht_path,
             overwrite=args.overwrite)

    # export to file if true
    if args.write_to_file:
        (ht
         .export(f'{output_ht_path}.tsv.bgz')
         )

    hl.stop()

    print("Done!")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--exome_cohort', help="One of <chd_ukbb> or <chd_ddd>",
                        type=str, default=None)
    parser.add_argument('--normalization_contig', help='Autosomal chr used to coverage normalization. Default chr20',
                        type=str, default='chr20')
    parser.add_argument('--chr_x', help='Chromosome X contig name',
                        type=str, default='chrX')
    parser.add_argument('--chr_y', help='Chromosome Y contig name',
                        type=str, default='chrY')
    parser.add_argument('--interval_to_include', help='Optional table of intervals to use in the computation. If defined,'
                        'variant falling in the interval are included when computing chromosome mean coverage.',
                        type=str, default=None)
    parser.add_argument('--interval_to_exclude', help='Optional table of intervals to use in the computation. If defined,'
                        'variant falling in the interval are excluded when computing chromosome mean coverage.',
                        type=str, default=None)
    parser.add_argument('--default_reference', help='One of GRCh37 and GRCh38', type=str, default='GRCh38')
    parser.add_argument('--overwrite', help='Overwrite results in HT output path if already exists...',
                        action='store_true')
    parser.add_argument('--write_to_file', help='Optionally, write HT to TSV file.', action='store_true')

    args = parser.parse_args()

    main(args)
