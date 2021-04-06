# eam
# 2021-03-30

"""
This pipeline compute coverage summary stats per site/variants across samples.
Optionality, the computation can be stratified by phenotype status (i.e. case/control).
If specified, a binomial test is run to access independence between the case/control status
and coverage.


usage: compute_coverage.py [-h] [--mt_input_path MT_INPUT_PATH]
                                [--ht_output_path HT_OUTPUT_PATH]
                                [--write_to_file] [--overwrite]
                                [--default_ref_genome DEFAULT_REF_GENOME]
                                [--compute_overall_coverage]
                                [--compute_phe_coverage]
                                [--run_binomial_test] [--phe_field PHE_FIELD]
                                [--pvalue_threshold PVALUE_THRESHOLD]
                                [--min_sample_proportion MIN_SAMPLE_PROPORTION]

optional arguments:
  -h, --help            show this help message and exit
  --mt_input_path MT_INPUT_PATH
                        Path to input Hail MatrixTable
  --ht_output_path HT_OUTPUT_PATH
                        Path to output HailTable with coverage stats
  --write_to_file       Write output to BGZ-compressed file
  --overwrite           Overwrite results in HT output path if already
                        exists...
  --default_ref_genome DEFAULT_REF_GENOME
                        Default reference genome to start Hail
  --compute_overall_coverage
                        Compute coverage stats over all samples.
  --compute_phe_coverage
                        Compute coverage stats stratified by phenotype.
                        Expected binary phenotype.
  --run_binomial_test   If compute_phe_coverage specified, run a binomial test
                        to access independence between the case/control status
                        and coverage.
  --phe_field PHE_FIELD
                        The phenotype field in the input MT. Expected to by
                        boolean (false: control, true: case)
  --pvalue_threshold PVALUE_THRESHOLD
                        The significant p-value threshold to access
                        independence between the case/control status and
                        coverage for every site/variant in the input MT
                        (binomial test).
  --min_sample_proportion MIN_SAMPLE_PROPORTION
                        The (minimal) proportion of samples with coverage
                        equal or higher than certain level (default 10X) to
                        consider a site well covered.

"""

import argparse
import logging
import sys
from typing import List

import hail as hl

from utils.generic import current_date
from utils.data_utils import get_sample_meta_data

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

nfs_dir = 'file:///home/ubuntu/data'


# NOTE: This function was modified from gnomad_methods module to work on dense MT, rather than sparse MT.
# original version: https://github.com/broadinstitute/gnomad_methods/blob/
# a5b107aa4e3d2715475a70c6c2ee8482b30f32df/gnomad/utils/sparse_mt.py#L651
def compute_coverage_stats(
        mt: hl.MatrixTable,
        reference_ht: hl.Table,
        coverage_over_x_bins: List[int] = [1, 5, 10, 15, 20, 25, 30, 50, 100],
) -> hl.Table:
    """
    Computes the following coverage statistics for every base of the `reference_ht` provided:
        - mean
        - median
        - total DP
        - fraction of samples with coverage above X, for each x in `coverage_over_x_bins`
    The `reference_ht` is a table that contains row for each locus coverage should be computed on.
    It needs to be keyed with the same keys as `mt`, typically either `locus` or `locus, alleles`.
    The `reference_ht` can e.g. be created using `get_reference_ht`
    :param mt: Input sparse MT
    :param reference_ht: Input reference HT
    :param coverage_over_x_bins: List of boundaries for computing samples over X
    :return: Table with per-base coverage stats
    """

    n_samples = mt.count_cols()
    print(f"Computing coverage stats on {n_samples} samples.")

    # Create an outer join with the reference Table
    mt = mt.select_entries("DP").select_cols().select_rows()
    col_key_fields = list(mt.col_key)
    t = mt._localize_entries("__entries", "__cols")
    t = t.join(reference_ht.key_by(*mt.row_key).select(_in_ref=True), how="outer")
    t = t.annotate(
        __entries=hl.or_else(
            t.__entries,
            hl.range(n_samples).map(lambda x: hl.null(t.__entries.dtype.element_type)),
        )
    )
    mt = t._unlocalize_entries("__entries", "__cols", col_key_fields)

    # Densify
    # mt = hl.experimental.densify(mt)

    # Filter rows where the reference is missing
    mt = mt.filter_rows(mt._in_ref)

    # Unfilter entries so that entries with no ref block overlap aren't null
    mt = mt.unfilter_entries()

    # Compute coverage stats
    coverage_over_x_bins = sorted(coverage_over_x_bins)
    max_coverage_bin = coverage_over_x_bins[-1]
    hl_coverage_over_x_bins = hl.array(coverage_over_x_bins)

    # This expression creates a counter DP -> number of samples for DP between 0 and max_coverage_bin
    coverage_counter_expr = hl.agg.counter(
        hl.min(max_coverage_bin, hl.or_else(mt.DP, 0))
    )

    # This expression aggregates the DP counter in reverse order of the coverage_over_x_bins
    # and computes the cumulative sum over them.
    #  It needs to be in reverse order because we want the sum over samples covered by > X.
    count_array_expr = hl.cumulative_sum(
        hl.array(
            [
                hl.int32(coverage_counter_expr.get(max_coverage_bin, 0))
            ]
            # The coverage was already floored to the max_coverage_bin, so no more aggregation is needed for the max bin
        ).extend(  # For each of the other bins, coverage needs to be summed between the boundaries
            hl.range(hl.len(hl_coverage_over_x_bins) - 1, 0, step=-1).map(
                lambda i: hl.sum(
                    hl.range(
                        hl_coverage_over_x_bins[i - 1], hl_coverage_over_x_bins[i]
                    ).map(lambda j: hl.int32(coverage_counter_expr.get(j, 0)))
                )
            )
        )
    )
    mean_expr = hl.agg.mean(hl.or_else(mt.DP, 0))

    # Annotate rows now
    return mt.select_rows(
        n_samples=n_samples,  # added @eam
        mean=hl.cond(hl.is_nan(mean_expr), 0, mean_expr),
        median_approx=hl.or_else(hl.agg.approx_median(hl.or_else(mt.DP, 0)), 0),
        total_DP=hl.agg.sum(mt.DP),
        **{
            f"over_{x}": count_array_expr[i] / n_samples
            for i, x in zip(
            range(
                len(coverage_over_x_bins) - 1, -1, -1
            ),  # Reverse the bin index as count_array_expr has the reverse order
            coverage_over_x_bins,
        )
        },
    ).rows()


def main(args):
    # init hail
    hl.init(default_reference=args.default_ref_genome)

    # import MT
    mt = hl.read_matrix_table(args.mt_input_path)

    n_variants, n_samples = mt.count()

    # Getting variant table. Basically, a table keyed by <locus> or <locus, alleles>
    # with all variants in the dataset and no extra fields (a.k.a reference table).
    tb_variants = (mt
                   .select_rows()
                   .rows()
                   )

    # compute overall coverage
    if args.compute_overall_coverage:
        logger.info(f"Computing coverage stats for {n_variants} variant over {n_samples} samples...")
        ht_cov_overall = compute_coverage_stats(mt=mt,
                                                reference_ht=tb_variants)

        tb_variants = (tb_variants
                       .annotate(overall=ht_cov_overall[tb_variants.key])
                       )

    # compute coverage stratified by phenotype status (expected binary)
    # force the input MT to have a case_control bool filed (is_case)
    # ***
    if args.compute_phe_coverage:
        logger.info(f"Computing coverage stats stratified by phenotype status...")

        # Annotate sample meta info
        # Note: Temporal solution, better to import annotated MT
        mt = (mt.annotate_cols(**get_sample_meta_data()[mt.col_key]))

        mt = (mt
              .annotate_cols(case_control=hl.if_else(mt[args.phe_field],
                                                     'case',
                                                     'control'))
              )

        strata = (mt
                  .aggregate_cols(hl.agg.collect_as_set(mt['case_control']))
                  )

        dict_strata_ht = {s: compute_coverage_stats(mt=mt.filter_cols(mt['case_control'] == s),
                                                    reference_ht=tb_variants)
                          for s in strata}

        for k in dict_strata_ht.keys():
            _tb = dict_strata_ht.get(k)
            tb_variants = tb_variants.annotate(**{k: _tb[tb_variants.key]})

        if args.run_binomial_test:
            logger.info(f"Running binomial test...")
            # perform a binomial test on coverage and case/control status
            # DOI: https://doi.org/10.1002/acn3.582
            tb_binomial = (tb_variants
                           .annotate(n_cases_over_10=hl.int(tb_variants.case.over_10 * 100),
                                     n_controls_over_10=hl.int(tb_variants.control.over_10 * 100),
                                     total_cases=tb_variants.case.n_samples,
                                     total_controls=tb_variants.control.n_samples,
                                     )
                           .select('n_cases_over_10',
                                   'n_controls_over_10',
                                   'total_cases',
                                   'total_controls')
                           )

            binomial_expr = {'p_value':
                hl.binom_test(
                    x=tb_binomial.n_cases_over_10,
                    n=tb_binomial.n_cases_over_10 + tb_binomial.n_controls_over_10,
                    p=tb_binomial.total_cases / (tb_binomial.total_cases + tb_binomial.total_controls),
                    alternative='two.sided')
            }

            tb_binomial = (tb_binomial
                           .annotate(**binomial_expr)
                           )

            tb_variants = (
                tb_variants
                    .annotate(binomial_stats=tb_binomial[tb_variants.key])
            )

    # make coverage filter expressions
    # Note: the default number of reads is set to 10X
    logger.info(f"Assigning per site coverage filters...")

    significant_level = args.pvalue_threshold
    min_sample_prop = args.min_sample_proportion

    coverage_filter_dict_expr = {}

    if args.compute_overall_coverage:
        coverage_filter_dict_expr.update(
            {'overall_hard_cutoff':
                 hl.if_else((tb_variants.overall.over_10 >= min_sample_prop),
                            "pass",
                            "fail")}
        )
    if args.compute_phe_coverage:
        # DOI: https://doi.org/10.1016/j.ajhg.2018.08.016
        coverage_filter_dict_expr.update(
            {'phe_hard_cutoff':
                 hl.if_else((tb_variants.case.over_10 >= min_sample_prop) &
                            (tb_variants.control.over_10 >= min_sample_prop),
                            "concordant",
                            "discordant")}
        )
    if args.run_binomial_test:
        coverage_filter_dict_expr.update(
            {'phe_binomial':
                 hl.if_else(tb_variants.binomial_stats.p_value < significant_level,
                            'dependent',
                            'independent')}

        )

    # annotate coverage filters
    tb_variants = (tb_variants
                   .annotate(coverage_filter=hl.struct(**coverage_filter_dict_expr))
                   )

    # add useful global annotations to final coverage stats ht
    # as well as affected/non-affected summary counts per filters
    global_ann_dict_expr = {
        'date': current_date(),
        'mt_path': args.mt_input_path,
        'min_sample_prop': min_sample_prop}
    if args.compute_overall_coverage:
        global_ann_dict_expr.update({'overall_hard_cutoff':
            tb_variants.aggregate(
                hl.agg.counter(
                    tb_variants.coverage_filter.overall_hard_cutoff))}
        )
    if args.compute_phe_coverage:
        global_ann_dict_expr.update({'phe_hard_cutoff':
            tb_variants.aggregate(
                hl.agg.counter(
                    tb_variants.coverage_filter.phe_hard_cutoff))}
        )
    if args.run_binomial_test:
        global_ann_dict_expr.update({'phe_binomial':
            tb_variants.aggregate(
                hl.agg.counter(
                    tb_variants.coverage_filter.phe_binomial)),
            'binomial_pvalue_cutoff': significant_level if args.run_binomial_test else hl.float('')}
        )

    tb_variants = (tb_variants
                   .annotate_globals(**global_ann_dict_expr)
                   )

    # check
    tb_variants.globals.show()
    tb_variants.describe()

    # write HT
    tb_variants = tb_variants.checkpoint(
        output=args.ht_output_path,
        overwrite=args.overwrite)

    # export to file if true
    if args.write_to_file:
        (tb_variants
         .export(f'{args.ht_output_path}.tsv.bgz')
         )

    hl.stop()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--mt_input_path',
                        help='Path to input Hail MatrixTable')

    parser.add_argument('--ht_output_path', help='Path to output HailTable with coverage stats')

    parser.add_argument('--write_to_file', help='Write output to BGZ-compressed file',
                        action='store_true')

    parser.add_argument('--overwrite', help='Overwrite results in HT output path if already exists...',
                        action='store_true')

    parser.add_argument('--default_ref_genome', help='Default reference genome to start Hail',
                        type=str, default='GRCh38')

    parser.add_argument('--compute_overall_coverage',
                        help='Compute coverage stats over all samples.',
                        action='store_true')

    parser.add_argument('--compute_phe_coverage',
                        help='Compute coverage stats stratified by phenotype. Expected binary phenotype.',
                        action='store_true')

    parser.add_argument('--run_binomial_test',
                        help='If compute_phe_coverage specified, run a binomial test to access independence between \
                              the case/control status and coverage.',
                        action='store_true')

    parser.add_argument('--phe_field',
                        help='The phenotype field in the input MT. Expected to by boolean (false: control, true: case)',
                        type=str,
                        default='phe.is_case')

    parser.add_argument('--pvalue_threshold',
                        help='The significant p-value threshold to access independence between the case/control status \
                              and coverage for every site/variant in the input MT (binomial test).',
                        type=float,
                        default=0.001)

    parser.add_argument('--min_sample_proportion',
                        help='The (minimal) proportion of samples with coverage equal or higher \
                              than certain level (default 10X) to consider a site well covered.',
                        type=float,
                        default=0.9)

    args = parser.parse_args()

    if args.compute_overall_coverage + args.compute_phe_coverage == 0:
        sys.exit('One of --compute_overall_coverage or --compute_phe_coverage must be specified...')

    main(args)
