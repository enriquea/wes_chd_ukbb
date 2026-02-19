# eam
# 2021-03-30

"""
This pipeline compute coverage summary stats per site/variants across samples.
Optionality, the computation can be stratified by phenotype status (i.e. case/control).
If specified, a binomial test is run to access independence between the case/control status
and coverage.


usage: compute_coverage.py [-h] [-i EXOME_COHORT] [--write_to_file]
                           [--overwrite]
                           [--default_ref_genome DEFAULT_REF_GENOME]
                           [--skip_sample_qc_filtering]
                           [--compute_overall_coverage]
                           [--compute_phe_coverage] [--run_binomial_test]
                           [--phe_field PHE_FIELD]
                           [--pvalue_threshold PVALUE_THRESHOLD]
                           [--min_sample_proportion MIN_SAMPLE_PROPORTION]
                           [--run_test_mode]

optional arguments:
  -h, --help            show this help message and exit
  -i EXOME_COHORT, --exome_cohort EXOME_COHORT
                        One of <chd_ukbb> or <chd_ddd>
  --write_to_file       Write output to BGZ-compressed file
  --overwrite           Overwrite results in HT output path if already
                        exists...
  --default_ref_genome DEFAULT_REF_GENOME
                        Default reference genome to start Hail
  --skip_sample_qc_filtering
                        Skip the sample QC filtering step
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
  --run_test_mode       Run pipeline on smaller chunk of data (chr20) for
                        testing propose


"""

import argparse
import logging
import sys
from typing import List

import hail as hl

from utils.generic import current_date
from utils.data_utils import (get_sample_meta_data,
                              get_mt_data,
                              get_qc_mt_path, get_variant_qc_ht_path)
from utils.qc import apply_sample_qc_filtering
from utils.config import NFS_DIR

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

nfs_dir = NFS_DIR


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


def load_mt(exome_cohort, run_test_mode, skip_sample_qc_filtering):
    """Load raw MatrixTable (test mode: chr20 sample) and optionally apply sample QC filtering."""
    if run_test_mode:
        logger.info('Running pipeline on test data...')
        mt = (get_mt_data(part='raw_chr20')
              .sample_rows(0.1)
              )
    else:
        logger.info('Running pipeline on raw MatrixTable...')
        mt = get_mt_data(part='raw')

    if not skip_sample_qc_filtering:
        logger.info('Applying per sample QC filtering...')
        mt = apply_sample_qc_filtering(mt)

    return mt


def compute_overall_coverage(mt, tb_variants):
    """Compute overall coverage stats and annotate the variant table."""
    n_variants, n_samples = mt.count()
    logger.info(f"Computing coverage stats for {n_variants} variant over {n_samples} samples...")
    ht_cov_overall = compute_coverage_stats(mt=mt, reference_ht=tb_variants)

    tb_variants = (tb_variants
                   .annotate(overall=ht_cov_overall[tb_variants.key])
                   )

    return tb_variants


def compute_phenotype_stratified_coverage(mt, tb_variants, phe_field):
    """Annotate cols with case/control label and compute per-stratum coverage stats."""
    mt = (mt.annotate_cols(**get_sample_meta_data()[mt.col_key]))

    mt = (mt
          .annotate_cols(case_control=hl.if_else(mt[phe_field],
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

    return tb_variants


def run_binomial_test(tb_variants):
    """Run a binomial test on coverage vs case/control status and annotate the variant table."""
    logger.info(f"Running binomial test...")
    tb_binomial = (tb_variants
                   .annotate(n_cases_over_10=hl.int(tb_variants.case.over_10 *
                                                    tb_variants.case.n_samples),
                             n_controls_over_10=hl.int(tb_variants.control.over_10 *
                                                       tb_variants.control.n_samples),
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

    return tb_variants


def build_coverage_filters(tb_variants, compute_overall, compute_phe, run_binom,
                           min_sample_prop, pvalue_threshold):
    """Build and annotate per-site coverage filter expressions."""
    logger.info(f"Assigning per site coverage filters...")

    coverage_filter_dict_expr = {}

    if compute_overall:
        coverage_filter_dict_expr.update(
            {'overall_hard_cutoff':
                 hl.if_else((tb_variants.overall.over_10 >= min_sample_prop),
                            "pass",
                            "fail")}
        )
    if compute_phe:
        # DOI: https://doi.org/10.1016/j.ajhg.2018.08.016
        coverage_filter_dict_expr.update(
            {'phe_hard_cutoff':
                 hl.if_else((tb_variants.case.over_10 >= min_sample_prop) &
                            (tb_variants.control.over_10 >= min_sample_prop),
                            "concordant",
                            "discordant")}
        )
    if run_binom:
        coverage_filter_dict_expr.update(
            {'phe_binomial':
                 hl.if_else(tb_variants.binomial_stats.p_value < pvalue_threshold,
                            'dependent',
                            'independent')}

        )

    tb_variants = (tb_variants
                   .annotate(coverage_filter=hl.struct(**coverage_filter_dict_expr))
                   )

    return tb_variants


def annotate_global_stats(tb_variants, exome_cohort, min_sample_prop, pvalue_threshold,
                          compute_overall, compute_phe, run_binom):
    """Annotate global fields with run metadata and per-filter summary counts."""
    global_ann_dict_expr = {
        'date': current_date(),
        'cohort': exome_cohort,
        'min_sample_prop': min_sample_prop}
    if compute_overall:
        global_ann_dict_expr.update({'overall_hard_cutoff':
            tb_variants.aggregate(
                hl.agg.counter(
                    tb_variants.coverage_filter.overall_hard_cutoff))}
        )
    if compute_phe:
        global_ann_dict_expr.update({'phe_hard_cutoff':
            tb_variants.aggregate(
                hl.agg.counter(
                    tb_variants.coverage_filter.phe_hard_cutoff))}
        )
    if run_binom:
        global_ann_dict_expr.update({'phe_binomial':
            tb_variants.aggregate(
                hl.agg.counter(
                    tb_variants.coverage_filter.phe_binomial)),
            'binomial_pvalue_cutoff': pvalue_threshold if run_binom else hl.float('')}
        )

    tb_variants = (tb_variants
                   .annotate_globals(**global_ann_dict_expr)
                   )

    return tb_variants


def export_results(tb_variants, exome_cohort, overwrite, write_to_file):
    """Checkpoint the coverage table and optionally export as flattened BGZ-compressed TSV."""
    ht_output_path = get_variant_qc_ht_path(dataset=exome_cohort, part='coverage_stats')
    tb_variants = tb_variants.checkpoint(output=ht_output_path, overwrite=overwrite)

    if write_to_file:
        (tb_variants
         .flatten()
         .export(f'{ht_output_path}.tsv.bgz')
         )

    return tb_variants


def main(args):
    # init hail
    hl.init(default_reference=args.default_ref_genome)

    ds = args.exome_cohort

    mt = load_mt(ds, args.run_test_mode, args.skip_sample_qc_filtering)

    tb_variants = (mt
                   .select_rows()
                   .rows()
                   )

    if args.compute_overall_coverage:
        tb_variants = compute_overall_coverage(mt, tb_variants)

    if args.compute_phe_coverage:
        tb_variants = compute_phenotype_stratified_coverage(mt, tb_variants, args.phe_field)

        if args.run_binomial_test:
            tb_variants = run_binomial_test(tb_variants)

    tb_variants = build_coverage_filters(
        tb_variants,
        args.compute_overall_coverage,
        args.compute_phe_coverage,
        args.run_binomial_test,
        args.min_sample_proportion,
        args.pvalue_threshold,
    )

    tb_variants = annotate_global_stats(
        tb_variants, ds,
        args.min_sample_proportion,
        args.pvalue_threshold,
        args.compute_overall_coverage,
        args.compute_phe_coverage,
        args.run_binomial_test,
    )

    # check
    tb_variants.globals.show()
    tb_variants.describe()

    export_results(tb_variants, ds, args.overwrite, args.write_to_file)

    hl.stop()


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--exome_cohort', help="One of <chd_ukbb> or <chd_ddd>",
                        type=str, default='chd_ukbb')

    parser.add_argument('--write_to_file', help='Write output to BGZ-compressed file',
                        action='store_true')

    parser.add_argument('--overwrite', help='Overwrite results in HT output path if already exists...',
                        action='store_true')

    parser.add_argument('--default_ref_genome', help='Default reference genome to start Hail',
                        type=str, default='GRCh38')

    parser.add_argument('--skip_sample_qc_filtering', help='Skip the sample QC filtering step',
                        action='store_true')

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

    parser.add_argument('--run_test_mode', help='Run pipeline on smaller chunk of data (chr20) for testing propose',
                        action='store_true')

    args = parser.parse_args()

    if args.compute_overall_coverage + args.compute_phe_coverage == 0:
        sys.exit('One of --compute_overall_coverage or --compute_phe_coverage must be specified...')

    main(args)
