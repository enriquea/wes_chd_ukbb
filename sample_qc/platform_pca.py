# Adapted from gnomad.methods

# Perform sample platform QC base on principal component analysis

# usage: platform_pca.py [-h] [--mt_input_path MT_INPUT_PATH]
#                        [--ht_output_path HT_OUTPUT_PATH]
#                        [--ht_intervals HT_INTERVALS] [--write_to_file]
#                        [--overwrite] [--default_ref_genome DEFAULT_REF_GENOME]
#                        [--binarization_threshold BINARIZATION_THRESHOLD]
#                        [--hdbscan_min_cluster_size HDBSCAN_MIN_CLUSTER_SIZE]
#                        [--hdbscan_min_samples HDBSCAN_MIN_SAMPLES]
#
# optional arguments:
#   -h, --help            show this help message and exit
#   --mt_input_path MT_INPUT_PATH
#                         Path to input Hail MatrixTable
#   --ht_output_path HT_OUTPUT_PATH
#                         Path to output HailTable with platform PCs
#   --ht_intervals HT_INTERVALS
#                         HailTable of intervals to use in the computation
#   --write_to_file       Write output to BGZ-compressed file
#   --overwrite           Overwrite results in HT output path if already
#                         exists...
#   --default_ref_genome DEFAULT_REF_GENOME
#                         Default reference genome to start Hail
#   --binarization_threshold BINARIZATION_THRESHOLD
#                         If set, the callrate MT is transformed to a 0/1 value
#                         MT based on the threshold. E.g. with the default
#                         threshold of 0.25, all entries with a callrate < 0.25
#                         are considered as 0s, others as 1s. Set to None if no
#                         threshold desired
#   --hdbscan_min_cluster_size HDBSCAN_MIN_CLUSTER_SIZE
#                         HDBSCAN `min_cluster_size` parameter. If not specified
#                         the smallest of 500 and 0.1*n_samples will be used.
#   --hdbscan_min_samples HDBSCAN_MIN_SAMPLES
#                         HDBSCAN `min_samples` parameter.
#

import logging
from typing import List, Optional, Tuple

import hail as hl
import numpy as np
import hdbscan
import argparse

from utils.expressions import bi_allelic_expr
from utils.filter import filter_to_autosomes, filter_genotypes_ab

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def compute_callrate_mt(
        mt: hl.MatrixTable,
        intervals_ht: hl.Table,
        bi_allelic_only: bool = True,
        autosomes_only: bool = True,
        match: bool = True,
) -> hl.MatrixTable:
    """
    Computes a sample/interval MT with each entry containing the call rate for that sample/interval.
    This can be used as input for imputing exome sequencing platforms.
    .. note::
        The input interval HT should have a key of type Interval.
        The resulting table will have a key of the same type as the `intervals_ht` table and
        contain an `interval_info` field containing all non-key fields of the `intervals_ht`.
    :param match:
    :param mt: Input MT
    :param intervals_ht: Table containing the intervals. This table has to be keyed by locus.
    :param bi_allelic_only: If set, only bi-allelic sites are used for the computation
    :param autosomes_only: If set, only autosomal intervals are used.
    :param matches: If set, returns all intervals in intervals_ht that overlap the locus in the input MT.
    :return: Callrate MT
    """
    logger.info("Computing call rate MatrixTable")

    if len(intervals_ht.key) != 1 or not isinstance(
            intervals_ht.key[0], hl.expr.IntervalExpression
    ):
        logger.warning(
            f"Call rate matrix computation expects `intervals_ht` with a key of type Interval. "
            f"Found: {intervals_ht.key}"
        )

    if autosomes_only:
        callrate_mt = filter_to_autosomes(mt)

    if bi_allelic_only:
        callrate_mt = callrate_mt.filter_rows(bi_allelic_expr(callrate_mt))

    intervals_ht = intervals_ht.annotate(_interval_key=intervals_ht.key)
    callrate_mt = callrate_mt.annotate_rows(
        _interval_key=intervals_ht.index(
            callrate_mt.locus, all_matches=match
        )._interval_key
    )

    if match:
        callrate_mt = callrate_mt.explode_rows("_interval_key")

    callrate_mt = callrate_mt.filter_rows(
        hl.is_defined(callrate_mt._interval_key.interval)
    )
    callrate_mt = callrate_mt.select_entries(
        GT=hl.or_missing(hl.is_defined(callrate_mt.GT), hl.struct())
    )
    callrate_mt = callrate_mt.group_rows_by(**callrate_mt._interval_key).aggregate(
        callrate=hl.agg.fraction(hl.is_defined(callrate_mt.GT))
    )
    intervals_ht = intervals_ht.drop("_interval_key")
    callrate_mt = callrate_mt.annotate_rows(
        interval_info=hl.struct(**intervals_ht[callrate_mt.row_key])
    )
    return callrate_mt


def run_platform_pca(
        callrate_mt: hl.MatrixTable, binarization_threshold: Optional[float] = 0.25
) -> Tuple[List[float], hl.Table, hl.Table]:
    """
    Runs a PCA on a sample/interval MT with each entry containing the call rate.
    When `binzarization_threshold` is set, the callrate is transformed to a 0/1 value based on the threshold.
    E.g. with the default threshold of 0.25, all entries with a callrate < 0.25 are considered as 0s, others as 1s.
    :param callrate_mt: Input callrate MT
    :param binarization_threshold: binzarization_threshold. None is no threshold desired
    :return: eigenvalues, scores_ht, loadings_ht
    """
    logger.info("Running platform PCA")

    if binarization_threshold is not None:
        callrate_mt = callrate_mt.annotate_entries(
            callrate=hl.int(callrate_mt.callrate > binarization_threshold)
        )
    # Center until Hail's PCA does it for you
    callrate_mt = callrate_mt.annotate_rows(
        mean_callrate=hl.agg.mean(callrate_mt.callrate)
    )
    callrate_mt = callrate_mt.annotate_entries(
        callrate=callrate_mt.callrate - callrate_mt.mean_callrate
    )
    eigenvalues, scores, loadings = hl.pca(
        callrate_mt.callrate, compute_loadings=True
    )  # TODO:  Evaluate whether computing loadings is a good / worthy thing
    logger.info("Platform PCA eigenvalues: {}".format(eigenvalues))

    return eigenvalues, scores, loadings


def assign_platform_from_pcs(
        platform_pca_scores_ht: hl.Table,
        pc_scores_ann: str = "scores",
        hdbscan_min_cluster_size: Optional[int] = None,
        hdbscan_min_samples: int = None,
) -> hl.Table:
    """
    Assigns platforms using HBDSCAN on the results of call rate PCA.
    :param platform_pca_scores_ht: Input table with the PCA score for each sample
    :param pc_scores_ann: Field containing the scores
    :param hdbscan_min_cluster_size: HDBSCAN `min_cluster_size` parameter. If not specified the smallest of 500 and 0.1*n_samples will be used.
    :param hdbscan_min_samples: HDBSCAN `min_samples` parameter
    :return: A Table with a `qc_platform` annotation containing the platform based on HDBSCAN clustering
    """

    logger.info("Assigning platforms based on platform PCA clustering")

    # Read and format data for clustering
    data = platform_pca_scores_ht.to_pandas()
    callrate_data = np.matrix(data[pc_scores_ann].tolist())
    logger.info("Assigning platforms to {} samples.".format(len(callrate_data)))

    # Cluster data
    if hdbscan_min_cluster_size is None:
        hdbscan_min_cluster_size = min(500, 0.1 * data.shape[0])
    clusterer = hdbscan.HDBSCAN(
        min_cluster_size=hdbscan_min_cluster_size, min_samples=hdbscan_min_samples
    )
    cluster_labels = clusterer.fit_predict(callrate_data)
    n_clusters = len(set(cluster_labels)) - (
            -1 in cluster_labels
    )  # NOTE: -1 is the label for noisy (un-classifiable) data points
    logger.info(
        "Found {} unique platforms during platform imputation.".format(n_clusters)
    )

    data["qc_platform"] = cluster_labels

    # Note: write pandas dataframe to disk and re-import as HailTable.
    # This a temporary solution until sort the hail's issue with the function 'hl.Table.from_pandas'
    # and different python versions between driver/executors.
    (data
     .drop(axis=1, labels=pc_scores_ann)
     .to_csv('data_tmp_hdbscan.tsv', index=False, sep='\t')
     )
    ht_tmp = (hl.import_table('data_tmp_hdbscan.tsv', impute=True)
              .key_by(*platform_pca_scores_ht.key)
              )

    ht = platform_pca_scores_ht.join(ht_tmp)

    # original/elegant solution (TODO: sort issue with 'from_pandas' function)
    # ht = hl.Table.from_pandas(data, key=[*platform_pca_scores_ht.key])

    # expand array structure and annotate scores (PCs) as individual fields.
    # drop array scores field before to export the results.
    n_pcs = len(ht[pc_scores_ann].take(1)[0])
    ht = (ht
          .annotate(**{f'platform_PC{i + 1}': ht[pc_scores_ann][i] for i in range(n_pcs)})
          .drop(pc_scores_ann)
          )

    ht = ht.annotate(qc_platform="platform_" + hl.str(ht.qc_platform))
    return ht


def main(args):

    # init hail
    hl.init(default_reference=args.default_ref_genome)

    # input MT
    mt = hl.read_matrix_table(args.mt_input_path)

    # filter high-quality genotype
    mt = filter_genotypes_ab(mt)

    # import capture interval table (intersect)
    intervals = hl.read_table(args.ht_intervals)

    # generate an interval x sample MT by computing per intervals callrate
    mt_callrate = compute_callrate_mt(mt=mt,
                                      intervals_ht=intervals)

    # run pca
    _, ht_pca, _ = run_platform_pca(callrate_mt=mt_callrate,
                                    binarization_threshold=args.binarization_threshold)

    # apply unsupervised clustering on PCs to infer samples platform
    ht_platform = assign_platform_from_pcs(platform_pca_scores_ht=ht_pca,
                                           pc_scores_ann='scores',
                                           hdbscan_min_cluster_size=args.hdbscan_min_cluster_size,
                                           hdbscan_min_samples=args.hdbscan_min_cluster_size)

    ht_platform.show()

    # write HT
    ht_platform.write(output=args.ht_output_path,
                      overwrite=args.overwrite)

    # export to file if true
    if args.write_to_file:
        (ht_platform
         .export(f'{args.ht_output_path}.tsv.bgz')
         )

    hl.stop()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--mt_input_path', help='Path to input Hail MatrixTable')
    parser.add_argument('--ht_output_path', help='Path to output HailTable with platform PCs')
    parser.add_argument('--ht_intervals', help='HailTable of intervals to use in the computation')
    parser.add_argument('--write_to_file', help='Write output to BGZ-compressed file',
                        action='store_true')
    parser.add_argument('--overwrite', help='Overwrite results in HT output path if already exists...',
                        action='store_true')
    parser.add_argument('--default_ref_genome', help='Default reference genome to start Hail',
                        type=str, default='GRCh38')
    parser.add_argument('--binarization_threshold',
                        help='If set, the callrate MT is transformed to a 0/1 value MT based on the threshold. '
                             'E.g. with the default threshold of 0.25, all entries with a callrate < 0.25 '
                             'are considered as 0s, others as 1s. Set to None if no threshold desired',
                        type=float, default=0.25)

    # HDBSCAN parameters
    parser.add_argument('--hdbscan_min_cluster_size',
                        help='HDBSCAN `min_cluster_size` parameter. If not specified the smallest of 500 and '
                             '0.1*n_samples will be used.', type=int, default=500)
    parser.add_argument('--hdbscan_min_samples', help='HDBSCAN `min_samples` parameter.',
                        type=int, default=None)

    args = parser.parse_args()

    main(args)
