"""
Compute relatedness estimates between individuals using a variant of the PC-Relate method.

usage: relatedness_inference.py [-h] [--mt_input_path MT_INPUT_PATH]
                                [--default_reference DEFAULT_REFERENCE]
                                [--write_to_file] [--overwrite]
                                [--skip_filter_data]
                                [--maf_threshold MAF_THRESHOLD]
                                [--skip_prune_ld] [--r2 R2] [--n_pcs N_PCS]
                                [--skip_compute_pc_relate]
                                [--min_individual_maf MIN_INDIVIDUAL_MAF]
                                [--min_kinship MIN_KINSHIP]

optional arguments:
  -h, --help            show this help message and exit
  --mt_input_path MT_INPUT_PATH
                        Path to MatrixTable with computed QC metrics
  --default_reference DEFAULT_REFERENCE
                        One of GRCh37 and GRCh38
  --write_to_file       Write results to TSV file
  --overwrite           Overwrite pre-existing data
  --skip_filter_data    Skip filtering the input MT. Load pre-existing
                        filtered MT to bi-allelic, high-callrate and common
                        SNPs...
  --maf_threshold MAF_THRESHOLD
                        Exclude variants with maf lower than maf_threshold
  --skip_prune_ld       Skip pruning variants in LD. It is recommended to
                        prune LD-variants before running PCA
  --r2 R2               Squared correlation threshold for LD pruning
  --n_pcs N_PCS         Number of PCs to be computed
  --skip_compute_pc_relate
                        Skip computing relatedness with the pc_relate method.
                        It will force to extract related samples from pre-
                        existing table with kinship stats.
  --min_individual_maf MIN_INDIVIDUAL_MAF
                        Individual specif MAF cutoff used to run pc_relate
                        method
  --min_kinship MIN_KINSHIP
                        Exclude pairs of samples with kinship lower than
                        min_kinship

"""

import argparse
import logging

import hail as hl

from utils.data_utils import (get_sample_meta_data,
                              import_fam_ht)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("relatedness_inference")
logger.setLevel(logging.INFO)

hdfs_dir = 'hdfs://spark-master:9820'
nfs_dir = 'file:///home/ubuntu/data'

MIN_KINSHIP = 0.125


def make_sample_rank_table(phe_ht: hl.Table) -> hl.Table:
    """
    Make table with rank of sample sorted by retention priority
    (lower rank has higher priority).
    It mainly uses two bits of information:
      - cases are prioritised over controls
      - samples are preferred based on the cohort info as follow: chd > ddd > ukbb
    :param phe_ht: Table with sample meta-data annotations (e.g. phenotype, cohort info...)
    :return: Hail Table
    """

    phe_ht = (phe_ht
              .annotate(case_control_rank=hl.int(phe_ht['phe.is_case']),  # 0: control, 1: cases
                        cohort_rank=hl.case()
                        .when(phe_ht.is_ukbb, 10)
                        .when(phe_ht.is_ddd, 100)
                        .when(phe_ht.is_chd, 1000)
                        .or_missing()
                        )
              .key_by()
              )

    phe_ht = (phe_ht
              .select('ega_id', 'case_control_rank', 'cohort_rank')
              )

    # sort table (descending)
    tb_rank = (phe_ht
               .order_by(hl.desc(phe_ht.case_control_rank),
                         hl.desc(phe_ht.cohort_rank))
               )

    tb_rank = (tb_rank
               .add_index(name='rank')
               .key_by('ega_id')
               )

    tb_rank = tb_rank.annotate(rank=tb_rank.rank + 1)

    return tb_rank


def get_related_samples_to_drop(rank_table: hl.Table, relatedness_ht: hl.Table) -> hl.Table:
    """
    Use the maximal independence function in Hail to intelligently prune clusters of related individuals, removing
    less desirable samples while maximizing the number of unrelated individuals kept in the sample set
    :param Table rank_table: Table with ranking annotations across exomes and genomes, computed via make_rank_file()
    :param Table relatedness_ht: Table with kinship coefficient annotations computed via pc_relate()
    :return: Table containing sample IDs ('s') to be pruned from the combined exome and genome sample set
    :rtype: Table
    """
    # Define maximal independent set, using rank list
    related_pairs = relatedness_ht.filter(
        relatedness_ht.kin > 0.125).select('i', 'j')
    n_related_samples = hl.eval(hl.len(
        related_pairs.aggregate(
            hl.agg.explode(
                lambda x: hl.agg.collect_as_set(x),
                [related_pairs.i, related_pairs.j]
            ),
            _localize=False)
    ))
    logger.info(
        '{} samples with at least 2nd-degree relatedness found in callset'.format(n_related_samples))
    max_rank = rank_table.count()
    related_pairs = related_pairs.annotate(
        id1_rank=hl.struct(
            id=related_pairs.i, rank=rank_table[related_pairs.i].rank),
        id2_rank=hl.struct(
            id=related_pairs.j, rank=rank_table[related_pairs.j].rank)
    ).select('id1_rank', 'id2_rank')

    def tie_breaker(l, r):
        return hl.or_else(l.rank, max_rank + 1) - hl.or_else(r.rank, max_rank + 1)

    related_samples_to_drop_ranked = hl.maximal_independent_set(related_pairs.id1_rank, related_pairs.id2_rank,
                                                                keep=False, tie_breaker=tie_breaker)
    return related_samples_to_drop_ranked


def main(args):

    # Init Hail
    hl.init(default_reference=args.default_reference)

    if not args.skip_compute_pc_relate:

        if not args.skip_filter_data:
            # Read MatrixTable
            mt = hl.read_matrix_table(args.mt_input_path)

            # filter variants (bi-allelic, high-callrate, common SNPs)
            logger.info(f"Filtering to bi-allelic, high-callrate, common SNPs ({args.maf_threshold}) for pc_relate...")

            mt = (mt
                  .filter_rows((hl.len(mt.alleles) == 2) & hl.is_snp(mt.alleles[0], mt.alleles[1]) &
                               (hl.agg.mean(mt.GT.n_alt_alleles()) / 2 > args.maf_threshold) &
                               (hl.agg.fraction(hl.is_defined(mt.GT)) > 0.99) &
                               ~mt.was_split)
                  .repartition(500, shuffle=False)
                  )

            # keep only GT entry field and force to evaluate expression
            (mt
             .select_entries(mt.GT)
             .write(f'{nfs_dir}/hail_data/sample_qc/chd_ukbb.filtered_high_confidence_variants.mt',
                    overwrite=args.overwrite)
             )

        mt = hl.read_matrix_table(f'{nfs_dir}/hail_data/sample_qc/chd_ukbb.filtered_high_confidence_variants.mt')

        if not args.skip_prune_ld:
            # LD pruning
            # Avoid filtering / missingness entries (genotypes) before run LP pruning
            # Zulip Hail support issue -> "BlockMatrix trouble when running pc_relate"
            # mt = mt.unfilter_entries()

            # Prune variants in linkage disequilibrium.
            # Return a table with nearly uncorrelated variants

            logger.info(f'Pruning variants in LD from MT with {mt.count_rows()} variants...')

            pruned_variant_table = hl.ld_prune(mt.GT, r2=args.r2)

            # Keep LD-pruned variants
            pruned_mt = (mt
                         .filter_rows(hl.is_defined(pruned_variant_table[mt.row_key]), keep=True)
                         )
            pruned_mt.write(f'{nfs_dir}/hail_data/sample_qc/chd_ukbb.ld_pruned.mt', overwrite=args.overwrite)

        pruned_mt = hl.read_matrix_table(f'{nfs_dir}/hail_data/sample_qc/chd_ukbb.ld_pruned.mt')
        v, s = pruned_mt.count()
        logger.info(f'{s} samples, {v} variants found in LD-pruned MT')

        pruned_mt = pruned_mt.select_entries(
            GT=hl.unphased_diploid_gt_index_call(pruned_mt.GT.n_alt_alleles()))

        # run pc_relate method...compute all stats
        logger.info('Running PCA for PC-Relate...')
        eig, scores, _ = hl.hwe_normalized_pca(pruned_mt.GT, k=10, compute_loadings=False)
        scores.write(f'{nfs_dir}/hail_data/sample_qc/chd_ukbb.pruned.pca_scores_for_pc_relate.ht',
                     overwrite=args.overwrite)

        logger.info(f'Running PC-Relate...')
        scores = hl.read_table(f'{nfs_dir}/hail_data/sample_qc/chd_ukbb.pruned.pca_scores_for_pc_relate.ht')
        relatedness_ht = hl.pc_relate(call_expr=pruned_mt.GT,
                                      min_individual_maf=args.min_individual_maf,
                                      scores_expr=scores[pruned_mt.col_key].scores,
                                      block_size=4096,
                                      min_kinship=args.min_kinship,
                                      statistics='all')

        logger.info(f'Writing relatedness table...')
        # Write/export table to file
        relatedness_ht.write(
            output=f'{nfs_dir}/hail_data/sample_qc/chd_ukbb.relatedness_kinship.ht',
            overwrite=args.overwrite)

        # Write PCs table to file (if specified)
        # if args.write_to_file:
        #    # Export table to file
        #    relatedness_ht.export(output=f'{args.ht_output_path}.tsv.bgz')

    # retrieve maximal independent set of related samples
    logger.info('Getting optimal set of related samples to prune...')

    relatedness_ht = hl.read_table(f'{nfs_dir}/hail_data/sample_qc/chd_ukbb.relatedness_kinship.ht')
    
    relatedness_ht = (relatedness_ht
                      .flatten()
                      .rename({'i.s':'i', 'j.s':'j'})
                      .repartition(100)
                     )

    # import trios info
    fam = import_fam_ht()
    mat_ids = hl.set(fam.mat_id.collect())
    fat_ids = hl.set(fam.pat_id.collect())

    # rank samples by retention priority (e.g. cases over controls)
    tb_rank = make_sample_rank_table(
        get_sample_meta_data()
    )

    # apply min kinship to consider related pairs
    relatedness_ht = (relatedness_ht.
                      filter(relatedness_ht.kin > MIN_KINSHIP))

    # run maximal_independent_set stratified by groups
    # Note: This method fails when considering all pairs together (e.g. it removes most of the index in trios, we want
    # keep them (index) since they are mostly affected individuals rather than parents).

    # defining pairs group
    # TODO: check groups with updated fam file
    relatedness_ht = (relatedness_ht
                      .annotate(pairs_group=hl.case()
                                .when(relatedness_ht.kin > 0.40, 'twins_or_dups')
                                .when(mat_ids.contains(relatedness_ht.i) | mat_ids.contains(relatedness_ht.j),
                                      'pairs_child_mat')
                                .when(fat_ids.contains(relatedness_ht.i) | fat_ids.contains(relatedness_ht.j),
                                      'pairs_child_fat')
                                .default('pairs_others')
                                )
                      )

    groups = (relatedness_ht
              .aggregate(hl.agg.collect_as_set(relatedness_ht['pairs_group']))
              )
    tbs = []
    for pair_group in groups:
        pair_ht = relatedness_ht.filter(relatedness_ht.pairs_group == pair_group)
        tb = get_related_samples_to_drop(rank_table=tb_rank,
                                         relatedness_ht=pair_ht)
        tbs.append(tb)

    related_samples_to_remove = hl.Table.union(*tbs)

    related_samples_to_remove.describe()

    related_samples_to_remove = related_samples_to_remove.checkpoint(
        f'{nfs_dir}/hail_data/sample_qc/chd_ukbb.related_samples_to_remove.ht',
        overwrite=args.overwrite)
    
    if args.write_to_file:
        (related_samples_to_remove
         .flatten()
         .export(f'{nfs_dir}/hail_data/sample_qc/chd_ukbb.related_samples_to_remove.tsv')
         )


    hl.stop()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--mt_input_path', help='Path to MatrixTable with computed QC metrics', type=str, default=None)
    parser.add_argument('--default_reference', help='One of GRCh37 and GRCh38', type=str, default='GRCh38')
    parser.add_argument('--write_to_file', help='Write results to TSV file', action='store_true')
    parser.add_argument('--overwrite', help='Overwrite pre-existing data', action='store_true')

    # actions for controling filtering step
    parser.add_argument('--skip_filter_data', help='Skip filtering the input MT. Load pre-existing filtered MT to \
                                                    bi-allelic, high-callrate and common SNPs...', action='store_true')
    parser.add_argument('--maf_threshold', help='Exclude variants with maf lower than maf_threshold',
                        type=float, default=0.001)
    # actions for ld_prune
    parser.add_argument('--skip_prune_ld', help='Skip pruning variants in LD. It is recommended to prune LD-variants \
                                                 before running PCA', action='store_true')
    parser.add_argument('--r2', help='Squared correlation threshold for LD pruning', type=float, default=0.1)

    # actions for pca
    parser.add_argument('--n_pcs', help='Number of PCs to be computed', type=int, default=10)

    # actions for pc_relate
    parser.add_argument('--skip_compute_pc_relate', help='Skip computing relatedness with the pc_relate method. \
                         It will force to extract related samples from pre-existing table with kinship stats.',
                        action='store_true')
    parser.add_argument('--min_individual_maf', help='Individual specif MAF cutoff used to run pc_relate method',
                        type=float, default=0.05)
    parser.add_argument('--min_kinship', help='Exclude pairs of samples with kinship lower than min_kinship',
                        type=float, default=0.1)

    args = parser.parse_args()

    main(args)
