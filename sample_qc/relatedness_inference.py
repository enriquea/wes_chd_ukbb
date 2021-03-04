"""
Compute relatedness estimates between individuals using a variant of the PC-Relate method.

Usage example:

python relatedness_inference.py --mt_path 'path/to/mt' \
                                --ht_output_path 'path/to/hts' \
                                --prune_ld \
                                --write_to_file
"""

import argparse
import logging

import hail as hl

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("relatedness_inference")
logger.setLevel(logging.INFO)

hdfs_dir = 'hdfs://spark-master:9820'


def main(args):

    # Init Hail
    hl.init(default_reference=args.default_reference)

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
    mt = (mt
          .select_entries(mt.GT)
          .checkpoint(f'{hdfs_dir}/tmp/chd_ukbb.filtered_high_confidence_variants.mt')
          )

    if args.prune_ld:
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
        pruned_mt.write(f'{hdfs_dir}/tmp/chd_ukbb.ld_pruned.mt', overwrite=args.overwrite)

    pruned_mt = hl.read_matrix_table(f'{hdfs_dir}/tmp/chd_ukbb.ld_pruned.mt')
    v, s = pruned_mt.count()
    logger.info(f'{s} samples, {v} variants found in LD-pruned MT')

    # run pc_relate method...compute all stats
    logger.info('Running PCA for PC-Relate...')
    eig, scores, _ = hl.hwe_normalized_pca(pruned_mt.GT, k=10, compute_loadings=False)
    scores.write(f'{hdfs_dir}/tmp/chd_ukbb.pruned.pca_scores_for_pc_relate.ht',
                 overwrite=args.overwrite)

    logger.info(f'Running PC-_Relate...')
    scores = hl.read_table(f'{hdfs_dir}/tmp/chd_ukbb.pruned.pca_scores.ht')
    relatedness_ht = hl.pc_relate(call_expr=pruned_mt.GT,
                                  min_individual_maf=args.min_individual_maf,
                                  scores_expr=scores[pruned_mt.col_key].scores,
                                  block_size=4096,
                                  min_kinship=args.min_kinship,
                                  statistics='kin2')

    # TODO: retrieve maximal independent sample set

    logger.info(f'Writing relatedness table...')
    # Write/export table to file
    tb = (relatedness_ht
          .flatten()
          .key_by('i.s', 'j.s')
          )

    tb.write(output=args.ht_output_path)

    # Write PCs table to file (if specified)
    if args.write_to_file:
        # Export table to file
        tb.export(output=f'{args.ht_output_path}.tsv.bgz')

    hl.stop()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--mt_input_path', help='Path to MatrixTable with computed QC metrics', type=str, default=None)
    parser.add_argument('--ht_output_path', help='Output HailTable path with Kinship stats', type=str, default=None)
    parser.add_argument('--default_reference', help='One of GRCh37 and GRCh38', type=str, default='GRCh38')
    parser.add_argument('--maf_threshold', help='Exclude variants with maf lower than maf_threshold',
                        type=float, default=0.001)
    parser.add_argument('--write_to_file', help='Write results to TSV file', action='store_true')

    # parameters for ld_prune
    parser.add_argument('--prune_ld', help='Perform LD pruning before PCA (recommended)', action='store_true')
    parser.add_argument('--r2', help='Squared correlation threshold for LD pruning', type=float, default=0.1)

    # parameters for pca
    parser.add_argument('--n_pcs', help='Number of PCs to be computed', type=int, default=10)

    # parameters for pc_relate
    parser.add_argument('--min_individual_maf', help='Individual specif MAF cutoff used to run pc_relate method',
                        type=float, default=0.05)
    parser.add_argument('--min_kinship', help='Exclude pairs of samples with kinship lower than min_kinship',
                        type=float, default=0.05)

    args = parser.parse_args()

    main(args)
