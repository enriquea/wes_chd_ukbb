"""
Compute relatedness estimates between individuals using a variant of the PC-Relate method.

usage: relatedness_inference.py [-h] [--mt_input_path MT_INPUT_PATH]
                                [--ht_output_path HT_OUTPUT_PATH]
                                [--default_reference DEFAULT_REFERENCE]
                                [--write_to_file] [--overwrite]
                                [--skip_filter_data]
                                [--maf_threshold MAF_THRESHOLD]
                                [--skip_prune_ld] [--r2 R2] [--n_pcs N_PCS]
                                [--min_individual_maf MIN_INDIVIDUAL_MAF]
                                [--min_kinship MIN_KINSHIP]

optional arguments:
  -h, --help            show this help message and exit
  --mt_input_path MT_INPUT_PATH
                        Path to MatrixTable with computed QC metrics
  --ht_output_path HT_OUTPUT_PATH
                        Output HailTable path with Kinship stats
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

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("relatedness_inference")
logger.setLevel(logging.INFO)

hdfs_dir = 'hdfs://spark-master:9820'


def main(args):
    # Init Hail
    hl.init(default_reference=args.default_reference)

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
         .write(f'{hdfs_dir}/tmp/chd_ukbb.filtered_high_confidence_variants.mt',
                overwrite=args.overwrite)
         )

    mt = hl.read_matrix_table(f'{hdfs_dir}/tmp/chd_ukbb.filtered_high_confidence_variants.mt')

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
        pruned_mt.write(f'{hdfs_dir}/tmp/chd_ukbb.ld_pruned.mt', overwrite=args.overwrite)

    pruned_mt = hl.read_matrix_table(f'{hdfs_dir}/tmp/chd_ukbb.ld_pruned.mt')
    v, s = pruned_mt.count()
    logger.info(f'{s} samples, {v} variants found in LD-pruned MT')

    pruned_mt = pruned_mt.select_entries(
        GT=hl.unphased_diploid_gt_index_call(pruned_mt.GT.n_alt_alleles()))

    # run pc_relate method...compute all stats
    logger.info('Running PCA for PC-Relate...')
    eig, scores, _ = hl.hwe_normalized_pca(pruned_mt.GT, k=10, compute_loadings=False)
    scores.write(f'{hdfs_dir}/tmp/chd_ukbb.pruned.pca_scores_for_pc_relate.ht',
                 overwrite=args.overwrite)

    logger.info(f'Running PC-Relate...')
    scores = hl.read_table(f'{hdfs_dir}/tmp/chd_ukbb.pruned.pca_scores_for_pc_relate.ht')
    relatedness_ht = hl.pc_relate(call_expr=pruned_mt.GT,
                                  min_individual_maf=args.min_individual_maf,
                                  scores_expr=scores[pruned_mt.col_key].scores,
                                  block_size=4096,
                                  min_kinship=args.min_kinship,
                                  statistics='all')

    # TODO: retrieve maximal independent sample set

    logger.info(f'Writing relatedness table...')
    # Write/export table to file
    relatedness_ht = relatedness_ht.checkpoint(output=args.ht_output_path,
                                               overwrite=args.overwrite)

    # Write PCs table to file (if specified)
    if args.write_to_file:
        # Export table to file
        relatedness_ht.export(output=f'{args.ht_output_path}.tsv.bgz')

    hl.stop()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--mt_input_path', help='Path to MatrixTable with computed QC metrics', type=str, default=None)
    parser.add_argument('--ht_output_path', help='Output HailTable path with Kinship stats', type=str, default=None)
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
    parser.add_argument('--min_individual_maf', help='Individual specif MAF cutoff used to run pc_relate method',
                        type=float, default=0.05)
    parser.add_argument('--min_kinship', help='Exclude pairs of samples with kinship lower than min_kinship',
                        type=float, default=0.1)

    args = parser.parse_args()

    main(args)