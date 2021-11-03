"""

Run principal components (PCs) analysis on a set of
high confidence polymorphic SNPs.

"""

import argparse
import logging

import hail as hl

from utils.data_utils import (get_sample_qc_ht_path,
                              get_qc_mt_path,
                              get_mt_checkpoint_path)

from utils.expressions import bi_allelic_expr

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def main(args):

    # Start Hail
    hl.init(default_reference=args.default_reference)

    if not args.skip_filter_step:
        logger.info("Importing data...")

        # import unfiltered MT
        mt = hl.read_matrix_table(get_qc_mt_path(dataset=args.exome_cohort,
                                                 part='unphase_adj_genotypes',
                                                 split=True))

        # filter to samples passing QC filters
        logger.info("Filtering MT to samples passing QC filters (hard filters, relatedness, european ancestries)...")
        sample_qc_ht = hl.read_table(get_sample_qc_ht_path(part='final_qc'))
        sample_qc_ht = (sample_qc_ht.filter(sample_qc_ht.pass_filters))
        mt = (mt
              .filter_cols(hl.is_defined(sample_qc_ht[mt.col_key]))
              )

        logger.info("Filtering joint MT to bi-allelic, high-callrate, common SNPs...")
        maf = args.maf_threshold
        mt = (mt
              .filter_rows(bi_allelic_expr(mt) &
                           hl.is_snp(mt.alleles[0], mt.alleles[1]) &
                           (hl.agg.mean(mt.GT.n_alt_alleles()) / 2 > maf) &
                           (hl.agg.fraction(hl.is_defined(mt.GT)) > 0.99))
              .naive_coalesce(500)
              )

        logger.info("Checkpoint: writing joint filtered MT before LD pruning...")
        mt = mt.checkpoint(get_mt_checkpoint_path(dataset=args.exome_cohort,
                                                  part='high_callrate_common_snp_biallelic'),
                           overwrite=args.overwrite)

        logger.info(f"Running ld_prune with r2 = {args.ld_prune_r2} on MT with {mt.count_rows()} variants...")
        # remove correlated variants
        pruned_variant_table = hl.ld_prune(mt.GT,
                                           r2=args.ld_prune_r2,
                                           bp_window_size=500000,
                                           memory_per_core=512)
        mt = (mt
              .filter_rows(hl.is_defined(pruned_variant_table[mt.row_key]))
              )

        logger.info("Writing filtered joint MT with variants in LD pruned...")
        (mt
         .write(get_qc_mt_path(dataset=args.exome_cohort,
                               part='high_callrate_common_snp_biallelic',
                               split=True,
                               ld_pruned=True),
                overwrite=args.overwrite)
         )

    logger.info("Importing filtered joint MT...")
    mt = hl.read_matrix_table(get_qc_mt_path(dataset=args.exome_cohort,
                                             part='high_callrate_common_snp_biallelic',
                                             split=True,
                                             ld_pruned=True))

    logger.info(f"Running PCA with {mt.count_rows()} variants...")
    # run pca on merged dataset
    eigenvalues, pc_scores, _ = hl.hwe_normalized_pca(mt.GT,
                                                      k=args.n_pcs)

    logger.info(f"Eigenvalues: {eigenvalues}")

    # Annotate eigenvalues as global field
    pc_scores = (pc_scores
                 .annotate_globals(**{'eigenvalues': eigenvalues})
                 )

    # Annotate PC array as independent fields.
    pca_table = (pc_scores
                 .annotate(**{'PC' + str(k + 1): pc_scores.scores[k] for k in range(0, args.n_pcs)})
                 .drop('scores')
                 )

    logger.info(f"Writing HT with PCA results...")
    # write as HT
    output_ht_path = args.output_ht
    pca_table = (pca_table
                 .checkpoint(output=output_ht_path,
                             overwrite=args.overwrite)
                 )

    if args.write_to_file:
        (pca_table
         .export(f'{output_ht_path}.tsv.bgz')
         )

    # Stop Hail
    hl.stop()

    print("PCA pipeline finalised...")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--exome_cohort', help="One of <chd_ukbb> or <chd_ddd>",
                        type=str, default=None)

    parser.add_argument('-o', '--output_ht', help="Path to output HailTable",
                        type=str, default=None)

    parser.add_argument('--skip_filter_step',
                        help='Skip filtering the MT. will load already existing filtered MT with \n'
                             'high-quality variants and samples passing qc filters.',
                        action='store_true')

    parser.add_argument('--overwrite', help='Overwrite pre-existing data',
                        action='store_true')

    parser.add_argument('--n_pcs', help='Number of PCs to compute',
                        type=int, default=10)

    parser.add_argument('--maf_threshold',
                        help='Minimal allele frequency (MAF) to filter common variants used to compute PCs',
                        type=float, default=0.01)

    parser.add_argument('--ld_prune_r2', help='Squared correlation threshold for ld_prune method',
                        type=float, default=0.2)

    parser.add_argument('--write_to_file', help='Write output to BGZ-compressed file',
                        action='store_true')

    parser.add_argument('--default_reference', help='One of GRCh37 and GRCh38',
                        type=str, default='GRCh38')

    args = parser.parse_args()

    main(args)
