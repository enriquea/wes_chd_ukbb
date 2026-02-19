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


def filter_and_prune_mt(exome_cohort, maf_threshold, ld_prune_r2, overwrite):
    """Filter MT to QC-passing samples and common bi-allelic SNPs, then LD-prune and write."""
    logger.info("Importing data...")
    mt = hl.read_matrix_table(get_qc_mt_path(dataset=exome_cohort,
                                             part='unphase_adj_genotypes',
                                             split=True))

    logger.info("Filtering MT to samples passing QC filters (hard filters, relatedness, european ancestries)...")
    sample_qc_ht = hl.read_table(get_sample_qc_ht_path(part='final_qc'))
    sample_qc_ht = (sample_qc_ht.filter(sample_qc_ht.pass_filters))
    mt = (mt
          .filter_cols(hl.is_defined(sample_qc_ht[mt.col_key]))
          )

    logger.info("Filtering joint MT to bi-allelic, high-callrate, common SNPs...")
    mt = (mt
          .filter_rows(bi_allelic_expr(mt) &
                       hl.is_snp(mt.alleles[0], mt.alleles[1]) &
                       (hl.agg.mean(mt.GT.n_alt_alleles()) / 2 > maf_threshold) &
                       (hl.agg.fraction(hl.is_defined(mt.GT)) > 0.99))
          .naive_coalesce(500)
          )

    logger.info("Checkpoint: writing filtered MT before LD pruning...")
    mt = mt.checkpoint(get_mt_checkpoint_path(dataset=exome_cohort,
                                              part='high_callrate_common_snp_biallelic'),
                       overwrite=overwrite)

    logger.info(f"Running ld_prune with r2 = {ld_prune_r2} on MT with {mt.count_rows()} variants...")
    pruned_variant_table = hl.ld_prune(mt.GT,
                                       r2=ld_prune_r2,
                                       bp_window_size=500000,
                                       memory_per_core=512)
    mt = (mt
          .filter_rows(hl.is_defined(pruned_variant_table[mt.row_key]))
          )

    logger.info("Writing filtered MT with ld-pruned variants...")
    (mt
     .write(get_qc_mt_path(dataset=exome_cohort,
                           part='high_callrate_common_snp_biallelic',
                           split=True,
                           ld_pruned=True),
            overwrite=overwrite)
     )


def run_pca(exome_cohort, n_pcs):
    """Load LD-pruned MT, run HWE-normalised PCA, and return annotated scores table."""
    logger.info("Importing filtered ld-pruned MT...")
    mt = hl.read_matrix_table(get_qc_mt_path(dataset=exome_cohort,
                                             part='high_callrate_common_snp_biallelic',
                                             split=True,
                                             ld_pruned=True))

    logger.info(f"Running PCA on {mt.count_rows()} variants...")
    eigenvalues, pc_scores, _ = hl.hwe_normalized_pca(mt.GT, k=n_pcs)

    logger.info(f"Eigenvalues: {eigenvalues}")

    pc_scores = (pc_scores
                 .annotate_globals(**{'eigenvalues': eigenvalues})
                 )

    pca_table = (pc_scores
                 .annotate(**{'PC' + str(k + 1): pc_scores.scores[k] for k in range(0, n_pcs)})
                 .drop('scores')
                 )

    return pca_table


def export_results(pca_table, output_ht_path, overwrite, write_to_file):
    """Checkpoint PCA table and optionally export as BGZ-compressed TSV."""
    logger.info(f"Writing HT with PCA results...")
    pca_table = (pca_table
                 .checkpoint(output=output_ht_path,
                             overwrite=overwrite)
                 )

    if write_to_file:
        (pca_table
         .export(f'{output_ht_path}.tsv.bgz')
         )

    return pca_table


def main(args):

    # Start Hail
    hl.init(default_reference=args.default_reference)

    if not args.skip_filter_step:
        filter_and_prune_mt(args.exome_cohort, args.maf_threshold, args.ld_prune_r2, args.overwrite)

    pca_table = run_pca(args.exome_cohort, args.n_pcs)

    export_results(pca_table, args.output_ht, args.overwrite, args.write_to_file)

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
