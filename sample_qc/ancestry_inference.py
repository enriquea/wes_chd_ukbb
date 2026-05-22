"""
Compute principal components (PCs) on the set of overlapping variants
between input MT and exomes in 1000 Genomes Phase 3 dataset.

Workflow:

1- Join (inner) MT and 1KGenomes dataset
2- Define/filter high-quality variants (bi-allelic, common SNP (MAF 0.1%), call_rate > 0.99)
3- Filter correlated variants (LD pruning)
4- Run PCA on the filtered variants set
5- Export results

"""

import argparse
import logging

import hail as hl

from utils.data_utils import (get_mt_data,
                              get_sample_qc_ht_path,
                              get_1kg_mt,
                              get_qc_mt_path,
                              get_mt_checkpoint_path)

from utils.expressions import bi_allelic_expr

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def join_with_1kg(exome_cohort: str, default_reference: str) -> hl.MatrixTable:
    """Import cohort MT and 1KG MT, then inner-join on GT entry field."""
    logger.info("Importing data...")

    # import unfiltered MT
    mt = get_mt_data(dataset=exome_cohort, part='raw')

    # Read MT from 1kgenome and keep only locus defined in interval
    mt_1kg = get_1kg_mt(default_reference)

    # Joining dataset (inner join). Keep only 'GT' entry field
    mt_joint = (mt
                .select_entries('GT')
                .union_cols(mt_1kg
                            .select_entries('GT'),
                            row_join_type='inner')
                )

    return mt_joint


def filter_joint_mt(mt_joint: hl.MatrixTable) -> hl.MatrixTable:
    """Filter joint MT to bi-allelic, high-callrate, common SNPs and coalesce."""
    logger.info("Filtering joint MT to bi-allelic, high-callrate, common SNPs...")
    mt_joint = (mt_joint
                .filter_rows(bi_allelic_expr(mt_joint) &
                             hl.is_snp(mt_joint.alleles[0], mt_joint.alleles[1]) &
                             (hl.agg.mean(mt_joint.GT.n_alt_alleles()) / 2 > 0.001) &
                             (hl.agg.fraction(hl.is_defined(mt_joint.GT)) > 0.99))
                .naive_coalesce(1000)
                )
    return mt_joint


def checkpoint_joint_mt(mt_joint: hl.MatrixTable, exome_cohort: str) -> hl.MatrixTable:
    """Checkpoint filtered joint MT to disk before LD pruning."""
    logger.info("Checkpoint: writing joint filtered MT before LD pruning...")
    mt_joint = mt_joint.checkpoint(get_mt_checkpoint_path(dataset=exome_cohort,
                                                          part='joint_1kg_high_callrate_common_snp_biallelic'),
                                   overwrite=True)
    return mt_joint


def ld_prune_and_write(mt_joint: hl.MatrixTable, exome_cohort: str,
                       ld_prune_r2: float, overwrite: bool) -> None:
    """LD-prune variants and write the final filtered joint MT."""
    logger.info(f"Running ld_prune with r2 = {ld_prune_r2} on MT with {mt_joint.count_rows()} variants...")
    # remove correlated variants
    pruned_variant_table = hl.ld_prune(mt_joint.GT,
                                       r2=ld_prune_r2,
                                       bp_window_size=500000,
                                       memory_per_core=512)
    mt_joint = (mt_joint
                .filter_rows(hl.is_defined(pruned_variant_table[mt_joint.row_key]))
                )

    logger.info("Writing filtered joint MT with variants in LD pruned...")
    (mt_joint
     .write(get_qc_mt_path(dataset=exome_cohort + '_1kg',
                           part='joint_high_callrate_common_snp_biallelic',
                           split=True,
                           ld_pruned=True),
            overwrite=overwrite)
     )


def read_filtered_joint_mt(exome_cohort: str) -> hl.MatrixTable:
    """Read the pre-computed filtered joint MT from disk."""
    logger.info("Importing filtered joint MT...")
    mt_joint = hl.read_matrix_table(get_qc_mt_path(dataset=exome_cohort + '_1kg',
                                                   part='joint_high_callrate_common_snp_biallelic',
                                                   split=True,
                                                   ld_pruned=True))
    return mt_joint


def run_pca(mt_joint: hl.MatrixTable, n_pcs: int) -> hl.Table:
    """Run HWE-normalised PCA and return a Table with per-sample PC scores."""
    logger.info(f"Running PCA with {mt_joint.count_rows()} variants...")
    # run pca on merged dataset
    eigenvalues, pc_scores, _ = hl.hwe_normalized_pca(mt_joint.GT,
                                                      k=n_pcs)

    logger.info(f"Eigenvalues: {eigenvalues}")  # TODO: save eigenvalues?

    # Annotate PC array as independent fields.
    pca_table = (pc_scores
                 .annotate(**{'PC' + str(k + 1): pc_scores.scores[k] for k in range(0, n_pcs)})
                 .drop('scores')
                 )
    return pca_table


def write_pca_results(pca_table: hl.Table, exome_cohort: str, write_to_file: bool) -> None:
    """Write PCA HailTable to disk and optionally export as BGZ-compressed TSV."""
    logger.info(f"Writing HT with PCA results...")
    # write as HT
    output_ht_path = get_sample_qc_ht_path(dataset=exome_cohort, part='joint_pca_1kg')
    pca_table.write(output=output_ht_path)

    if write_to_file:
        (pca_table
         .export(f'{output_ht_path}.tsv.bgz')
         )


def main(args):
    # Start Hail
    hl.init(default_reference=args.default_reference)

    if not args.skip_filter_step:
        mt_joint = join_with_1kg(args.exome_cohort, args.default_reference)
        mt_joint = filter_joint_mt(mt_joint)
        mt_joint = checkpoint_joint_mt(mt_joint, args.exome_cohort)
        ld_prune_and_write(mt_joint, args.exome_cohort, args.ld_prune_r2, args.overwrite)

    mt_joint = read_filtered_joint_mt(args.exome_cohort)
    pca_table = run_pca(mt_joint, args.n_pcs)
    write_pca_results(pca_table, args.exome_cohort, args.write_to_file)

    # Stop Hail
    hl.stop()

    print("Done!")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exome_cohort', help="One of <chd_ukbb> or <chd_ddd>",
                        type=str, default=None)
    parser.add_argument('--skip_filter_step',
                        help='Skip filtering the joint MT. will load already existing filtered joint MT with '
                             'high-quality variants',
                        action='store_true')
    parser.add_argument('--overwrite', help='Overwrite pre-existing data', action='store_true')
    parser.add_argument('--n_pcs', help='Number of PCs to compute', type=int, default=20)
    parser.add_argument('--ld_prune_r2', help='Squared correlation threshold for ld_prune method',
                        type=float, default=0.2)
    parser.add_argument('--write_to_file', help='Write output to BGZ-compressed file',
                        action='store_true')
    parser.add_argument('--default_reference', help='One of GRCh37 and GRCh38', type=str, default='GRCh38')

    args = parser.parse_args()

    main(args)
