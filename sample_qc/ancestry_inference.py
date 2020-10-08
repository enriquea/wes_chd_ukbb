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

from resources.data_utils import (get_mt_data,
                                  get_sample_qc_ht_path,
                                  get_1kg_mt,
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
        mt = get_mt_data(dataset=args.exome_cohort, part='unfiltered')

        # Read MT from 1kgenome and keep only locus defined in interval
        mt_1kg = get_1kg_mt(args.default_reference)

        # Joining dataset (inner join). Keep only 'GT' entry field
        mt_joint = (mt
                    .select_entries('GT')
                    .union_cols(mt_1kg
                                .select_entries('GT'),
                                row_join_type='inner')
                    )

        logger.info("Filtering joint MT to bi-allelic, high-callrate, common SNPs...")
        mt_joint = (mt_joint
                    .filter_rows(bi_allelic_expr(mt_joint) &
                                 hl.is_snp(mt_joint.alleles[0], mt_joint.alleles[1]) &
                                 (hl.agg.mean(mt_joint.GT.n_alt_alleles()) / 2 > 0.001) &
                                 (hl.agg.fraction(hl.is_defined(mt_joint.GT)) > 0.99))
                    .naive_coalesce(1000)
                    )

        logger.info("Checkpoint: writing joint filtered MT before LD pruning...")
        mt_joint = mt_joint.checkpoint(get_mt_checkpoint_path(dataset=args.exome_cohort,
                                                              part='joint_1kg_high_callrate_common_snp_biallelic'),
                                       overwrite=True)

        logger.info(f"Running ld_prune with r2 = {args.ld_prune_r2} on MT with {mt_joint.count_rows()} variants...")
        # remove correlated variants
        pruned_variant_table = hl.ld_prune(mt_joint.GT,
                                           r2=args.ld_prune_r2,
                                           bp_window_size=500000,
                                           memory_per_core=512)
        mt_joint = (mt_joint
                    .filter_rows(hl.is_defined(pruned_variant_table[mt_joint.row_key]))
                    )

        logger.info("Writing filtered joint MT with variants in LD pruned...")
        (mt_joint
         .write(get_qc_mt_path(dataset=args.exome_cohort + '_1kg',
                               part='joint_high_callrate_common_snp_biallelic',
                               split=True,
                               ld_pruned=True),
                overwrite=args.overwrite)
         )

    logger.info("Importing filtered joint MT...")
    mt_joint = hl.read_matrix_table(get_qc_mt_path(dataset=args.exome_cohort + '_1kg',
                                                   part='joint_high_callrate_common_snp_biallelic',
                                                   split=True,
                                                   ld_pruned=True))

    logger.info(f"Running PCA with {mt_joint.count_rows()} variants...")
    # run pca on merged dataset
    eigenvalues, pc_scores, _ = hl.hwe_normalized_pca(mt_joint.GT,
                                                      k=args.n_pcs)

    logger.info(f"Eigenvalues: {eigenvalues}")  # TODO: save eigenvalues?

    # Annotate PC array as independent fields.
    pca_table = (pc_scores
                 .annotate(**{'PC' + str(k + 1): pc_scores.scores[k] for k in range(0, args.n_pcs)})
                 .drop('scores')
                 )

    logger.info(f"Writing HT with PCA results...")
    # write as HT
    output_ht_path = get_sample_qc_ht_path(dataset=args.exome_cohort, part='joint_pca_1kg')
    pca_table.write(output=output_ht_path)

    if args.write_to_file:
        (pca_table
         .export(f'{output_ht_path}.tsv.bgz')
         )

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
