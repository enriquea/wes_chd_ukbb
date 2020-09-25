"""
Compute relatedness estimates between individuals using a variant of the PC-Relate method.

Usage example:

python sample_relatedness --mt_path 'path/to/mt' \
                          --ht_output_path 'path/to/hts' \
                          --maf_threshold 0.05 \
                          --ld_pruning \
                          --write_to_file
"""

import argparse
import logging

import hail as hl

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("relatedness_inference")
logger.setLevel(logging.INFO)


def main(args):
    # Initializing Hail on cluster mode
    hl.init()

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

    if args.sample_to_keep is not None:
        sample_table = hl.import_table(paths=args.sample_to_keep,
                                       no_header=True)
        sample_set = hl.literal(sample_table.aggregate(hl.agg.collect_as_set(sample_table.f0)))
        mt = mt.filter_cols(sample_set.contains(mt.s), keep=True)

    if args.ld_pruning:
        # LD pruning
        # Avoid filtering / missingness entries (genotypes) before run LP pruning
        # Zulip Hail support issue -> "BlockMatrix trouble when running pc_relate"
        # mt = mt.unfilter_entries()

        # Prune variants in linkage disequilibrium.
        # Return a table with nearly uncorrelated variants

        logger.info(f'Pruning variants in LD from MT with {mt.count_rows()} variants...')

        pruned_variant_table = hl.ld_prune(mt.GT,
                                           r2=args.r2,
                                           bp_window_size=args.bp_window_size,
                                           memory_per_core=1024)

        # Keep LD-pruned variants
        mt = (mt
              .filter_rows(hl.is_defined(pruned_variant_table[mt.locus, mt.alleles]), keep=True)
              )

    # run pc_relate method...compute all stats
    logger.info(f'Running pc_relate method with {args.n_pcs} PCs and {mt.count_rows()} variants...')

    relatedness = hl.pc_relate(call_expr=mt.GT,
                               min_individual_maf=args.min_individual_maf,
                               k=args.n_pcs,
                               min_kinship=args.min_kinship,
                               statistics='kin',
                               block_size=1024)

    # TODO: retrieve maximal independent sample set

    # Write/export table to file
    tb = (relatedness
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
    parser.add_argument('--sample_to_keep', help='Text file (one-column, no header) listing the samples to keep',
                        type=str, default=None)
    parser.add_argument('--n_pcs', help='Number of PCs to be computed', type=int, default=10)
    parser.add_argument('--min_individual_maf', help='Individual specif MAF cutoff used to run pc_relate method',
                        type=float, default=0.01)
    parser.add_argument('--min_kinship', help='Exclude pairs of samples with kinship lower than min_kinship',
                        type=float, default=0.05)
    parser.add_argument('--maf_threshold', help='Exclude variants with maf lower than maf_threshold',
                        type=float, default=0.05)
    parser.add_argument('--ld_pruning', help='Perform LD pruning before PCA (recommended)', action='store_true')
    parser.add_argument('--r2', help='Squared correlation threshold for LD pruning', type=float, default=0.2)
    parser.add_argument('--bp_window_size', help='Window size in bps for ld_prune', type=int, default=500000)
    parser.add_argument('--write_to_file', help='Write results to TSV file', action='store_true')

    args = parser.parse_args()

    main(args)
