import argparse
import hail as hl
import os

"""
Compute PCs on the set of SNPs defined in the Purcell interval.
The 1000 Genomes Phase 3 dataset aligned against hg38 is used 
as reference to accesss known sample population.

Usage example:

python ancestry_inference_pca.py \
        --mt_input_path 'file:///home/ubuntu/data/hail_data/mts/chd_ukbb_split_v2_09092020.mt' \
        --mt_input_1kg 'file:///home/ubuntu/data/hail_data/mts/phase3_1kg_snp_biallelic_hg38.mt' \
        --ht_interval 'file:///home/ubuntu/data/resources/intervals/purcell_5k_intervals/purcell5k.ht' \
        --ht_output_path 'hdfs://spark-master:9820/dir/chd_ukb_purcell5k_1kg_merged_pca_maf0.05_14092020.ht' \
        --write_to_file 
"""


def main(args):
    # Start Hail
    hl.init(default_reference=args.default_ref_genome)

    # Read Table with defined intervals
    # intervals = hl.read_table(args.ht_interval)

    # collect intervals as set for filtering (expected small set)
    # interval_set = hl.set(intervals.aggregate(hl.agg.collect_as_set(intervals.locus)))

    # Read Hail MatrixTable
    mt = hl.read_matrix_table(args.mt_input_path)

    # drop entries fields (except GT) to harmonaze entry schemes between dataset
    entry_list = list(mt.entry)
    entry_fields_to_drop = [x for x in entry_list if x != 'GT']
    mt = mt.drop(*entry_fields_to_drop)

    # Filter by intervals and join
    # mt = (mt.
    #     filter_rows(interval_set.contains(mt.locus))
    #     )

    # Read MT from 1kgenome and keep only locus defined in interval
    mt_1kg = hl.read_matrix_table(args.mt_input_1kg)
    # mt_1kg = (mt_1kg.
    #         filter_rows(interval_set.contains(mt_1kg.locus))
    #         )

    # Join dataset (inner join)
    mt_merged = mt.union_cols(mt_1kg)

    print(mt_merged.count())

    # compute variant qc
    mt_merged = hl.variant_qc(mt_merged)

    # keep high quality variants
    mt_merged = (mt_merged
                 .filter_rows(hl.is_snp(mt_merged.alleles[0], mt_merged.alleles[1]) &
                              (mt_merged.variant_qc.AF[1] > 0.05) &
                              (mt_merged.variant_qc.call_rate >= 0.99),
                              keep=True)
                 )

    # remove correlated variants
    pruned_variant_table = hl.ld_prune(mt_merged.GT,
                                       r2=0.2,
                                       bp_window_size=500000)
    mt_merged = (mt_merged
                 .filter_rows(hl.is_defined(pruned_variant_table[mt_merged.row_key]))
                 )

    print(mt_merged.count())

    # run pca on merged dataset
    eigenvalues, pc_scores, _ = hl.hwe_normalized_pca(mt_merged.GT,
                                                      k=args.n_pcs)

    # getting PCs
    pca_table = (pc_scores
                 .annotate(**{'PC' + str(k + 1):
                                  pc_scores.scores[k] for k in range(0, args.n_pcs)})
                 .drop('scores')
                 )

    # write as HT
    output_ht_path = args.ht_output_path
    pca_table.write(output=output_ht_path)

    if args.write_to_file:
        (pca_table
         .export(f'{output_ht_path}.tsv.bgz')
         )

    # Stop Hail
    hl.stop()

    print("Finished!")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--mt_input_path', help='Path to Hail MatrixTable')
    parser.add_argument('--mt_input_1kg', help='Path to 1K Genome MatrixTable')
    parser.add_argument('--ht_interval', help='Path to Hail Table with defined intervals')
    parser.add_argument('--ht_output_path', help='Output path to HailTable with computed PCs')
    parser.add_argument('--n_pcs', help='Number of PCs to compute', type=int, default=10)
    parser.add_argument('--write_to_file', help='Write output to BGZ-compressed file',
                        action='store_true')
    parser.add_argument('--default_ref_genome', help='Default reference genome to start Hail',
                        type=str, default='GRCh38')

    args = parser.parse_args()

    main(args)
