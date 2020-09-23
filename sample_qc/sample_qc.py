import argparse
import hail as hl
import os

"""
Given a Hail MatrixTable, compute sample QC metrics.
Export result as separated HailTables. Optionally, export results as 
GBZ-compresed TSV files.

Usage example:

nohup python sample_qc.py \
             --mt_input_path /path/to/matrixtable \
             --ht_output_path /path/to/output/hailtable \
             --write_to_file \
             > nohup.log 2>&1 &
"""


def main(args):
    # Start Hail
    hl.init(default_reference=args.default_ref_genome)

    # Read Hail MatrixTable
    mt = hl.read_matrix_table(args.mt_input_path)

    # compute sample qc (on autosomes)
    mt = hl.sample_qc(mt
                      .filter_rows(mt.locus.in_autosome(), keep=True)
                      )

    # write sample qc hailtable
    tb_sample_qc = (mt
                    .select_cols('sample_qc')
                    .cols()
                    .flatten()
                    .key_by('s')
                    )

    output_ht_path = f'{args.ht_output_path}_sample_qc_autosomes.ht'
    tb_sample_qc.write(output=output_ht_path)

    if args.write_to_file:
        (hl.read_table(output_ht_path)
         .export(f'{output_ht_path}.tsv.bgz')
         )

    # Stop Hail
    hl.stop()

    print("Finished!")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--mt_input_path', help='Path to Hail MatrixTable')
    parser.add_argument('--ht_output_path', help='Path to sample/variant qced HailTable')
    parser.add_argument('--write_to_file', help='Write output to BGZ-compressed file',
                        action='store_true')
    parser.add_argument('--default_ref_genome', help='Default reference genome to start Hail',
                        type=str, default='GRCh38')

    args = parser.parse_args()

    main(args)
