"""

Annotate chromosomal sex from inbreeding coefficient (F)

"""

import hail as hl
import argparse

from utils.generic import unphase_mt


def annotate_sex(mt: hl.MatrixTable,
                 male_threshold: float = 0.6,
                 female_threshold: float = 0.4) -> hl.MatrixTable:
    """
    Imputes sex, exports data, and annotates mt with this data
    NOTE:
    :param female_threshold:
    :param male_threshold:
    :param MatrixTable mt: MT containing samples to be ascertained for sex
able
    """
    mt = mt.filter_rows(mt.locus.contig == 'chrX',
                        keep=True)
    mt = unphase_mt(mt)

    sex_ht = hl.impute_sex(mt.GT,
                           aaf_threshold=0.05,
                           female_threshold=female_threshold,
                           male_threshold=male_threshold,
                           include_par=False)

    sex_colnames = ['f_stat', 'is_female']
    sex_ht = sex_ht.select(*sex_colnames)
    mt = mt.annotate_cols(**sex_ht[mt.col_key])
    return mt


def main(args):
    # Start Hail
    hl.init(default_reference=args.default_ref_genome)

    # Read Hail MatrixTable
    mt = hl.read_matrix_table(args.mt_input_path)

    # impute sex
    ht_sex = (annotate_sex(mt)
              .cols()
              .persist()
              )

    # write as HT
    output_ht_path = args.ht_output_path
    ht_sex.write(output=output_ht_path)

    # export to file if true
    if args.write_to_file:
        (ht_sex
         .export(f'{output_ht_path}.tsv.bgz')
         )

    hl.stop()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--mt_input_path', help='Path to Hail MatrixTable', type=str, default=None)
    parser.add_argument('--ht_output_path', help='Output path to HailTable with sex annotation')
    parser.add_argument('--write_to_file', help='Optionality, write table output to BGZ-compressed file',
                        action='store_true')
    parser.add_argument('--default_ref_genome', help='Default reference genome to start Hail',
                        type=str, default='GRCh38')

    args = parser.parse_args()

    main(args)
