# eam
# 2021-03-25

"""
Annotate variant Hail Table with prediction scores from dbNSFP database.

usage: annotate_dbnsfp_scores.py [-h] [--ht_dbnsfp HT_DBNSFP]
                                 [--ht_variant HT_VARIANT]
                                 [--ht_output HT_OUTPUT]
                                 [--transcript_field TRANSCRIPT_FIELD]
                                 [--write_to_file] [--overwrite]
                                 [--add_clinvar] [--add_gnomad] [--add_exac]
                                 [--add_1000Gp3] [--add_ESP6500]

optional arguments:
  -h, --help            show this help message and exit
  --ht_dbnsfp HT_DBNSFP
                        Path to dbNSFP Hail Table. Generated with
                        <parse_dbnsfp.py>
  --ht_variant HT_VARIANT
                        Path to HailTable with variants to be annotated with
                        dbNSFP. Basically keyed by <locus> and <alleles>
  --ht_output HT_OUTPUT
                        Path to output HailTable with annotated scores from
                        dbNSFP
  --transcript_field TRANSCRIPT_FIELD
                        Transcript field in the input <ht_variants>. By
                        default, it is expected to be annotated with VEP. It
                        will be used to annotate transcript-specific scores
                        from dbNSFP
  --write_to_file       Write output to BGZ-compressed file
  --overwrite           Overwrite results in HT output path if already
                        exists...
  --add_clinvar         Add Clinvar annotations from dbNSFP.
  --add_gnomad          Add gnomAD annotations from dbNSFP.
  --add_exac            Add ExAC annotations from dbNSFP.
  --add_1000Gp3         Add 1000Genomes phase 3 annotations from dbNSFP.
  --add_ESP6500         Add ESP6500 annotations from dbNSFP.

"""

import argparse

import hail as hl

nfs_dir = 'file:///home/ubuntu/data'


def main(args):
    hl.init(default_reference='GRCh38')

    transcript_field = args.transcript_field

    # import variant HT
    ht_variants_path = args.ht_variant
    ht_variants = hl.read_table(
        ht_variants_path
    ).select(transcript_field)

    # import dbNSFP HT
    ht_dbnsfp_path = args.ht_dbnsfp
    ht_dbnsfp = hl.read_table(
        ht_dbnsfp_path
    )

    # annotate scores from dbNSFP
    # prediction scores fields to annotate
    score_fields = [f for f in ht_dbnsfp.row if f.endswith('_score') or f == 'CADD_phred']

    ht_variants = (ht_variants
                   .annotate(**ht_dbnsfp.select(*score_fields)[ht_variants.key])
                   )

    #  Match score with specific transcript
    ht_variants = (ht_variants
                   .annotate(**{f: ht_variants[f].get(ht_variants[transcript_field])
                               for f in score_fields}
                            )
                   )

    # Annotate extra info from dbNSFP
    # Note: Expected extra fields (as struct) from dbNSFP: ['gnomAD', 'ExAC', '1000Gp3', 'ESP6500', 'clinvar']

    ann_expr_dict = {}

    if args.add_clinvar:
        ann_expr_dict.update(
            {'clinvar': ht_dbnsfp[ht_variants.key]['clinvar']}
        )

    if args.add_gnomad:
        ann_expr_dict.update(
            {'gnomAD': ht_dbnsfp[ht_variants.key]['gnomAD']}
        )

    if args.add_exac:
        ann_expr_dict.update(
            {'ExAC': ht_dbnsfp[ht_variants.key]['ExAC']}
        )

    if args.add_1000Gp3:
        ann_expr_dict.update(
            {'1000Gp3': ht_dbnsfp[ht_variants.key]['1000Gp3']}
        )

    if args.add_ESP6500:
        ann_expr_dict.update(
            {'ESP6500': ht_dbnsfp[ht_variants.key]['ESP6500']}
        )

    if len(ann_expr_dict) > 0:
        ht_variants = (ht_variants
                       .annotate(**ann_expr_dict)
                       )

    # write annotated table
    # write HT
    ht_variants = ht_variants.checkpoint(
        output=args.ht_output,
        overwrite=args.overwrite)

    # export to file if true
    if args.write_to_file:
        (ht_variants
         .export(f'{args.ht_output}.tsv.bgz')
         )

    hl.stop()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--ht_dbnsfp',
                        help='Path to dbNSFP Hail Table. Generated with <parse_dbnsfp.py>')

    parser.add_argument('--ht_variant', help='Path to HailTable with variants to be annotated with dbNSFP. \
                                              Basically keyed by <locus> and <alleles>')

    parser.add_argument('--ht_output', help='Path to output HailTable with annotated scores from dbNSFP')

    parser.add_argument('--transcript_field',
                        help='Transcript field in the input <ht_variants>. By default, it is expected to be annotated \
                              with VEP. It will be used to annotate transcript-specific scores from dbNSFP',
                        type=str, default='Feature')

    parser.add_argument('--write_to_file', help='Write output to BGZ-compressed file',
                        action='store_true')

    parser.add_argument('--overwrite', help='Overwrite results in HT output path if already exists...',
                        action='store_true')

    parser.add_argument('--add_clinvar',
                        help='Add Clinvar annotations from dbNSFP.',
                        action='store_true')

    parser.add_argument('--add_gnomad',
                        help='Add gnomAD annotations from dbNSFP.',
                        action='store_true')

    parser.add_argument('--add_exac',
                        help='Add ExAC annotations from dbNSFP.',
                        action='store_true')

    parser.add_argument('--add_1000Gp3',
                        help='Add 1000Genomes phase 3 annotations from dbNSFP.',
                        action='store_true')

    parser.add_argument('--add_ESP6500',
                        help='Add ESP6500 annotations from dbNSFP.',
                        action='store_true')

    args = parser.parse_args()

    main(args)
