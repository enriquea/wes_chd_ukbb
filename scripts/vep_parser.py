import argparse
import hail as hl
import sys

"""

Script to parse a VEP-annotated VCF file. It returns a HailTable with a selected transcript per variant/consequence
(most severe). Multi-allelic variants, if any, are split into biallelic variants. Optionally, multi-allelic variant can
be filtered out.

The VEP-annotated VCF file is expected to be generated with the following vep setting:
[vep run setting here]

Example usage:
python vep_parser.py --vcf_vep_path path/to/vcf /
                     --tb_output_path path/to/output/dir /
                     --force_bgz /
                     --write_to_file /
                     --keep_info_fields VQSLOD

"""


# Define useful parser functions
def get_vep_fields(vcf_path: str,
                   vep_csq_field: str) -> list:
    """
    Parse VCF meta information (header) and extract
    annotated VEP field names.

    :param vep_csq_field: VEP consequence field name
    :param vcf_path: Path to VEP annotated VCF file
    :return: List of annotated fields with VEP
    """

    # Extract VCF meta-info
    meta_info = hl.get_vcf_metadata(vcf_path)

    # parse names in CSQ field
    csq = (meta_info
           .get('info')
           .get(vep_csq_field))
    fields = (csq
              .get('Description')
              .split('Format:')[1]
              .strip()
              .split('|')
              )

    return fields


def annotate_from_array(ht: hl.Table,
                        array_field: str,
                        field_names: list) -> hl.Table:
    """
    Expand an array structure and add new fields.

    :param ht: HailTable
    :param array_field: The array field to be expanded
    :param field_names: The pre-defined fields names (ordered). Number of fields must match with array length
    :return: Annotated HailTable
    """

    # number of fields to be annotated
    n_fields = field_names.__len__()

    # get array field length
    array_len = hl.len(ht[array_field]).take(1)[0]

    if array_len == n_fields:
        tb_expanded = ht.transmute(**{field_names[i]: ht[array_field][i] for i in range(n_fields)})
    else:
        print("Number of fields don't match with array length...")
        sys.exit()
    return tb_expanded


# select one a transcript from array using pre-defined rules/conditions.
def pick_transcript(ht: hl.Table,
                    csq_array: str) -> hl.Table:
    # TODO: This function could be improved by scanning the array (just once) and sorting it as suggested here:
    # TODO: https://hail.zulipchat.com/#narrow/stream/123010-Hail-0.2E2.20support/topic/pick.20transcript.20from.20array
    # TODO: /near/190400193

    """
    Annotate an extra field (tx) with the selected transcript.
    This function will pick one transcript per variant/consequence based on the impact of the variant in the transcript
    (from more severe to less severe).

    :param ht: Hail table with VEP annotations
    :param csq_array: Parsed CSQ field name. Expected to be an array of dict(s). Transcript expected to be as a dict.
    :return: Hail table with an annotated extra field (tx). The transcript selected from the array based on a set of
    pre-defined criteria.
    """

    # Set transcript (tx) field initially to 'NA' and update it sequentially based on a set of pre-defined criteria
    # (order matters)
    ht = (ht
          .annotate(tx=ht[csq_array].find(lambda x: False))
          )

    # getting current keys from dict
    keys = ht[csq_array].take(1)[0][0]

    # select tx if LoF == 'HC'
    if 'LoF' in keys:
        ht = (ht
              .annotate(tx=hl.cond(hl.is_missing(ht.tx),
                                   ht[csq_array].find(lambda x:
                                                      x['LoF'] == 'HC'),
                                   ht.tx)
                        )
              )

    # select transcript based on the consequence impact (high -> moderate -> low)
    if 'IMPACT' in keys:
        # select tx if IMPACT == HIGH
        ht = (ht
              .annotate(tx=hl.cond(hl.is_missing(ht.tx),
                                   ht[csq_array].find(lambda x: x['IMPACT'] == 'HIGH'),
                                   ht.tx)
                        )
              )
        # select tx if IMPACT == MODERATE
        ht = (ht
              .annotate(tx=hl.cond(hl.is_missing(ht.tx),
                                   ht[csq_array].find(lambda x: x['IMPACT'] == 'MODERATE'),
                                   ht.tx)
                        )
              )
        # select tx if IMPACT == LOW
        ht = (ht
              .annotate(tx=hl.cond(hl.is_missing(ht.tx),
                                   ht[csq_array].find(lambda x: x['IMPACT'] == 'LOW'),
                                   ht.tx)
                        )
              )

    # select tx if CANONICAL
    ht = (ht
          .annotate(tx=hl.cond(hl.is_missing(ht.tx),
                               ht[csq_array].find(lambda x: x['CANONICAL'] == 'YES'),
                               ht.tx)
                    )
          )

    # if tx is still missing, set tx as the first annotated transcript
    ht = (ht
          .annotate(tx=hl.cond(hl.is_missing(ht.tx),
                               ht[csq_array][0],
                               ht.tx)
                    )
          )
    return ht


def annotate_from_dict(ht: hl.Table,
                       dict_field: str) -> hl.Table:
    """
    Expand an dict field and add new fields.

    :param ht: HailTable
    :param dict_field: The dict field to be expanded
    :return: Annotated HailTable
    """

    # number of fields to be annotated
    dict_keys = ht[dict_field].keys().take(1)[0]

    # print(dict_keys)

    ht = (ht
          .annotate(**{dict_keys[i]: ht[dict_field].get(dict_keys[i]) for i in range(len(dict_keys))})
          )

    return ht


def cast_str(ht: hl.Table,
             field_names: list,
             output_type: str) -> hl.Table:
    """

    :param ht:
    :param field_names:
    :param output_type:
    :return:
    """
    if output_type == 'float':
        ht = (ht
              .transmute(**{k: hl.cond(~ht[k].matches('\p{Digit}'),
                                       hl.float(0),
                                       hl.float(ht[k])) for k in field_names})
              )
    if output_type == 'int':
        ht = (ht
              .transmute(**{k: hl.cond(~ht[k].matches('\p{Digit}'),
                                       hl.int(0),
                                       hl.int(ht[k])) for k in field_names})
              )

    return ht


def filter_biallelic(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Remove multi-allelic variant. Keep bi-allelic only.

    :param mt: Hail MatrixTable
    :return: Filtered MatrixTable
    """
    n_total = mt.count_rows()
    mt = (mt
          .filter_rows(hl.len(mt.alleles) == 2,
                       keep=True)

          )
    n_filtered = n_total - mt.count_rows()

    pct = (n_filtered / n_total) * 100

    print(f'Filtered {n_filtered} ({round(pct, 2)}%) multi-allelic variants out of {n_total}.')

    return mt


def main(args):
    # Init Hail with hg38 genome build as default
    hl.init(default_reference=args.default_ref_genome)

    # Import VEPed VCF file as MatrixTable and get VCF file meta-data
    vcf_path = args.vcf_vep_path
    mt = hl.import_vcf(path=vcf_path,
                       force_bgz=args.force_bgz)

    # getting annotated VEP fields names from VCF-header
    vep_fields = get_vep_fields(vcf_path=vcf_path,
                                vep_csq_field=args.csq_field)

    if args.exclude_multi_allelic:
        # TODO: This option should skip the split_multi step...
        # Filter out multi-allelic variants. Keep only bi-allelic
        mt = filter_biallelic(mt)

    # split multi-allelic variants
    mt = hl.split_multi_hts(mt)

    # flatten nested structure (e.g. 'info') and get a HailTable with all rows fields
    tb_csq = (mt
              .rows()
              .flatten()
              .key_by('locus', 'alleles')
              )

    # select locus/alleles, info fields and CSQ field.
    if len(args.keep_info_fields) > 0:
        info_fields_to_keep = ['info.' + x for x in args.keep_info_fields]
    else:
        info_fields_to_keep = []

    tb_csq = (tb_csq
              .annotate(csq_array=tb_csq['info.' + args.csq_field])
              .select('a_index', 'was_split', 'csq_array', *info_fields_to_keep)
              )

    # Convert/annotate all transcripts per variants with a structure of type array<dict<str, str>>.
    # The transcript(s) are represented as a dict<k,v>, the keys are the field names extracted from the VCF header, the
    # values are the current annotated values in the CSQ field.
    tb_csq = (tb_csq
              .annotate(csq_array=tb_csq.csq_array.map(lambda x:
                                                       hl.dict(hl.zip(vep_fields, x.split('[|]')))
                                                       ))
              )

    # Keep transcript(s) matching with the allele index.
    # It requires having the flag "ALLELE_NUM" annotated by VEP
    # Apply only were the alleles were split.
    # TODO: Handle exception when the flag "ALLELE_NUM" is not present
    tb_csq = (tb_csq
              .annotate(csq_array=hl.cond(tb_csq.was_split,
                                          tb_csq.csq_array.filter(lambda x:
                                                                  (hl.int(x["ALLELE_NUM"]) == tb_csq.a_index)
                                                                  ),
                                          tb_csq.csq_array
                                          )
                        )
              )

    # select and annotate one transcript per variant based on pre-defined rules
    tb_csq = pick_transcript(ht=tb_csq,
                             csq_array='csq_array')

    # Expand selected transcript (dict) annotations adding independent fields.
    tb_csq = annotate_from_dict(ht=tb_csq, dict_field='tx')

    # Parse the "Consequence" field. Keep only the more severe consequence.
    # Avoid the notation "consequence_1&consequence_2"
    tb_csq = (tb_csq
              .transmute(Consequence=tb_csq.Consequence.split('&')[0])
              )

    # print fields overview
    tb_csq.describe()

    # write table as HailTable to disk
    (tb_csq
     .drop('csq_array', 'tx')
     .write(output=args.tb_output_path)
     )

    if args.write_to_file:
        # write table to disk as a BGZ-compressed TSV file
        (tb_csq
         .drop('csq_array', 'tx')
         .export(args.tb_output_path + '.tsv.bgz')
         )

    # Stop Hail
    hl.stop()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--vcf_vep_path', help='VEP-annotated VCF file path',
                        type=str, default=None)
    parser.add_argument('--tb_output_path', help='Output HailTable path with VEP annotations parsed',
                        type=str, default=None)
    parser.add_argument('--write_to_file', help='Write also the HailTable as a BGZ-compressed TSV file',
                        action='store_true')
    parser.add_argument('--force_bgz', help='Forces BGZ decoding for VCF files with .gz extension.',
                        action='store_true')
    parser.add_argument('--default_ref_genome', help='Default reference genome to start Hail',
                        type=str, default='GRCh38')
    parser.add_argument('--csq_field', help='Consequence field name in the VCF file',
                        type=str, default='CSQ')
    parser.add_argument('--keep_info_fields', help='List of non-VEP fields in the VCF files to keep in the output',
                        nargs="*", type=str, default=[])
    parser.add_argument('--exclude_multi_allelic', help='Remove multi allelic variants if any',
                        action='store_true')

    args = parser.parse_args()

    main(args)
