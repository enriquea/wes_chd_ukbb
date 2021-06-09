"""

Script to parse a VEP-annotated VCF file. It returns a HailTable with a selected transcript per variant/consequence
(most severe). Multi-allelic variants, if any, are split into biallelic variants. Optionally, multi-allelic variant can
be filtered out.

The VEP-annotated VCF file is expected to be generated with the following vep setting:
[vep run setting here]

usage: vep_parser.py [-h] [--vcf_vep_path VCF_VEP_PATH]
                     [--tb_output_path TB_OUTPUT_PATH] [--csq_field CSQ_FIELD]
                     [--force_bgz] [--write_to_file] [--overwrite]
                     [--default_ref_genome DEFAULT_REF_GENOME]
                     [--split_multi_allelic]

optional arguments:
  -h, --help            show this help message and exit
  --vcf_vep_path VCF_VEP_PATH
                        Path to VEP-annotated VCF file
  --tb_output_path TB_OUTPUT_PATH
                        Path to output HailTable with VEP annotations parsed
  --csq_field CSQ_FIELD
                        Consequence field name in the VCF file
  --force_bgz           If True, load .vcf.gz files as blocked gzip files,
                        assuming that they were actually compressed using the
                        BGZ codec.
  --write_to_file       Optionally, export the HailTable as a BGZ-compressed
                        TSV file
  --overwrite           Overwrite pre-existing data
  --default_ref_genome DEFAULT_REF_GENOME
                        Default reference genome to start Hail
  --split_multi_allelic
                        Split multi allelic variants if any



"""

#  Note: Hail will import the info field nested into a structure.
#        This fields need to be split after apply split_multi_allelic.
#
#  Example info structure imported:
#
# 'info': struct
# {
#     AC: array < int32 >,
#     AF: array < float64 >,
#     AN: int32,
#     AS_BaseQRankSum: array < float64 >,
#     AS_FS: array < float64 >,
#     AS_InbreedingCoeff: array < float64 >,
#     AS_MQ: array < float64 >,
#     AS_MQRankSum: array < float64 >,
#     AS_QD: array < float64 >,
#     AS_RAW_BaseQRankSum: str,
#     AS_RAW_MQ: str,
#     AS_RAW_MQRankSum: str,
#     AS_RAW_ReadPosRankSum: str,
#     AS_ReadPosRankSum: array < float64 >,
#     AS_SB_TABLE: str,
#     AS_SOR: array < float64 >,
#     BaseQRankSum: float64,
#     CSQ: array < str >,
#     DB: bool,
#     DP: int32,
#     DS: bool,
#     END: int32,
#     ExcessHet: float64,
#     FS: float64,
#     InbreedingCoeff: float64,
#     MLEAC: array < int32 >,
#     MLEAF: array < float64 >,
#     MQ: float64,
#     MQRankSum: float64,
#     NEGATIVE_TRAIN_SITE: bool,
#     POSITIVE_TRAIN_SITE: bool,
#     QD: float64,
#     RAW_MQandDP: array < int32 >,
#     ReadPosRankSum: float64,
#     SOR: float64,
#     VQSLOD: float64,
#     culprit: str
# }

import argparse
import hail as hl
import sys

from utils.data_utils import (get_variant_qc_ht_path,
                              get_vep_vqsr_vcf_path)

# info fields which need to be split after apply multi-allelic split.
INFO_FIELDS = ['AC',
               'AF',
               'AS_BaseQRankSum',
               'AS_FS',
               'AS_InbreedingCoeff',
               'AS_MQ',
               'AS_MQRankSum',
               'AS_QD',
               'AS_ReadPosRankSum',
               'AS_SOR',
               'MLEAC',
               'MLEAF',
               'RAW_MQandDP']


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
                       dict_field: str,
                       output_filed: str) -> hl.Table:
    """
    Expand an dict field and add new fields.

    :param ht: HailTable
    :param dict_field: The dict field to be expanded
    :param output_filed: The output filed name (annotated as structure)
    :return: Annotated HailTable
    """

    # retrieve dict keys to be annotated as fields
    dict_keys = ht[dict_field].keys().take(1)[0]

    # structure annotation expression
    struct_expr = hl.struct(**{dict_keys[i]: ht[dict_field].get(dict_keys[i]) for i in range(len(dict_keys))})

    ht = (ht
          .annotate(_tmp_field_=struct_expr)
          )
    ht = ht.rename({'_tmp_field_': output_filed})

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
    # Init Hail
    hl.init(default_reference=args.default_ref_genome)

    # Import VEPed VCF file as MatrixTable and get VCF file meta-data
    # vcf_path = args.vcf_vep_path
    mt = hl.import_vcf(path=get_vep_vqsr_vcf_path(),
                       force_bgz=args.force_bgz)

    # getting annotated VEP fields names from VCF-header
    vep_fields = get_vep_fields(vcf_path=get_vep_vqsr_vcf_path(),
                                vep_csq_field=args.csq_field)

    if args.split_multi_allelic:
        # split multi-allelic variants
        mt = hl.split_multi_hts(mt)

        # split/annotate fields in the info field (use allele index )
        mt = mt.annotate_rows(info=mt.info.annotate(**{field: mt.info[field][mt.a_index - 1]
                                                       for field in INFO_FIELDS}))

    # parse/annotate the CSQ field in a different structure
    tb_csq = mt.rows()
    tb_csq = (tb_csq
              .annotate(csq_raw=tb_csq.info[args.csq_field])
              )

    # Convert/annotate all transcripts per variants with a structure of type array<dict<str, str>>.
    # The transcript(s) are represented as a dict<k,v>, where keys are the field names extracted from the VCF header and
    # the values are the current annotated values in the CSQ field.
    tb_csq = (tb_csq
              .annotate(csq_raw=tb_csq.csq_raw.map(lambda x:
                                                   hl.dict(hl.zip(vep_fields, x.split('[|]')))
                                                   ))
              )

    # Keep transcript(s) matching with the allele index (only used if variant were split with split_multi_hts)
    # It requires having the flag "ALLELE_NUM" annotated by VEP
    # Apply only were the alleles were split.
    # TODO: Handle exception when the flag "ALLELE_NUM" is not present
    if all([x in list(tb_csq._fields.keys()) for x in ['was_split', 'a_index']]):
        tb_csq = (tb_csq
                  .annotate(csq_raw=hl.cond(tb_csq.was_split,
                                            tb_csq.csq_raw.filter(lambda x:
                                                                  (hl.int(x["ALLELE_NUM"]) == tb_csq.a_index)
                                                                  ),
                                            tb_csq.csq_raw
                                            )
                            )
                  )

    # select and annotate one transcript per variant based on pre-defined rules
    tb_csq = pick_transcript(ht=tb_csq,
                             csq_array='csq_raw',
                             )

    # Expand selected transcript (dict) annotations adding independent fields.
    tb_csq = annotate_from_dict(ht=tb_csq,
                                dict_field='tx',
                                output_filed='vep')

    # Parse the "Consequence" field. Keep only the more severe consequence.
    # Avoid the notation "consequence_1&consequence_2"
    tb_csq = (tb_csq
              .annotate(vep=tb_csq.vep.annotate(Consequence=tb_csq.vep.Consequence.split('&')[0]))
              )

    # drop redundant/temp fields
    tb_csq = (tb_csq
              .drop('csq_raw', 'tx')
              .repartition(500)
              )

    # print fields overview
    tb_csq.describe()

    # write table as HailTable to disk
    # (tb_csq
    # .write(output=args.tb_output_path,
    #        overwrite=args.overwrite)
    # )

    output_path = get_variant_qc_ht_path(part='vep_vqsr', split=args.split_multi_allelic)
    tb_csq = (tb_csq
              .checkpoint(output=output_path,
                          overwrite=args.overwrite)
              )

    if args.write_to_file:
        # write table to disk as a BGZ-compressed TSV file
        (tb_csq
         .export(f'{output_path}.tsv.bgz')
         )

    # Stop Hail
    hl.stop()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    # parser.add_argument('--vcf_vep_path', help='Path to VEP-annotated VCF file',
    #                    type=str, default=None)
    # parser.add_argument('--tb_output_path', help='Path to output HailTable with VEP annotations parsed',
    #                    type=str, default=None)
    parser.add_argument('--csq_field', help='Consequence field name in the VCF file',
                        type=str, default='CSQ')
    parser.add_argument('--force_bgz',
                        help='If True, load .vcf.gz files as blocked gzip files, assuming that they were actually'
                             ' compressed using the BGZ codec.', action='store_true')
    parser.add_argument('--write_to_file', help='Optionally, export the HailTable as a BGZ-compressed TSV file',
                        action='store_true')
    parser.add_argument('--overwrite', help='Overwrite pre-existing data', action='store_true')
    parser.add_argument('--default_ref_genome', help='Default reference genome to start Hail',
                        type=str, default='GRCh38')
    parser.add_argument('--split_multi_allelic', help='Split multi allelic variants if any',
                        action='store_true')

    args = parser.parse_args()

    main(args)
