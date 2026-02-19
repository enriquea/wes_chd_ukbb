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

from utils.vep import (get_vep_fields,
                       pick_transcript,
                       vep_protein_domain_ann_expr)

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


def import_vcf_and_split(vcf_path, force_bgz, split_multi_allelic):
    """Import VEP-annotated VCF as MatrixTable and optionally split multi-allelic sites."""
    mt = hl.import_vcf(path=vcf_path, force_bgz=force_bgz)

    if split_multi_allelic:
        mt = hl.split_multi_hts(mt)
        mt = mt.annotate_rows(info=mt.info.annotate(**{field: mt.info[field][mt.a_index - 1]
                                                       for field in INFO_FIELDS}))

    return mt


def parse_csq_transcripts(mt, csq_field, vep_fields):
    """Extract rows, parse the CSQ field into per-transcript dicts, and filter by allele index."""
    tb_csq = mt.rows()
    tb_csq = (tb_csq
              .annotate(csq_raw=tb_csq.info[csq_field])
              )

    tb_csq = (tb_csq
              .annotate(csq_raw=tb_csq.csq_raw.map(lambda x:
                                                   hl.dict(hl.zip(vep_fields, x.split('[|]')))
                                                   ))
              )

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

    return tb_csq


def pick_and_expand_transcript(tb_csq, vep_fields):
    """Pick one transcript per variant, expand fields, parse Consequence and DOMAINS, drop temp fields."""
    tb_csq = pick_transcript(ht=tb_csq, csq_array='csq_raw')

    tb_csq = annotate_from_dict(ht=tb_csq, dict_field='tx', output_filed='vep')

    tb_csq = (tb_csq
              .annotate(vep=tb_csq.vep.annotate(Consequence=tb_csq.vep.Consequence.split('&')[0]))
              )

    if 'DOMAINS' in vep_fields:
        tb_csq = (tb_csq
                  .annotate(vep=tb_csq.vep.annotate(DOMAINS=vep_protein_domain_ann_expr(tb_csq.vep['DOMAINS'])))
                  )

    tb_csq = (tb_csq
              .drop('csq_raw', 'tx')
              .repartition(500)
              )

    return tb_csq


def checkpoint_and_export(ht, output_path, overwrite, write_to_file):
    """Checkpoint table and optionally export as BGZ-compressed TSV."""
    ht = (ht
          .checkpoint(output=output_path,
                      overwrite=overwrite)
          )

    if write_to_file:
        ht.export(f'{output_path}.tsv.bgz')

    return ht


def main(args):
    # Init Hail
    hl.init(default_reference=args.default_ref_genome)

    vep_fields = get_vep_fields(vcf_path=get_vep_vqsr_vcf_path(),
                                vep_csq_field=args.csq_field)

    mt = import_vcf_and_split(get_vep_vqsr_vcf_path(), args.force_bgz, args.split_multi_allelic)

    tb_csq = parse_csq_transcripts(mt, args.csq_field, vep_fields)

    tb_csq = pick_and_expand_transcript(tb_csq, vep_fields)

    # print fields overview
    tb_csq.describe()

    output_path = get_variant_qc_ht_path(part='vep_vqsr', split=args.split_multi_allelic)
    checkpoint_and_export(tb_csq, output_path, args.overwrite, args.write_to_file)

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
