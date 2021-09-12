# eam
# 2021-09-12

"""

Utils functions to parse and process VEP annotations

"""

import hail as hl

# Main pre-defined protein domain source of annotations
PROTEIN_DOMAIN_DB = ['CDD',
                     'Pfam',
                     'SMART',
                     'PANTHER']


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


# select one a transcript from array using pre-defined rules/conditions.
def pick_transcript(ht: hl.Table,
                    csq_array: str,
                    pick_protein_coding: bool = True,
                    pick_canonical: bool = True) -> hl.Table:
    # TODO: This function could be improved by scanning the array (just once) and sorting it as suggested here:
    # TODO: https://hail.zulipchat.com/#narrow/stream/123010-Hail-0.2E2.20support/topic/pick.20transcript.20from.20array
    # TODO: /near/190400193

    """
    Annotate an extra field (tx) with the selected transcript.
    This function will pick one transcript per variant/consequence based on the impact of the variant in the transcript
    (from more severe to less severe).

    :param pick_canonical: keep only canonical transcript(s)
    :param pick_protein_coding: keep only protein-coding transcript(s)
    :param ht: Hail table with VEP annotations
    :param csq_array: Parsed CSQ field name. Expected to be an array of dict(s). Transcript expected to be as a dict.
    :return: Hail table with an annotated extra field (tx). The transcript selected from the array based on a set of
    pre-defined criteria.
    """

    # getting current keys from dict
    keys = ht[csq_array].take(1)[0][0]

    # filter to protein-coding transcript(s)
    if pick_protein_coding and 'BIOTYPE' in keys:
        ht = (ht
              .annotate(**{csq_array:
                               ht[csq_array].filter(lambda x: x['BIOTYPE'] == 'protein_coding')
                           }
                        )
              )

    # filter to canonical transcript(s)
    if pick_canonical and 'CANONICAL' in keys:
        ht = (ht
              .annotate(**{csq_array:
                               ht[csq_array].filter(lambda x: x['CANONICAL'] == 'YES')
                           }
                        )
              )

    # Set transcript (tx) field initially to 'NA' and update it sequentially based on a set of pre-defined criteria
    # (order matters)
    ht = (ht
          .annotate(tx=ht[csq_array].find(lambda x: False))
          )

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

    # if tx is still missing, set tx as the first annotated transcript
    ht = (ht
          .annotate(tx=hl.cond(hl.is_missing(ht.tx) & (hl.len(ht[csq_array]) > 0),
                               ht[csq_array][0],
                               ht.tx)
                    )
          )
    return ht


def vep_protein_domain_ann_expr(s: hl.expr.StringExpression) -> hl.expr.DictExpression:
    """
    Parse and annotate protein domain(s) from VEP annotation.
    Expected StringExpression as input (e.g. 'Pfam:PF13853&Prints:PR00237&PROSITE_profiles:PS50262')
    It will generate a dict<k,v> where keys (k) represent source/database and values (v) the annotated domain_id.

    :param s: hl.expr.StringExpression
    :return: hl.expr.DictExpression
    """
    a1 = s.split(delim="&")

    # keep only well-annotated domain(s) (i.e. <source:domain_id>)
    a2 = a1.map(lambda x: x.split(delim=":"))
    a2 = a2.filter(lambda x: x.length() == 2)

    d = (hl.case()
         .when(hl.len(a2) > 0,
               hl.dict(hl.zip(a2.map(lambda x: x[0]),  # TODO: Optimize by scanning array just one.
                              a2.map(lambda x: x[1]))))
         .or_missing()
         )

    return d


def vep_protein_domain_filter_expr(d: hl.expr.DictExpression) -> hl.expr.BooleanExpression:
    """

    Return True of False if any protein domain source(s) are contained within pre-defined protein domain sources.
    Expected as input dict<k,v> where keys (k) represent source/database and values (v) the annotated domain_id.

    :param d: hl.DictExpression
    :return: hl.BoolExpression
    """
    domain_dbs = hl.set(PROTEIN_DOMAIN_DB)
    return (d.key_set()
             .intersection(domain_dbs)
             .length() >= 1)

