from typing import Union

import hail as hl
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
    # unphase MT
    mt = unphase_mt(mt)

    # impute data
    sex_ht = hl.impute_sex(mt.GT,
                           aaf_threshold=0.05,
                           female_threshold=female_threshold,
                           male_threshold=male_threshold,
                           include_par=False)

    sex_colnames = ['f_stat', 'is_female']
    sex_ht = sex_ht.select(*sex_colnames)
    mt = mt.annotate_cols(**sex_ht[mt.col_key])
    return mt


def annotate_variant_id(t: Union[hl.Table, hl.MatrixTable],
                        field_name: str = 'vid') -> Union[hl.Table, hl.MatrixTable]:
    """
    Expected input dataset with bi-allelic variant, and fields `locus` and `alleles`.
    Annotate variant ids as follows 'chr:position:ref:alt'.

    :param field_name: variant id field name
    :param t: dataset
    :return: HailTable or MatrixTable
    """

    variant_id_ann_exp = {
        field_name: hl.delimit([hl.str(t.locus.contig),
                                hl.str(t.locus.position),
                                hl.str(t.alleles[0]),
                                hl.str(t.alleles[1])],
                               delimiter=":")
    }

    if isinstance(t, hl.Table):
        return t.annotate(**variant_id_ann_exp)
    else:
        return t.annotate_rows(**variant_id_ann_exp)


def annotate_sample_blind_id(t: hl.Table,
                             field_name: str = 'biid'):
    """
    Add generic unique IDs to rows (e.g., non-duplicated samples) in input Table.

    :param t: Table
    :param field_name: Output field name

    :return: Hail Table
    """

    t = t.add_index()
    ann_expr = {field_name: 'BIID_' + hl.str(t.idx + 1_000_000)}
    t = t.annotate(**ann_expr).drop('idx')

    return t


def annotate_af(t: Union[hl.Table, hl.MatrixTable],
                field_name: str = 'AF') -> Union[hl.Table, hl.MatrixTable]:
    """
    Annotate AF field in Table or MatrixTable.

    :param t: Table or MatrixTable
    :param field_name: Output field name
    :return: Table or MatrixTable
    """
    if isinstance(t, hl.MatrixTable):
        t = t.annotate_rows(**{field_name: hl.agg.mean(t.GT.n_alt_alleles())})
    else:
        t = t.annotate(**{field_name: hl.agg.mean(t.GT.n_alt_alleles())})

    return t


def annotate_platform(t: Union[hl.Table, hl.MatrixTable],
                      field_name: str = 'platform') -> Union[hl.Table, hl.MatrixTable]:
    """
    Annotate platform field in Table or MatrixTable.

    :param t: Table or MatrixTable
    :param field_name: Output field name
    :return: Table or MatrixTable
    """
    if isinstance(t, hl.MatrixTable):
        t = t.annotate_cols(**{field_name: hl.if_else(hl.is_defined(t.platform),
                                                      t.platform,
                                                      'unknown')})
    else:
        t = t.annotate(**{field_name: hl.if_else(hl.is_defined(t.platform),
                                                 t.platform,
                                                 'unknown')})

    return t


def annotate_ancestry(t: Union[hl.Table, hl.MatrixTable],
                      field_name: str = 'ancestry') -> Union[hl.Table, hl.MatrixTable]:
    """
    Annotate ancestry field in Table or MatrixTable.

    :param t: Table or MatrixTable
    :param field_name: Output field name
    :return: Table or MatrixTable
    """
    if isinstance(t, hl.MatrixTable):
        t = t.annotate_cols(**{field_name: hl.if_else(hl.is_defined(t.ancestry),
                                                      t.ancestry,
                                                      'unknown')})
    else:
        t = t.annotate(**{field_name: hl.if_else(hl.is_defined(t.ancestry),
                                                 t.ancestry,
                                                 'unknown')})

    return t
