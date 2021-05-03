"""
Utils general-propose Hail expression
"""

import logging
from typing import List, Optional, Tuple, Union

import hail as hl


def bi_allelic_expr(t: Union[hl.Table, hl.MatrixTable]) -> hl.expr.BooleanExpression:
    """
    Returns a boolean expression selecting bi-allelic sites only,
    accounting for whether the input MT/HT was split.
    :param t: Input HT/MT
    :return: Boolean expression selecting only bi-allelic sites
    """
    return ~t.was_split if "was_split" in t.row else (hl.len(t.alleles) == 2)


def af_filter_expr(t: Union[hl.Table, hl.MatrixTable],
                   af_field: str,
                   af_cutoff: float = 0.001) -> hl.expr.BooleanExpression:
    """
    Returns a boolean expression filtering (rare) sites by allelic frequency.
    Note: keep sites with missing maf.

    :param t: Input HT or MT
    :param af_field: Name of the allelic frequency field
    :param af_cutoff: Cutoff to apply
    :return: Boolean expression selecting sites with AF <= cutoff
    """
    return (t[af_field] <= af_cutoff) | hl.is_missing(t[af_field])
