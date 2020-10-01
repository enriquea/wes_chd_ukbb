"""

general-propose filter functions

"""

from typing import List, Optional, Tuple, Union
from utils.reference_genome import get_reference_genome

import hail as hl


def filter_to_autosomes(
        t: Union[hl.MatrixTable, hl.Table]
) -> Union[hl.MatrixTable, hl.Table]:
    """
    Filters the Table or MatrixTable to autosomes only.
    This assumes that the input contains a field named `locus` of type Locus
    :param t: Input MT/HT
    :return:  MT/HT autosomes
    """
    reference = get_reference_genome(t.locus)
    autosomes = hl.parse_locus_interval(
        f"{reference.contigs[0]}-{reference.contigs[21]}", reference_genome=reference
    )
    return hl.filter_intervals(t, [autosomes])


def filter_genotypes_ab(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Filter high-quality genotypes based on allelic-balance expression.
    Expected AD and GT in entries fields.
    Rules:
      hom_ref: ab <= 0.1
      hets: 0.2 >= ab <= 0.8
      hom_var: ab >= 0.9

    :param mt: Input MT
    :return: Genotype-filtered MT
    """
    ab = mt.AD[1] / hl.sum(mt.AD)
    filter_condition_ab = ((mt.GT.is_hom_ref() & (ab <= 0.1)) |
                           (mt.GT.is_het() & (ab >= 0.2) & (ab <= 0.8)) |
                           (mt.GT.is_hom_var() & (ab >= 0.9)))
    return mt.filter_entries(filter_condition_ab, keep=True)


