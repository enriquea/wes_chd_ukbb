"""

general-propose filter functions

"""

from typing import List, Optional, Tuple, Union
from utils.reference_genome import get_reference_genome

from utils.reference_data import (get_lcr_ht,
                                  get_segdups_ht,
                                  get_telomeres_and_centromeres_ht,
                                  import_cds_intervals_from_gtf)

from utils.data_utils import get_capture_interval_ht

import functools
import operator

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


def filter_to_cds_regions(t: Union[hl.Table, hl.MatrixTable]) -> Union[hl.Table, hl.MatrixTable]:
    """
    Keep sites within CDS regions. CDS regions were extracted from gencode v29 GTF file
    (currently support only hg38).

    :param t:  MatrixTable or Table to filter
    :return:  MatrixTable or Table
    """
    cds_intervals = import_cds_intervals_from_gtf(overwrite=False)

    if isinstance(t, hl.MatrixTable):
        t = t.filter_rows(hl.is_defined(cds_intervals[t.locus]))
    else:
        t = t.filter(hl.is_defined(cds_intervals[t.locus]))
    return t


def remove_telomeres_centromes(t: Union[hl.Table, hl.MatrixTable]) -> Union[hl.Table, hl.MatrixTable]:
    """
    Remove sites overlapping telomeres and centromeres regions

    :param t:  MatrixTable or Table to filter
    :return:  MatrixTable or Table
    """
    tc_intervals = get_telomeres_and_centromeres_ht(overwrite=False)

    if isinstance(t, hl.MatrixTable):
        t = t.filter_rows(hl.is_missing(tc_intervals[t.locus]))
    else:
        t = t.filter(hl.is_missing(tc_intervals[t.locus]))
    return t


# TODO: filter_low_conf_region
def filter_low_conf_regions(
        mt: Union[hl.MatrixTable, hl.Table],
        filter_lcr: bool = True,
        filter_segdup: bool = True) -> Union[hl.MatrixTable, hl.Table]:
    """
    Filters low-confidence regions

    :param mt: MatrixTable or Table to filter
    :param filter_lcr: Whether to filter LCR regions
    # :param filter_decoy: Whether to filter decoy regions
    :param filter_segdup: Whether to filter Segdup regions
    # :param filter_exome_low_coverage_regions: Whether to filter exome low confidence regions
    # :param high_conf_regions: Paths to set of high confidence regions to restrict to (union of regions)
    :return: MatrixTable or Table with low confidence regions removed
    """

    criteria = []
    if filter_lcr:
        lcr = get_lcr_ht()
        criteria.append(hl.is_missing(lcr[mt.locus]))

    if filter_segdup:
        segdup = get_segdups_ht()
        criteria.append(hl.is_missing(segdup[mt.locus]))

    # if filter_decoy:
    #    decoy = resources.decoy_intervals.ht()
    #    criteria.append(hl.is_missing(decoy[mt.locus]))

    # if filter_exome_low_coverage_regions:
    #    high_cov = resources.high_coverage_intervals.ht()
    #    criteria.append(hl.is_missing(high_cov[mt.locus]))

    # if high_conf_regions is not None:
    #    for region in high_conf_regions:
    #        region = hl.import_locus_intervals(region)
    #        criteria.append(hl.is_defined(region[mt.locus]))

    if criteria:
        filter_criteria = functools.reduce(operator.iand, criteria)
        if isinstance(mt, hl.MatrixTable):
            mt = mt.filter_rows(filter_criteria)
        else:
            mt = mt.filter(filter_criteria)

    return mt


def filter_capture_intervals(t: Union[hl.Table, hl.MatrixTable],
                             capture_intervals: tuple = ('ssv2', 'ssv3', 'ssv4', 'ssv5_idt_intersect'),
                             reference: str = "GRCh38") -> Union[hl.Table, hl.MatrixTable]:
    """
    Filter sites defined in capture intervals. Sites are kept if they are defined in all input captures intervals.
    Note: Both SureSelect V5 and IDT (UKBB) have been intersected and the known 'bad' regions from UKBB removed.
          Thus, they are provided here as a joint capture interval HT.

    :param t: MatrixTable or Table to filter
    :param capture_intervals: List of capture intervals names. Defined in the module `intervals`
    :param reference: Genome reference
    :return: MatrixTable or Table
    """

    # evaluate if input capture interval names are valid...
    defined_captures_interval = ('ssv2', 'ssv3', 'ssv4', 'ssv5_idt_intersect')
    if any([x not in defined_captures_interval for x in capture_intervals]):
        raise ValueError(
            f"Invalid input capture interval name(s). Expected capture interval(s): {defined_captures_interval}"
        )

    # list of intervals HT
    capture_intervals_hts = [get_capture_interval_ht(name=x, reference=reference)
                             for x in capture_intervals]

    # build filter expression
    filter_interval_expr = [hl.is_defined(ht_interval[t.locus]) for
                            ht_interval in capture_intervals_hts]

    # filtering
    # Note: Keep locus that are defined in all (input) intervals (intersect)
    if filter_interval_expr:
        filter_criteria = functools.reduce(operator.iand, filter_interval_expr)
        if isinstance(t, hl.MatrixTable):
            t = t.filter_rows(filter_criteria)
        else:
            t = t.filter(filter_criteria)

    return t
