# eam
# 16.03.22


import hail as hl

from utils.data_utils import (get_variant_qc_ht_path,
                              get_sample_qc_ht_path)
from utils.expressions import (af_filter_expr)
from utils.filter import (filter_low_conf_regions,
                          remove_telomeres_centromes)


def apply_sample_qc_filtering(mt: hl.MatrixTable,
                              keep_rare_variants: bool = True,
                              maf_threshold: float = 0.01) -> hl.MatrixTable:
    """
    Apply sample QC filtering, compute internal allelic frequencies on samples passing qc and
    adjusted phenotypes. Optionally, return MT filtered to rare variants.

    :param mt: hl.MatrixTable
    :param keep_rare_variants: Filter MT to rare variants
    :param maf_threshold: allelic frequency cutoff
    :return: hl.MatrixTable
    """
    # import variant qc final table
    sample_qc_ht = hl.read_table(
        get_sample_qc_ht_path(part='final_qc')
    )
    sample_qc_ht = (sample_qc_ht
                    .filter(sample_qc_ht.pass_filters)
                    )
    mt = (mt
          .filter_cols(hl.is_defined(sample_qc_ht[mt.col_key]))
          )
    # compute cohort-specific (internal) allelic frequencies on samples passing qc
    mt = (mt
          .annotate_rows(gt_stats=hl.agg.call_stats(mt.GT, mt.alleles))
          )
    mt = (mt
          .annotate_rows(internal_af=mt.gt_stats.AF[1],
                         internal_ac=mt.gt_stats.AC[1])
          )
    # filter out common variants base don internal af
    if keep_rare_variants:
        mt = (mt
              .filter_rows(af_filter_expr(mt, 'internal_af', maf_threshold))
              )

    return mt


def apply_variant_qc_filtering(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Apply variant QC filtering

    :param mt: hl.MatrixTable
    :return: hl.MatrixTable
    """
    # import variant qc final table
    variant_qc_ht = hl.read_table(
        get_variant_qc_ht_path(part='final_qc')
    )
    mt = (mt
          .annotate_rows(**variant_qc_ht[mt.row_key])
          )
    mt = (mt
          .filter_rows(~mt.fail_inbreeding_coeff &
                       ~mt.AC0 &
                       ~mt.fail_vqsr &
                       ~mt.fail_rf &
                       mt.is_coveraged_gnomad_genomes &
                       mt.is_defined_capture_intervals)
          )
    # filter low conf regions
    mt = filter_low_conf_regions(
        mt,
        filter_lcr=True,  # TODO: include also decoy and low coverage exome regions
        filter_segdup=True
    )
    # filter telomeres/centromes
    mt = remove_telomeres_centromes(mt)

    return mt

