# eam
# 16.03.22


import hail as hl

from utils.data_utils import (get_variant_qc_ht_path,
                              get_sample_qc_ht_path)
from utils.expressions import (af_filter_expr)
from utils.filter import (filter_low_conf_regions,
                          remove_telomeres_centromes)


def apply_sample_qc_filtering(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Apply sample QC filtering

    :param mt: hl.MatrixTable
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

