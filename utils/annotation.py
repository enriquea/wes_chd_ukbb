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
