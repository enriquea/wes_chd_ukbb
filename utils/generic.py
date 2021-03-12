"""

Useful generic functions

"""

import hail as hl
import time


# unphase MatrixTable (from gnomad_methods)
def unphase_mt(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Generate unphased version of MatrixTable (assumes call is in mt.GT and is diploid or haploid only)
    """
    return mt.annotate_entries(GT=hl.case()
                               .when(mt.GT.is_diploid(), hl.call(mt.GT[0], mt.GT[1], phased=False))
                               .when(mt.GT.is_haploid(), hl.call(mt.GT[0], phased=False))
                               .default(hl.null(hl.tcall))
                               )


# print current date
def current_date() -> str:
    return time.strftime("%d%m%Y")


# get Hail version
def hail_version() -> str:
    return hl.__version__
