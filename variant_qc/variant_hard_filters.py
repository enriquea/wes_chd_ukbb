"""

Annotate different filter flags (e.g. hard filters) useful for
downstream analysis.

Expected a bi-allelic variant Hail Table (keyed) by <locus, alleles>,
annotated with VEP and QCed by 'variant_qc'.

"""

import hail as hl
import logging
import functools
import operator

from utils.constants import *

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# start hail on cluster
hl.init(default_reference='GRCh38')

# work dir
dir = 'file:///home/ubuntu'

# read HT
tb = hl.read_table(f'{dir}/data/hail_data/hts/chd_ukbb_variants_veped_raw_qced_v2_16092020.ht')

logger.info(
    "Processing input dataset with {} variants".format(tb.count())
)

# read table interval SureSelect V5/UKBB overlap
agilent_ukbb = hl.read_table(f'{dir}/data/resources/intervals/agilent_ukbb/agilent_ukbb_overlap.no_alt_affected.bed.ht')

# define hard filters expressions
filter_exprs = {'AC0': (tb.raw_qc['variant_qc.dp_stats.max'] < 10) & 
                        (tb.raw_qc['variant_qc.gq_stats.max'] < 20),
                'is_coding_csq': hl.set(CSQ_CODING).contains(tb.Consequence),
                'is_defined_ssv5_ukbb': hl.is_defined(agilent_ukbb[tb.locus]),
                'high_callrate_99': (tb.raw_qc['variant_qc.call_rate'] > 0.99),
                'is_rare': (tb.raw_qc['variant_qc.AF'][1] < 0.001),
                'is_ultrarare_ac5': (tb.raw_qc['variant_qc.AC'][1] <= 5)}

# annotate hard filters
tb = (tb
     .annotate(hard_filters=hl.struct(**filter_exprs))
     )

tb.describe()
                             
logger.info(
    "Flagged {} variants out of {} by AC0 filter..."
    .format(tb.filter(tb.hard_filters.AC0).count(), tb.count())
)
                             
                             
logger.info(
    "Flagged {} variants out of {} by is_defined_ssv5_ukbb filter..."
    .format(tb.filter(tb.hard_filters.is_defined_ssv5_ukbb).count(), tb.count())
)
                             
logger.info(
    "Flagged {} variants out of {} by is_coding_csq filter..."
    .format(tb.filter(tb.hard_filters.is_coding_csq).count(), tb.count())
)

                             
combined_filter_expr =  [~tb.hard_filters.AC0,
                         tb.hard_filters.is_coding_csq,
                         tb.hard_filters.is_defined_ssv5_ukbb]
                             
logger.info(
    "Flagged {} variants out of {} by AC0, is_defined_ssv5_ukbb, is_coding_csq filters"
    .format(tb.filter(functools.reduce(operator.iand, combined_filter_expr)).count(), tb.count())
)
                   
# write annotated HT
tb.write(f'{dir}/data/hail_data/hts/chd_ukbb_variants_veped_rawQCed_hardFilters_v3_21092020.ht')


