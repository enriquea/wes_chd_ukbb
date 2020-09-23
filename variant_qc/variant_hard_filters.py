'''

Annotate different filter flags (e.g. hard filters) useful for
downstream anlysis.

Expected a biallelic variant Hail Table (keyed) by <locus, alleles>,
annotated with VEP and QCed by 'variant_qc'.

'''

import hail as hl
import logging
import functools
import operator

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


# Updated from https://github.com/broadinstitute/gnomad_methods.git (vep.py)
# Note that this is the current as of v81 with some included for backwards compatibility (VEP <= 75)
# Still valid for this version (v97) (checked!)

CSQ_CODING_HIGH_IMPACT = [
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
]

CSQ_CODING_MEDIUM_IMPACT = [
    "start_lost",  # new in v81
    "initiator_codon_variant",  # deprecated
    "transcript_amplification",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",  # new in v79
    "splice_region_variant",
]

CSQ_CODING_LOW_IMPACT = [
    "incomplete_terminal_codon_variant",
    "start_retained_variant",  # new in v92
    "stop_retained_variant",
    "synonymous_variant",
    "coding_sequence_variant",
]

CSQ_CODING = (
        CSQ_CODING_HIGH_IMPACT
        + CSQ_CODING_MEDIUM_IMPACT
        + CSQ_CODING_LOW_IMPACT
)

CSQ_NON_CODING = [
    "mature_miRNA_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant",
    "non_coding_exon_variant",  # deprecated
    "intron_variant",
    "NMD_transcript_variant",
    "non_coding_transcript_variant",
    "nc_transcript_variant",  # deprecated
    "upstream_gene_variant",
    "downstream_gene_variant",
    "TFBS_ablation",
    "TFBS_amplification",
    "TF_binding_site_variant",
    "regulatory_region_ablation",
    "regulatory_region_amplification",
    "feature_elongation",
    "regulatory_region_variant",
    "feature_truncation",
    "intergenic_variant",
]

CSQ_ORDER = (
        CSQ_CODING_HIGH_IMPACT
        + CSQ_CODING_MEDIUM_IMPACT
        + CSQ_CODING_LOW_IMPACT
        + CSQ_NON_CODING
)


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


