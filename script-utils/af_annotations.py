# eam
# 2021-09-08

"""

Annotate variant HailTable with allelic frequencies from different (external) sources
(e.g., gnomad exomes and genomes)

"""

import hail as hl

from utils.data_utils import (get_gnomad_genomes_v3_af_ht,
                              get_bonn_af_ht,
                              get_germ_af_ht,
                              get_rum_af_ht,
                              get_vep_annotation_ht)

from utils.generic import current_date

hl.init(default_reference='GRCh38')

nfs_dir = 'file:///home/ubuntu/data'
nfs_tmp = 'file:///home/ubuntu/data/tmp'
hdfs_dir = 'hdfs://spark-master:9820/dir/hail_data'

## import variant table
variant_ht = get_vep_annotation_ht()

## import af tables
# In-hause german allelic frequencies (Tuebingen)
ht_ger_af = get_germ_af_ht()

# In-hause german allelic frequencies (Bonn)
bonn_af = get_bonn_af_ht()

# RUMC cohort allelic frequencies
ht_rumc_af = get_rum_af_ht()

# Gnomad genomes (version 3.0.0) allelic frequencies
gnomad_af_ht = get_gnomad_genomes_v3_af_ht()

## build af annotation expression
af_ann_expr = {
    'ger_af': ht_ger_af[variant_ht.key].GERPOP_AF,
    'rumc_af': ht_rumc_af[variant_ht.key].RUMC_AF,
    'bonn_af': bonn_af[variant_ht.key].AF,
    'gnomad_genomes_af': gnomad_af_ht[variant_ht.key].AF,
    'gnomad_exomes_af': hl.parse_float(variant_ht.vep.gnomAD_AF),
    'gnomad_exomes_ac': hl.parse_int(variant_ht.vep.gnomAD_AC)  # NOTE: AF/AC gnomad fields are already VEP-annotated
}

## annotate afs
variant_ht = (variant_ht
              .annotate(**af_ann_expr)
              )

af_fields = list(af_ann_expr.keys())
variant_ht = (variant_ht
              .select(af_fields)
              )

## add global annotation
date = current_date()
global_ann_expr = {'date': date,
                   'af_fields': af_fields}
variant_ht = (variant_ht
              .annotate_globals(**global_ann_expr)
              )

## export af table
variant_ht.write(
    f'{nfs_dir}/hail_data/hts/chd_ukbb.variants.af.annotations.external.{date}.ht'
)

