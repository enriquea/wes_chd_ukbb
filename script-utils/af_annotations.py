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

from utils.config import NFS_DIR, HDFS_DIR

nfs_dir = NFS_DIR
nfs_tmp = f'{NFS_DIR}/tmp'
hdfs_dir = f'{HDFS_DIR}/dir/hail_data'


def load_af_tables() -> tuple:
    """Load variant table and all external allelic-frequency tables, return them as a tuple."""
    variant_ht = get_vep_annotation_ht()
    ht_ger_af = get_germ_af_ht()
    bonn_af = get_bonn_af_ht()
    ht_rumc_af = get_rum_af_ht()
    gnomad_af_ht = get_gnomad_genomes_v3_af_ht()
    return variant_ht, ht_ger_af, bonn_af, ht_rumc_af, gnomad_af_ht


def annotate_af(
    variant_ht: hl.Table,
    ht_ger_af: hl.Table,
    bonn_af: hl.Table,
    ht_rumc_af: hl.Table,
    gnomad_af_ht: hl.Table,
    date: str,
) -> hl.Table:
    """Build AF annotation expressions, annotate variant table, select AF fields, and add globals."""
    ## build af annotation expression
    af_ann_expr = {
        'ger_af': ht_ger_af[variant_ht.key].GERPOP_AF,
        'rumc_af': ht_rumc_af[variant_ht.key].RUMC_AF,
        'bonn_af': bonn_af[variant_ht.key].AF,
        'gnomad_genomes_af': gnomad_af_ht[variant_ht.key].AF
    }

    # NOTE: AF gnomad exome fields are already VEP-annotated
    gnomad_exomes_af_fields = ['gnomAD_AF',
                               'gnomAD_SAS_AF',
                               'gnomAD_OTH_AF',
                               'gnomAD_NFE_AF',
                               'gnomAD_FIN_AF']

    gnomad_exomes_af_expr = {f: hl.parse_float(variant_ht.vep[f]) for f in gnomad_exomes_af_fields}

    # add gnomad exomes AF expression to annotation dict
    af_ann_expr.update(gnomad_exomes_af_expr)

    ## annotate afs
    variant_ht = (variant_ht
                  .annotate(**af_ann_expr)
                  )

    af_fields = list(af_ann_expr.keys())
    variant_ht = (variant_ht
                  .select(*af_fields)
                  )

    ## add global annotation
    global_ann_expr = {'date': date,
                       'af_fields': af_fields}
    variant_ht = (variant_ht
                  .annotate_globals(**global_ann_expr)
                  )

    return variant_ht


def export_af_table(variant_ht: hl.Table, output_path_ht: str) -> None:
    """Checkpoint the AF-annotated table to disk and export a BGZ-compressed TSV."""
    ## export af table

    # write to Hail table
    variant_ht = (variant_ht
                  .checkpoint(output_path_ht, overwrite=True)
                  )

    # write to TSV file
    (variant_ht
     .export(f'{output_path_ht}.tsv.bgz')
     )


def main() -> None:
    """Orchestrate AF annotation: init Hail, load tables, annotate, export, stop Hail."""
    hl.init(default_reference='GRCh38')

    ## import variant table and af tables
    variant_ht, ht_ger_af, bonn_af, ht_rumc_af, gnomad_af_ht = load_af_tables()

    ## add global annotation
    date = current_date()

    variant_ht = annotate_af(variant_ht, ht_ger_af, bonn_af, ht_rumc_af, gnomad_af_ht, date)

    ## export af table
    output_path_ht = f'{nfs_dir}/hail_data/hts/chd_ukbb.variants.af.annotations.external.{date}.ht'

    export_af_table(variant_ht, output_path_ht)

    hl.stop()


if __name__ == '__main__':
    main()
