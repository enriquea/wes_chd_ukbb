# eam
# 2021-02-17

"""
This script generates a Hail Table from the dbNSFP v4.1a.
Used for annotations and downstream analysis.

Source: http://database.liulab.science/dbNSFP

It takes as input an unique block-compressed gz file.
It can be generated as follow:
zcat *.gz | bgzip -c > dbNSFP4.1a_variant.bgz
"""

import hail as hl
import os

from utils.config import NFS_DIR


def import_dbnsfp_table(path_file_in: str) -> hl.Table:
    """Import the dbNSFP BGZ-compressed flat file as a Hail Table."""
    ht = hl.import_table(paths=path_file_in,
                         min_partitions=1000,
                         impute=False,
                         missing='.',
                         force_bgz=True)
    return ht


def parse_chromosome_and_variant_key(ht: hl.Table) -> hl.Table:
    """Rename the chr field, build a variant key, and parse it into locus/alleles."""
    # parse chromosome field
    ht = ht.rename({'#chr': 'chr'})
    ht = ht.annotate(chr='chr' + hl.str(ht['chr']))

    # annotate variant field (chr:pos:ref:alt)
    variant_key_expr = hl.array([ht.chr,
                                 hl.str(ht['pos(1-based)']),
                                 ht.ref,
                                 ht.alt])

    ht = ht.annotate(variant_key=hl.delimit(variant_key_expr, ':'))

    # parse variant field (-> locus, alleles)
    ht = (ht
          .annotate(**hl.parse_variant(ht.variant_key, reference_genome='GRCh38'))
          )
    return ht


def rekey_table(ht: hl.Table) -> hl.Table:
    """Rearrange table fields and generate a table keyed by locus and alleles."""
    # rearrange tables fields and generate keyed table
    tb_fields = ht.row

    ht = (ht.select('locus',
                    'alleles',
                    *[f for f in tb_fields if f not in ['locus', 'alleles']])
          .key_by('locus', 'alleles')
          )
    return ht


def map_scores_to_transcripts(ht: hl.Table) -> hl.Table:
    """Split Ensembl transcript IDs and map score fields to per-transcript dicts."""
    # Map scores to transcript. If a score is site-specific rather than
    # transcript-specific, map the same value to all transcripts.
    # TODO: split genename, protein name, VEP_canonical ect...
    ht = (ht
          .annotate(Ensembl_transcriptid=ht.Ensembl_transcriptid.split(";"))
          )

    scores_fields = [f for f in ht.row if f.endswith('_score') or f == 'CADD_phred']

    ht = (ht
          .annotate(**{f: hl.if_else(ht[f].contains(";"),
                                     hl.dict(hl.zip(ht.Ensembl_transcriptid,
                                                    hl.map(lambda x:
                                                           hl.parse_float(x),
                                                           ht[f].split(";")))),
                                     hl.dict(hl.zip(ht.Ensembl_transcriptid,
                                                    hl.map(lambda x:
                                                           hl.parse_float(ht[f]),
                                                           ht.Ensembl_transcriptid)))
                                     )
                       for f in scores_fields
                       })
          )
    return ht


def nest_field_groups(ht: hl.Table, tb_fields) -> hl.Table:
    """Transmute related fields into nested structs and print the table schema."""
    # nest related fields into structures (easier to analyse/filter later)
    field_groups = ['gnomAD', 'ExAC', '1000Gp3', 'ESP6500', 'clinvar']  # fold these groups into structures
    # TODO: split parse int, float fields
    ht = (ht
          .transmute(**{f_root: hl.struct(**{field: ht[field]
                                             for field in tb_fields if field.startswith(f_root)}
                                          )
                        for f_root in field_groups}
                     )
          )

    ht.describe()
    return ht


def main() -> None:
    """Orchestrate dbNSFP table import, parsing, annotation, and export."""
    hl.init(default_reference='GRCh38')

    nfs_dir = NFS_DIR

    path_file_in = f'{nfs_dir}/resources/dbNSFP/variants/dbNSFP4.1a_variant.bgz'
    path_ht_out = f'{os.path.splitext(path_file_in)[0]}.ht'

    ht = import_dbnsfp_table(path_file_in)

    ht = parse_chromosome_and_variant_key(ht)

    ht = rekey_table(ht)

    # capture tb_fields before score mapping changes the schema
    tb_fields = ht.row

    ht = map_scores_to_transcripts(ht)

    ht = nest_field_groups(ht, tb_fields)

    # export table
    ht.write(path_ht_out,
             overwrite=True)

    hl.stop()


if __name__ == '__main__':
    main()
