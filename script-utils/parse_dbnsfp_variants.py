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

hl.init(default_reference='GRCh38')

nfs_dir = 'file:///home/ubuntu/data'

path_file_in = f'{nfs_dir}/resources/dbNSFP/variants/dbNSFP4.1a_variant.bgz'
path_ht_out = f'{os.path.splitext(path_file_in)[0]}.ht'

ht = hl.import_table(paths=path_file_in,
                     min_partitions=1000,
                     impute=False,
                     missing='.',
                     force_bgz=True)

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

# rearrange tables fields and generate keyed table
tb_fields = ht.row

ht = (ht.select('locus',
                'alleles',
                *[f for f in tb_fields if f not in ['locus', 'alleles']])
      .key_by('locus', 'alleles')
      )

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

# export table
ht.write(path_ht_out,
         overwrite=True)

hl.stop()
