# eam
# 2021-02-16

from gnomad.utils.liftover import *

import hail as hl

project_dir = '~/projects/github/wes_chd_ukbb'


hl.init()


# prepare genome references
rg37 = hl.get_reference('GRCh37')
rg38 = hl.get_reference('GRCh38')
rg37.add_liftover(f'{project_dir}/data/resources/grch37_to_grch38.over.chain.gz', rg38)


# lift over german AFs (hg37 -> hg38)
path_ht_ger_af_hg37 = f'{project_dir}/data/resources/german_pop_af.ht'
output_path = f'{project_dir}/data/resources/german_af_hg38'

lift_data(t=hl.read_table(path_ht_ger_af_hg37),
          gnomad=False,
          data_type='exomes',
          rg=rg38,
          path=output_path,
          overwrite=False)


# lift over RUMC exomes AFs (hg37 -> hg38)
path_ht_rumc_af_hg37 = f'{project_dir}/data/resources/rumc_af_18032020.ht'
output_path = f'{project_dir}/data/resources/rumc_af_hg38'

lift_data(t=hl.read_table(path_ht_rumc_af_hg37),
          gnomad=False,
          data_type='exomes',
          rg=rg38,
          path=output_path,
          overwrite=False)


