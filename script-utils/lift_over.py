# eam
# 2021-02-16

from gnomad.utils.liftover import *

from utils.reference_genome import rg38

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

hl.init()

project_dir = 'file:///home/ubuntu/data/projects/wes_chd_ukbb'


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

hl.stop()
