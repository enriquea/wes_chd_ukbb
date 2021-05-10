# eam
# 2021-02-16

"""

Use the gnomad liftover utilities to translate
some data from one genome coordinate to other.

"""

from gnomad.utils.liftover import *

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

hl.init()

nfs_dir = 'file:///home/ubuntu/data'
project_dir = 'file:///home/ubuntu/data/projects/wes_chd_ukbb'

rg37 = hl.get_reference("GRCh37")
rg38 = hl.get_reference("GRCh38")

if not rg37.has_liftover("GRCh38"):
    rg37.add_liftover(
        f'{nfs_dir}/resources/liftover/grch37_to_grch38.over.chain.gz', rg38
    )

# lift over german AFs (hg37 -> hg38)
path_ht_ger_af_hg37 = f'{nfs_dir}/resources/german_pop_af.ht'
output_path = f'{nfs_dir}/resources/german_af_hg38'

lift_data(t=hl.read_table(path_ht_ger_af_hg37),
          gnomad=False,
          data_type='exomes',
          rg=rg38,
          path=output_path,
          overwrite=False)


# lift over RUMC exomes AFs (hg37 -> hg38)
path_ht_rumc_af_hg37 = f'{nfs_dir}/resources/rumc_af_18032020.ht'
output_path = f'{nfs_dir}/resources/rumc_af_hg38'

lift_data(t=hl.read_table(path_ht_rumc_af_hg37),
          gnomad=False,
          data_type='exomes',
          rg=rg38,
          path=output_path,
          overwrite=False)


# lift over Bonn exomes AFs (hg37 -> hg38)
path_ht_bonn_af_hg37 = f'{nfs_dir}/resources/Cohort_Bonn_AF_0521.ht'
output_path = f'{nfs_dir}/resources/Cohort_Bonn_AF_0521'

lift_data(t=hl.read_table(path_ht_bonn_af_hg37),
          gnomad=False,
          data_type='exomes',
          rg=rg38,
          path=output_path,
          overwrite=False)


# lift over denovo data (hg37 -> hg38)
path_denovo_ht = f'{nfs_dir}/resources/denovo/DNM_Jin2017_Sifrim2016.ht'
output_path = f'{nfs_dir}/resources/denovo/DNM_Jin2017_Sifrim2016_GRCh38'

lift_data(t=hl.read_table(path_denovo_ht),
          gnomad=False,
          data_type='exomes',
          rg=rg38,
          path=output_path,
          overwrite=False)

hl.stop()
