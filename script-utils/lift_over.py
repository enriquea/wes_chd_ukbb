# eam
# 2021-02-16

"""
Use the gnomad liftover utilities to translate
some data from one genome coordinate to other.
"""

from gnomad.utils.liftover import *

from utils.config import NFS_DIR

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def setup_liftover_references(nfs_dir: str) -> tuple:
    """Set up GRCh37 and GRCh38 references and add liftover chain."""
    rg37 = hl.get_reference("GRCh37")
    rg38 = hl.get_reference("GRCh38")

    if not rg37.has_liftover("GRCh38"):
        rg37.add_liftover(
            f'{nfs_dir}/resources/liftover/grch37_to_grch38.over.chain.gz', rg38
        )

    return rg37, rg38


def lift_german_af(nfs_dir: str, rg38: hl.genetics.ReferenceGenome) -> None:
    """Lift over German population AFs from hg37 to hg38."""
    path_ht_ger_af_hg37 = f'{nfs_dir}/resources/annotation/german_pop_af.ht'
    output_path = f'{nfs_dir}/resources/annotation/german_af_hg38'

    lift_data(t=hl.read_table(path_ht_ger_af_hg37),
              gnomad=False,
              data_type='exomes',
              rg=rg38,
              path=output_path,
              overwrite=True)


def lift_rumc_af(nfs_dir: str, rg38: hl.genetics.ReferenceGenome) -> None:
    """Lift over RUMC exomes AFs from hg37 to hg38."""
    path_ht_rumc_af_hg37 = f'{nfs_dir}/resources/annotation/rumc_af_18032020.ht'
    output_path = f'{nfs_dir}/resources/annotation/rumc_af_hg38'

    lift_data(t=hl.read_table(path_ht_rumc_af_hg37),
              gnomad=False,
              data_type='exomes',
              rg=rg38,
              path=output_path,
              overwrite=True)


def lift_bonn_af(nfs_dir: str, rg38: hl.genetics.ReferenceGenome) -> None:
    """Lift over Bonn exomes AFs from hg37 to hg38."""
    path_ht_bonn_af_hg37 = f'{nfs_dir}/resources/annotation/Cohort_Bonn_AF_0521.ht'
    output_path = f'{nfs_dir}/resources/annotation/Cohort_Bonn_AF_0521'

    lift_data(t=hl.read_table(path_ht_bonn_af_hg37),
              gnomad=False,
              data_type='exomes',
              rg=rg38,
              path=output_path,
              overwrite=True)


def lift_denovo(nfs_dir: str, rg38: hl.genetics.ReferenceGenome) -> None:
    """Lift over de novo variant data from hg37 to hg38."""
    path_denovo_ht = f'{nfs_dir}/resources/denovo/DNM_Jin2017_Sifrim2016.ht'
    output_path = f'{nfs_dir}/resources/denovo/DNM_Jin2017_Sifrim2016_GRCh38'

    lift_data(t=hl.read_table(path_denovo_ht),
              gnomad=False,
              data_type='exomes',
              rg=rg38,
              path=output_path,
              overwrite=True)


def main() -> None:
    """Run all liftover operations from GRCh37 to GRCh38."""
    hl.init()

    nfs_dir = NFS_DIR
    project_dir = f'{NFS_DIR}/projects/wes_chd_ukbb'  # noqa: F841

    rg37, rg38 = setup_liftover_references(nfs_dir)

    # lift over german AFs (hg37 -> hg38)
    lift_german_af(nfs_dir, rg38)

    # lift over RUMC exomes AFs (hg37 -> hg38)
    lift_rumc_af(nfs_dir, rg38)

    # lift over Bonn exomes AFs (hg37 -> hg38)
    lift_bonn_af(nfs_dir, rg38)

    # lift over denovo data (hg37 -> hg38)
    lift_denovo(nfs_dir, rg38)

    hl.stop()


if __name__ == '__main__':
    main()
