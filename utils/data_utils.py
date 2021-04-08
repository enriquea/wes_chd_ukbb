"""
Set of function to retrieve paths from jointly called data as MatrixTable / HailTable
"""

import hail as hl

nfs_dir = 'file:///home/ubuntu/data'
hdfs_dir = 'hdfs://spark-master:9820/dir/hail_data'
hdfs_checkpoint_dir = 'hdfs://spark-master:9820/checkpoint'


def get_mt_data(dataset: str = 'chd_ukbb', part: str = None) -> hl.MatrixTable:
    """
    Return MatrixTable for different version of CHD/UKBB cohort. Multi-allelic sites were split.

    :param dataset: dataset name (Refers just to one ds at the moment, but included to handle multiple ds...)
    :type part: object
    :param part: One of:
                unfiltered: unfiltered MatrixTable with split multi-allelic sites.
                filtered_high_callrate: filtered MT for common SNPs (MAF 0.1%), bi-allelic, callrate > 0.99
    :return MT path

    """
    # TODO: add different MT versions here...
    if part == 'unfiltered':
        return hl.read_matrix_table(f'{nfs_dir}/hail_data/mts/{dataset}_split_v2_09092020.mt')
    elif part == 'filtered_high_callrate':
        return hl.read_matrix_table(f'{nfs_dir}/hail_data/mts/{dataset}.high_callrate.common_snp.biallelic.mt')
    else:
        raise DataException("Select one valid version of the dataset: unfiltered or filtered_high_callrate...")


def get_qc_mt_path(dataset: str = 'chd_ukbb', part: str = None, split=False, ld_pruned=False) -> str:
    """
    Generate path to MT where the results of a QC process should be read/written from/to.

    :param dataset: Name of the exome cohort to be processed (e.g. <chd_ukbb> or <chd_ddd>)
    :param part: Short description of the current QC process
    :param split: Is the MT split?
    :param ld_pruned: Is the MT LD pruned?

    :return: Path to MT
    """
    split = '.split' if split else ''
    ld_pruned = '.ld_pruned' if ld_pruned else ''
    return f'{nfs_dir}/hail_data/mt_qc/{dataset}.qc.{part}{split}{ld_pruned}.mt'


def get_sample_qc_ht_path(dataset: str = 'chd_ukbb', part: str = None) -> str:
    if part not in ('sex_chrom_coverage', 'hard_filters', 'joint_pca_1kg'):
        raise DataException('Expected part one of: sex_chrom_coverage, hard_filters')

    return f'{nfs_dir}/hail_data/sample_qc/{dataset}.sample_qc.{part}.ht'


def get_interval_ht(name: str, reference: str) -> hl.Table:
    """
    Return interval HT from different capture products

    :param name: name/description of the interval HT to be returned.
                 - ssv2: Agilent SureSelect Clinical Exomes V2
                 - ssv3: Agilent SureSelect V3
                 - ssv4: Agilent SureSelect V4
                 - ssv5: Agilent SureSelect V5
                 - idt_xgen: IDT XGEN panel V1
                 - ssv5_idt_intersect: Interval intersect SSV5/IDT
    :param reference: genome reference. One of GRCh37 and GRCh38.
    :return: Interval HT
    """
    if name not in ('ssv2', 'ssv3', 'ssv4', 'ssv5', 'idt_xgen', 'ssv5_idt_intersect'):
        raise DataException('Invalid interval name...')
    if reference not in ('GRCh37', 'GRCh38'):
        raise DataException('Invalid genome reference...must be one of GRCh37 and GRCh38')

    return hl.read_table(f'{nfs_dir}/resources/intervals/{name}.intervals.{reference}.ht')


def get_1kg_mt(reference: str = 'GRCh38') -> hl.MatrixTable:
    """
    Return MT 1K genome phase 3 dataset.

    :param reference: genome reference. One of GRCh37 and GRCh38.
    :return: MatrixTable
    """
    return hl.read_matrix_table(f'{nfs_dir}/resources/1kgenome/phase3_1kg.snp_biallelic.{reference}.mt')


def get_mt_checkpoint_path(dataset: str = 'dataset', part: str = '') -> str:
    return f'{hdfs_checkpoint_dir}/{dataset}.{part}.checkpoint.mt'


def get_sample_meta_data() -> hl.Table:
    """
    Return HT with sample meta information (e.g. phenotypes, cohort, release permission, etc...)
    :return: Hail Table
    """
    return hl.import_table(
        f"{nfs_dir}/projects/wes_chd_ukbb/data/annotation/samples/sample.annotation.wes50k.final.u01032021.tsv",
        min_partitions=50,
        impute=True,
        key='ega_id'
    )


def get_chd_denovo_ht() -> hl.Table():
    """
    Return a list of de novo mutations called from CHD trios.
    Curated from two studies, Jin 2017 and Sifrim-Hitz 2016.

    :return: Hail Table
    """
    return hl.read_table(f'{nfs_dir}/resources/denovo/DNM_Jin2017_Sifrim2016_GRCh38_lift.ht')


##### gnomad resources #####

def get_gnomad_genomes_coverage_ht() -> hl.Table:
    return hl.import_table(
        f"{nfs_dir}/resources/gnomad/gnomad.genomes.r3.0.1.coverage.summary.tsv.bgz",
        min_partitions=1000,
        impute=True,
        key='locus'
    )


def get_gene_lof_metrics_ht() -> hl.Table:
    return hl.import_table(
        f"{nfs_dir}/resources/gnomad/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz",
        min_partitions=100,
        impute=True,
        key='gene'
    )


def get_transcript_lof_metrics_ht() -> hl.Table:
    return hl.import_table(
        f"{nfs_dir}/resources/gnomad/gnomad.v2.1.1.lof_metrics.by_transcript.txt.bgz",
        min_partitions=100,
        impute=True,
        key='transcript'
    )


class DataException(Exception):
    pass
