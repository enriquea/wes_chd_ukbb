"""
Set of function to retrieve paths from jointly called data as MatrixTable / HailTable
"""

import hail as hl

nfs_dir = 'file:///home/ubuntu/data'
hdfs_dir = 'hdfs://spark-master:9820/dir'
hdfs_tmp = 'hdfs://spark-master:9820/tmp'


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
    if part == 'unfiltered':
        return hl.read_matrix_table(f'{nfs_dir}/hail_data/mts/{dataset}_split_v2_09092020.mt')
    elif part == 'filtered_high_callrate':
        return hl.read_matrix_table(f'{nfs_dir}/hail_data/mts/{dataset}.high_callrate.common_snp.biallelic.mt')
    else:
        raise DataException("Select one valid version of the dataset: unfiltered or filtered_high_callrate...")


def get_qc_mt_path(dataset: str = 'chd_ukbb', part: str = None, split=False) -> str:
    split = '.split' if split else ''
    return f'{nfs_dir}/hail_data/mts/{dataset}.qc.{part}{split}.mt'


def get_sample_qc_ht_path(dataset: str = 'chd_ukbb', part: str = None) -> str:
    return f'{nfs_dir}/hail_data/sample_qc/{dataset}.sample_qc.{part}.ht'


def get_interval_ht(name: str = None, reference: str = None) -> hl.Table:
    """
    Return interval HT. Used mostly for filtering...

    :param name: name/description of the interval to be returned.
                 - ssv5: Agilent SureSelect V5
                 - idt_xgen: IDT XGEN panel V1
                 - ssv5_idt_intersect: Interval intersect SSV5/IDT
    :param reference: genome reference. One the One of GRCh37 and GRCh38.
    :return: Interval HT
    """
    if name or reference is None:
        raise DataException('Both name and reference must be specified...')

    return hl.read_table(f'{nfs_dir}/resources/intervals/{name}.intervals.{reference}.ht')


class DataException(Exception):
    pass
