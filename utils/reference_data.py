# eam
# 2021-05-01

import hail as hl


nfs_dir = 'file:///home/ubuntu/data'


def get_lcr_ht(overwrite: bool = False) -> hl.Table:
    lcr_interval = hl.import_locus_intervals(
        f'{nfs_dir}/resources/grch38/LCRFromHengHg38.txt.bgz',
        skip_invalid_intervals=True,
        min_partitions=50,
        reference_genome='GRCh38'
    )
    return lcr_interval.checkpoint(
        f'{nfs_dir}/resources/grch38/LCRFromHengHg38.ht',
        overwrite=overwrite,
        _read_if_exists=not overwrite
    )


def get_telomeres_and_centromeres_ht(overwrite: bool = False) -> hl.Table:
    tc_interval = hl.import_bed(
        f'{nfs_dir}/resources/grch38/hg38.telomeresAndMergedCentromeres.bed',
        skip_invalid_intervals=True,
        min_partitions=10,
        reference_genome='GRCh38'
    )
    return tc_interval.checkpoint(
        f'{nfs_dir}/resources/grch38/hg38.telomeresAndMergedCentromeres.ht',
        overwrite=overwrite,
        _read_if_exists=not overwrite
    )


def get_segdups_ht(overwrite: bool = False) -> hl.Table:
    segdup_interval = hl.import_bed(
        f'{nfs_dir}/resources/grch38/GRCh38_segdups.bed',
        skip_invalid_intervals=True,
        min_partitions=50,
        reference_genome='GRCh38'
    )
    return segdup_interval.checkpoint(
        f'{nfs_dir}/resources/grch38/GRCh38_segdups.ht',
        overwrite=overwrite,
        _read_if_exists=not overwrite
    )


def import_cds_from_gtf(overwrite: bool = False) -> hl.Table:
    """
    Creates a HT with a row for each base / gene that is in a CDS in gencode v29
    :param bool overwrite: Overwrite existing table
    :return: HT
    :rtype: Table
    """
    gtf = hl.experimental.import_gtf(
        f'{nfs_dir}/resources/hail-common/references/gencode/gencode.v29.annotation.gtf.bgz',
        reference_genome='GRCh38',
        skip_invalid_contigs=True,
        min_partitions=200
    )
    gtf = gtf.filter((gtf.feature == 'CDS') & (gtf.transcript_type == 'protein_coding') & (gtf.tag == 'basic'))
    gtf = gtf.annotate(locus=hl.range(gtf.interval.start.position,
                                      gtf.interval.end.position).map(lambda x:
                                                                     hl.locus(gtf.interval.start.contig, x, 'GRCh38')))
    gtf = gtf.key_by().select('gene_id', 'locus', 'gene_name').explode('locus')
    gtf = gtf.key_by('locus', 'gene_id').distinct()
    return gtf.checkpoint(
        f'{nfs_dir}/resources/grch38/gencode_grch38.gene_by_base.cds.ht',
        overwrite=overwrite,
        _read_if_exists=not overwrite
    )


def import_cds_intervals_from_gtf(overwrite: bool = False) -> hl.Table:
    """
    Creates a HT with a interval row for each CDS in gencode v29
    :param bool overwrite: Overwrite existing table
    :return: HT
    :rtype: Table
    """
    gtf = hl.experimental.import_gtf(
        f'{nfs_dir}/resources/hail-common/references/gencode/gencode.v29.annotation.gtf.bgz',
        reference_genome='GRCh38',
        skip_invalid_contigs=True,
        min_partitions=200
    )
    gtf = gtf.filter((gtf.feature == 'CDS') & (gtf.transcript_type == 'protein_coding') & (gtf.tag == 'basic'))

    return gtf.checkpoint(
        f'{nfs_dir}/resources/grch38/gencode.cds_intervals.GRCh38.ht',
        overwrite=overwrite,
        _read_if_exists=not overwrite
    )
