"""
Set of functions to handle import/export of sample metadata, MT, intervals and resources.

"""

import hail as hl

nfs_dir = 'file:///home/ubuntu/data'
hdfs_dir = 'hdfs://spark-master:9820/dir/hail_data'
hdfs_checkpoint_dir = 'hdfs://spark-master:9820/checkpoint'


##### MatrixTable utils #####

def get_mt_data(dataset: str = 'chd_ukbb', part: str = None, split: bool = True) -> hl.MatrixTable:
    """
    Return MatrixTable for different version of CHD/UKBB cohort. Multi-allelic sites were split.

    :param split: Multi-allelic variants were split
    :param dataset: dataset name (Refers just to one ds at the moment, but should handle multiple ds)
    :param part: One of:
                raw: unfiltered MatrixTable with split multi-allelic sites.
                filtered_high_callrate: filtered MT for common SNPs (MAF 0.1%), bi-allelic, callrate > 0.99
    :return MT

    """
    # TODO: add different MT versions here...
    parts = ['raw', 'raw_chr20', 'filtered_high_callrate']
    split = '.split' if split else ''

    if part == 'raw':
        return hl.read_matrix_table(f'{nfs_dir}/hail_data/mts/{dataset}.{part}{split}.mt')
    elif part == 'raw_chr20':
        return hl.read_matrix_table(f'{nfs_dir}/hail_data/mts/{dataset}.{part}{split}.mt')
    elif part == 'filtered_high_callrate':
        return hl.read_matrix_table(f'{nfs_dir}/hail_data/mts/{dataset}.high_callrate.common_snp.biallelic.mt')
    else:
        raise DataException(f"Select one valid version of the dataset: {parts}...")


def get_qc_mt_path(dataset: str = 'chd_ukbb', part: str = None, split=True, ld_pruned=False) -> str:
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


def get_mt_checkpoint_path(dataset: str = 'dataset', part: str = '') -> str:
    return f'{hdfs_checkpoint_dir}/{dataset}.{part}.checkpoint.mt'


def get_1kg_mt(reference: str = 'GRCh38') -> hl.MatrixTable:
    """
    Return MT 1K genome phase 3 dataset.

    :param reference: genome reference. One of GRCh37 and GRCh38.
    :return: MatrixTable
    """
    return hl.read_matrix_table(f'{nfs_dir}/resources/1kgenome/phase3_1kg.snp_biallelic.{reference}.mt')


##### Intervals utils #####

def get_capture_interval_ht(name: str, reference: str) -> hl.Table:
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


def get_chd_denovo_ht() -> hl.Table:
    """
    Return a list of de novo mutations called from CHD trios.
    Curated from two studies, Jin 2017 and Sifrim-Hitz 2016.

    :return: Hail Table
    """
    return hl.read_table(f'{nfs_dir}/resources/denovo/DNM_Jin2017_Sifrim2016_GRCh38_lift.ht')


##### sample-metadata annotations #####

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


def get_fam_path() -> str:
    """
    Return fam file path
    :return: Path to fam file
    """
    return (
        f"{nfs_dir}/projects/wes_chd_ukbb/data/annotation/samples/sample.complete_trios.wes50k.u01032021.noheader.fam"
    )


def import_fam_ht() -> hl.Table:
    """
    Return fam table (trio information)
    :return: Hail Table
    """
    return hl.import_fam(
        get_fam_path()
    )


##### sample-qc #####

def get_sample_pop_qc() -> hl.Table:
    return hl.import_table(
        f"{nfs_dir}/projects/wes_chd_ukbb/data/annotation/samples/chd_ukbb_population_predicted_pca_rf_09102020.txt",
        min_partitions=50,
        no_header=False,
        impute=True,
        quote='"',
    ).key_by('s')


def get_sample_qc_ht_path(dataset: str = 'chd_ukbb', part: str = None) -> str:
    qc_parts = ['sex_chrom_coverage',
                'hard_filters',
                'joint_pca_1kg',
                'platform_pca',
                'population_qc',
                'high_conf_autosomes',
                'stratified_metrics_filter',
                'final_qc']

    if part not in qc_parts:
        raise DataException(f'Expected part one of: {qc_parts}')

    return f'{nfs_dir}/hail_data/sample_qc/{dataset}.sample_qc.{part}.ht'


##### variant-qc #####

def get_variant_qc_ht_path(dataset: str = 'chd_ukbb', part: str = None, split=True) -> str:
    qc_parts = ['vep_vqsr',
                'hard_filters',
                'rf_result',
                'coverage_stats',
                'final_qc']

    if part not in qc_parts:
        raise DataException(f'Expected part one of: {qc_parts}')

    split = '.split' if split else ''
    return f'{nfs_dir}/hail_data/variant_qc/{dataset}.variant_qc.{part}{split}.ht'


def get_variant_af_pops_ht_path(dataset: str = 'chd_ukbb',
                                split=True) -> str:
    split = '.split' if split else ''
    return f'{nfs_dir}/hail_data/variant_qc/{dataset}.af_pops{split}.ht'


def get_vep_vqsr_vcf_path() -> str:
    return f'{nfs_dir}/chd_ukbb_vep/recal_snp_recal_indelgenome.sorted.vcf.gz'


def get_vep_annotation_ht() -> hl.Table:
    return hl.read_table(
        get_variant_qc_ht_path(part='vep_vqsr')
    )


#### pathogenic prediction scores ####

def get_vep_scores_ht() -> hl.Table:
    """
    Return HT with pathogenic prediction scores annotated with VEP (e.g. CADD)

    :return: HailTable
    """
    return hl.read_table(
        f'{nfs_dir}/hail_data/scores/chd_ukbb.pathogenic_scores.vep.split.ht'
    ).key_by('locus', 'alleles')


def get_ccr_ht() -> hl.Table:
    """
    Return HT with coding-constraint regions (CCRs) and percentile lifted over hg38

    :return: HailTable
    """
    return hl.read_table(
        f"{nfs_dir}/resources/ccr/ccr_hg38.ht"
    )


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


def get_gnomad_genomes_v3_af_ht() -> hl.Table:
    return hl.read_table(
        f"{nfs_dir}/resources/gnomad/gnomad_v3/gnomad_3.0_sites_AF.ht",
    )


##### af annotation tables #####

def get_bonn_af_ht() -> hl.Table:
    """
    In-hause german allelic frequencies from Bonn (May 2021)
    Liftover from hg37 -> hg38

    :return: HailTable
    """
    return hl.read_table(
        f'{nfs_dir}/resources/annotation/Cohort_Bonn_AF_0521_lift.ht',
    )


def get_germ_af_ht() -> hl.Table:
    """
    In-hause german allelic frequencies from Tuebingen
    Liftover from hg37 -> hg38

    :return: HailTable
    """
    return hl.read_table(
        f'{nfs_dir}/resources/annotation/german_af_hg38_lift.ht',
    )


def get_rum_af_ht() -> hl.Table:
    """
    Allelic frequencies from RUMC (Dutch) cohort.
    Liftover from hg37 -> hg38
    DOI: 10.1038/s41586-020-2832-5

    :return: HailTable
    """
    return hl.read_table(
        f'{nfs_dir}/resources/annotation/rumc_af_hg38_lift.ht',
    )


def get_af_annotation_ht() -> hl.Table:
    """
    CHD/UKBB variants (multi-allelic split) HT with annotated allelic frequencies
    from external sources (e.g. gnomad exomes and genomes). Genome build hg38.

    :return: HailTable
    """
    return hl.read_table(
        f'{nfs_dir}/hail_data/hts/chd_ukbb.variants.af.annotations.external.08092021.ht'
    )


class DataException(Exception):
    pass
