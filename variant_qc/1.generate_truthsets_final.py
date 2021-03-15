# eam
# 2021-02-08

# Modified from @pavlos-pa10 
# Generate files required for RF part 1
import os
import hail as hl
import pandas as pd
import pyspark
import json
import sys
import re
from pathlib import Path
import logging
from typing import Any, Counter, List, Optional, Tuple, Union
from bokeh.plotting import output_file, save, show
from gnomad.resources.grch38 import gnomad
from gnomad.utils.annotations import annotate_adj
from gnomad.utils.annotations import unphase_call_expr, add_variant_type
from gnomad.utils.annotations import annotate_adj, bi_allelic_expr, bi_allelic_site_inbreeding_expr
from gnomad.utils.filtering import filter_to_autosomes
from gnomad.sample_qc.relatedness import (
    SIBLINGS,
    generate_sib_stats_expr,
    generate_trio_stats_expr,
)
from hail import Table

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

LABEL_COL = "rf_label"
TRAIN_COL = "rf_train"
PREDICTION_COL = "rf_prediction"

TRUTH_DATA = ["hapmap", "omni", "mills", "kgp_phase1_hc"]
INBREEDING_COEFF_HARD_CUTOFF = -0.3

################################
# Define the hail persistent storage directory
nfs_dir = 'file:///home/ubuntu/data'
hdfs_dir = 'hdfs://spark-master:9820/dir/hail_data'
# hdfs_checkpoint_dir = 'hdfs://spark-master:9820/checkpoint'

project_dir = 'file:///home/ubuntu/data/projects/wes_chd_ukbb'

######################################

omni = f'{nfs_dir}/resources/gnomad-training-sets/1000G_omni2.5.hg38.ht'
omni_ht = hl.read_table(omni)
mills = f'{nfs_dir}/resources/gnomad-training-sets/Mills_and_1000G_gold_standard.indels.hg38.ht'
mills_ht = hl.read_table(mills)
thousand_genomes = f'{nfs_dir}/resources/gnomad-training-sets/1000G_phase1.snps.high_confidence.hg38.ht'
thousand_genomes_ht = hl.read_table(thousand_genomes)
hapmap = f'{nfs_dir}/resources/gnomad-training-sets/hapmap_3.3.hg38.ht'
hapmap_ht = hl.read_table(hapmap)

######################################


def generate_trio_stats(
    mt: hl.MatrixTable, autosomes_only: bool = True, bi_allelic_only: bool = True
) -> hl.Table:
    """
    Default function to run `generate_trio_stats_expr` to get trio stats stratified by raw and adj
    .. note::
        Expects that `mt` is it a trio matrix table that was annotated with adj and if dealing with
        a sparse MT `hl.experimental.densify` must be run first.
        By default this pipeline function will filter `mt` to only autosomes and bi-allelic sites.
    :param mt: A Trio Matrix Table returned from `hl.trio_matrix`. Must be dense
    :param autosomes_only: If set, only autosomal intervals are used.
    :param bi_allelic_only: If set, only bi-allelic sites are used for the computation
    :return: Table with trio stats
    """
    if autosomes_only:
        mt = filter_to_autosomes(mt)
    if bi_allelic_only:
        mt = mt.filter_rows(bi_allelic_expr(mt))

    logger.info(f"Generating trio stats using {mt.count_cols()} trios.")
    trio_adj = mt.proband_entry.adj & mt.father_entry.adj & mt.mother_entry.adj

    ht = mt.select_rows(
        **generate_trio_stats_expr(
            mt,
            transmitted_strata={"raw": True, "adj": trio_adj},
            de_novo_strata={"raw": True, "adj": trio_adj},
            ac_strata={"raw": True, "adj": trio_adj},
        )
    ).rows()

    return ht


def generate_allele_data(mt: hl.MatrixTable) -> hl.Table:
    """
    Writes bi-allelic sites MT with the following annotations:
     - allele_data (nonsplit_alleles, has_star, variant_type, and n_alt_alleles)
    :param MatrixTable mt: Full unsplit MT
    :return: Table with allele data annotations
    :rtype: Table
    """
    ht = mt.rows().select()
    allele_data = hl.struct(nonsplit_alleles=ht.alleles,
                            has_star=hl.any(lambda a: a == '*', ht.alleles))
    ht = ht.annotate(allele_data=allele_data.annotate(
        **add_variant_type(ht.alleles)))

    ht = hl.split_multi_hts(ht)
    allele_type = (hl.case()
                   .when(hl.is_snp(ht.alleles[0], ht.alleles[1]), 'snv')
                   .when(hl.is_insertion(ht.alleles[0], ht.alleles[1]), 'ins')
                   .when(hl.is_deletion(ht.alleles[0], ht.alleles[1]), 'del')
                   .default('complex')
                   )
    ht = ht.annotate(allele_data=ht.allele_data.annotate(allele_type=allele_type,
                                                         was_mixed=ht.allele_data.variant_type == 'mixed'))
    return ht


def generate_ac(mt: hl.MatrixTable, fam_file: str) -> hl.Table:
    """
    Creates Table with QC samples, QC samples removing children and release samples raw and adj ACs.
    """
    # mt = mt.filter_cols(mt.meta.high_quality)
    fam_ht = hl.import_fam(fam_file, delimiter="\t")
    mt = mt.annotate_cols(unrelated_sample=hl.is_missing(fam_ht[mt.s]))
    mt = mt.filter_rows(hl.len(mt.alleles) > 1)
    mt = annotate_adj(mt)
    mt = mt.annotate_rows(
        ac_qc_samples_raw=hl.agg.sum(mt.GT.n_alt_alleles()),
        # ac_qc_samples_unrelated_raw=hl.agg.filter(~mt.meta.all_samples_related, hl.agg.sum(mt.GT.n_alt_alleles())),
        # ac_release_samples_raw=hl.agg.filter(mt.meta.release, hl.agg.sum(mt.GT.n_alt_alleles())),
        ac_qc_samples_adj=hl.agg.filter(
            mt.adj, hl.agg.sum(mt.GT.n_alt_alleles())),
        # ac_qc_samples_unrelated_adj=hl.agg.filter(~mt.meta.all_samples_related & mt.adj, hl.agg.sum(mt.GT.n_alt_alleles())),
        # ac_release_samples_adj=hl.agg.filter(mt.meta.release & mt.adj, hl.agg.sum(mt.GT.n_alt_alleles())),
    )
    return mt.rows()


def get_truth_ht() -> Table:
    """
    Returns a table with the following annotations from the latest version of the corresponding truth data:
    - hapmap
    - kgp_omni (1000 Genomes intersection Onni 2.5M array)
    - kgp_phase_1_hc (high confidence sites in 1000 genonmes)
    - mills (Mills & Devine indels)
    :return: A table with the latest version of popular truth data annotations
    """
    omni_ht = hl.read_table(omni)
    mills_ht = hl.read_table(mills)
    thousand_genomes_ht = hl.read_table(thousand_genomes)
    hapmap_ht = hl.read_table(hapmap)
    return (
        hapmap_ht
        .select(hapmap=True)
        .join(omni_ht.select(omni=True), how="outer")
        .join(thousand_genomes_ht.select(kgp_phase1_hc=True), how="outer")
        .join(mills_ht.select(mills=True), how="outer")
        .repartition(200, shuffle=False)
        .persist()
    )


if __name__ == "__main__":
    # need to create spark cluster first before intiialising hail
    # sc = pyspark.SparkContext()
    
    hl.stop()
    hl.init(default_reference="GRCh38")
    # s3 credentials required for user to access the datasets in farm flexible compute s3 environment
    # you may use your own here from your .s3fg file in your home directory

    group = "raw"

    mt = hl.read_matrix_table(
        f'{nfs_dir}/hail_data/mts/chd_ukbb_split_v2_09092020.mt')

    # Truthset
    truthset_ht = get_truth_ht()
    truthset_ht.write(f'{nfs_dir}/hail_data/hts/truthset.ht', overwrite=True)
    truthset_ht = hl.read_table(f'{nfs_dir}/hail_data/hts/truthset.ht')

    # Trio data
    # trio annotation:
    mt_adj = annotate_adj(mt)
    fam = f"{project_dir}/data/annotation/samples/sample.complete_trios.wes50k.02022021.noheader.fam"
    pedigree = hl.Pedigree.read(fam)
    trio_dataset = hl.trio_matrix(mt_adj, pedigree, complete_trios=True)
    trio_dataset.checkpoint(
        f'{hdfs_dir}/chd_ukbb.trios.adj.mt', overwrite=True)
    trio_stats_ht = generate_trio_stats(
        trio_dataset, autosomes_only=True, bi_allelic_only=True)
    trio_stats_ht.write(
        f'{hdfs_dir}/chd_ukbb.trios.stats.ht', overwrite=True)

    # inbreeding ht
    mt_inbreeding = mt.annotate_rows(
        InbreedingCoeff=bi_allelic_site_inbreeding_expr(mt.GT))

    mt = mt.key_rows_by('locus').distinct_by_row(
    ).key_rows_by('locus', 'alleles')

    ht_inbreeding = mt_inbreeding.rows()

    # allele data and qc_ac ht
    allele_data_ht = generate_allele_data(mt)

    qc_ac_ht = generate_ac(mt, fam)

    ht_inbreeding.write(
        f'{hdfs_dir}/chd_ukbb.inbreeding.ht', overwrite=True)
    qc_ac_ht.write(
        f'{hdfs_dir}/chd_ukbb.qc_ac.ht', overwrite=True)
    allele_data_ht.write(
        f'{hdfs_dir}/chd_ukbb.allele_data.ht', overwrite=True)
