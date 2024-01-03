# wes_chd_ukbb


## Meta-analysis of a large-scale wes dataset using Hail

This repository contains a series of pipelines (mostly command-line tools) to analyse the 
CHD case-control exome cohort linked to the manuscript: 
*"Assessing the contribution of rare variants to congenital heart disease through 
a large-scale case-control exome study"* ([medRxiv link](https://doi.org/10.1101/2023.12.23.23300495)) .


Most pipelines were adapted from the gnomad repository [1]. Pipelines were built using the python-like library
Hail (https://hail.is).


*Exome dataset*: Congenital Heart Disease (CHD) cases were mainly sequenced as part of an initiative from the German 
Competence Network for Congenital Heart Defects, the Deciphering Developmental Disorder (DDD) project and the 
University of Nottingham (UK); controls were sequenced as part of the UK Biobank (UKBB). 


### Sample QC
1. *Hard filters*: Mark samples with unspecific chromosomal sex, low call rate and/or low coverage.
    (`sample_qc/apply_hard_filters.py`)
    
2. *Population ancestries inferring*: Impute sample ancestries using the the 1000 Genomes Phase 3 sequence dataset. 
    (`sample_qc/ancestry_inference.py`)
    
3. *Inferring sample relatedness*: Identify twins/duplicated samples as well as first- and second-relatives. 
    (`sample_qc/relatedness_inference.py`)
    
4. *Platform inference*: Assign capture platform to samples using unsupervised clustering. 
    (`sample_qc/platform_pca.py`)
    
5. *Platform- and population-specific outliers filtering*: Detect sample outliers stratified by population/platforms. 
    (`sample_qc/sample_qc.py`)
    
6. *Final sample QC*: Mark samples failing QC and generate HailMatrix with high-quality genotypes (dept of coverage >= 10, 
    genotype quality >= 20 and genotype allele balance of heterozygotes > 0.20).
    (`sample_qc/finalise_sample_qc.py`)

### Variant QC
1. *Hard filters*: Mark variants failing hard filters if it a) showed an excess of heterozygotes
    (inbreeding coefficient < -0.3) and b) contains absence of at least one sample with a high-quality genotype
                      
2. *RF model*: Application of a random forest (RF) model to distinguish true variations from potential false positives.
    (`variant_qc/3.train_apply_finalise_RF.py`)
 
3. *VQSR filter*: Application of the GATK Variant Quality Score Recalibration (VQSR) tool.

4. *Coverage*: Mark variants as covered if it a) is defined in the major capture platforms intervals used in the 
               assembled cohort and b) showed a coverage of 10X or more in at least the 90% of the samples in 
               the gnomAD genome dataset (version 3.1.0).
               
5. *Final variant QC*: Mark variants failing hard, random forest, VQSR and/or coverage filters.
    (`variant_qc/finalise_variant_qc.py`)                

               
### Burden testing
*Gene burden test*: Run gene-based case-control burden test (Fisher Exact) stratified by variant functional category 
and proband syndromic status. (`pipelines/gene_burden_fet.py`)

*Gene-set burden test*: Run gene set-based case-control burden test (logistic regression) stratified by variant functional category 
and proband syndromic status. (`pipelines/geneset_burden_logreg.py`)


### Util scripts
*VEP parser*: Parse a VCF file annotated with the Variant Effect Predictor (VEP) tool
 (`pipelines/vep_parser.py`)

*dbNSFP parser*: Generate a HailTable from the dbNSFP (version 4.1a) for annotations and downstream analysis
 (`script-utils/parse_dbnsfp_variants.py`)


### References
[1] Karczewski, KJ et al. The mutational constraint spectrum quantified from variation in 141,456 humans.
    Nature 581, 434–443 (2020).
         


