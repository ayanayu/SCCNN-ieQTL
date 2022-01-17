# SCCNN-ieQTL

SCCNN-ieQTL
SCCNN-ieQTL is an analytical framework for decoding cell type-dependent genetic regulation. It integrates single-cell-based genetic interaction QTLs and CNN predictions to confidently infer genetic regulatory effects as well as the functional variants and effector cell types.

Input formats
Four inputs are required for the analyses with SCCNN-ieQTL: genotype and phenotype data of bulk tissue samples, cell type expression profiles and cell type composition files. The pipeline is originally developed based on bulk tissue samples in GTEx and single-cell transcriptome data in HCL. You can download all preprocessed files except genotype data (which are protected and available via dbGaP) under the folder named ‘data’. However, you can also apply the pipeline to any other dataset consistent with the following requirements:
1)	Genotype data must be in PLINK format generated from a VCF using plink2 --make-bed.
2)	Phenotype data must be normalized between samples using TMM (Robinson & Oshlack, Genome Biology, 2010) and provided in BED format with a single header line starting with # and first four columns corresponding to chr, start, end, phenotype_id. You can use the function eqtl_prepare_expression.py available in a GTEx eQTL discovery pipeline to generate the BED template from a gene annotation in GTF format.
3)	Cell type expression profiles should be generated from single-cell transcriptome data by aggregating the average data from multiple cells in the same cell cluster. We suggest performing the Markov affinity-based graph imputation of cells (MAGIC) algorithm for imputation of the normalized single-cell expression matrix before cellular aggregation.
4)	Cell type compositions across all the bulk tissue samples should be estimated with signature matrices from the same single-cell transcriptome dataset, which can be implemented with an in silico deconvolution tool CIBERSORTx.

Requirements
Required Tools
plink
liftOver
Required R packages
peer
argparser
dplyr
data.table
readr
ggplot2
mgcv
tidymv
splines
Required Python packages
h5py>=2.7.0
numpy>=1.14.2
pandas==0.22.0
scipy>=0.19.1
six>=1.11.0
xgboost==0.7.post4
pyfasta>=0.5.2
torch>=1.0.0
argparse
tensorQTL

Running SCCNN-ieQTL
