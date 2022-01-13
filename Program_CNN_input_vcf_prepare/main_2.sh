#!/bin/bash

#### The script is used to get all variants in LD to top ieQTL 

# define paths
vcf_directory=/media/ggj/My\ Passport1/Xiaoyanyu/CNN/vcf/OtherTissue/ 
bfile_directory=/tmpfile/79471/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v8.p2.c1.GRU/GenotypeFiles/phg001219.v1.GTEx_v8_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze

cd $vcf_directory
/home/ggj/plink --bfile ${bfile_directory}\
      --r2\
      --ld-snp-list topsnp.txt\
      --ld-window-kb 1000\
      --ld-window 99999\
      --ld-window-r2 0.6
