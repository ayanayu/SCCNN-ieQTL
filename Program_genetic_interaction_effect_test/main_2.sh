#!/bin/bash

#### The script is used to (1)estimate PEER factors (2)estimate genotype pcs (3)combine all covariates

# define paths
genotype_pc_script=/home/ggj/QTL/nohup_script/genotype_pc.r # genotype_pc.r scirpt pathway
run_PEER_script=/home/ggj/QTL/script/run_PEER.R # run_PEER.R scirpt pathway
combine_covariates_script=/home/ggj/QTL/script/combine_covariates.py # combine_covariates.py scirpt pathway
PEER_directory=/home/ggj/QTL/PEER/ # PEER directory
bfile_directory=/tmpfile/79471/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v8.p2.c1.GRU/GenotypeFiles/phg001219.v1.GTEx_v8_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze ## genotype directory
ieqtl_output_directory=/media/ggj/Guo-4T-MS2/XiaoYanyu/qtl_result # ieqtl output directory
all_celltype_list='AdultArtery3 AdultArtery4 AdultArtery5 AdultArtery6 AdultArtery7' ## all cell types used for following ieQTL test (have been used for deconvolution) 


##
cd $PEER_directory
tissuename=$(ls)
#echo $tissuename

for i in $tissuename;do
echo $i start;
echo $i PEER run;
cd $PEER_directory;
cd $i;
/usr/local/bin/Rscript $run_PEER_script choosesample_exp.txt GTEx_V8 60 --max_iter 200;
echo $i PEER finish;
echo $i get genotype pcs start;
cd $PEER_directory;
cd $i;
testout=./;
/home/ggj/plink2 --allow-extra-chr --threads 20 -bfile ${bfile_directory} --pca 5 --out ${testout} --keep genotype_pc_analyse_sample.txt
/usr/local/bin/Rscript $genotype_pc_script --tissue $i;
echo $i get genotype pcs finish;
echo $i combine covariates start;
source /home/ggj/anaconda2/bin/activate py3.7;
prefix=GTEx_V8;
genotype_pcs=genotype_pcs.txt;
cd $PEER_directory;
cd $i;
python $combine_covariates_script ${prefix}.PEER_covariates.txt ${prefix} \
    --genotype_pcs ${genotype_pcs}
echo $i all finish	;
done


for i in $all_celltype_list;do
cd $ieqtl_output_directory;
mkdir $i;
done





