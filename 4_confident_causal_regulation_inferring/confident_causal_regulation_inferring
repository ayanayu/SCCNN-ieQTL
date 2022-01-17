
#### The script is used for (1)preprocessing CNN prediction results 
####                        (2)adding CNN prediction information to ieQTL results for each tissue
####                        (3)inferring confident causal genetic regulation

# define paths
vcf_directory='/media/ggj/My Passport1/Xiaoyanyu/CNN/vcf/OtherTissue/' 
CNN_single_cell_model_result_directory = '/media/ggj/My Passport1/Xiaoyanyu/CNN/predict_result/OtherTissue/' # predictions of CNN (single-cell based)
output_directory='/media/ggj/My Passport/YanyuX/'
CNN_200kb_result = '/home/ggj/QTL/Script_collating/resource/CNN_200kb_result.txt'
ieqtl_output_directory='/media/ggj/Guo-4T-MS2/XiaoYanyu/qtl_result/' # ieqtl output directory


#
library(dplyr)
library(data.table)
library(readr)

library(ggplot2)
theme_set(theme_bw())
library(dplyr)
library(mgcv)
library(tidymv)
library(splines)

file_list <- list.files(path = ieqtl_output_directory)
file_list

allcelltype_top_assoc_file <- data.frame()
for (f in file_list){
  if (file.exists(paste0(ieqtl_output_directory,f,'/GTEx_V8_samples.cis_qtl_top_assoc.txt.gz'))==TRUE){
   top_assoc <- read.csv(paste0(f,'/GTEx_V8_samples.cis_qtl_top_assoc.txt.gz'),sep='\t' )
   allcelltype_top_assoc_file <- rbind(allcelltype_top_assoc_file,top_assoc)
 }
}
head(allcelltype_top_assoc_file)

##(1)preprocessing CNN prediction results 

match <- fread( paste0(vcf_directory,'h19/unique_pos.txt'),stringsAsFactors = F,sep="\t")
match <-as.data.frame(match)
match$index <- paste0(match$V1,'_',match$V2,'_',match$V4,'_',match$V5 )
match$V7 <- gsub('-.*','',match$V6)

#####single model

file_list <- list.files(path = CNN_single_cell_model_result_directory)

setwd(CNN_single_cell_model_result_directory)
single <- data.frame()
for (i in file_list){
   single_tmp <- fread(i,sep=',',stringsAsFactors = F)
   single_tmp <- as.data.frame(single_tmp)
   single <- rbind(single,single_tmp)
  }

single$variant_ID <- paste0(single$'0','_',single$'1','_',single$'3','_',single$'4')
single$index <- paste0(single$variant_ID,'_',single$gene)
single <- single[!duplicated(single$index),]
single <- single[nchar(as.character(single$'3'))%in%1&nchar(as.character(single$'4'))%in%1, ]

single <- single[,c(8:ncol(single))]
single$strand <- NULL

match1 <- match[match$index%in%as.character(single$variant_ID), ]
match1 <- match1[!duplicated(match1$index), ]
rownames(match1) <- match1$index

single$variant_ID_h38 <- match1[as.character(single$variant_ID),]$V7

write.table(single,file='single.txt',quote=F,row.names=F,col.names=T,sep='\t' )



##(2)adding CNN prediction information to ieQTL results for each tissue

all_top_assoc_bei <- fread(allcelltype_top_assoc_file,sep='\t',stringsAsFactors = F)
all_top_assoc_bei <- as.data.frame(all_top_assoc_bei)

single_alltissue <- fread(paste0(CNN_single_cell_model_result_directory,'single.txt'),sep='\t',stringsAsFactors = F)
single_alltissue <- as.data.frame(single_alltissue)

load(paste0(vcf_directory,'LD.RData'))

tissue_list <- unique(gsub('\\d','',colnames(single_alltissue)[grep(colnames(single_alltissue),pattern='Adult')]))
for (tissue in tissue_list){

 celltype_list <- colnames(single_alltissue)[grep(colnames(single_alltissue),pattern=tissue)]

 ##predict_result
 result_single_all <- data.frame()
 for (i in celltype_list){
   single_tmp <- data.frame(single_alltissue[,c('variant_ID_h38','gene','dist',i)])
   colnames(single_tmp) <- c('variant_ID_h38','gene','dist','single')
   single_tmp$celltype <- rep(i,nrow(single_tmp))
   result_single_all <- rbind(result_single_all,single_tmp)
 }

 #####ieQTL result
 all_top_assoc <- all_top_assoc_bei[all_top_assoc_bei$celltype%in%celltype_list, ]

 all_top_assoc$gene <- gsub('\\..*','',all_top_assoc$phenotype_id)
 all_top_assoc$index <- paste0(all_top_assoc$celltype,'_',all_top_assoc$gene)
 rownames(all_top_assoc) <- all_top_assoc$index

 all_top_assoc_20kb <- all_top_assoc[abs(all_top_assoc$tss_distance)<20000,]

 result_single_all$index <- paste0(result_single_all$celltype,'_',result_single_all$gene )

 ###
 result_single_all <- result_single_all[result_single_all$index%in%rownames(all_top_assoc_20kb), ]

 result_single_all$top_ieqtl <- all_top_assoc_20kb[as.character(result_single_all$index), ]$variant_id
 result_single_all$top_ieqtl_pval_gi_emt <- all_top_assoc_20kb[as.character(result_single_all$index), ]$pval_emt
 result_single_all$top_ieqtl_pval_gi <- all_top_assoc_20kb[as.character(result_single_all$index), ]$pval_gi
 result_single_all$top_ieqtl_pval_g <- all_top_assoc_20kb[as.character(result_single_all$index), ]$pval_g
 result_single_all$maf <- all_top_assoc_20kb[as.character(result_single_all$index), ]$maf
 result_single_all$b_g <- all_top_assoc_20kb[as.character(result_single_all$index), ]$b_g
 result_single_all$b_gi <- all_top_assoc_20kb[as.character(result_single_all$index), ]$b_gi

 result_single_all$tests_emt <- all_top_assoc_20kb[as.character(result_single_all$index), ]$tests_emt
 result_single_all$top_ieqtl_pval_g_emt <- result_single_all$top_ieqtl_pval_g * result_single_all$tests_emt
 result_single_all$top_ieqtl_pval_g_emt[result_single_all$top_ieqtl_pval_g_emt>1] <- 1

 result_single_all <- na.omit(result_single_all)

 result_single_all$index1 <- paste0(result_single_all$top_ieqtl,'_',result_single_all$variant_ID_h38 )

 LD1 <- LD[LD$SNP_A%in%result_single_all$top_ieqtl&LD$SNP_B%in%result_single_all$variant_ID_h38, ]
 LD1$index1 <- paste0(LD1$SNP_A,'_',LD1$SNP_B)
 LD1 <- LD1[LD1$index1%in%as.character(result_single_all$index1),]

 result_single_all<-merge(result_single_all,LD1,all=T,by="index1")
 result_single_all <- na.omit(result_single_all)

 result_single_all <- result_single_all[result_single_all$R2 > 0.6, ]

 result_single_all_bei <- result_single_all
 save(result_single_all_bei,file=paste0(output_directory,tissue,'_result_single_all.RData'))
 }



##(3)inferring confident causal genetic regulation


### estimate the 95% prediction interval bounds

ss <- fread(CNN_200kb_result,sep='\t',stringsAsFactors = F)##
ss  <- as.data.frame(ss)

ss$dist <- abs(ss$dist)

model <- lm(log(abs(single)) ~ dist, data = ss)
pred.int <- predict(model, interval = "prediction")
mydata <- cbind(ss, pred.int)


###infer putative regulatory variant-gene pairs

all_top_assoc_bei$pval_emt_g <- all_top_assoc_bei$pval_g*all_top_assoc_bei$tests_emt
all_top_assoc_bei$pval_emt_g[all_top_assoc_bei$pval_emt_g >1] <-1

LD <- LD[LD$R2>0.6, ]

##
file_list <- list.files(path = output_directory , pattern = 'RData')
file_list

all_celltype_cnn_identified <- data.frame()
for  (i in file_list){   
   tissue <- gsub('_result_single_all.RData','',i)
   load( paste0(output_directory,i) )
   result_single_all <- result_single_all_bei

   #####significant in silico mutagenesis
   fit_upr <- exp(predict(model,data.frame(dist=abs(result_single_all$dist)),interval = "prediction"))
   result_single_all$fit_upr <- fit_upr[,3]
   result_single_all_cnn_identified <- result_single_all[abs(result_single_all$single) > result_single_all$fit_upr, ]

   ##get top_ieqtl in LD
   all_top_assoc_tissue <- all_top_assoc_bei[gsub('\\d','',all_top_assoc_bei$celltype)%in%tissue, ]
   LD1 <- LD[LD$SNP_A%in%as.character(all_top_assoc_tissue$variant_id),]
   LD2 <- LD1[LD1$SNP_B%in%as.character(result_single_all_cnn_identified$variant_ID_h38),c('SNP_A','SNP_B') ];
   colnames(LD2) <- c('top_ieqtl_redifined','variant_ID_h38')
   result_single_all_cnn_identified <- merge(result_single_all_cnn_identified,LD2,all=T,by="variant_ID_h38")
   result_single_all_cnn_identified <- na.omit(result_single_all_cnn_identified)
   
   ##significance of top_ieqtl
   result_single_all_cnn_identified$index3 <- paste0(result_single_all_cnn_identified$top_ieqtl_redifined,'-',result_single_all_cnn_identified$gene )
   all_top_assoc_tissue$index3 <- paste0(all_top_assoc_tissue$variant_id,'-',gsub('\\..*','',all_top_assoc_tissue$phenotype_id) )
   all_top_assoc_tissue <- all_top_assoc_tissue[order(all_top_assoc_tissue$index3,all_top_assoc_tissue$pval_emt),]
   all_top_assoc_tissue <- all_top_assoc_tissue[!duplicated(all_top_assoc_tissue$index3),]
   rownames(all_top_assoc_tissue) <- as.character(all_top_assoc_tissue$index3)
   result_single_all_cnn_identified$top_ieqtl_redifined_min_p <- all_top_assoc_tissue[as.character(result_single_all_cnn_identified$index3),]$pval_emt
   result_single_all_cnn_identified <- na.omit(result_single_all_cnn_identified)
   
   ##pval_emt < 0.05 reserved
   result_single_all_cnn_identified <- result_single_all_cnn_identified[result_single_all_cnn_identified$top_ieqtl_redifined_min_p < 0.05, ]
   result_single_all_cnn_identified$index4 <- paste0(result_single_all_cnn_identified$variant_ID_h38,'-',
                                             result_single_all_cnn_identified$gene,'-',
											 result_single_all_cnn_identified$celltype)
   result_single_all_cnn_identified <- result_single_all_cnn_identified[order(result_single_all_cnn_identified$index4, result_single_all_cnn_identified$top_ieqtl_redifined_min_p),]
   result_single_all_cnn_identified <- result_single_all_cnn_identified[!duplicated(result_single_all_cnn_identified$index4),]
   #print(table(result_single_all_cnn_identified$celltype))
   all_celltype_cnn_identified <- rbind(all_celltype_cnn_identified,result_single_all_cnn_identified)
   print(i)
 } 

table(all_celltype_cnn_identified$celltype)

write.table(all_celltype_cnn_identified,sep='\t',
            file=paste0(output_directory,'All_putative_causal_variant_gene_pairs.txt'),quote=F,
			row.names=F,col.names=T) ##putative causal variant-gene pairs
























