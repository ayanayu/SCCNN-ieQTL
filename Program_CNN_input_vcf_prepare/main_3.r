
#### The script is used to write h38 bed of variants

# define paths
vcf_directory='/media/ggj/My Passport1/Xiaoyanyu/CNN/vcf/OtherTissue/' 

LD <- read.table( paste0(vcf_directory,'plink.ld') ,header=T)
save(LD,file=paste0(vcf_directory,'LD.RData'))

LD1 <- LD[LD$R2>0.6, ]
LD_tmp <- data.frame(top_ieqtl=LD1$SNP_A , variant_id=LD1$SNP_B )

all_top_assoc <- read.csv(paste0(vcf_directory,'all_top_assoc.txt'),sep='\t')
all_top_assoc_tmp <- data.frame(gene=all_top_assoc$gene,top_ieqtl=all_top_assoc$variant_id)
all_top_assoc_tmp$index <- paste0(all_top_assoc_tmp$gene,'-',all_top_assoc_tmp$top_ieqtl )
all_top_assoc_tmp <- all_top_assoc_tmp[!duplicated(all_top_assoc_tmp$index),]

hb<-merge(LD_tmp,all_top_assoc_tmp,all=T,by="top_ieqtl")
hb <- na.omit(hb)

hb$index <- paste0(hb$variant_id,'-',hb$gene)
hb <- hb[!duplicated(hb$index), ]

aa<-matrix(unlist(strsplit(as.character(hb$variant_id),"_")),ncol = 5,byrow = T)
all_predict_vcf <- data.frame(chr=aa[,1],start=aa[,2],end=as.numeric(aa[,2])+1,ref=aa[,3],mut=aa[,4],index=hb$index)

options(scipen = 200)

all_predict_vcf_bed <- all_predict_vcf
all_predict_vcf_bed$ref <- NULL
all_predict_vcf_bed$mut <- NULL

write.table(all_predict_vcf_bed,file=paste0(vcf_directory,'unique_pos_bed.txt'),sep='\t',quote=F,row.names=F,col.names=F)

