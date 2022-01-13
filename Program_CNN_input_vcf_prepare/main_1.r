
#### The script is used to get all top ieQTL (after genetic interaction effect test)

# define paths
ieqtl_output_directory='/media/ggj/Guo-4T-MS2/XiaoYanyu/qtl_result/' # ieqtl output directory
vcf_directory='/media/ggj/My Passport1/Xiaoyanyu/CNN/vcf/OtherTissue/' # vcf output directory

file_list <- list.files(path = ieqtl_output_directory, pattern = 'Adult')
all_top_assoc <- data.frame() 
for (i in file_list){
   if (file.exists(paste0(ieqtl_output_directory,i,"/GTEx_V8_samples.cis_qtl_top_assoc.txt.gz"))==TRUE){
     top_assoc <- read.csv(paste0(ieqtl_output_directory,i,'/GTEx_V8_samples.cis_qtl_top_assoc.txt.gz'),sep='\t')
     top_assoc$celltype <- rep(i,nrow(top_assoc))
     all_top_assoc <- rbind(all_top_assoc,top_assoc)
	 print(i)
 }
}
all_top_assoc <- all_top_assoc[abs(all_top_assoc$tss_distance) < 200000,] # define the range of CNN predictions

all_top_assoc$gene <- gsub('\\..*','',all_top_assoc$phenotype_id)
all_top_assoc$index <- paste0(all_top_assoc$celltype,'_',all_top_assoc$gene)
rownames(all_top_assoc) <- all_top_assoc$index

write.table(unique(all_top_assoc$variant_id),file=paste0(vcf_directory,'topsnp.txt'),row.names=F,col.names=F,quote=F)
write.table(all_top_assoc,file=paste0(vcf_directory,'all_top_assoc.txt'),row.names=F,col.names=T,quote=F,sep='\t')
