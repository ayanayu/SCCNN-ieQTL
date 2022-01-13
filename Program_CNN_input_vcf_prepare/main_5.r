
###### The script is used to merge h19 bed with h38 index

# define paths
vcf_directory='/media/ggj/My Passport1/Xiaoyanyu/CNN/vcf/OtherTissue/' 


h19_bed <- read.csv(paste0(vcf_directory,'h19/unique_pos_bed.txt'),sep='\t',header=F)
aa<-matrix(unlist(strsplit(as.character(h19_bed$V4),"_")),ncol = 5,byrow = T)
h19_bed$ref <- aa[,3]
h19_bed$mut <- aa[,4]
h19_bed <- h19_bed[,c('V1','V2','V3','ref','mut','V4')]
write.table(h19_bed,file=paste0(vcf_directory,'h19/unique_pos.txt'),sep='\t',quote=F,row.names=F,col.names=F)
