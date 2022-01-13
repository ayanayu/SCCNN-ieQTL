
#### The script is used to extract 
#### (1)normalized expression matrix of deconvolution samples for each tissue to estimate PEER factors
#### (2)matched individial id for genotype pc calculation

# define paths
deconvolution_result_directory = '/home/ggj/QTL/deconvolution_result/absolute_fraction/Transformed/' # transformed deconvolution result directory(文件以'组织_fraction.txt'命名)   
PEER_directory = "/home/ggj/QTL/PEER/" # PEER directory
expression_bed_file = '/home/ggj/QTL/GTEx_V8.expression.bed.gz' # TMM normalized expression data in BED format


filename <- list.files(path = deconvolution_result_directory)
tissuename <- gsub('_fraction.txt','',filename)

for (i in tissuename){
dir.create(paste0(PEER_directory,i))
}

nrows <- as.integer(system(paste0("zcat ", expression_bed_file, " | wc -l | cut -d' ' -f1 "), intern=TRUE, wait=TRUE))
df <- read.table(expression_bed_file, sep="\t", nrows= nrows, header=TRUE, check.names=FALSE, comment.char="")
row.names(df) <- df[, 4]
df <- df[, 5:ncol(df)]

for (i in tissuename) {
   setwd(paste0(PEER_directory,'',i))     
   cell_fraction <- read.table(paste0(deconvolution_result_directory,i,'_fraction.txt') ,header = 1 , row.names=1,sep='\t')
   choosesample <- rownames(cell_fraction)
   data <- df[,choosesample]
   write.table(data, file='choosesample_exp.txt' ,row.names=T,col.names=T , sep='\t',quote=F) # expression matrix of deconvolution samples for each tissue
} 

for (i in tissuename) {
  cell_fraction <- read.table(paste0(deconvolution_result_directory,i,'_fraction.txt') ,header = 1 , row.names=1,sep='\t')
  choosesample <- rownames(cell_fraction)
  all <- unlist(strsplit(choosesample,'-'))
  a <- grep(all, pattern='GTEX')
  b <- a+1
  choosesample_participant <- paste0(all[a],'-',all[b])
  analyse_sample <- data.frame(name=choosesample_participant)
  write.table(analyse_sample,file=paste0(PEER_directory,i,'/genotype_pc_analyse_sample.txt'),sep='\t',col.names=F,row.names=F,quote=F) # matched individial id
}


