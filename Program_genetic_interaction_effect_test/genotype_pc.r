library(argparser, quietly=TRUE)
p <- arg_parser("")
p <- add_argument(p, "--tissue", short="-c", help="tissue name")
argv <- parse_args(p)

tissue = as.character(argv$tissue)

setwd(paste0('/home/ggj/QTL/PEER/',tissue))

pca <- read.table("./.eigenvec",sep = "\t",header = F)
rownames(pca) <- pca$V2
pca <- pca[,3:ncol(pca)]
colnames(pca) <- paste0('pca','',c(1:5))

cell_fraction <- read.table(paste0('/home/ggj/QTL/deconvolution_result/absolute_fraction/Transformed/',tissue,'_fraction.txt') ,header = 1 , row.names=1,sep='\t')
choosesample <- rownames(cell_fraction)
   
all <- unlist(strsplit(choosesample,'-'))
a <- grep(all, pattern='GTEX')
b <- a+1
choosesample_participant <- paste0(all[a],'-',all[b])

pca_sample <- pca[choosesample_participant,]
rownames(pca_sample) <- choosesample

write.table(t(pca_sample),file='genotype_pcs.txt',sep='\t',quote=F,col.names=T,row.names=T)
