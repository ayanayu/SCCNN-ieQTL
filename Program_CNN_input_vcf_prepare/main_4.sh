
###### The script is used to transform hg38 to hg19 bed 

# define paths
vcf_directory=/media/ggj/My\ Passport1/Xiaoyanyu/CNN/vcf/OtherTissue/ 
coordinate_file=/home/ggj/hg38ToHg19.over.chain.gz

mkdir ${vcf_directory}h19
liftOver ${vcf_directory}unique_pos_bed.txt $coordinate_file ${vcf_directory}h19/unique_pos_bed.txt ${vcf_directory}h19/'NotMap'
