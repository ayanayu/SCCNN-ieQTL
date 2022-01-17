####### define user input parameter ####
run_script_folder=$(dirname $0)
bfile_directory=${run_script_folder}"/data/xxx"
expr_bed_file=${run_script_folder}"/home/ggj/QTL/GTEx_V8.expression.bed.gz"
tissue_20_impute_dir=${run_script_folder}"/data/xxxxx"
deconv_result_dir=${run_script_folder}"/home/ggj/QTL/deconvolution_result/absolute_fraction/Transformed/"
output_dir=${run_script_folder}"/media/ggj/Guo-4T-MS2/XiaoYanyu/result/"

celltype_tissue_anno=${run_script_folder}'/home/ggj/QTL/Celltype-GTExSampleNumber-Tissue-Annotation.txt'
sample_participant_anno=${run_script_folder}'/home/ggj/QTL/sample_participant_lookup.txt'

function print_help_info() {
	echo "**********this page showes how to use SCCNN-ieQTL*********"
	echo "		-b,  --bfile_dir,	xxxxxintroduction,    default=${bfile_directory}"
	echo "          -expr, --expr_bed_file,		xxxintroduction,	default=${expr_bed_file}"
	echo "          -tissue, --tissue_20_impute_dir, 	xxxintroduction,    default=${tissue_20_impute_dir}"
	echo "          -deconv, --deconv_result_dir, 	xxxintroduction,    default=${deconv_result_dir}"
	echo "          -output, --output_dir, 	xxxintroduction,    default=${output_dir}"
	echo "          -h, --help,	help info"
}

if [ $# -ge 1 ]; then
	if [ $1 = "-h" ] || [ $1 = "--help" ]; then
		print_help_info; exit 0;
	fi
fi

while [ $# -ge 2 ]; do
	case "$1" in
		-b|--bfile_dir)
			bfile_directory=$2; shift 2;;
		-expr|--expr_bed_file_dir)
			expr_bed_file=$2; shift2;;
		-tissue|--tissue_20_impute_dir)
			tissue_20_impute_dir=$2; shift2;;
		-deconv|--deconv_result_dir)
			deconv_result_dir=$2; shift 2;;
		*)
			echo "unknown input parameter!"; exit 1;;
	esac
done

########## define internal path ###############
PEER_directory=${run_script_folder}"/home/ggj/QTL/PEER/"
coordinate_file=${run_script_folder}"resource/xxx"


alias Rscript='/usr/local/bin/Rscript'
alias Plink2='/home/ggj/plink2'

##### part I: Program_genetic_interaction_effect_test #######
Rscript ./Program_genetic_interaction_effect_test/main_1.r --deconv_result_dir ${deconv_result_dir} \
							   --peer_dir ${PEER_directory} \
							   --expr_bed_file ${expr_bed_file}
tissuename=$(ls ${PEER_directory})

for i in ${tissuename}; do
	# peer run
	echo "${i} start";
	echo "${i} PEER run";
	cd ${PEER_directory};
	cd ${i};
	Rscript $run_PEER_script choosesample_exp.txt GTEx_V8 60 --max_iter 200;
	echo "${i} PEER finish";

	# get genotype
	echo "${i} get genotype pcs start";
	cd ${PEER_directory};
	cd ${i};
	local testout="./";
	Plink2 --allow-extra-chr --threads 20 -bfile ${bfile_directory} --pca 5 --out ${testout} --keep \
								        genotype_pc_analyse_sample.txt
	Rscript ./Program_genetic_interaction_effect_test/genotype_pc.r --deconv_result_dir ${deconv_result_dir} \
									--peer_dir ${PEER_directory} \
									--tissue ${i};
	echo "${i} get genotype pcs finish";

	# combine covariates
	echo "${i} combine covariates start";
	local prefix="GTEx_V8";
	local genotype_pcs="genotype_pcs.txt";
	cd ${PEER_directory};
	cd ${i};
	python3 ./Program_genetic_interaction_effect_test/combine_covariates.py ${prefix}.PEER_covariates.txt ${prefix} \
									--genotype_pcs ${genotype_pcs}
	echo "${i} combine covariates finish";

	echo "${i} all finish";
done

all_celltype_list=("AdultArtery3" "AdultArtery4" "AdultArtery5" "AdultArtery6" "AdultArtery7")
for i in ${all_celltype_list}; do
	mkdir -p ${output_dir}"qtl_result/"${i}
done

python3 ./Program_genetic_interaction_effect_test/main_3.py --peer_dir ${PEER_directory} --bfile_dir ${bfile_directory} \
							    --expr_bed_file ${expr_bed_file} --deconv_result_dir ${deconv_result_dir} \
							    --celltype_tissue_anno ${celltype_tissue_anno}
							    --sample_participant_anno ${sample_participant_anno}
							    --qtl_output_dir ${output_dir}"qtl_result/"

##### part II: Program_CNN_input_vcf_prepare #######
mkdir -p ${output_dir}"cnn_input/"

Rscript ./Program_CNN_input_vcf_prepare/cnn_input_vcf_prepare.r --ieqtl_output_dir ${output_dir}"qtl_result/" \
								--vcf_dir ${output_dir}"cnn_input/" \
								--bfile_dir ${bfile_directory} \
								--coord_file ${coordinate_file}

##### part III: Program_CNN_predict #######



