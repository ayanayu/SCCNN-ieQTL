#shell
git clone https://github.com/FunctionLab/ExPecto.git
cd ExPecto
sh download_resources.sh 
tar xf resources_20190807.tar.gz

conda create -n Predictor python=3.7
conda activate Predictor
pip3 install -r requirements.txt

## model traning
cd ../Model_Single/
conda activate Predictor
INDEX=1
for i in $( head -n 1 ../dge_SingleCell/Tissue20_imputation_celltype_group_exp.csv  | sed 's/\"//g' | tr ',' '\n')
do

python ./train.py --expFile ../dge_SingleCell/Tissue20_imputation_celltype_group_exp.csv --targetIndex $INDEX --output ../Model_Single_0804/models/model.$i --evalFile ../Model_Single_0804/pred_target/model.$i.csv 

INDEX=$((INDEX + 1))
done

## computes the chromatin effects of the variants
cd ../Expecto
python chromatin.py ./example_20210728_20kb/AdultArtery_20kb.vcf --cuda &

## model prediction
python predict.py --coorFile ./example_20210728_20kb/AdultArtery_20kb.vcf --geneFile ./example_20210728_20kb/AdultArtery.vcf.bed.sorted.bed_20kb.closestgene --snpEffectFilePattern ./example_20210728_20kb/AdultArtery_20kb.vcf.shift_SHIFT.diff.h5 --modelList ../Model_Single/modellist --output ../Model_Single/Result_AdultArtery_20kb_output.csv &