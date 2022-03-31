#!/bin/bash -l

temp_path=$(pwd)
export PATH=$PATH:$temp_path

. "/scicomp/home/ycm6/miniconda3/etc/profile.d/conda.sh"
conda activate strep_lab
export PERL5LIB=""

###Load Modules###
#. /usr/share/Modules/init/bash
#module load perl/5.30.1-mt

readPair_1=${PARAM[0]}
readPair_2=${PARAM[1]}

mType=${1}
month=${2}
mPath=${3}
echo "emm type: $mType || output path: $mPath || month interval: $month"
#perl /scicomp/home/ycm6/PROJECTS_StrepLab/2018/GAS_clustrTool-v2_7-31-2018/mashFilter.pl -t "${mType}" -o "${mPath}"
if [[ "$month" -ne "18" ]]
then
    perl strepClust_v1.0.pl -t "${mType}" -m "${month}" -o "${mPath}"
else
    perl strepClust_v1.0.pl -t "${mType}" -o "${mPath}"
fi

conda deactivate
#module unload perl/5.30.1-mt
