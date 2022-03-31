#!/bin/bash -l

#. /usr/share/Modules/init/bash
###This wrapper script validates the input arguments and sends out a qsub cluster job for each new emm_type observed since last run.###

temp_path=$(pwd)
export PATH=$PATH:$temp_path

while getopts :p:m:o: option
do
    case $option in
        p) prevRun=$OPTARG;;
	m) moInput=$OPTARG;;
        o) outPath=$OPTARG;;
    esac
done


if [[ ! -z "$prevRun" ]]
then
    if [[ "$prevRun" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2} ]]
    then
	echo "The previous run occurred on: $prevRun"
    else 
	echo "The previous date is not in the right format. Should be: YYYY-MM-DD"
	exit 1
    fi
else
    echo "No previous run input argument given. Argument format should be: YYYY-MM-DD"
    exit 1
fi

month="18"
if [[ ! -z "$moInput" ]]
then
    re='^[0-9]+$'
    if [[ $moInput =~ $re ]] && [[ $moInput -gt 0 ]]
    then
	month="$moInput"
	echo "The cluster search interval is $month months"
    else
	echo "The search interval (in number of months) has to be a number greater than 0"
	exit 1
    fi
else
    echo "The cluster search interval is 18 months (default)"
fi

curDate=$(date +'%Y'-'%m'-'%d')
echo "Current Date: $curDate"
outDefalt="/scicomp/groups/OID/NCIRD/DBD/RDB/Strep_Lab/StrepLab_Clustr_Detection/GAS_Clusters"
if [[ ! -z "$outPath" ]]
then
    outPath=$(echo $outPath | sed 's/\/$//g')
    if [[ -d "$outPath" ]]
    then
	echo "Output directory exists: $outPath"
	outDir="$outPath/strepClust_$curDate"
    else
	echo "Output directory was created: $outPath"
	mkdir "$outPath"
	outDir="$outPath/strepClust_$curDate"
    fi
else
    echo "Default output path will be used:"
    echo "/scicomp/groups/OID/NCIRD/DBD/RDB/Strep_Lab/StrepLab_Clustr_Detection/GAS_Clusters"
    outDir="$outDefalt/strepClust_$curDate"
fi


##Start Doing Stuff###
echo "output directory path: $outDir"
newTypes=$(mysql -u streplab -h xxxx -pxxxx -e "Use LabDatabase; select DISTINCT L_WGS_emm_subtype from Main_GAS_2015thru2018_Ben where CULT2 IS NULL OR CULT2 between '$prevRun' and now();" | grep -v '[a-zA-Z]' | sed 's/\([0-9]\+\)\..*/\1/g' | sort | uniq)
echo "new types: $newTypes"
# Save current IFS
SAVEIFS=$IFS
# Change IFS to new line. 
IFS=$'\n'
newTypes=($newTypes)
# Restore IFS
IFS=$SAVEIFS
#newTypes=(28 11)
#newTypes=(4)

for m in "${newTypes[@]}"  
do  
    emmDir="$outDir/emm${m}_${curDate}"
    emm_qsub="$emmDir/qsub_output"
    mkdir -p "$emm_qsub"
    echo  "new emm type: $m || emm type dir: $emmDir || search interval: $month"
    #Run qsub command here#
    qsub -q all.q -o "$emm_qsub" -e "$emm_qsub" -cwd strepClust_qsub.sh "$m" "$month" "$emmDir"
done

