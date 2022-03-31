#!/bin/bash -l

temp_path=$(pwd)
export PATH=$PATH:$temp_path

## -- begin embedded SGE options --
#read -a PARAM <<< $(/bin/sed -n ${SGE_TASK_ID}p $1)
## -- end embedded SGE options --

###Load Modules###
module load MUMmer/3.9.4

###Start Doing Stuff###
pathFile=${1}
fullPath=${2}
index=${SGE_TASK_ID}
file_out=nucmer-"$index"_dist.txt
declare -a nucmer_arr
#file_out2=nucmer-"$index"_allDist.txt
index2=$(echo "$index + 1" | bc)
#echo "$index || $pathFile"

genome1=$(sed -n "$index"p $pathFile | cut -d"," -f1)
sample1=$(sed -n "$index"p $pathFile | cut -d"," -f2)
count=0
while IFS= read -r genome; do
    genome2=$(echo "$genome" | cut -d"," -f1)
    sample2=$(echo "$genome" | cut -d"," -f2)
    echo "$genome1 || $genome2"
    nucmer --mum --prefix=TEMP_nucDist_"$index" $genome1 $genome2
    delta-filter -1 TEMP_nucDist_"$index".delta > TEMP_nucDist-Fg_"$index".delta
    #delta-filter -g TEMP_nucDist_"$index".delta > TEMP_nucDist-Fg_"$index".delta
    #/scicomp/home/ycm6/tool_install_testing/test_mummer-3.23/MUMmer3.23/nucmer --mum --prefix=TEMP_nucDist_"$index" $genome1 $genome2
    #/scicomp/home/ycm6/tool_install_testing/test_mummer-3.23/MUMmer3.23/delta-filter -g TEMP_nucDist_"$index".delta > TEMP_nucDist-Fg_"$index".delta
    if [[ -s TEMP_nucDist-Fg_"$index".delta ]]
    then
	nuc_dist=$(show-snps -CHI TEMP_nucDist-Fg_"$index".delta | awk '$6 >= 100 {print $0}' | wc -l) #Changed from 1000
	#nuc_dist=$(show-snps -CIH TEMP_nucDist-Fg_"$index".delta | wc -l) #OR show-snps -HI TEMP_nucDist-Fg_395.delta | awk '$6 >=5 {print $0}'

	#nuc_dist=$(/scicomp/home/ycm6/tool_install_testing/test_mummer-3.23/MUMmer3.23/show-snps -CIH TEMP_nucDist-Fg_"$index".delta | wc -l)
	nuc_out=$(echo "e(-($nuc_dist/2.3-1))" | bc -l) #inverse: (-ln($nuc_exp)+1)*2.3
	#echo "$line || $nuc_out" >> $file_out
	#sample1=$(echo $genome1 | awk -F"/" '{print $(NF-2)}')
	#sample2=$(echo $genome2 | awk -F"/" '{print $(NF-2)}')
	#printf "$sample1 $sample2 $nuc_out $nuc_dist\n" >> $file_out2
	printf "$sample1 $sample2 $nuc_out\n" >> $file_out
	#nucmer_arr[$count]=$nuc_out
    else
	printf "$sample1 $sample2 NA ERROR\n"
        #nucmer_arr[$count]="NA"
    fi
    rm TEMP_nucDist*_"$index".*
    count=$(( $count + 1 ))
done < $fullPath
#done < <(tail -n "+$index2" $pathFile)

rm TEMP_nucDist*_"$index".*
#echo ${nucmer_arr[*]} > TEMP_nucmer-"$index"_arr.txt
