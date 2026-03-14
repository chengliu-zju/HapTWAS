#!/bin/bash
Tissues=("Lymphoblastoid")
cd ./HapTWAS/PredictDB/Haplotype

###Cycle submit task
for tissue in ${Tissues[*]}
do
	for ttt in "1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22"
	do
		awk 'NR>1{print $0}' ./covariances/${tissue}/${tissue}_Model_training_chr${ttt}_covariances.txt >> ./covariances/${tissue}/chr.tmp #Remove column names and merge files
	done
	sed -i '1i\GENE RSID1 RSID2 VALUE' ./covariances/${tissue}/chr.tmp #Add column name
	gzip -c ./covariances/${tissue}/chr.tmp > ./covariances/${tissue}/${tissue}_Model_training_covariances.txt.gz
	rm ./covariances/${tissue}/chr.tmp #Delete temporary files
done

