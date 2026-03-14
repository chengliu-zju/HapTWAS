#!/bin/bash
tissue=$1
work_path=./HapTWAS/PredictDB/Haplotype


cd ${work_path}/summary/${tissue}
for tt in "1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22"
do
	cd ${tt}
	cat ${tissue}_Model_training_chr${tt}.*.model_summaries.txt > ../${tissue}_Model_training_chr${tt}_model_summaries.txt
	sed -i '1i\gene_id\tgene_name\tgene_type\talpha\tn_snps_in_window\tn_snps_in_model\tlambda_min_mse\ttest_R2_avg\ttest_R2_sd\trho_avg\trho_se\trho_avg_squared\tzscore_est\tzscore_pval\tnested_cv_fisher_pval\tcv_R2_avg\tcv_R2_sd\tcv_rho_avg\tcv_rho_se\tcv_rho_avg_squared\tcv_zscore_est\tcv_zscore_pval\tcv_pval_est\tin_sample_R2\tin_sample_rho\tin_sample_rho_squared' ../${tissue}_Model_training_chr${tt}_model_summaries.txt
	rsync -av ${tissue}_Model_training_chr${tt}_summary.txt ../
	cd ..
done

cd ${work_path}/covariances/${tissue}
for ttt in "1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22"
do
	cd ${ttt}
	cat ${tissue}_Model_training_chr${ttt}.*.covariances.txt > ../${tissue}_Model_training_chr${ttt}_covariances.txt
	sed -i '1i\GENE RSID1 RSID2 VALUE' ../${tissue}_Model_training_chr${ttt}_covariances.txt
	cd ..
done

cd ${work_path}/weights/${tissue}
for tttt in "1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22"
do
	cd ${tttt}
	cat ${tissue}_Model_training_chr${tttt}.*.weights.txt > ../${tissue}_Model_training_chr${tttt}_weights.txt
	sed -i '1i\gene_id\trsid\tvarID\tref\talt\tbeta' ../${tissue}_Model_training_chr${tttt}_weights.txt
	cd ..
done

