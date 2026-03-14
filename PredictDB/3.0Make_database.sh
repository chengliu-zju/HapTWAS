#!/bin/bash
Tissues=("Lymphoblastoid")

# Make the database
for tissue in ${Tissues[*]}
do
# Merge the result files of 2.Build_models_parallel_RUN.R------------
  bash ./HapTWAS/PredictDB/Haplotype/code/3.1combined.sh ${tissue}
# Build the model DB files for each tissue------------
  Rscript ./HapTWAS/PredictDB/Haplotype/code/3.2create_filter_db.R ${tissue}
done

# Combine and compress covariate files for all chromosomes of each tissue model----------------
bash ./HapTWAS/PredictDB/Haplotype/code/3.3make_allchr_covariances_gz.sh

