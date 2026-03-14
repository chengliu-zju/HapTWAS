### ====================================================================================================================
### ====================================================================================================================
### ====== Prepare the input files, with the contents and formats being consistent with that of MetaXcan modeling ======
### ========================== (gene_annot/snp_annot/genotype_chr/expression_chr/covariates) ===========================
### ================================= https://github.com/hakyimlab/PredictDB-Tutorial ==================================
### ====================================================================================================================
### ====================================================================================================================

# # files
# tissue <- "Lymphoblastoid"
# setwd("/disk201/liucheng/TWASmethod.haplo_8.30/1kGP/PredictDB/SNP/summary")
# dir.create(tissue)
# setwd(paste0("/disk201/liucheng/TWASmethod.haplo_8.30/1kGP/PredictDB/SNP/summary/",tissue))
# for (chr in 1:22) {
#   dir.create(as.character(chr))
# }
# 
# # gene annot
# tissue <- "Lymphoblastoid"
# annot <- data.frame(fread("/disk201/liucheng/TWASmethod.haplo_8.30/1kGP/data/expression_filteredGenes.MAGE.v1.0.txt"))
# annot <- annot[which(annot$geneType==c("protein_coding","lncRNA")),]
# annot <- annot[which(!annot$chrom=="chrX"),]
# gene_annotation <- annot[,c(1,5,6,2,3,7)]
# colnames(gene_annotation) <- c("chr","gene_id","gene_name","start","end","gene_type")
# gene_annotation$chr <- as.numeric(gsub("chr","",gene_annotation$chr))
# gene_annotation <- gene_annotation[order(gene_annotation$chr, gene_annotation$start),]
# fwrite(gene_annotation, paste0("/disk201/liucheng/TWASmethod.haplo_8.30/1kGP/PredictDB/SNP/tissue_gene_annotation/",tissue,".gene_annot.txt"),sep="\t")
# 
# # expression
# expression <- data.frame(fread("/disk201/liucheng/TWASmethod.haplo_8.30/1kGP/data/inverse_normal_TMM.filtered.TSS.MAGE.v1.0.bed"))
# expression <- expression[match(gene_annotation$gene_id, expression$ID),]
# expression <- expression[,-(1:3)]
# colnames(expression)[1] <- "Gene_Name"
# n <- expression$Gene_Name
# expression_transpose <- data.frame(t(expression[,-1]))
# colnames(expression_transpose) <- n
# write.table(expression_transpose, paste0("/disk201/liucheng/TWASmethod.haplo_8.30/1kGP/PredictDB/SNP/tissue_expression/",tissue,".transformed_expression.txt"), sep = "\t", row.names = T)
# 
# # covariate
# covariates <- read.table("/disk201/liucheng/TWASmethod.haplo_8.30/1kGP/data/eQTL_covariates.tab")
# x <- covariates[1,][-1]
# y <- covariates$V1[-1]
# covariates <- covariates[-1,-1]
# colnames(covariates) <- x
# rownames(covariates) <- y
# write.table(covariates,paste0("/disk201/liucheng/TWASmethod.haplo_8.30/1kGP/PredictDB/SNP/tissue_covar/",tissue,".covariates.txt"),sep = "\t")
# # 共有部分复制
# cd /disk201/liucheng/TWASmethod.haplo_8.30/1kGP/PredictDB/Haplotype
# find . -maxdepth 1 -type d ! -name "." -exec cp -r /disk201/liucheng/TWASmethod.haplo_8.30/1kGP/PredictDB/SNP/* {} \;
# 
# 
# # 基因型文件
# # 下载
# for chrom in {1..22} X; do
# aria2c -c -x 16 -s 16 \
# -d "/disk201/liucheng/TWASmethod.haplo_8.30/1kGP/data/20220422_3202_phased_SNV_INDEL_SV" \
# "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr${chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
# done
# # 核对
# base_path="/disk201/liucheng/TWASmethod.haplo_8.30/1kGP/data/20220422_3202_phased_SNV_INDEL_SV/"
# for chrom in {1..22} X; do
# filename="1kGP_high_coverage_Illumina.chr${chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
# filepath="${base_path}${filename}"
# md5=$(md5sum "$filepath" | awk '{print $1}')
# echo "$filename: $md5"
# done
# # 提SNP
# max_jobs=10
# for chr in {1..22} X; do
# cd /disk201/liucheng/TWASmethod.haplo_8.30/1kGP/data/20220422_3202_phased_SNV_INDEL_SV
# /disk191/zzy/software/bcftools-1.12/bcftools view -v snps -m2 -M2 ./1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz -Oz -o ./1kGP_high_coverage_Illumina.chr${chr}.filtered.SNP_phased_panel.vcf.gz > chr${chr}.log 2>&1 &
# while (( $(jobs -p | wc -l) >= max_jobs )); do
# sleep 1
# done
# done
# # 提sample
# max_jobs=10
# for chr in {1..22} X; do
# cd /disk201/liucheng/TWASmethod.haplo_8.30/1kGP/data/20250724_731_phased_SNP
# /disk191/zzy/software/bcftools-1.12/bcftools view \
# --threads 8 \
# --samples-file ./samples.list \
# -Oz -o ./1kGP.chr${chr}.filtered.SNP_phased_panel.731.vcf.gz \
# /disk201/liucheng/TWASmethod.haplo_8.30/1kGP/data/20250723_3202_phased_SNP/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNP_phased_panel.vcf.gz \
# > chr${chr}.log 2>&1 &
# while (( $(jobs -p | wc -l) >= max_jobs )); do
# sleep 1
# done
# done
# 
# # 过滤位点
# INPUT_DIR="/disk201/liucheng/TWASmethod.haplo_8.30/1kGP/data/20250724_731_phased_SNP"
# OUTPUT_DIR="/disk201/liucheng/TWASmethod.haplo_8.30/1kGP/data/20250724_731_phased_SNP/filtered_results"
# for CHR in {1..22} X; do
# INPUT_VCF="${INPUT_DIR}/1kGP.chr${CHR}.filtered.SNP_phased_panel.731.vcf.gz"
# PREFIX="${OUTPUT_DIR}/chr${CHR}"
# # vcf转bedbimfam
# /disk191/chenzt/software/plink --vcf $INPUT_VCF \
# --double-id \
# --make-bed \
# --out ${PREFIX}
# # 缺失率和maf过滤
# /disk191/chenzt/software/plink --bfile ${PREFIX} \
# --geno 0.05 \
# --maf 0.05 \
# --make-bed \
# --out ${PREFIX}_geno_maf
# # # LD过滤
# # /disk191/chenzt/software/plink --bfile ${PREFIX}_geno_maf \
# # --indep-pairwise 50 5 0.99 \
# # --out ${PREFIX}_ld
# # /disk191/chenzt/software/plink --bfile ${PREFIX}_geno_maf \
# # --extract ${PREFIX}_ld.prune.in \
# # --make-bed \
# # --out ${PREFIX}_final
# # BEDBIMFAM转vcf.gz
# /disk191/chenzt/software/plink2 --bfile ${PREFIX}_geno_maf \
# --export vcf bgz \
# --out ${PREFIX}_geno_maf
# done
# 
# # PHASE
# cd /disk201/liucheng/TWASmethod.haplo_8.30/1kGP/data/SNP
# for CHR in {15..22}; do
# java -Xmx220g -Djava.io.tmpdir=./tmp -jar /disk191/chenzt/software/beagle.21Apr21.304.jar ne=100 \
# gt=/disk201/liucheng/TWASmethod.haplo_8.30/1kGP/data/20250724_731_phased_SNP/filtered_results/chr${CHR}_geno_maf.vcf.gz \
# seed=4137 nthreads=10 \
# out=chr${CHR}
# /disk191/chenzt/software/plink --vcf chr${CHR}.vcf.gz \
# --double-id \
# --make-bed \
# --out chr${CHR}
# /disk191/chenzt/software/plink --bfile chr${CHR} \
# --recode A \
# --out chr${CHR}
# done
# 
# 
# # geno_annot
# library(data.table)
# tissue <- "Lymphoblastoid"
# for (haplo in c("LD_0.51","LD_0.6","LD_0.7","LD_0.8","SW_2","SW_5","SW_10","SW_20")) {
#   for (i in 1:22) {
#     snp_annotation <- data.frame(fread(paste0("/disk201/liucheng/TWASmethod.haplo_8.30/1kGP/data/Haplotype/",haplo,"/chr",i,".bim")))
#     snp_annotation <- snp_annotation[,c(1,4,2,6,5,2)]
#     colnames(snp_annotation) <- c("chromosome","pos","varID","ref_vcf","alt_vcf","rsid")
#     fwrite(snp_annotation, paste0("/disk201/liucheng/TWASmethod.haplo_8.30/1kGP/PredictDB/Haplotype/",haplo,"/tissue_snp_annotation/",tissue,"/",tissue,".snp_annot.chr",i,".txt"),sep="\t")
#   }
# }
# 
# rsid <- function(i) {
#   setwd('/disk201/liucheng/TWASmethod.haplo_8.30/1kGP/data/GRCh38')
#   bcftools <- '/disk191/zzy/software/bcftools-1.12/bcftools'
#   system(paste0(bcftools,' view -Oz -r ',i,' GCF_000001405.40.gz -o chr_',i,'.gz'))
# }
# snp_annot <- function(i) {
#   library(data.table)
#   tissue <- "Lymphoblastoid"
#   snp_annotation <- data.frame(fread(paste0("/disk201/liucheng/TWASmethod.haplo_8.30/1kGP/data/SNP/chr",i,".bim")))
#   snp_annotation <- snp_annotation[,c(1,4,2,6,5,2)]
#   colnames(snp_annotation) <- c("chromosome","pos","varID","ref_vcf","alt_vcf","rsid")
# 
#   rsid <- data.frame(fread(paste0("/disk201/liucheng/TWASmethod.haplo_8.30/1kGP/data/GRCh38/chr",i,".gz"),skip = 680))
#   rsid <- rsid[nchar(rsid$ALT)==1 & nchar(rsid$REF)==1,]
#   match_row <- snp_annotation$pos%in%rsid$POS
#   snp_annotation$rsid[match_row] <- rsid$ID[match(snp_annotation$pos[match_row],rsid$POS)]
#   fwrite(snp_annotation, paste0("/disk201/liucheng/TWASmethod.haplo_8.30/1kGP/PredictDB/SNP/tissue_snp_annotation/",tissue,"/",tissue,".snp_annot.chr",i,".txt"),sep="\t")
# }
# 
# 
# # genotype
# geno_haplo <- function(haplo){
#   library(data.table)
#   tissue <- "Lymphoblastoid"
#   for (i in 1:22) {
#     genotype <- fread(paste0("/disk201/liucheng/TWASmethod.haplo_8.30/1kGP/data/Haplotype/",haplo,"/chr",i,".raw"))
#     genotype <- t(genotype)
#     colnames(genotype) <- genotype[2,]
#     genotype <- genotype[-(1:6),]
#     rownames(genotype) <- as.character(substr(rownames(genotype),1,nchar(rownames(genotype))-2))
#     genotype <- cbind(rownames(genotype),genotype)
#     colnames(genotype)[1] <- "varID"
#     fwrite(genotype,paste0("/disk201/liucheng/TWASmethod.haplo_8.30/1kGP/PredictDB/Haplotype/",haplo,"/tissue_genotype_matrix/",tissue,"/",tissue,".genotype.chr",i,".txt"), row.names = F, sep = "\t")
#   }
# }
# 
# library(data.table)
# tissue <- "Lymphoblastoid"
# for (i in 1:22) {
#   genotype <- fread(paste0("/disk201/liucheng/TWASmethod.haplo_8.30/1kGP/data/SNP/chr",i,".raw"))
#   genotype <- t(genotype)
#   colnames(genotype) <- genotype[1,]
#   genotype <- genotype[-(1:6),]
#   rownames(genotype) <- as.character(substr(rownames(genotype),1,nchar(rownames(genotype))-2))
#   genotype <- cbind(rownames(genotype),genotype)
#   colnames(genotype)[1] <- "varID"
#   fwrite(genotype,paste0("/disk201/liucheng/TWASmethod.haplo_8.30/1kGP/PredictDB/SNP/tissue_genotype_matrix/",tissue,"/",tissue,".genotype.chr",i,".txt"), row.names = F, sep = "\t")
# }

