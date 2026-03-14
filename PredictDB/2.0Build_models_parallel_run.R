# Building models——Parallel submission of tasks to the server by tissues and genes============

for (chrom in 1:22) {
  library(parallel)
  suppressMessages(library(dplyr))
  suppressMessages(library(glmnet))
  suppressMessages((library(reshape2)))
  suppressMessages(library(methods))
  "%&%" <- function(a,b) paste(a,b, sep='')
  
  tissue <- "Lymphoblastoid"
  work_path <- "./HapTWAS/PredictDB/Haplotype"
  
  setwd(work_path)
  gene_annot_file <- "./tissue_gene_annotation/" %&% tissue %&% ".gene_annot.txt"
  expression_file <- "./tissue_expression/" %&% tissue %&% ".transformed_expression.txt"
  covariates_file <- "./tissue_covar/" %&% tissue %&% ".covariates.txt"
  snp_annot_file <- "./tissue_snp_annotation/" %&% tissue %&% "/" %&% tissue %&% ".snp_annot.chr" %&% chrom %&% ".txt"
  genotype_file <- "./tissue_genotype_matrix/" %&% tissue %&% "/" %&% tissue %&% ".genotype.chr" %&% chrom %&% ".txt"
  prefix <- "Model_training"
  
  get_gene_annotation <- function(gene_annot_file, chrom, gene_types=c('protein_coding', 'lncRNA')){
    gene_df <- read.table(gene_annot_file, header = TRUE, stringsAsFactors = FALSE) %>%
      filter((chr == chrom) & gene_type %in% gene_types)
    gene_df
  }
  get_gene_expression <- function(expression_file, gene_annot) {
    expr_df <- as.data.frame(read.table(expression_file, header = T, stringsAsFactors = F, row.names = 1))
    expr_df <- expr_df %>% select(one_of(intersect(gene_annot$gene_id, colnames(expr_df))))
    expr_df
  }
  
  gene_annot <- get_gene_annotation(gene_annot_file, chrom)
  expr_df <- get_gene_expression(expression_file, gene_annot)
  n_gene <- length(expr_df)
  
  parallel_main <- function(i) {
    source("./code/2.1parallel_nested_cv_elnet.R")
    main(snp_annot_file, gene_annot, genotype_file, expr_df, covariates_file, as.numeric(chrom), prefix, tissue, as.numeric(i), null_testing=FALSE)
  }
  
  cl <- makeCluster(10, type = "FORK")
  clusterEvalQ(cl, c("dplyr", "glmnet", "reshape2", "methods"))
  clusterExport(cl, c("snp_annot_file", "gene_annot", "genotype_file", "expr_df", "covariates_file", "prefix", "tissue", "chrom"))
  parLapply(cl, 1:n_gene, parallel_main)
  stopCluster(cl)
}

