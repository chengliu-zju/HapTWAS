cov <- function(tissue) {
  options(stringsAsFactors = FALSE)
  library(limma)
  library(edgeR)
  library(peer)
  library(preprocessCore)
  library(RNOmni)
  library(data.table)
  library(R.utils)
  library(SNPRelate)
  "%&%" = function(a, b) { paste0(a, b) }
  
  # Genotype PCA
  vcf.fn <- "./genotype/VCF_tissues/" %&% tissue %&% "_maf.recode.vcf.gz"
  snpgdsVCF2GDS(vcf.fn, "./covariates/pca_peer/" %&% tissue %&% ".ccm.gds", method = "biallelic.only")
  genofile <- openfn.gds("./covariates/pca_peer/" %&% tissue %&% ".ccm.gds")
  ccm_pca <- snpgdsPCA(genofile,num.thread=23)
  
  pca_genotype <- ccm_pca$eigenvect[, 1:10]
  colnames(pca_genotype) <- paste0("pc", 1:10)
  rownames(pca_genotype) <- ccm_pca$sample.id
  pca_genotype0 = data.frame(SampleID=ccm_pca$sample.id,pca_genotype)
  pca_var0 = data.frame(pc=1:10,eigenval=ccm_pca$eigenval[1:10],varprop=ccm_pca$varprop[1:10])
  pca_vect = pca_genotype0
  
  n_samples = nrow(pca_vect)
  if(n_samples<200){
    cov_pc = pca_vect[,2:6]
  } else if(n_samples>=200){
    cov_pc = pca_vect[,2:11]
  }
  
  fwrite(pca_genotype0, "./covariates/pca_peer/" %&% tissue %&% ".PCA_eigenvect.txt", sep = "\t", row.names = F, quote = FALSE)
  fwrite(pca_var0, "./covariates/pca_peer/" %&% tissue %&% ".PCA_var.txt", sep = "\t", row.names = F, quote = FALSE)
  
  
  # Peer factor estimation
  bed <- data.frame(fread("./expression_tpm/" %&% tissue %&% ".expr_tpm.bed.gz"))
  rownames(bed) = bed$gene_id
  bed$gene_id = NULL
  bed$X.Chr = NULL
  bed$start = NULL
  bed$end = NULL
  expr_peer = bed
  
  k = 10
  model = PEER()
  PEER_setPhenoMean(model, as.matrix(t(expr_peer))) #NULL response means no err# #N rows (samples); G columns (Genes) #have been sorted based on genotype samples
  dim(PEER_getPhenoMean(model))
  PEER_setNk(model, k)
  PEER_setNmax_iterations(model, 1000)
  PEER_getNk(model)
  PEER_update(model)
  factors = PEER_getX(model)
  
  rownames(factors) = colnames(expr_peer)
  colnames(factors) = paste0("peer",c(1:k))
  
  Alpha = PEER_getAlpha(model)
  
  alpha0 = data.frame(Peer=c(1:k), Alpha=Alpha, Relevance = 1.0 / Alpha)
  residuals0 = t(PEER_getResiduals(model))
  residuals0 = data.frame(GENE_ID=rownames(expr_peer),residuals0)
  colnames(residuals0) = colnames(expr_peer)
  covariates0 <- data.frame(SampleID = colnames(expr_peer), factors)
  
  
  write.table(covariates0, "./covariates/pca_peer/" %&% tissue %&% ".PEER_covariates.Nk" %&% k %&% ".txt", sep = "\t", row.names = F, quote = FALSE)
  write.table(alpha0, "./covariates/pca_peer/" %&% tissue %&% ".PEER_alpha.Nk" %&% k %&% ".txt", sep = "\t", row.names = F, quote = FALSE)
  write.table(residuals0, "./covariates/pca_peer/" %&% tissue %&% ".PEER_residuals.Nk" %&% k %&% ".txt", sep = "\t", row.names = F, quote = FALSE)
  
  
  # merge all covariates
  peer_vect = covariates0
  cov_peer = peer_vect[,2:11]
  cov0 = as.matrix(t(cbind(pca_vect[,1],cov_pc,cov_peer)))
  colnames(cov0) = cov0[1,]
  cov0 = cov0[-1,]
  cov_output = data.frame(rownames(cov0),cov0)
  colnames(cov_output)[1] = ""
  fwrite(cov_output,"./covariates/" %&% tissue %&% ".covariates.txt",sep="\t",quote=F)
}


library(parallel)
cl <- makeCluster(7)
system.time(res <- parLapply(cl,c('Adipose','Blood','Brain','Bursa','Cecum','Cerebellum','Duodenum'),cov))
stopCluster(cl)

