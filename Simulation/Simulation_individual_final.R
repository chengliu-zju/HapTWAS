library(MASS)
library(gtools)
library(glmnet)
library(pROC)
library(PRROC)
library(parallel)
library(doParallel)
library(foreach)
library(purrr)


# -----------------------------
# Step 1: SNPs and Haplotypes
simulate_snp_hap <- function(n, n_genes, snps_per_gene, n_blocks_per_gene, snps_per_block, haps_per_block, hap_freq_conc) {
  snp_list <- list()
  hap_list <- list()
  
  for (g in 1:n_genes) {
    snp_mat_gene <- NULL
    hap_mat_gene <- NULL
    
    for (b in 1:n_blocks_per_gene) {
      p <- snps_per_block
      K <- haps_per_block
      
      # 定义 block 内的 haplotype 模式，确保多样性和MAF约束
      unique_haps <- FALSE
      while(!unique_haps) {
        pat <- matrix(rbinom(K * p, size = 1, prob = 0.3), nrow = K, ncol = p)      # 0.3合理，0.5SNP完全随机不真实、0.1大部分是0多样性差
        # 确保每个 SNP 至少有一个但不全为 1 （不然SNP没变异）且每个单倍型独特
        if(all(colSums(pat) > 0 & colSums(pat) < K) & nrow(unique(pat)) == K) {
          unique_haps <- TRUE
        }
      }
      
      variant_all <- FALSE
      while(!variant_all) {
        freqs <- as.numeric(gtools::rdirichlet(1, rep(hap_freq_conc, K)))      # 单倍型频率(Dirichlet)
        
        # 每个个体随机抽两条 hap → SNP 基因型
        hap_index <- replicate(2, sample(K, n, replace=TRUE, prob=freqs))
        snp_block <- pat[hap_index[,1], ] + pat[hap_index[,2], ]
        
        # hap 剂量矩阵
        hap_block <- matrix(0, n, K)
        for (i in 1:n) {
          hap_block[i, hap_index[i,1]] <- hap_block[i, hap_index[i,1]] + 1
          hap_block[i, hap_index[i,2]] <- hap_block[i, hap_index[i,2]] + 1
        }
        
        # 确保每个 SNP 和 hap 都有变异
        all_same_cols <- apply(cbind(snp_block, hap_block), 2, function(col) all(col == col[1]))
        if(all(!all_same_cols)) {
          variant_all <- TRUE
        }
      }
      
      # 拼接到基因整体矩阵
      snp_mat_gene <- cbind(snp_mat_gene, snp_block)
      hap_mat_gene <- cbind(hap_mat_gene, hap_block)
    }
    
    snp_list[[g]] <- snp_mat_gene
    hap_list[[g]] <- hap_mat_gene
  }
  
  list(snp=snp_list, hap=hap_list)
}

# -----------------------------
# Step 2: gene expression
simulate_expression <- function(geno_snp, geno_hap, n_causal_snp, n_causal_hap, prop_hap_driven, h2_G) {
  n_genes <- length(geno_snp)
  n <- nrow(geno_snp[[1]])
  G_mat <- matrix(NA, n, n_genes)
  causal_info <- list()
  
  for (g in 1:n_genes) {
    geno_s <- scale(geno_snp[[g]], center=TRUE, scale=TRUE)
    geno_h <- scale(geno_hap[[g]], center=TRUE, scale=TRUE)
    
    causal_s <- sample.int(ncol(geno_s), n_causal_snp)
    causal_h <- sample.int(ncol(geno_h), n_causal_hap)
    w_s <- rnorm(length(causal_s), 0, 1)
    w_h <- rnorm(length(causal_h), 0, 1.5)
    
    signal_s <- as.numeric(geno_s[, causal_s, drop=FALSE] %*% w_s)
    signal_h <- as.numeric(geno_h[, causal_h, drop=FALSE] %*% w_h)
    
    signal <- signal_s*(1-prop_hap_driven) + signal_h*prop_hap_driven
    var_signal <- var(signal)
    sigma_g2 <- var_signal * (1 - h2_G) / h2_G
    
    G <- signal + rnorm(n, 0, sqrt(sigma_g2))
    G_mat[, g] <- G
    
    causal_info[[g]] <- list(causal_snp=causal_s, causal_hap=causal_h, weights_snp=w_s, weights_hap=w_h, sigma_g2=sigma_g2)
  }
  
  list(G = G_mat, causal_info = causal_info)
}

# -----------------------------
# Step 3: phenotype
simulate_phenotype <- function(G_mat, causal_genes, h2_Y) {
  n <- nrow(G_mat)
  alpha <- rnorm(length(causal_genes), 0, 1)
  
  signal <- as.numeric(G_mat[, causal_genes] %*% alpha)
  var_signal <- var(signal)
  sigma_eta2 <- var_signal * (1 - h2_Y) / h2_Y
  
  Y <- signal + rnorm(n, 0, sqrt(sigma_eta2))
  
  list(Y=Y, causal_genes=causal_genes, alpha=alpha, sigma_eta2=sigma_eta2)
}

# -----------------------------
# Step 4: simulate dataset
simulate_dataset <- function(n, n_genes, snps_per_gene, n_blocks_per_gene, snps_per_block, haps_per_block, hap_freq_conc,
                             n_causal_snp, n_causal_hap, prop_hap_driven, h2_G, n_causal_gene, h2_Y) {
  
  geno <- simulate_snp_hap(n, n_genes, snps_per_gene, n_blocks_per_gene, snps_per_block, haps_per_block, hap_freq_conc)
  
  expr <- simulate_expression(geno$snp, geno$hap, n_causal_snp, n_causal_hap, prop_hap_driven, h2_G)
  
  if (prop_hap_driven == 1) {
    driver_type <- rep("haplo", n_genes)
  } else if (prop_hap_driven == 0) {
    driver_type <- rep("snp", n_genes)
  } else {
    driver_type <- rep("mixed", n_genes)
  }
  
  causal_genes <- sample(1:n_genes, n_causal_gene)
  pheno <- simulate_phenotype(expr[["G"]], causal_genes, h2_Y)
  
  list(
    snp_geno = geno$snp,
    hap_geno = geno$hap,
    expression = expr[["G"]],
    causal_var = expr$causal_info,
    driver_type = driver_type,
    causal_genes = causal_genes,
    pheno = pheno$Y,
    alpha = pheno$alpha,
    sigma_eta2 = pheno$sigma_eta2
  )
}

# -----------------------------
# Step 5: SNP-TWAS and haplo-TWAS
run_twas <- function(sim_data, n_train, n_test) {
  n_all <- nrow(sim_data$expression)
  n_genes <- ncol(sim_data$expression)
  
  idx_all <- sample(1:n_all)
  idx_train <- idx_all[1:n_train]
  idx_test <- idx_all[(n_train+1):n_all]
  
  Y_test <- as.numeric(scale(sim_data$pheno[idx_test]))
  
  pvals_snp <- rep(NA, n_genes)  # significance p in pred_expr associated pheno
  pvals_hap <- rep(NA, n_genes)
  r2_pred_snp <- rep(NA, n_genes)  # prediction R2 in pred vs true expr
  r2_pred_hap <- rep(NA, n_genes)
  
  # For each gene train elastic net on reference expression
  for(g in seq_len(n_genes)) {
    X_snp_all <- sim_data$snp_geno[[g]]
    X_hap_all <- sim_data$hap_geno[[g]]
    # reference and test splits
    Xs_train <- X_snp_all[idx_train, , drop=FALSE]
    Xs_test <- X_snp_all[idx_test, , drop=FALSE]
    Xh_train <- X_hap_all[idx_train, , drop=FALSE]
    Xh_test <- X_hap_all[idx_test, , drop=FALSE]
    
    y_train <- as.numeric(scale(sim_data$expression[idx_train, g]))
    
    # ---- SNP elastic-net (cv) ----
    try({
      cv_snp <- cv.glmnet(Xs_train, y_train, alpha = 0.5, nfolds = 5, standardize = TRUE)
      # beta_snp <- coef(cv_snp, s = "lambda.min")
      # predicted expression in test
      pred_snp <- predict(cv_snp, Xs_test, s = "lambda.min")[,1]
      # regress Y_test ~ pred_snp
      fit_snp <- lm(Y_test ~ pred_snp)
      pvals_snp[g] <- summary(fit_snp)$coefficients["pred_snp","Pr(>|t|)"]
      # compute prediction R2
      true_G_test <- as.numeric(scale(sim_data$expression[idx_test, g]))
      r2_pred_snp[g] <- cor(pred_snp, true_G_test)^2
    }, silent = TRUE)
    
    # ---- haplotype elastic-net (cv) ----
    try({
      cv_hap <- cv.glmnet(Xh_train, y_train, alpha = 0.5, nfolds = 5, standardize = TRUE)
      # beta_hap <- coef(cv_hap, s = "lambda.min")[-1]
      pred_hap <- predict(cv_hap, Xh_test, s = "lambda.min")[,1]
      fit_hap <- lm(Y_test ~ pred_hap)
      pvals_hap[g] <- summary(fit_hap)$coefficients["pred_hap","Pr(>|t|)"]
      true_G_test <- as.numeric(scale(sim_data$expression[idx_test, g]))
      r2_pred_hap[g] <- cor(pred_hap, true_G_test)^2
    }, silent = TRUE)
  }
  
  pvals_snp[is.na(pvals_snp)] <- 1
  pvals_hap[is.na(pvals_hap)] <- 1
  r2_pred_snp[is.na(r2_pred_snp)] <- 0
  r2_pred_hap[is.na(r2_pred_hap)] <- 0
  
  list(pvals_snp = pvals_snp, pvals_hap = pvals_hap, r2_pred_snp = r2_pred_snp, r2_pred_hap = r2_pred_hap)
}

# -----------------------------
# Step 6: metrics computation helpers
compute_metrics <- function(pvals, causal, alpha = 0.05) {
  qvals <- p.adjust(pvals, method="BH")
  discovered <- which(qvals <= alpha)
  TP <- sum(causal[discovered] == 1)
  FP <- length(discovered) - TP
  FN <- sum(causal == 1) - TP
  TN <- length(causal) - sum(causal == 1) - FP
  
  power <- ifelse(sum(causal==1)==0, NA, TP / sum(causal==1))  # power即recall
  precision <- ifelse((TP+FP)==0, NA, TP / (TP+FP))
  FDR <- ifelse((TP+FP)==0, NA, FP / (TP+FP))
  
  F1 <- ifelse(is.na(precision) | is.na(power) | (precision+power)==0, NA, 2*precision*power/(precision+power))
  num <- TP*TN - FP*FN
  den <- sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  MCC <- ifelse(den==0, NA, num/den)
  
  roc_obj <- tryCatch(roc(response = causal, predictor = -log10(pvals + 1e-300), quiet=TRUE), error = function(e) NULL)
  auc_val <- ifelse(is.null(roc_obj), NA, as.numeric(roc_obj$auc))
  
  pr_obj <- tryCatch(pr.curve(scores.class0 = -log10(pvals + 1e-300)[causal == 1],
                              scores.class1 = -log10(pvals + 1e-300)[causal == 0], curve = TRUE), error = function(e) NULL)
  aupr_val <- ifelse(is.null(pr_obj), NA, as.numeric(pr_obj$auc.davis.goadrich))
  
  list(power = power, precision = precision, FDR = FDR, F1 = F1, MCC = MCC, auc_val = auc_val, aupr_val = aupr_val)
}

# -----------------------------
# Step 7: run one replicate simulation
run_one_replicate <- function(params) {
  n <- params$n
  n_genes <- params$n_genes
  snps_per_gene <- params$snps_per_gene
  n_blocks_per_gene <- params$n_blocks_per_gene
  snps_per_block <- params$snps_per_block
  haps_per_block <- params$haps_per_block
  hap_freq_conc <- params$hap_freq_conc
  
  n_causal_snp <- params$n_causal_snp
  n_causal_hap <- params$n_causal_hap
  prop_hap_driven <- params$prop_hap_driven
  h2_G <- params$h2_G
  n_causal_gene <- params$n_causal_gene
  h2_Y <- params$h2_Y
  
  n_train <- params$n_train
  n_test <- params$n_test
  
  
  sim_data <- simulate_dataset(n, n_genes, snps_per_gene, n_blocks_per_gene, snps_per_block, haps_per_block, hap_freq_conc,
                               n_causal_snp, n_causal_hap, prop_hap_driven, h2_G, n_causal_gene, h2_Y)
  
  twas_result <- run_twas(sim_data, n_train, n_test)
  
  causal_genes <- rep(0, n_genes); causal_genes[sim_data$causal_genes] <- 1
  metrics_snp <- compute_metrics(twas_result$pvals_snp, causal_genes)
  metrics_hap <- compute_metrics(twas_result$pvals_hap, causal_genes)
  
  res_tab <- data.frame(
    gene = seq_len(n_genes),
    driver = sim_data$driver_type,
    causal = causal_genes,
    p_snp = twas_result$pvals_snp,
    p_hap = twas_result$pvals_hap,
    r2_pred_snp = twas_result$r2_pred_snp,
    r2_pred_hap = twas_result$r2_pred_hap,
    q_snp = p.adjust(twas_result$pvals_snp, method="BH"),
    q_hap = p.adjust(twas_result$pvals_hap, method="BH"),
    stringsAsFactors = FALSE
  )
  
  list(metrics_snp = metrics_snp, metrics_hap = metrics_hap, res_tab = res_tab)
}

# -----------------------------
# Step 8: multiple replicates
run_replicates <- function(params, parallel = TRUE) {
  nrep <- params$nrep
  n_cores <- params$n_cores
  
  seeds <- sample.int(1e6, nrep)
  res_list <- vector("list", nrep)
  
  if(parallel) {
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    res_list <- foreach::foreach(i = 1 : nrep,
                                 .packages = c("MASS","gtools","glmnet","pROC","PRROC"),
                                 .export = c("simulate_snp_hap","simulate_expression","simulate_phenotype","simulate_dataset",
                                             "run_twas","compute_metrics","run_one_replicate")) %dopar% {
                                               set.seed(seeds[i])
                                               res <- run_one_replicate(params)
                                               list(result = res, seed = seeds[i])
                                             }
    parallel::stopCluster(cl)
  } else {
    for(i in seq_len(nrep)) {
      set.seed(seeds[i])
      res_list[[i]] <- list(result = run_one_replicate(params), seed = seeds[i])
      if(i %% 20 == 0) cat("Completed replicate", i, "/", nrep, "\n")
    }
  }
  
  return(res_list)
}

# -----------------------------
# Step 9: summarize replicate results
summarize_replicates <- function(res_list, output) {
  seeds <- do.call(rbind, lapply(seq_along(res_list), function(i) {
    data.frame(replicate = i, seed = res_list[[i]]$seed)
  }))
  metrics_snp <- do.call(rbind, lapply(seq_along(res_list), function(i) {
    df <- data.frame(res_list[[i]]$result$metrics_snp)
    df$replicate <- i
    df
  }))
  metrics_hap <- do.call(rbind, lapply(seq_along(res_list), function(i) {
    df <- data.frame(res_list[[i]]$result$metrics_hap)
    df$replicate <- i
    df
  }))
  results_gene <- do.call(rbind, lapply(seq_along(res_list), function(i) {
    df <- data.frame(res_list[[i]]$result$res_tab)
    df$replicate <- i
    df
  }))
  
  agg_metrics <- function(m, r, type) {
    data.frame(
      power_mean = mean(m$power, na.rm=TRUE), power_sd = sd(m$power, na.rm=TRUE),
      precision_mean = mean(m$precision, na.rm=TRUE), precision_sd = sd(m$precision, na.rm=TRUE),
      FDR_mean = mean(m$FDR, na.rm=TRUE), FDR_sd = sd(m$FDR, na.rm=TRUE),
      F1_mean = mean(m$F1, na.rm=TRUE), F1_sd = sd(m$F1, na.rm=TRUE),
      MCC_mean = mean(m$MCC, na.rm=TRUE), MCC_sd = sd(m$MCC, na.rm=TRUE),
      AUC_mean = mean(m$auc_val, na.rm=TRUE), AUC_sd = sd(m$auc_val, na.rm=TRUE),
      AUPR_mean = mean(m$aupr_val, na.rm=TRUE), AUPR_sd = sd(m$aupr_val, na.rm=TRUE),
      R2_mean = mean(r[[paste0("r2_pred_",type)]], na.rm=TRUE), R2_sd = sd(r[[paste0("r2_pred_",type)]], na.rm=TRUE)
    )
  }
  
  metrics_snp_summary <- agg_metrics(metrics_snp, results_gene, "snp")
  metrics_hap_summary <- agg_metrics(metrics_hap, results_gene, "hap")
  metrics_summary <- rbind(cbind(Method="SNP-TWAS", metrics_snp_summary), cbind(Method="Haplo-TWAS", metrics_hap_summary))
  
  write.csv(seeds, paste0(output, ".seeds.csv"), row.names = FALSE)
  write.csv(metrics_snp, paste0(output, ".metrics_snp.csv"), row.names = FALSE)
  write.csv(metrics_hap, paste0(output, ".metrics_hap.csv"), row.names = FALSE)
  write.csv(results_gene, paste0(output, ".results_gene.csv"), row.names = FALSE)
  write.csv(metrics_summary, paste0(output, ".metrics_summary.csv"), row.names = FALSE)
}

# -----------------------------
# Step 10: full pipeline wrapper
run_full_simulation <- function(params, output) {
  res_list <- run_replicates(params)
  
  summarize_replicates(res_list, output)
}


# -----------------------------
# Running the parameter combination
base_params <- list(
  n = 5500,                   # n = n_train + n_test
  n_train = 500,              # 100 500
  n_test = 5000,              # 1000 5000
  
  h2_Y = 0.1,                 # 0.1 0.3 0.6
  h2_G = 0.1,                 # 0.1 0.3
  n_genes = 20,               # 20 50
  n_causal_gene = 2,          # 2 4 6
  n_causal_snp = 2,           # 2 10 20
  n_causal_hap = 2,           # n_causal_hap = n_causal_snp
  
  snps_per_gene = 200,
  snps_per_block = 5,
  n_blocks_per_gene = 40,     # n_blocks_per_gene = snps_per_gene / snps_per_block
  haps_per_block = 5,
  hap_freq_conc = 1,          # concentration for Dirichlet 越小越稀疏（一两个主导），越大越均匀
  prop_hap_driven = 0.5,
  nrep = 200,
  n_cores = 40
)

param_grid <- expand.grid(
  h2_Y = c(0.1, 0.3, 0.6),
  h2_G = c(0.1, 0.3),
  n_genes = c(20, 50),
  n_causal_gene = c(2, 4, 6),
  n_causal_var = c(2, 10, 20),
  stringsAsFactors = FALSE
)
param_grid$sim <- seq_len(nrow(param_grid))
write.csv(param_grid, "./simulation/result1/param_grid.csv", row.names = FALSE)

modify_params <- function(base, new_values) {
  modifyList(base, as.list(new_values))
}

for (i in seq_len(nrow(param_grid))) {
  current_params <- modify_params(base_params, list(h2_Y = param_grid$h2_Y[i],
                                                    h2_G = param_grid$h2_G[i],
                                                    n_genes = param_grid$n_genes[i],
                                                    n_causal_gene = param_grid$n_causal_gene[i],
                                                    n_causal_snp = param_grid$n_causal_var[i],
                                                    n_causal_hap = param_grid$n_causal_var[i]))
  output_prefix <- paste0("./simulation/result1/500train_5000test/sim",i)
  
  message("=== Start running the parameter combination  ", i, " / ", nrow(param_grid), " ===")
  set.seed(71717)
  run_full_simulation(current_params, output_prefix)
}

