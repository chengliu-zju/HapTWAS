library(MASS)
library(gtools)
library(glmnet)
library(psych)
library(parallel)
library(doParallel)
library(foreach)
library(dplyr)


# -----------------------------
# Step 1: SNPs and Haplotypes
simulate_snp_hap <- function(n, n_genes, snps_per_gene, n_blocks_per_gene, snps_per_block, haps_per_block, hap_freq_conc) {
  snp_list <- list()
  hap_list <- list()
  pat_list <- list()
  
  for (g in 1:n_genes) {
    snp_mat_gene <- NULL
    hap_mat_gene <- NULL
    pat_lis_gene <- list()
    
    for (b in 1:n_blocks_per_gene) {
      p <- snps_per_block
      K <- haps_per_block
      
      # 定义 block 内的 haplotype 模式，确保多样性和MAF约束
      unique_haps <- FALSE
      while(!unique_haps) {
        pat <- matrix(rbinom(K * p, size = 1, prob = 0.3), nrow = K, ncol = p)      # 0.3合理，0.5SNP完全随机不真实、0.1大部分是0多样性差
        # 确保每个 SNP 至少有一个但不全为 1 （不然SNP没变异）且每个单倍型独特
        if(all(colSums(pat) > 0 & colSums(pat) < K) & nrow(unique(pat)) == K & all(rowSums(pat) > 0 & rowSums(pat) < p)) {
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
      pat_lis_gene[[b]] <- pat
    }
    
    snp_list[[g]] <- snp_mat_gene
    hap_list[[g]] <- hap_mat_gene
    pat_list[[g]] <- pat_lis_gene
  }
  
  list(snp=snp_list, hap=hap_list, pat=pat_list)
}

# -----------------------------
# Step 2: gene expression
simulate_expression <- function(geno_snp, geno_hap, n_causal_snp, n_causal_hap, prop_hap_driven, h2_G) {
  n_genes <- length(geno_snp)
  n <- nrow(geno_snp[[1]])
  G_mat <- matrix(NA, n, n_genes)
  causal_info <- list()
  
  for (g in 1:n_genes) {
    geno_s <- scale(geno_snp[[g]])
    geno_h <- scale(geno_hap[[g]])
    
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
  G_mat <- scale(G_mat)
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
  pheno <- simulate_phenotype(expr$G, causal_genes, h2_Y)
  
  list(
    snp_geno = geno$snp,
    hap_geno = geno$hap,
    pat = geno$pat,
    expression = expr$G,
    causal_var = expr$causal_info,
    driver_type = driver_type,
    causal_genes = causal_genes,
    pheno = pheno$Y,
    alpha = pheno$alpha,
    sigma_eta2 = pheno$sigma_eta2
  )
}

# -----------------------------
# Step 5: train haplotype-based gene expression prediction model
train_pred_model <- function(n_genes, sim_data, idx_train) {
  cv_fits <- list()
  weights <- list()
  covariances <- list()
  
  for(g in seq_len(n_genes)) {
    x_train <- sim_data$hap_geno[[g]][idx_train, , drop=FALSE]
    y_train <- as.numeric(scale(sim_data$expression[idx_train, g]))
    
    cv_hap <- NULL
    beta_hap <- rep(0, ncol(x_train))
    try({
      set.seed(111 + g)
      cv_hap <- cv.glmnet(x_train, y_train, alpha = 0.5, nfolds = 5, standardize = TRUE)
      beta_hap <- as.numeric(coef(cv_hap, s = "lambda.min")[-1])
    }, silent = TRUE)
    
    cv_fits[[g]] <- cv_hap
    weights[[g]] <- beta_hap
    covariances[[g]] <- cor(sim_data$hap_geno[[g]])
  }
  
  list(cv_fits = cv_fits, weights = weights, covariances = covariances)
}

# -----------------------------
# Step 6: GWAS summary statistic
GWAS_summary <- function(geno, idx_test, Y_test) {
  results_list <- lapply(seq_along(geno), function(g) {
    geno_block <- geno[[g]][idx_test, , drop=FALSE]
    
    res <- data.frame(
      rs = paste0(g, "_", seq_len(ncol(geno_block))),
      beta = NA_real_,
      se = NA_real_,
      z = NA_real_,
      p = NA_real_,
      stringsAsFactors = FALSE
    )
    
    for (i in seq_len(ncol(geno_block))) {
      vi <- as.numeric(scale(geno_block[, i]))
      if (is.na(var(vi)) || var(vi) < 1e-8) next
      try({
        fit <- lm(Y_test ~ vi)
        s <- summary(fit)$coefficients
        if (nrow(s) > 1) {
          res[i, 2:5] <- s[2, ]
        }
      }, silent = TRUE)
    }
    
    return(res)
  })
  return(results_list)
}

# -----------------------------
# Step 7: individual & summary haplo-TWAS
individual_haplo_TWAS <- function(n_genes, sim_data, idx_test, model_train, Y_test) {
  results <- data.frame(
    gene = paste0("gene_", seq_len(n_genes)),
    beta = NA_real_,
    se = NA_real_,
    zscore = NA_real_,
    pvalue = NA_real_,
    stringsAsFactors = FALSE
  )
  
  for (g in seq_len(n_genes)) {
    if (is.null(model_train$cv_fits[[g]])) next
    Xh_test <- sim_data$hap_geno[[g]][idx_test, , drop=FALSE]
    
    try({
      pred_hap <- predict(model_train$cv_fits[[g]], Xh_test, s = "lambda.min")[,1]
      if (is.na(var(pred_hap)) || var(pred_hap) < 1e-8) next
      fit_hap <- lm(Y_test ~ scale(pred_hap))
      s <- summary(fit_hap)$coefficients
      if (nrow(s) > 1) {
        results[g, 2:5] <- s[2, ]
      }
    }, silent = TRUE)
  }
  
  return(results)
}

summary_haplo_TWAS <- function(n_genes, model_train, GWAS_res, n_test) {
  results <- data.frame(
    gene = paste0("gene_", seq_len(n_genes)),
    beta = NA_real_,
    se = NA_real_,
    zscore = NA_real_,
    pvalue = NA_real_,
    stringsAsFactors = FALSE
  )
  
  for (g in seq_len(n_genes)) {
    w <- model_train$weights[[g]]
    if (is.null(w) || all(w == 0) || all(is.na(w))) next
    GWAS <- GWAS_res[[g]]
    
    valid_idx <- !is.na(GWAS$z)
    if (sum(valid_idx) == 0) next
    w <- w[valid_idx]
    if (all(w == 0)) next
    GWAS <- GWAS[valid_idx, , drop = FALSE]
    if (!is.data.frame(GWAS) || ncol(GWAS) < 5) next
    sigma <- model_train$covariances[[g]][valid_idx, valid_idx, drop = FALSE]
    
    # lambda <- mean(diag(sigma)) * 1e-4
    # sigma_inv <- solve(sigma + diag(lambda, ncol(sigma)))
    numerator <- t(w) %*% GWAS[["z"]]
    denominator <- sqrt(t(w) %*% sigma %*% w)
    if (is.na(denominator) || denominator == 0) next
    zscore_g <- numerator / denominator
    pvalue_g <- 2 * pnorm(-abs(zscore_g))
    
    R2 <- zscore_g^2 / (zscore_g^2 + n_test -2)
    se_g <- sqrt((1 - R2) / (n_test - 2))
    beta_g <- zscore_g * se_g
    
    results[g, 2:5] <- c(beta_g, se_g, zscore_g, pvalue_g)
  }
  return(results)
}

# -----------------------------
# Step 8: summary statistic transform from SNP to haplotype
summary_transform <- function(n_genes, GWAS_snp, sim_data, n_test) {
  GWAS_hap_list <- list()
  
  for (g in seq_len(n_genes)) {
    GWAS_snp_gene <- GWAS_snp[[g]]
    pat_gene <- sim_data$pat[[g]]
    ld_gene <- cov(sim_data$snp_geno[[g]])
    GWAS_hap_gene <- list()
    
    for (b in seq_len(length(pat_gene))) {
      pat <- pat_gene[[b]]
      p <- ncol(pat)      # snps_per_block
      K <- nrow(pat)      # haps_per_block
      
      start <- (b - 1) * p + 1
      end <- b * p
      gwas <- GWAS_snp_gene[start : end, ]
      ld <- ld_gene[start : end, start : end]
      
      res <- data.frame(
        rs = paste0(g, "_", ((b - 1) * K + 1) : (b * K)),
        beta = NA_real_,
        se = NA_real_,
        z = NA_real_,
        p = NA_real_,
        stringsAsFactors = FALSE
      )
      
      valid_idx <- !is.na(gwas$z)
      if (sum(valid_idx) < p / 2) {
        GWAS_hap_gene[[b]] <- res
      } else {
        gwas <- gwas[valid_idx, ]
        H <- t(pat[ , valid_idx])
        L <- ld[valid_idx, valid_idx]
        
        tryCatch({
          V <- diag(gwas$se) %*% L %*% diag(gwas$se)
          W <- solve(V + diag(1e-6, ncol(V)))
          Ht_W_H <- t(H) %*% W %*% H
          Ht_W_H_inv <- solve(Ht_W_H + diag(1e-4, ncol(Ht_W_H)))
          
          numerator <- Ht_W_H_inv %*% t(H) %*% W %*% gwas$z
          denominator <- sqrt(diag(Ht_W_H_inv %*% t(H) %*% W %*% t(W) %*% H %*% Ht_W_H_inv))
          z_hap <- ifelse(denominator == 0, NA, numerator / denominator)
          p_hap <- 2 * pnorm(-abs(z_hap))
          
          R2 <- z_hap^2 / (z_hap^2 + n_test - 2)
          se_hap <- sqrt((1 - R2) / (n_test - 2))
          beta_hap <- z_hap * se_hap
          
          res <- data.frame(
            rs = paste0(g, "_", ((b - 1) * K + 1) : (b * K)),
            beta = as.numeric(beta_hap),
            se = as.numeric(se_hap),
            z = as.numeric(z_hap),
            p = as.numeric(p_hap),
            stringsAsFactors = FALSE
          )
          
          GWAS_hap_gene[[b]] <- res
        }, error = function(e) {GWAS_hap_gene[[b]] <- res})
      }
    }
    
    GWAS_hap_list[[g]] <- do.call(rbind, GWAS_hap_gene)
  }
  return(GWAS_hap_list)
}

# -----------------------------
# Step 9: metrics computation helpers
cor_lm <- function(x1, x2, x3, x4, y1, y2, y3, y4) {
  cor_beta <- cor.test(x1, y1)
  cor_se <- cor.test(x2, y2)
  cor_z <- cor.test(x3, y3)
  cor_p <- cor.test(x4, y4)
  
  res <- data.frame(
    comparison = c("beta", "se", "z", "p"),
    correlation = c(cor_beta$estimate, cor_se$estimate, cor_z$estimate, cor_p$estimate),
    p_value = c(cor_beta$p.value, cor_se$p.value, cor_z$p.value, cor_p$p.value),
    ci_lower = c(cor_beta$conf.int[1], cor_se$conf.int[1], cor_z$conf.int[1], cor_p$conf.int[1]),
    ci_upper = c(cor_beta$conf.int[2], cor_se$conf.int[2], cor_z$conf.int[2], cor_p$conf.int[2]),
    stringsAsFactors = FALSE
  )
  
  model_beta <- lm(x1 ~ y1)
  model_se <- lm(x2 ~ y2)
  model_z <- lm(x3 ~ y3)
  model_p <- lm(x4 ~ y4)
  
  res$slope <- c(coef(model_beta)[2], coef(model_se)[2], coef(model_z)[2], coef(model_p)[2])
  res$intercept <- c(coef(model_beta)[1], coef(model_se)[1], coef(model_z)[1], coef(model_p)[1])
  res$r2 <- c(summary(model_beta)$r.squared, summary(model_se)$r.squared, summary(model_z)$r.squared, summary(model_p)$r.squared)
  
  return(res)
}

ccc <- function(x, y) {
  mean_x <- mean(x)
  mean_y <- mean(y)
  var_x <- var(x)
  var_y <- var(y)
  cov_xy <- cov(x, y)
  
  rho <- cov_xy / sqrt(var_x * var_y)      # Pearson correlation
  ccc <- (2 * rho * sqrt(var_x) * sqrt(var_y)) / (var_x + var_y + (mean_x - mean_y)^2)
  
  return(ccc)
}

jac_kap_rho <- function(j, k) {
  q1 <- p.adjust(j, method="BH");      q2 <- p.adjust(k, method="BH")
  discovered1 <- q1 <= 0.05;           discovered2 <- q2 <= 0.05
  
  jaccard <- sum(discovered1 & discovered2) / sum(discovered1 | discovered2)
  if (all(!discovered1) || all(!discovered2)) {
    kappa <- 0
  } else {
    kappa <- cohen.kappa(table(discovered1, discovered2))$kappa
  }
  rho <- cor(-log10(q1), -log10(q2), method = "spearman")
  
  return(c(jaccard, kappa, rho))
}

compute_metrics <- function(TWAS_ind, TWAS_sum, TWAS_sum_trans, GWAS_hap, GWAS_hap_trans) {
  colnames(TWAS_ind)[-1] <- paste0(colnames(TWAS_ind)[-1], "_ind")
  colnames(TWAS_sum) <- paste0(colnames(TWAS_sum), "_sum")
  colnames(TWAS_sum_trans) <- paste0(colnames(TWAS_sum_trans), "_sum_trans")
  twas <- cbind(TWAS_ind, TWAS_sum[ ,-1], TWAS_sum_trans[ ,-1])
  t <- na.omit(twas)
  
  twas_res1 <- cor_lm(twas$beta_ind, twas$se_ind, twas$zscore_ind, twas$pvalue_ind, 
                      twas$beta_sum, twas$se_sum, twas$zscore_sum, twas$pvalue_sum)
  twas_res2 <- cor_lm(twas$beta_sum, twas$se_sum, twas$zscore_sum, twas$pvalue_sum, 
                      twas$beta_sum_trans, twas$se_sum_trans, twas$zscore_sum_trans, twas$pvalue_sum_trans)
  twas_res1$ccc <- c(ccc(t$beta_ind, t$beta_sum), ccc(t$se_ind, t$se_sum), 
                     ccc(t$zscore_ind, t$zscore_sum), ccc(t$pvalue_ind, t$pvalue_sum))
  twas_res2$ccc <- c(ccc(t$beta_sum, t$beta_sum_trans), ccc(t$se_sum, t$se_sum_trans), 
                     ccc(t$zscore_sum, t$zscore_sum_trans), ccc(t$pvalue_sum, t$pvalue_sum_trans))
  colnames(twas_res1)[-1] <- paste0(colnames(twas_res1)[-1], "_IS")
  colnames(twas_res2)[-1] <- paste0(colnames(twas_res2)[-1], "_SST")
  twas_res <- cbind(twas_res1, twas_res2[-1])
  
  jkr1 <- jac_kap_rho(t$pvalue_ind, t$pvalue_sum)
  jkr2 <- jac_kap_rho(t$pvalue_sum, t$pvalue_sum_trans)
  jkr <- c(jkr1, jkr2)
  
  
  gwas_hap <- do.call(rbind, GWAS_hap)
  gwas_hap_trans <- do.call(rbind, GWAS_hap_trans)
  colnames(gwas_hap_trans) <- paste0(colnames(gwas_hap_trans), "_trans")
  gwas <- cbind(gwas_hap, gwas_hap_trans[-1])
  g <- na.omit(gwas)
  
  gwas_res <- cor_lm(gwas$beta, gwas$se, gwas$z, gwas$p, gwas$beta_trans, gwas$se_trans, gwas$z_trans, gwas$p_trans)
  gwas_res$ccc <- c(ccc(g$beta, g$beta_trans), ccc(g$se, g$se_trans), ccc(g$z, g$z_trans), ccc(g$p, g$p_trans))
  gwas_res$rho <- c(cor(g$beta, g$beta_trans, method = "spearman"), 
                    cor(g$se, g$se_trans, method = "spearman"), 
                    cor(g$z, g$z_trans, method = "spearman"), 
                    cor(-log10(p.adjust(g$p, method = "BH")), -log10(p.adjust(g$p_trans, method = "BH")), method = "spearman"))
  
  
  list(TWAS = twas, TWAS_RES = twas_res, JKR = jkr, GWAS = gwas, GWAS_RES = gwas_res)
}

# -----------------------------
# Step 10: run one replicate simulation
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
  
  idx_all <- sample(1:n)
  idx_train <- idx_all[1:n_train]
  idx_test <- idx_all[(n_train+1):n]
  Y_test <- as.numeric(scale(sim_data$pheno[idx_test]))
  
  model_train <- train_pred_model(n_genes, sim_data, idx_train)
  
  TWAS_ind <- individual_haplo_TWAS(n_genes, sim_data, idx_test, model_train, Y_test)
  
  GWAS_hap <- GWAS_summary(sim_data$hap_geno, idx_test, Y_test)
  TWAS_sum <- summary_haplo_TWAS(n_genes, model_train, GWAS_hap, n_test)
  
  GWAS_snp <- GWAS_summary(sim_data$snp_geno, idx_test, Y_test)
  GWAS_hap_trans <- summary_transform(n_genes, GWAS_snp, sim_data, n_test)
  TWAS_sum_trans <- summary_haplo_TWAS(n_genes, model_train, GWAS_hap_trans, n_test)
  
  
  results <- compute_metrics(TWAS_ind, TWAS_sum, TWAS_sum_trans, GWAS_hap, GWAS_hap_trans)
  return(results)
}

# -----------------------------
# Step 11: multiple replicates
run_replicates <- function(params, parallel = TRUE) {
  nrep <- params$nrep
  n_cores <- params$n_cores
  
  seeds <- sample.int(1e6, nrep)
  res_list <- vector("list", nrep)
  
  if(parallel) {
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    res_list <- foreach::foreach(i = 1 : nrep,
                                 .packages = c("MASS","gtools","glmnet","psych"),
                                 .export = c("simulate_snp_hap","simulate_expression","simulate_phenotype","simulate_dataset",
                                             "train_pred_model","GWAS_summary","individual_haplo_TWAS","summary_haplo_TWAS",
                                             "summary_transform","cor_lm","ccc","jac_kap_rho","compute_metrics","run_one_replicate")) %dopar% {
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
# Step 12: summarize replicate results
summarize_replicates <- function(res_list, output) {
  seeds <- do.call(rbind, lapply(seq_along(res_list), function(i) {
    data.frame(replicate = i, seed = res_list[[i]]$seed)
  }))
  results_twas <- do.call(rbind, lapply(seq_along(res_list), function(i) {
    df <- data.frame(res_list[[i]]$result$TWAS)
    df$replicate <- i
    df
  }))
  results_gwas <- do.call(rbind, lapply(seq_along(res_list), function(i) {
    df <- data.frame(res_list[[i]]$result$GWAS)
    df$replicate <- i
    df
  }))
  
  metrics_twas <- do.call(rbind, lapply(seq_along(res_list), function(i) {
    df <- data.frame(res_list[[i]]$result$TWAS_RES)
    df$replicate <- i
    df
  }))
  metrics_gwas <- do.call(rbind, lapply(seq_along(res_list), function(i) {
    df <- data.frame(res_list[[i]]$result$GWAS_RES)
    df$replicate <- i
    df
  }))
  
  jkr_twas <- data.frame(do.call(rbind, lapply(res_list, function(i) i$result$JKR)))
  colnames(jkr_twas) <- c("jaccard_IS", "kappa_IS", "rho_IS", "jaccard_SST", "kappa_SST", "rho_SST")
  jkr_twas$replicate <- seq_along(res_list)
  
  
  metrics_twas$comparison <- factor(metrics_twas$comparison, levels = c("beta", "se", "z", "p"))
  metrics_summary_twas <- metrics_twas %>%
    group_by(comparison) %>%
    summarise(across(
      .cols = c(correlation_IS, p_value_IS, ci_lower_IS, ci_upper_IS, slope_IS, intercept_IS, r2_IS, ccc_IS, 
                correlation_SST, p_value_SST, ci_lower_SST, ci_upper_SST, slope_SST, intercept_SST, r2_SST, ccc_SST),
      .fns = list(mean = \(x) mean(x, na.rm = TRUE), sd = \(x) sd(x, na.rm = TRUE)),
      .names = "{.col}_{.fn}"
    ))
  
  metrics_gwas$comparison <- factor(metrics_gwas$comparison, levels = c("beta", "se", "z", "p"))
  metrics_summary_gwas <- metrics_gwas %>%
    group_by(comparison) %>%
    summarise(across(
      .cols = c(correlation, p_value, ci_lower, ci_upper, slope, intercept, r2, ccc, rho),
      .fns = list(mean = \(x) mean(x, na.rm = TRUE), sd = \(x) sd(x, na.rm = TRUE)),
      .names = "{.col}_{.fn}"
    ))
  
  mean_jkr <- jkr_twas %>%
    summarise(across(-replicate, \(x) mean(x, na.rm = TRUE))) %>%
    mutate(replicate = "mean")
  sd_jkr <- jkr_twas %>%
    summarise(across(-replicate, \(x) sd(x, na.rm = TRUE))) %>%
    mutate(replicate = "sd")
  jkr_summary_twas <- rbind(mean_jkr, sd_jkr, jkr_twas)
  
  
  write.csv(seeds, paste0(output, ".seeds.csv"), row.names = FALSE)
  write.csv(results_twas, paste0(output, ".results_twas.csv"), row.names = FALSE)
  write.csv(results_gwas, paste0(output, ".results_gwas.csv"), row.names = FALSE)
  
  write.csv(metrics_twas, paste0(output, ".metrics_twas.csv"), row.names = FALSE)
  write.csv(metrics_gwas, paste0(output, ".metrics_gwas.csv"), row.names = FALSE)
  
  write.csv(metrics_summary_twas, paste0(output, ".metrics_summary_twas.csv"), row.names = FALSE)
  write.csv(metrics_summary_gwas, paste0(output, ".metrics_summary_gwas.csv"), row.names = FALSE)
  write.csv(jkr_summary_twas, paste0(output, ".jkr_twas.csv"), row.names = FALSE)
}

# -----------------------------
# Step 13: full pipeline wrapper
run_full_simulation <- function(params, output) {
  res_list <- run_replicates(params)
  
  summarize_replicates(res_list, output)
}


# -----------------------------
# Running the parameter combination
base_params <- list(
  n = 1500,                   # n = n_train + n_test
  n_train = 500,
  n_test = 1000,              # 1000 5000 10000
  
  h2_Y = 0.6,                 # 0.1 0.3 0.6
  h2_G = 0.3,                 # 0.1 0.3
  n_genes = 50,
  n_causal_gene = 5,
  n_causal_snp = 10,
  n_causal_hap = 10,
  
  snps_per_gene = 200,        # 200 500
  snps_per_block = 5,         # 5 10 20
  n_blocks_per_gene = 40,     # n_blocks_per_gene = snps_per_gene / snps_per_block  40 20 10 / 100 50 25
  haps_per_block = 6,         # 6 / 10 / 16  即 1.2 / 1 / 0.8 ( haps_per_block / snps_per_block )
  
  hap_freq_conc = 1,
  prop_hap_driven = 0.5,
  nrep = 200,
  n_cores = 20
)

param_grid <- expand.grid(
  h2_Y = c(0.1, 0.3, 0.6),
  h2_G = c(0.1, 0.3),
  snps_per_gene = c(200, 500),
  snps_per_block = c(5, 10, 20),
  stringsAsFactors = FALSE
)
param_grid$n_blocks_per_gene <- param_grid$snps_per_gene / param_grid$snps_per_block
param_grid$haps_per_block <- ifelse(param_grid$snps_per_block == 5, 6,
                                    ifelse(param_grid$snps_per_block == 10, 10, 16))
param_grid$sim <- seq_len(nrow(param_grid))
write.csv(param_grid, "./simulation/result2/param_grid.csv", row.names = FALSE)

modify_params <- function(base, new_values) {
  modifyList(base, as.list(new_values))
}

for (i in seq_len(nrow(param_grid))) {
  current_params <- modify_params(base_params, list(h2_Y = param_grid$h2_Y[i],
                                                    h2_G = param_grid$h2_G[i],
                                                    snps_per_gene = param_grid$snps_per_gene[i],
                                                    snps_per_block = param_grid$snps_per_block[i],
                                                    n_blocks_per_gene = param_grid$n_blocks_per_gene[i],
                                                    haps_per_block = param_grid$haps_per_block[i]))
  output_prefix <- paste0("./simulation/result2/500train_1000test/sim",i)
  
  message("=== Start running the parameter combination  ", i, " / ", nrow(param_grid), " ===")
  set.seed(69410)
  run_full_simulation(current_params, output_prefix)
}

