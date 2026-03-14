# ---- library R package ----
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(future.apply)
  library(optparse)
  library(MASS)
})


# ---- Retrieve command-line parameters ----
option_list <- list(
  make_option(c("--snp_gwas"), type = "character", help = "Path to SNP GWAS summary statistic file"),
  make_option(c("--n_test"), type = "integer", help = "GWAS sample size"),
  make_option(c("--species"), type = "character", help = "Species (human, pig, cattle, or chicken)"),
  make_option(c("--tissue"), type = "character", help = "Tissue (capitalize the first letter)"),
  make_option(c("--chr"), type = "character", help = "Chromosome ID to process"),
  make_option(c("--sw"), type = "integer", default = 2, help = "Haplotype block size"),
  make_option(c("--threads"), type = "integer", default = 1, help = "Number of parallel threads invoked"),
  make_option(c("--out"), type = "character", help = "Path to output file"),
  
  make_option(c("--col_chr"), type = "character", help = "Column name for chromosome ID in GWAS file"),
  make_option(c("--col_rs"), type = "character", help = "Column name for SNP ID (rsid) in GWAS file"),
  make_option(c("--col_ps"), type = "character", help = "Column name for physical position (ps) in GWAS file"),
  make_option(c("--col_effect_allele"), type = "character", help = "Column name for effect allele in GWAS file"),
  make_option(c("--col_non_effect_allele"), type = "character", help = "Column name for non-effect allele in GWAS file"),
  make_option(c("--col_beta"), type = "character", default = "beta", help = "Column name for beta coefficient in GWAS file"),
  make_option(c("--col_se"), type = "character", default = "se", help = "Column name for standard error in GWAS file"),
  make_option(c("--col_z"), type = "character", default = "z", help = "Column name for Z-score in GWAS file"),
  make_option(c("--col_p"), type = "character", default = "p", help = "Column name for P-value in GWAS file")
)
opt <- parse_args(OptionParser(option_list = option_list))

n_test <- opt$n_test
species <- opt$species
tissue <- opt$tissue
chr <- opt$chr
sw <- opt$sw

col_map <- c(
  "chr" = opt$col_chr, "rs" = opt$col_rs, "ps" = opt$col_ps,
  "allele1" = opt$col_effect_allele, "allele0" = opt$col_non_effect_allele,
  "beta" = opt$col_beta, "se" = opt$col_se, "z" = opt$col_z, "p" = opt$col_p
)


# ---- Read the input files ----
gwas <- data.frame(fread(opt$snp_gwas))

essential_vars <- c("chr", "rs", "ps", "allele1", "allele0")
missing_essentials <- setdiff(col_map[essential_vars], names(gwas))
if(length(missing_essentials) > 0) stop("Error: Please provide all required columns exactly: --col_chr, --col_rs, --col_ps, --col_effect_allele, --col_non_effect_allele ! ! !")

gwas_chr <- gwas[gwas[[opt$col_chr]] == chr, ]
gwas_chr[[opt$col_chr]] <- as.numeric(gwas_chr[[opt$col_chr]])

actual_map <- col_map[col_map %in% names(gwas_chr)]
gwas_sub <- gwas_chr[, actual_map, drop = FALSE]
names(gwas_sub) <- names(actual_map)
gwas_chr <- gwas_sub

if (all(c("se", "z") %in% names(gwas_chr))) {
  print("Start to transform ! ! !")
} else if (all(c("beta", "se") %in% names(gwas_chr))) {
  gwas_chr$z <- gwas_chr$beta / gwas_chr$se
  print("Start to transform ! ! !")
} else if (all(c("beta", "z") %in% names(gwas_chr))) {
  gwas_chr$se <- gwas_chr$beta / gwas_chr$z
  print("Start to transform ! ! !")
} else if (all(c("beta", "p") %in% names(gwas_chr))) {
  safe_p <- pmax(gwas_chr$p, 1e-300)
  z_abs <- qnorm(1 - safe_p / 2)
  gwas_chr$z <- sign(gwas_chr$beta) * z_abs
  gwas_chr$se <- gwas_chr$beta / gwas_chr$z
  print("Start to transform ! ! !")
} else {
  stop("Error: GWAS file is missing necessary information (z se / beta se / beta z / beta p) ! ! !")
}

panel_chr <- data.frame(fread(paste0("./Pro1_liucheng_HapTWAS/", species, "GTEx/PredictDB/SNP/tissue_snp_annotation/", tissue, "/", tissue, ".snp_annot.chr", chr, ".txt")))
gwas_chr_full <- panel_chr %>%
  filter(!(pos %in% gwas_chr$ps)) %>%
  dplyr::select(chromosome, varID, pos, alt_vcf, ref_vcf) %>%
  dplyr::rename(chr = chromosome, rs = varID, ps = pos, allele1 = alt_vcf, allele0 = ref_vcf) %>%
  bind_rows(gwas_chr)
gwas_chr_full <- gwas_chr_full[match(panel_chr$pos, gwas_chr_full$ps), ]
gwas_chr_full$rs <- panel_chr$varID

ld_chr <- data.frame(fread(paste0("./Pro1_liucheng_HapTWAS/LD_panel/", species, "_", tissue, "/chr", chr, ".ld")))

hap_pat <- expand.grid(rep(list(c(1, 0)), sw)) %>% as.matrix()


# ---- Set up parallelism and logging ----
plan(multisession, workers = opt$threads)
options(future.globals.maxSize = 2 * 1024^3)
log_file <- sub("\\.txt$", ".log", opt$out)
writeLines(paste0("Haplotype summary statistic construction started at ", round(Sys.time())), log_file)


# ---- Transform summary statistic from SNP to HAP ----
construct_haplo_summary <- function(gwas_b, hap_pat, ld_mat, bid, n_test) {
  valid_idx <- !is.na(gwas_b$z) & !is.na(gwas_b$se)
  gwas_b <- gwas_b[valid_idx, ]
  H <- t(hap_pat[ , valid_idx])
  L <- ld_mat[valid_idx, valid_idx]
  
  tryCatch({
    V <- diag(gwas_b$se) %*% L %*% diag(gwas_b$se)
    W <- tryCatch(solve(V + diag(1e-6, ncol(V))), error = function(e) NULL)
    Ht_W_H <- t(H) %*% W %*% H
    Ht_W_H_inv <- tryCatch(solve(Ht_W_H + diag(1e-4, ncol(Ht_W_H))), error = function(e) NULL)
    
    numerator <- Ht_W_H_inv %*% t(H) %*% W %*% gwas_b$z
    denominator <- sqrt(diag(Ht_W_H_inv %*% t(H) %*% W %*% t(W) %*% H %*% Ht_W_H_inv))
    z_hap <- ifelse(denominator == 0, NA, numerator / denominator)
    p_hap <- 2 * pnorm(-abs(z_hap))
    
    R2 <- z_hap^2 / (z_hap^2 + n_test - 2)
    se_hap <- sqrt((1 - R2) / (n_test - 2))
    beta_hap <- z_hap * se_hap
  }, error = function(e) {return(NULL)})
  
  map_df(1:nrow(hap_pat), function(i) {
    hap <- hap_pat[i, ]
    hap_code <- paste0(ifelse(hap == 1, gwas_b$allele1, gwas_b$allele0), collapse = "")
    ps_mid <- floor((min(gwas_b$ps) + max(gwas_b$ps)) / 2)
    
    tibble(
      chr = as.numeric(unique(gwas_b$chr)),
      rs = paste0("CHR", unique(gwas_b$chr), "_B", bid, "_", min(gwas_b$ps), "_", max(gwas_b$ps), "_", hap_code),
      ps = ps_mid,
      allele1 = "H",
      allele0 = "N",
      beta = beta_hap[i],
      se = se_hap[i],
      z_score = z_hap[i],
      p_value = p_hap[i]
    )
  })
}


# ---- Process each sliding window block in parallel ----
block_indices <- seq(1, nrow(panel_chr), by = sw)
if (tail(block_indices, 1) + sw - 1 > nrow(panel_chr)) {
  block_indices <- block_indices[-length(block_indices)]
}

haplotype_summary_statistic <- future_lapply(seq_along(block_indices), function(bid) {
  result <- tryCatch({
    start <- block_indices[bid]
    end <- start + sw - 1
    
    gwas_b <- gwas_chr_full[start : end, ]
    valid_idx <- !is.na(gwas_b$z) & !is.na(gwas_b$se)
    if (sum(valid_idx) <= sw / 2) {
      log_line <- paste0(round(Sys.time()), " - Block ", bid, " / ", length(block_indices), " : ", gwas_b$ps[1], gwas_b$ps[sw], " - skipped (too few valid SNPs: ", sum(valid_idx), "/", sw, ")")
      write(log_line, log_file, append = TRUE)
      return(NULL)
    }
    
    ld <- ld_chr %>% filter(SNP_A %in% gwas_b$rs & SNP_B %in% gwas_b$rs)
    snp_ids <- unique(c(ld$SNP_A, ld$SNP_B))
    ld_mat <- matrix(0, nrow = sw, ncol = sw, dimnames = list(snp_ids, snp_ids))
    for (i in seq_len(nrow(ld))) {
      a <- ld$SNP_A[i]; b <- ld$SNP_B[i]; r <- ld$R[i]
      ld_mat[a, b] <- r
      ld_mat[b, a] <- r
    }
    diag(ld_mat) <- 1
    
    log_line <- paste0(round(Sys.time()), " - Block ", bid, " / ", length(block_indices), " : ", gwas_b$ps[1], gwas_b$ps[sw], " - processed successfully")
    write(log_line, log_file, append = TRUE)
    
    construct_haplo_summary(gwas_b, hap_pat, ld_mat, bid, n_test)
  }, error = function(e) {
    log_line <- paste0(round(Sys.time()), " - Block ", bid, " / ", length(block_indices), " : ", gwas_b$ps[1], gwas_b$ps[sw], " - ld_mat failed - ", conditionMessage(e))
    write(log_line, log_file, append = TRUE)
    return(NULL)
  })
  return(result)
}) %>% bind_rows() %>% filter(!is.na(p_value) & p_value < 1)


# ---- Output the transformation results ----
fwrite(haplotype_summary_statistic, opt$out, quote = FALSE, sep = "\t")
write(paste0("Finished all ", length(block_indices), " blocks on CHR ", chr, " at ", round(Sys.time())), log_file, append = TRUE)

