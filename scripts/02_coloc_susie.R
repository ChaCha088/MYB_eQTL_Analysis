# ==============================================================================
# Title: Automated Colocalization Analysis (Metadata-driven)
# Description: Performs fine-mapping and colocalization using a metadata file.
# Dependencies: coloc, susieR, dplyr, data.table
# ==============================================================================

suppressMessages({
  library(coloc)
  library(susieR)
  library(dplyr)
  library(data.table)
})

# 1.  Preferences ----------------------------------------------------------
setwd("/Filepath/coloc") 
metadata_path <- "metadata.csv"  # 메타데이터 파일 경로
output_file   <- "Coloc_SuSiE_Results_Final.csv"

# LD related files
ld_file  <- "chr*_ld.ld"
snp_file <- "chr*_ld.snplist"
bim_file <- "../GW.E-GEUV-3.EUR.MAF005.HWE1e-06.bim"

#2. Defining a helper function ------------------------------------------------------------

# [Maintaining existing functions] Data cleaning and formatting
clean_and_format <- function(file_path, type="GWAS") {
  if (!file.exists(file_path)) stop(paste("File not found:", file_path))
  dt <- fread(file_path, header=TRUE)
  
  if (type == "GWAS") {
    df <- dt %>% select(snp = rsID, beta = Beta, se = SE, maf = EAF, n = any_of("Neff"), p = Pval, A1 = EffectAllele, A2 = NonEffectAllele)
  } else {
    cols <- names(dt)
    get_col <- function(candidates) {
      match <- intersect(candidates, cols); if(length(match) > 0) return(match[1]) else return(NULL)
    }
    df <- dt %>% select(
      snp = all_of(get_col(c("rsid", "SNP", "snp", "variant_id"))),
      beta = all_of(get_col(c("beta", "b", "Beta"))),
      se = all_of(get_col(c("se", "SE", "standard_error"))),
      maf = all_of(get_col(c("Freq", "maf", "EAF"))),
      p = all_of(get_col(c("p", "P", "pval", "P-value"))),
      A1 = all_of(get_col(c("A1", "effect_allele", "ref"))),
      A2 = all_of(get_col(c("A2", "other_allele", "alt")))
    )
  }
  df %>% mutate(across(c(beta, se, maf, p), as.numeric)) %>% filter(!is.na(beta), !is.na(se), !is.na(maf))
}

# [Keep existing functions] Allele Harmonization
harmonize_data <- function(raw_data, ref_bim, N_input) {
  target <- raw_data %>% filter(snp %in% ref_bim$snp)
  ref    <- ref_bim  %>% filter(snp %in% target$snp)
  target <- target[match(ref$snp, target$snp), ]
  if(nrow(target) == 0) return(NULL)
  
  is_match <- (target$A1 == ref$A1_ld & target$A2 == ref$A2_ld) | (target$A1 == ref$A2_ld & target$A2 == ref$A1_ld)
  keep_idx <- which(is_match)
  if(length(keep_idx) == 0) return(NULL)
  
  target <- target[keep_idx, ]; ref <- ref[keep_idx, ]
  to_flip <- target$A1 == ref$A2_ld
  target$beta[to_flip] <- -target$beta[to_flip]
  
# N processing: If it is a column name, it is retrieved from the data, and if it is a number, it is used as is.
  n_scalar <- if (N_input %in% names(target)) median(target[[N_input]], na.rm=TRUE) else as.numeric(N_input)
  
  return(list(beta = target$beta, se = target$se, n = n_scalar, snp = target$snp))
}

# 3. main analysis process -------------------------------------------------------

# metadata load
if(!file.exists(metadata_path)) stop("The metadata file does not exist.")
meta <- fread(metadata_path)

# Load LD Matrix (only done once)
cat("Loading LD Matrix...\n")
ld_snps_all <- fread(snp_file, header=FALSE)$V1
ld_mat_all  <- as.matrix(fread(ld_file, header=FALSE))
colnames(ld_mat_all) <- ld_snps_all
rownames(ld_mat_all) <- ld_snps_all
bim_df <- fread(bim_file, col.names=c("chr", "snp", "cm", "bp", "A1_ld", "A2_ld")) %>% 
  select(snp, A1_ld, A2_ld) %>% filter(snp %in% ld_snps_all)

results_list <- list()

# Loop through each row of metadata
for (i in 1:nrow(meta)) {
  row <- meta[i, ]
  cat(sprintf("\n[%d/%d] Processing: %s vs %s\n", i, nrow(meta), row$trait_label, row$gene_label))
  
  # load data
  d1_raw <- clean_and_format(row$gwas_file, "GWAS")
  d2_raw <- clean_and_format(row$eqtl_file, "eQTL")
  
  # common SNP extraction
  common_snps <- intersect(ld_snps_all, intersect(d1_raw$snp, d2_raw$snp))
  if (length(common_snps) < 50) { cat("  [SKIP] Too few SNPs.\n"); next }
  ref_sub <- bim_df %>% filter(snp %in% common_snps)
  
  # Sorting and filtering
  D1 <- harmonize_data(d1_raw, ref_sub, row$gwas_n)
  D2 <- harmonize_data(d2_raw, ref_sub, row$eqtl_n)
  
  if (is.null(D1) || is.null(D2)) { cat("  [SKIP] Harmonization failed.\n"); next }
  
  final_snps <- intersect(D1$snp, D2$snp)
  idx1 <- match(final_snps, D1$snp); idx2 <- match(final_snps, D2$snp)
  
  # Prepare data for analysis
  LD_final <- ld_mat_all[final_snps, final_snps]
  diag(LD_final) <- diag(LD_final) + 1e-4 # LD Regularization
  
# SuSiE execution function
  run_susie <- function(b, s, n, R) {
    tryCatch(susie_rss(bhat = b, shat = s, R = R, n = n, L = 10, check_R = TRUE), error = function(e) NULL)
  }
  
  fit1 <- run_susie(D1$beta[idx1], D1$se[idx1], D1$n, LD_final)
  fit2 <- run_susie(D2$beta[idx2], D2$se[idx2], D2$n, LD_final)
  
  if (is.null(fit1) || is.null(fit2)) { cat("  [FAIL] SuSiE convergence failed.\n"); next }
  
  # Run coloc
  res <- coloc.susie(fit1, fit2)
  
  if (!is.null(res$summary) && nrow(res$summary) > 0) {
    summ <- res$summary
    summ$Trait_Label <- row$trait_label
    summ$Gene_Label  <- row$gene_label
    results_list[[length(results_list) + 1]] <- summ
    cat(sprintf("  Success! Max PP.H4: %.4f\n", max(summ$PP.H4.abf)))
  }
}

# 4. save results -----------------------------------------------------------------
if (length(results_list) > 0) {
  final_df <- do.call(rbind, results_list) %>% arrange(desc(PP.H4.abf))
  write.csv(final_df, output_file, row.names = FALSE)
  cat(sprintf("\nDone! Results saved to %s\n", output_file))
} else {
  cat("\nNo significant results found.\n")
}
