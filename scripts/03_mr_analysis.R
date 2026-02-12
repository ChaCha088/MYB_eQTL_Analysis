# ==============================================================================
# Title: Automated Robust Two-sample MR (Final Version)
# Features: Local Clumping, Custom Column Matching, Safe Visualization
# ==============================================================================

# 1. Load Packages and Set Environment -----------------------------------------
suppressMessages({
  library(TwoSampleMR)
  library(MRPRESSO)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ieugwasr) 
})

# Set working directory (User should adjust this path)
setwd("/File_path/mr_analysis")
metadata_path <- "metadata.csv"

# Create output directories if they do not exist
if(!dir.exists("mr_results")) dir.create("mr_results")
if(!dir.exists("mr_plots")) dir.create("mr_plots")

# 2. Load Metadata -------------------------------------------------------------
if(!file.exists(metadata_path)) stop("Metadata file not found.")
meta <- fread(metadata_path)

all_summary_results <- list()

# 3. Start Analysis Loop -------------------------------------------------------
for (i in 1:nrow(meta)) {
  row <- meta[i, ]
  prefix <- paste0(row$exp_label, "_vs_", row$out_label)
  
  cat(sprintf("\n[%d/%d] Working on: %s\n", i, nrow(meta), prefix))
  
  # --- A. Exposure Loading & F-stat ---
  exp_dat <- read_exposure_data(
    filename = row$exp_file, sep = "\t",
    snp_col = "SNP", beta_col = "b", se_col = "se",
    effect_allele_col = "A1", other_allele_col = "A2",
    eaf_col = "Freq", pval_col = "p"
  )
  
  # QC: P < 5e-8 & F-statistic
  exp_dat <- subset(exp_dat, pval.exposure < 5e-8)
  if(nrow(exp_dat) == 0) { cat("  [SKIP] No significant SNPs.\n"); next }
  
  exp_dat$F_stat <- (exp_dat$beta.exposure^2) / (exp_dat$se.exposure^2)
  avg_f <- mean(exp_dat$F_stat)
  
  # --- B. Local Clumping (using PLINK) ---
  cat("  Performing local clumping using PLINK...\n")
  
  # Backup original data
  exp_backup <- exp_dat 
  
  # Create dataframe for ld_clump input
  clump_input <- data.frame(
    rsid = exp_dat$SNP, 
    pval = exp_dat$pval.exposure, 
    id = exp_dat$id.exposure
  )
  
  # Execute local PLINK
  # [NOTE] Please update 'plink_bin' and 'bfile' paths according to your system environment.
  exp_dat_clumped <- ld_clump(
    dat = clump_input,
    clump_kb = 10000,
    clump_r2 = 0.001,
    plink_bin = "Filepath/bin/plink", 
    bfile = "GW.E-GEUV-3.EUR.MAF005.HWE1e-06"
  )
  
  # Keep only clumped SNPs
  exp_dat <- exp_backup[exp_backup$SNP %in% exp_dat_clumped$rsid, ]
  cat(paste("  Instruments after clumping:", nrow(exp_dat), "\n"))
  
  # --- C. Outcome Loading (Smart Auto-Matching) ---
  cat("  Reading Outcome header to detect column names...\n")
  
  # 1. Read only the header to identify column names
  header <- names(fread(row$out_file, nrows = 0))
  
  # 2. Smart matching function: Find existing columns from candidates
  find_col <- function(candidates, file_header) {
    # Check match ignoring case
    match <- intersect(tolower(candidates), tolower(file_header))
    if (length(match) > 0) {
      # Return the actual column name (first match)
      return(file_header[which(tolower(file_header) == match[1])])
    }
    return(NULL)
  }
  
  # 3. Define priority candidates (Left is higher priority)
  # Include both "hm_" prefixed and non-prefixed versions.
  col_snp  <- find_col(c("rsid", "hm_rsid", "snp", "markername", "variant_id"), header)
  col_beta <- find_col(c("beta", "hm_beta", "effect", "log_or", "b"), header)
  col_se   <- find_col(c("standard_error", "hm_standard_error", "se", "std_error"), header)
  col_a1   <- find_col(c("effect_allele", "hm_effect_allele", "a1", "allele1", "reference_allele"), header)
  col_a2   <- find_col(c("other_allele", "hm_other_allele", "a2", "allele2", "other_allele"), header)
  col_eaf  <- find_col(c("effect_allele_frequency", "hm_effect_allele_frequency", "eaf", "freq", "a1_freq"), header)
  col_pval <- find_col(c("p_value", "hm_p_value", "pval", "p", "p-value"), header)
  
  # 4. Skip if critical columns are missing (Error prevention)
  if (is.null(col_snp) || is.null(col_beta) || is.null(col_se) || is.null(col_pval)) {
    cat("  [SKIP] Critical columns missing in this file. Check header manually.\n")
    cat(paste("  Found headers:", paste(header, collapse=", "), "\n"))
    next
  }
  
  cat(sprintf("  Detected columns: SNP=%s, Beta=%s, P=%s\n", col_snp, col_beta, col_pval))
  
  # 5. Load Data
  out_dat <- read_outcome_data(
    snps = exp_dat$SNP, 
    filename = row$out_file, 
    sep = "\t",
    snp_col = col_snp,
    beta_col = col_beta,
    se_col = col_se,
    effect_allele_col = col_a1,
    other_allele_col = col_a2,
    eaf_col = col_eaf,
    pval_col = col_pval
  )
  
  if(is.null(out_dat) || nrow(out_dat) == 0) { cat("  [SKIP] No matching SNPs in Outcome.\n"); next }
  
  # --- D. Harmonization & Steiger Filtering ---
  dat <- harmonise_data(exp_dat, out_dat, action = 2)
  dat$samplesize.exposure <- as.numeric(row$n_exp)
  dat$samplesize.outcome  <- as.numeric(row$n_out)
  
  # Steiger Filtering
  dat_steiger <- steiger_filtering(dat)
  dat <- subset(dat_steiger, steiger_dir == TRUE & steiger_pval < 0.05)
  
  if(nrow(dat) < 1) { cat("  [SKIP] No SNPs passed Steiger Filtering.\n"); next }
  
  # --- E. MR Analysis & Sensitivity Analysis ---
  res <- mr(dat)
  n_snps <- nrow(dat)
  cat(sprintf("  Number of SNPs for MR: %d\n", n_snps))
  
  final_res <- res # Set default value
  
  if (n_snps >= 3) {
    # >= 3 SNPs: Perform all analyses including MR-PRESSO
    het <- mr_heterogeneity(dat)
    pleio <- mr_pleiotropy_test(dat)
    final_res <- combine_all_mrresults(res, het, pleio)
    
    try({
      preso <- run_mr_presso(
        BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
        SdOutcome = "se.outcome", SdExposure = "se.exposure", 
        OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, 
        NbDistribution = 1000, SignifThreshold = 0.05
      )
      capture.output(preso, file = paste0("mr_results/", prefix, "_MR_PRESSO.txt"))
    }, silent = TRUE)
    
  } else if (n_snps == 2) {
    # 2 SNPs: Perform heterogeneity test only
    het <- mr_heterogeneity(dat)
    # final_res <- combine_all_mrresults(res, het) 
    # Note: combine_all_mrresults might fail if pleio is missing, so using basic results + manual check recommended.
    
  } else {
    # 1 SNP: Perform Wald Ratio only
    cat("  [Info] Only 1 SNP. Using Wald Ratio.\n")
  }
  
  # Add common info and save
  final_res$Avg_F_stat <- avg_f
  final_res$Exposure_Label <- row$exp_label
  final_res$Outcome_Label  <- row$out_label
  write.csv(final_res, paste0("mr_results/", prefix, "_Main_Results.csv"), row.names = F)
  all_summary_results[[i]] <- final_res
  
  # --- F. Visualization (Safe Mode) ---
  cat("  Generating Plots...\n")
  
  # 1. Scatter Plot (Always available)
  p1 <- mr_scatter_plot(res, dat)[[1]]
  ggsave(paste0("mr_plots/", prefix, "_1_Scatter.pdf"), p1, width=7, height=7)
  
  # 2. Forest Plot (Always available)
  res_single <- mr_singlesnp(dat)
  p2 <- mr_forest_plot(res_single)[[1]]
  ggsave(paste0("mr_plots/", prefix, "_2_Forest.pdf"), p2, width=7, height=10)
  
  # 3. & 4. Leave-one-out & Funnel (Only if > 1 SNP)
  if (n_snps > 1) {
    res_loo <- mr_leaveoneout(dat)
    p3 <- mr_leaveoneout_plot(res_loo)[[1]]
    ggsave(paste0("mr_plots/", prefix, "_3_LeaveOneOut.pdf"), p3, width=7, height=10)
    
    p4 <- mr_funnel_plot(res_single)[[1]]
    ggsave(paste0("mr_plots/", prefix, "_4_Funnel.pdf"), p4, width=7, height=7)
  } else {
    cat("  [Info] Skipping Leave-one-out & Funnel plots (need > 1 SNP).\n")
  }
}

# 4. Save Summary of Results ---------------------------------------------------
if(length(all_summary_results) > 0) {
  total_summary <- do.call(rbind, all_summary_results)
  write.csv(total_summary, "Total_MR_Summary_Report.csv", row.names = F)
  cat("\nAll analyses finished successfully!\n")
}
