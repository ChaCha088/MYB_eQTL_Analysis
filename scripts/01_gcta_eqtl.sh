#!/bin/bash
# ==============================================================================
# Title: GCTA-MLMA cis-eQTL Analysis Pipeline
# Description: Primary and Conditional MLMA analysis for MYB family and IFT52.
# ==============================================================================

1. Setting global variables (common use files)
BFILE="GW.E-GEUV-3.EUR.MAF005.HWE1e-06"
GRM="GW.E-GEUV-3.EUR.MAF005.HWE1e-06"
THREADS=10

echo "Starting GCTA analysis workflow..."

#2. Primary cis-eQTL analysis
for GENE in MYBL1 MYBL2 MYB IFT52
do
    echo "Processing primary eQTL for $GENE..."
    gcta64 --bfile $BFILE --mlma --grm $GRM --pheno ${GENE}.phen --thread-num $THREADS --out ${GENE}_Primary_Result
done

#3. Conditional Analysis
echo "Processing conditional analyses..."

# MYBL2: SNP (rs142261424) correction and IFT52 expression level correction
gcta64 --bfile $BFILE --mlma --grm $GRM --pheno MYBL2.phen --qcovar rs142261424.covar --thread-num $THREADS --out MYBL2_cond_rs142261424
gcta64 --bfile $BFILE --mlma --grm $GRM --pheno MYBL2.phen --qcovar IFT52.phen --thread-num $THREADS --out MYBL2_cond_IFT52_exp

# IFT52: SNP (rs285205) correction and MYBL2 expression level correction
gcta64 --bfile $BFILE --mlma --grm $GRM --pheno IFT52.phen --qcovar rs285205.covar --thread-num $THREADS --out IFT52_cond_rs285205
gcta64 --bfile $BFILE --mlma --grm $GRM --pheno IFT52.phen --qcovar MYBL2.phen --thread-num $THREADS --out IFT52_cond_MYBL2_exp

echo "All GCTA analyses completed successfully."
