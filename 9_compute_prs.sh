#!/usr/bin/env bash
# =====================================================================
# 9_compute_prs.sh
#
# Purpose:
#   Apply the PRS-CS posterior SNP weights (step 6) to the UK Biobank
#   imputed genotypes to compute a per-individual insulin-resistance
#   polygenic risk score. PLINK2 --score is run per autosome (each
#   individual's dosages weighted by the PRS-CS posterior effect and
#   summed within the chromosome), then the 22 per-chromosome sums are
#   added together to yield a single IR PRS per individual. The total
#   number of scored variants (summed over autosomes) is also printed.
#
# Inputs:
#   - ${OUT_DIR}/prscs_weights/IR_pst_eff_a1_b0.5_phi1e-02_chr{1..22}.txt
#     PRS-CS weights (col 2 = SNP, col 4 = A1, col 6 = BETA_POST)
#   - ${UKB_BFILE_PREFIX}${chr}.{bed,bim,fam}  per-chromosome UKB genotypes
#
# Outputs:
#   - ${WORK_DIR}/prs_ir/PRS_result_chr{1..22}.sscore  per-chromosome scores
#   - ${OUT_DIR}/UKBB_PRS_IR.txt   FID  IID  PRS_IR (summed over autosomes)
#
# Requires:
#   - PLINK v2.0 (--score, cols=+scoresums)
#   - R (>= 4.x) with data.table, tidyverse  (via 9_compute_prs.R)
#
# Run after: 6_prscs_weights (produces the PRS-CS posterior weights)
# =====================================================================
#
# CONFIG VARS USED:
#   PLINK2            PLINK 2.0 executable
#   UKB_BFILE_PREFIX  per-chromosome UKB PLINK bfile prefix (chr appended)
#   OUT_DIR           output directory (weights input, final PRS output)
#   WORK_DIR          working directory (per-chromosome .sscore)
# =====================================================================

set -euo pipefail
source config.sh

WEIGHTS_DIR="${OUT_DIR}/prscs_weights"
SCORE_DIR="${WORK_DIR}/prs_ir"
mkdir -p "$SCORE_DIR"

# ---- PLINK2 --score per autosome ------------------------------------
# SCORE1_SUM (cols=+scoresums) is used downstream because it correctly
# accounts for mean-imputation of missing genotypes.
for chr in $(seq 1 22); do
  weights="${WEIGHTS_DIR}/IR_pst_eff_a1_b0.5_phi1e-02_chr${chr}.txt"
  "$PLINK2" \
    --bfile "${UKB_BFILE_PREFIX}${chr}" \
    --score "$weights" 2 4 6 cols=+scoresums \
    --out "${SCORE_DIR}/PRS_result_chr${chr}" &
done
wait

# ---- Total number of scored variants (summed across autosomes) ------
total_variants=0
for chr in $(seq 1 22); do
  weights="${WEIGHTS_DIR}/IR_pst_eff_a1_b0.5_phi1e-02_chr${chr}.txt"
  n=$(wc -l < "$weights")
  total_variants=$(( total_variants + n ))
done
echo "Total number of scored variants (22 autosomes): ${total_variants}"

# ---- Sum SCORE1_SUM across the 22 autosomes (R helper) --------------
Rscript 9_compute_prs.R "$SCORE_DIR" "${OUT_DIR}/UKBB_PRS_IR.txt"
