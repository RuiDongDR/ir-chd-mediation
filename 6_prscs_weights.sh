#!/usr/bin/env bash
# =====================================================================
# 6_prscs_weights.sh
#
# Purpose:
#   Run PRS-CS to derive posterior SNP effect sizes (weights) for the
#   insulin-resistance latent factor (F1) from the multivariate GWAS
#   summary statistics, under the Bayesian continuous-shrinkage prior.
#   A fixed global shrinkage parameter phi = 1e-2 is used, together with
#   the package default scale/shape hyperparameters (a = 1, b = 0.5).
#   The Gibbs sampler is run for 1000 iterations with 500 burn-in and a
#   thinning interval of 5, on each autosome (1..22) in parallel, using
#   the 1000 Genomes Phase 3 EUR LD reference. Only variants shared
#   between the F1 summary statistics and the target (UK Biobank)
#   imputation panel are retained (via --bim_prefix).
#
# Inputs:
#   - ${F1_SUMSTATS}                 F1 GWAS sumstats, PRS-CS format: SNP A1 A2 BETA SE
#                                    (F1_prscs_sst.txt from step 4)
#   - ${F1_N_EFF}                    F1 GWAS effective sample size (--n_gwas)
#   - ${PRSCS_LD_REF}                PRS-CS 1000G EUR LD reference directory (ldblk_1kg_eur)
#   - ${UKB_BFILE_PREFIX}${chr}.bim  per-chromosome target-panel .bim (variant intersection)
#
# Outputs:
#   - ${OUT_DIR}/prscs_weights/IR_pst_eff_a1_b0.5_phi1e-02_chr{1..22}.txt
#     columns: CHR  SNP  BP  A1  A2  BETA_POST
#
# Requires:
#   - PLINK 2.0 (--make-just-bim, to build the target .bim per chromosome)
#   - Python 3 with PRS-CS (PRScs.py); dependencies scipy, h5py, numpy
#
# Run after: 4_multivariate_gwas (produces ${F1_SUMSTATS}); F1_N_EFF is the
#            F1-GWAS effective sample size, which you supply in config.sh.
# =====================================================================
#
# CONFIG VARS USED:
#   PLINK2            PLINK 2.0 executable
#   PYTHON            python3 executable
#   PRSCS_DIR         PRS-CS install directory (contains PRScs.py)
#   PRSCS_LD_REF      PRS-CS 1000G EUR LD reference directory
#   F1_SUMSTATS       F1 multivariate GWAS summary statistics (sst_file)
#   F1_N_EFF          F1 GWAS effective sample size (--n_gwas; you set this)
#   UKB_GENO_DIR      UK Biobank imputed BGEN directory (per-chromosome)
#   UKB_BFILE_PREFIX  per-chromosome target .bim prefix (chr appended)
#   OUT_DIR           output directory
# =====================================================================

set -euo pipefail
source config.sh

OUTDIR="${OUT_DIR}/prscs_weights"
mkdir -p "$OUTDIR"

# --out_dir is used by PRS-CS as an output PREFIX; per-chromosome files
# are written as <prefix>_pst_eff_a1_b0.5_phi1e-02_chr<chr>.txt
OUT_PREFIX="${OUTDIR}/IR"

# ---- (0) Build per-chromosome target .bim files ---------------------
# PRS-CS restricts to variants shared with the target panel via
# --bim_prefix. Build a .bim per autosome from the UK Biobank imputed
# BGEN dosages (ref-first) with PLINK2 --make-just-bim.
for chr in $(seq 1 22); do
  "$PLINK2" \
    --bgen  "${UKB_GENO_DIR}/ukb22828_c${chr}_b0_v3.bgen" ref-first \
    --sample "${UKB_GENO_DIR}/ukb22828_c${chr}_b0_v3.sample" \
    --make-just-bim \
    --out "${UKB_BFILE_PREFIX}${chr}" &
done
wait

# ---- (1) PRS-CS per autosome, in parallel ---------------------------
for chr in $(seq 1 22); do
  "$PYTHON" "${PRSCS_DIR}/PRScs.py" \
    --ref_dir="$PRSCS_LD_REF" \
    --bim_prefix="${UKB_BFILE_PREFIX}${chr}" \
    --sst_file="$F1_SUMSTATS" \
    --n_gwas="$F1_N_EFF" \
    --chrom="$chr" \
    --phi=1e-2 \
    --a=1 \
    --b=0.5 \
    --n_iter=1000 \
    --n_burnin=500 \
    --thin=5 \
    --out_dir="$OUT_PREFIX" &
done
wait
echo "PRS-CS complete: posterior weights written to ${OUTDIR}"
