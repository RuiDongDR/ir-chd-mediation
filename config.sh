#!/usr/bin/env bash
# =====================================================================
# config.sh  —  paths for the shell steps of the IR -> CHD pipeline
# Edit the placeholders, then `source config.sh` at the top of each
# shell script. No absolute paths are shipped; set your own below.
# =====================================================================

# ---- external tools (must be on PATH or give full command) ----
export PLINK="plink"            # PLINK 1.9  (LD clumping)
export PLINK2="plink2"          # PLINK 2.0  (--score, .sscore output)
export PYTHON="python3"

# ---- reference panels ----
export REF_1000G_EUR="PATH/TO/1000G_phase3_EUR"   # PLINK bfile prefix, EUR
export PRSCS_DIR="PATH/TO/PRScs"                  # PRS-CS install (PRScs.py)
export PRSCS_LD_REF="PATH/TO/ldblk_1kg_eur"       # PRS-CS 1000G EUR LD reference

# ---- pipeline inputs / outputs ----
export F1_SUMSTATS="PATH/TO/F1_prscs_sst.txt"         # from step 4 (PRS-CS sst_file: SNP A1 A2 BETA SE)
export F1_N_EFF="SET_N_GWAS"                          # F1 (IR-factor) GWAS effective sample size -> PRS-CS --n_gwas
export UKB_GENO_DIR="PATH/TO/ukb_imputed"             # per-chromosome genotypes (directory)
export UKB_BFILE_PREFIX="PATH/TO/ukb_imputed/ukb_chr" # per-chromosome PLINK bfile prefix; chromosome number is appended (e.g. ukb_chr1)
export WORK_DIR="work"
export OUT_DIR="results"

mkdir -p "$WORK_DIR" "$OUT_DIR"
