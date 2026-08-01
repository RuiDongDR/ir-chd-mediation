#!/usr/bin/env bash
# =====================================================================
# 5_ld_clumping.sh
#
# Purpose:
#   LD-based clumping in PLINK to define independent genome-wide-
#   significant loci for BOTH (i) the IR latent-factor (F1) association
#   and (ii) the per-SNP heterogeneity test (Q_SNP) from the F1
#   multivariate GWAS. Clumping is performed against the 1000 Genomes
#   Phase 3 EUR reference panel using an index-SNP threshold of
#   P < 5e-8, an LD threshold of r2 < 0.1, and a 1 Mb window.
#
# Inputs:
#   - ${WORK_DIR}/f1.assoc    two columns (SNP  P), F1 association p-values
#   - ${WORK_DIR}/qsnp.assoc  two columns (SNP  P), Q_SNP p-values
#     (both derived from the step-4 multivariate GWAS results; keep only
#      SNPs with a finite p-value in (0, 1])
#   - ${REF_1000G_EUR}.{bed,bim,fam}  1000G Phase 3 EUR PLINK bfile
#
# Outputs:
#   - ${WORK_DIR}/f1_clump.clumped    one row per independent F1 locus
#   - ${WORK_DIR}/qsnp_clump.clumped  one row per independent Q_SNP locus
#   - matching *.log files (report "clumps formed")
#
# Requires:
#   - PLINK v1.9 (LD clumping)
#
# Run after: 4_multivariate_gwas (produces the F1 / Q_SNP summary stats)
# =====================================================================
#
# CONFIG VARS USED:
#   PLINK           PLINK 1.9 executable
#   REF_1000G_EUR   1000G Phase 3 EUR PLINK bfile prefix
#   WORK_DIR        working directory (assoc inputs and clumped outputs)
# =====================================================================

set -euo pipefail
source config.sh

F1_ASSOC="${WORK_DIR}/f1.assoc"
QSNP_ASSOC="${WORK_DIR}/qsnp.assoc"

# ---- LD clumping (index P<5e-8, r2<0.1, 1 Mb window) -----------------
clump() {
  local assoc="$1" out="$2"
  "$PLINK" \
    --bfile "$REF_1000G_EUR" \
    --clump "$assoc" \
    --clump-p1 5e-8 \
    --clump-r2 0.1 \
    --clump-kb 1000 \
    --clump-snp-field SNP \
    --clump-field P \
    --out "$out"
}

clump "$F1_ASSOC"   "${WORK_DIR}/f1_clump"
clump "$QSNP_ASSOC" "${WORK_DIR}/qsnp_clump"

echo "clumping done"
grep -h "clumps formed" "${WORK_DIR}/f1_clump.log" "${WORK_DIR}/qsnp_clump.log" || true
