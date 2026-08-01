## =====================================================================
## 7_ukb_cohort_qc.R
##
## Title:   UK Biobank sample-level quality control (European cohort)
##
## Purpose: Apply sequential sample-level QC to the UK Biobank phenotype
##          table to derive the European-ancestry analysis cohort used
##          throughout the IR -> CHD analyses. Filtering follows standard
##          UK Biobank practice:
##            (1) genetic-ancestry restriction (European descent),
##            (2) exclusion of heterozygosity / missingness outliers,
##            (3) exclusion of sex-chromosome aneuploidy,
##            (4) removal of individuals missing age or sex,
##            (5) reported vs genetic sex concordance.
##          Sample sizes are reported after every step.
##
## Inputs:  cfg$ukb_pheno        - extracted UK Biobank phenotype table (one
##                                 row per participant; wide "pFIELD[_iN/_aN]"
##                                 column naming).
##
## Outputs: cfg$sample_qc_ids      - one-column ("eid") text file listing the
##                                 QC-passing participants that define the
##                                 analysis cohort.
##
## Requires: R (>= 4.1)
##           data.table (>= 1.14)
##           tidyverse  (>= 1.3.0)   [dplyr, tibble]
##
## Run after: (none) - this is the first step of the analysis pipeline.
##
## ---------------------------------------------------------------------
## CONFIG KEYS USED:
##   ukb_pheno       extracted UKB phenotype table
##   sample_qc_ids   output: QC-passing EID list (defines the analysis cohort)
##   out_dir         final-output directory
##
## UK Biobank public fields used to define the cohort (kept for
## reproducibility; see Methods):
##   22006 Genetic ethnic grouping (Caucasian)
##   22027 Outliers for heterozygosity or missing rate
##   22019 Sex chromosome aneuploidy
##   21022 Age at recruitment
##   31    Sex (reported)
##   22001 Genetic sex
## =====================================================================

source("config.R")

suppressMessages({
  library(data.table)
  library(tidyverse)
})

## ---- Load phenotype table --------------------------------------------
cat("Loading UK Biobank phenotype table...\n")
if (!file.exists(cfg$ukb_pheno)) {
  stop(sprintf("Phenotype table not found at cfg$ukb_pheno: %s", cfg$ukb_pheno))
}
# na.strings: treat empty cells and the literal string "NA" as missing.
ukbb_raw <- fread(cfg$ukb_pheno, na.strings = c("", "NA"))
cat(sprintf("Initial sample size: %d individuals\n", nrow(ukbb_raw)))
cat(sprintf("Total variables:     %d\n\n", ncol(ukbb_raw)))

## Small helper to print retention after each step.
report_step <- function(label, n_before, n_after) {
  cat(sprintf("\n=== %s ===\n", label))
  cat(sprintf("Removed:   %d\n", n_before - n_after))
  cat(sprintf("Remaining: %d individuals\n", n_after))
  cat(sprintf("Retention: %.2f%%\n", 100 * n_after / n_before))
}

## ---- Step 1: Genetic ancestry (Field 22006) --------------------------
## Restrict to participants of European (Caucasian) genetic ancestry.
n0 <- nrow(ukbb_raw)
ukbb_qced <- ukbb_raw %>% filter(p22006 == "Caucasian")
report_step("Step 1: Genetic-ancestry restriction (Field 22006)", n0, nrow(ukbb_qced))

## ---- Step 2: Heterozygosity / missingness outliers (Field 22027) -----
## Field 22027 flags heterozygosity- or missing-rate outliers; keep NA (not
## flagged).
n1 <- nrow(ukbb_qced)
ukbb_qced <- ukbb_qced %>% filter(is.na(p22027))
report_step("Step 2: Heterozygosity / missing-rate outliers (Field 22027)", n1, nrow(ukbb_qced))

## ---- Step 3: Sex chromosome aneuploidy (Field 22019) -----------------
## Field 22019 flags putative sex-chromosome aneuploidy; keep NA (not
## flagged).
n2 <- nrow(ukbb_qced)
ukbb_qced <- ukbb_qced %>% filter(is.na(p22019))
report_step("Step 3: Sex-chromosome aneuploidy (Field 22019)", n2, nrow(ukbb_qced))

## ---- Step 4: Missing age / sex (Fields 21022, 31, 22001) -------------
## Require non-missing genetic sex (22001) and reported sex (31); age at
## recruitment (21022) is also checked.
n3 <- nrow(ukbb_qced)
ukbb_qced <- ukbb_qced %>% filter(!is.na(p22001) & !is.na(p31))
report_step("Step 4: Missing age / sex (Fields 21022, 31, 22001)", n3, nrow(ukbb_qced))

## ---- Step 5: Reported vs genetic sex concordance (31 vs 22001) -------
## Retain only participants whose reported sex (Field 31) matches their
## genetically determined sex (Field 22001).
n4 <- nrow(ukbb_qced)
sex_concordant <- (
  (ukbb_qced$p31 == "Female" & ukbb_qced$p22001 == "Female") |
  (ukbb_qced$p31 == "Male"   & ukbb_qced$p22001 == "Male")
) & !is.na(ukbb_qced$p31) & !is.na(ukbb_qced$p22001)
ukbb_qced <- ukbb_qced[sex_concordant, ]
report_step("Step 5: Reported vs genetic sex concordance (Fields 31, 22001)", n4, nrow(ukbb_qced))

## ---- Write QC-passing EID list ---------------------------------------
cat(sprintf("\nFinal QC-passing cohort: %d individuals\n", nrow(ukbb_qced)))
write.table(ukbb_qced[, "eid"],
            cfg$sample_qc_ids,
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
cat(sprintf("QC-passing EID list written to: %s\n", cfg$sample_qc_ids))
