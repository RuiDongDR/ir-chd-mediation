## =====================================================================
## 9_compute_prs.R
##
## Purpose:
##   Sum the PLINK2 per-chromosome SCORE1_SUM values across all 22
##   autosomes to obtain a single per-individual insulin-resistance
##   polygenic risk score (PRS_IR). SCORE1_SUM is used (rather than
##   SCORE1_AVG x ALLELE_CT) because it correctly accounts for the
##   mean-imputation of missing genotypes performed by PLINK --score.
##   Called by 9_compute_prs.sh.
##
## Inputs:
##   argv[1]  score_dir   directory with PRS_result_chr{1..22}.sscore
##   argv[2]  out_file    output path for FID IID PRS_IR
##
## Outputs:
##   <out_file>   tab-separated: FID  IID  PRS_IR
##
## Requires:
##   - R (>= 4.x) with data.table, tidyverse
##
## Run after: PLINK2 --score step in 9_compute_prs.sh
## =====================================================================
##
## CONFIG KEYS USED:
##   out_dir   (default output directory, used if argv[2] not supplied)
## =====================================================================

suppressMessages({
  library(tidyverse)
  library(data.table)
})
options(scipen = 6, digits = 4)

source("config.R")

args      <- commandArgs(trailingOnly = TRUE)
score_dir <- if (length(args) >= 1) args[1] else file.path(cfg$work_dir, "prs_ir")
out_file  <- if (length(args) >= 2) args[2] else file.path(cfg$out_dir, "UKBB_PRS_IR.txt")

if (!dir.exists(score_dir)) stop(paste("Directory not found:", score_dir))

## ---- Load per-chromosome scores ------------------------------------
prs_list <- list()
for (chr in 1:22) {
  chr_file <- file.path(score_dir, paste0("PRS_result_chr", chr, ".sscore"))
  if (!file.exists(chr_file)) {
    warning(paste("File not found:", chr_file)); next
  }
  chr_data <- fread(chr_file)
  colnames(chr_data)[1] <- "FID"
  chr_data <- chr_data %>%
    select(FID, IID, SCORE1_SUM) %>%
    rename(!!paste0("SCORE1_SUM_chr", chr) := SCORE1_SUM)
  prs_list[[chr]] <- chr_data
}
cat("Loaded per-chromosome scores\n")

## ---- Merge and sum across autosomes --------------------------------
prs_all <- prs_list[[1]]
for (chr in 2:22) {
  if (is.null(prs_list[[chr]])) next
  prs_all <- prs_all %>% left_join(prs_list[[chr]], by = c("FID", "IID"))
}

setDT(prs_all)
prs_ir <- prs_all[, .(
  FID = FID,
  IID = IID,
  PRS_IR = rowSums(.SD, na.rm = TRUE)
), .SDcols = patterns("^SCORE1_SUM")]

cat("Computed IR PRS for", nrow(prs_ir), "individuals\n")

## ---- Save ----------------------------------------------------------
out_dir <- dirname(out_file)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(prs_ir, out_file, sep = "\t")
cat("Results saved to:", out_file, "\n")
