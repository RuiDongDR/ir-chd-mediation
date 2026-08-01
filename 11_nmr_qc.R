## =====================================================================
## 11_nmr_qc.R
##
## Title:   Nightingale NMR quality control and normalisation with ukbnmr.
##
## Purpose: Process the UK Biobank Nightingale NMR metabolite biomarkers
##          to remove technical variation and flag outliers:
##            1. Restrict to White-European, genotype-QC'd participants.
##            2. Remove technical variation (effects of sample-prep time,
##               spectrometer batch/drift, and well position) using
##               ukbnmr::remove_technical_variation().
##            3. Keep the baseline assessment visit only (visit_index 0).
##            4. Keep the 168 absolute-concentration biomarkers (drop
##               derived ratios/percentages).
##            5. Flag/remove outliers (median +/- 4*IQR -> NA), following
##               Julkunen et al. (2023).
##            6. Apply log1p() transformation.
##
## Inputs:  cfg$ukb_nmr        Nightingale NMR biomarkers + QC flags
##                             (raw UK Biobank Nightingale NMR data)
##          cfg$sample_qc_ids  one-column table of QC-passed European ids
##
## Outputs: <out_dir>/metab_white_eur_qced.rds
##                                 full ukbnmr result list (all components)
##          <out_dir>/metab_white_eur_qced.csv
##                                 technical-variation-removed biomarkers
##          <out_dir>/metab_white_eur_qced_no_outliers_log1p.csv
##                                 baseline, 168 absolute biomarkers,
##                                 outliers removed, log1p-transformed
##
## Requires: R packages data.table, ukbnmr, dplyr.
##           (ukbnmr encodes the Nightingale processing steps and ships
##           the `nmr_info` biomarker metadata table.)
##
## Run after: 7_ukb_cohort_qc (for cfg$sample_qc_ids); requires the raw
##            UK Biobank Nightingale NMR data.
##
## ---------------------------------------------------------------------
## CONFIG KEYS USED: ukb_nmr, sample_qc_ids, out_dir
## =====================================================================

rm(list = ls())

suppressMessages(library(data.table))
suppressMessages(library(ukbnmr))
suppressMessages(library(dplyr))

source("config.R")

## ---------------------------------------------------------------------
## 1. Load data and restrict to White-European, QC'd participants
## ---------------------------------------------------------------------
metab_data <- fread(cfg$ukb_nmr)
cat("Raw metabolomics data:", nrow(metab_data), "rows x",
    ncol(metab_data), "columns\n")

white_eur_ids <- fread(cfg$sample_qc_ids)
cat("White-European individuals in QC list:", nrow(white_eur_ids), "\n")

metab_filtered <- metab_data[eid %in% white_eur_ids$eid]
cat("After filtering to White-European:", nrow(metab_filtered), "rows x",
    ncol(metab_filtered), "columns\n")
cat("Individuals removed:", nrow(metab_data) - nrow(metab_filtered), "\n")

## ---------------------------------------------------------------------
## 2. Remove technical variation with ukbnmr
##    Removes the effects of technical covariates (sample-preparation
##    time, spectrometer batch and drift, well position). The input file
##    already carries the per-sample QC flags merged with the biomarkers.
##    Returns a list: biomarkers, biomarker_qc_flags, sample_processing,
##    algorithm_version, etc.
## ---------------------------------------------------------------------
cat("\nApplying ukbnmr::remove_technical_variation() ...\n")
start_time <- Sys.time()

result <- remove_technical_variation(metab_filtered)

cat("Duration:",
    round(as.numeric(difftime(Sys.time(), start_time, units = "secs"))),
    "seconds\n")
cat("Result components:", paste(names(result), collapse = ", "), "\n")
cat("Algorithm version:", result$algorithm_version, "\n")

metab_cleaned <- result$biomarkers
cat("Biomarkers after technical-variation removal:",
    nrow(metab_cleaned), "rows x", ncol(metab_cleaned), "columns\n")

## Persist the full result list and the cleaned biomarker table.
saveRDS(result, file = file.path(cfg$out_dir, "metab_white_eur_qced.rds"))
write.csv(metab_cleaned,
          file = file.path(cfg$out_dir, "metab_white_eur_qced.csv"),
          row.names = FALSE)

## ---------------------------------------------------------------------
## 3. Keep baseline visit only (visit_index == 0)
## ---------------------------------------------------------------------
metab_i0 <- metab_cleaned[visit_index == 0]
cat("After filtering to baseline (visit_index 0):",
    nrow(metab_i0), "rows x", ncol(metab_i0), "columns\n")

## ---------------------------------------------------------------------
## 4. Keep the 168 absolute-concentration biomarkers
##    Drop derived biomarkers (ratios and percentages), identified by the
##    ukbnmr naming convention (_pct, _by_, _to_).
## ---------------------------------------------------------------------
all_cols   <- names(metab_i0)
metab_cols <- setdiff(all_cols, c("eid", "visit_index"))

derived_pattern <- "_pct_|_pct$|_by_|_to_"
abs_conc_cols   <- metab_cols[!grepl(derived_pattern, metab_cols)]
derived_cols    <- metab_cols[grepl(derived_pattern, metab_cols)]
cat("Absolute-concentration biomarkers kept:", length(abs_conc_cols), "\n")
cat("Derived (ratio/percentage) biomarkers dropped:", length(derived_cols), "\n")

cols_to_keep <- c("eid", abs_conc_cols)
metab_168    <- metab_i0[, ..cols_to_keep]

## Cross-check biomarker types against the ukbnmr metadata table.
nmr_info_in_168 <- nmr_info %>% filter(Biomarker %in% colnames(metab_168))
print(table(nmr_info_in_168$Type))

## ---------------------------------------------------------------------
## 5. Flag/remove outliers: values outside median +/- 4*IQR set to NA
##    (Julkunen et al., 2023).
## ---------------------------------------------------------------------
cat("\nRemoving outliers (median +/- 4*IQR) for",
    length(abs_conc_cols), "biomarkers ...\n")

metab_no_outliers <- copy(metab_168)
outlier_counts <- data.table(
  metabolite   = abs_conc_cols,
  n_outliers   = 0L,
  pct_outliers = 0.0
)

for (i in seq_along(abs_conc_cols)) {
  col    <- abs_conc_cols[i]
  values <- metab_no_outliers[[col]]

  med <- median(values, na.rm = TRUE)
  q1  <- quantile(values, 0.25, na.rm = TRUE)
  q3  <- quantile(values, 0.75, na.rm = TRUE)
  iqr <- q3 - q1

  lower_bound <- med - 4 * iqr
  upper_bound <- med + 4 * iqr
  is_outlier  <- !is.na(values) & (values < lower_bound | values > upper_bound)

  n_outliers_count   <- sum(is_outlier)
  n_valid_count      <- sum(!is.na(values))
  pct_outliers_value <- round(100 * n_outliers_count / n_valid_count, 3)

  set(metab_no_outliers, which(is_outlier), col, NA_real_)
  set(outlier_counts, i, "n_outliers", n_outliers_count)
  set(outlier_counts, i, "pct_outliers", pct_outliers_value)
}

total_outliers <- sum(outlier_counts$n_outliers)
total_values   <- nrow(metab_no_outliers) * length(abs_conc_cols)
cat("Total outliers removed:", total_outliers,
    sprintf("(%.3f%%)\n", 100 * total_outliers / total_values))
cat("Biomarkers with any outliers:", sum(outlier_counts$pct_outliers > 0), "\n")

## ---------------------------------------------------------------------
## 6. log1p transformation: log(1 + x)
##    Handles zeros gracefully and approximates log(x) for large values.
## ---------------------------------------------------------------------
cat("\nApplying log1p transformation to",
    length(abs_conc_cols), "biomarkers ...\n")

metab_log1p <- copy(metab_no_outliers)
metab_log1p[, (abs_conc_cols) := lapply(.SD, log1p), .SDcols = abs_conc_cols]

cat("Final data:", nrow(metab_log1p), "rows x", ncol(metab_log1p), "columns\n")

write.csv(
  metab_log1p,
  file = file.path(cfg$out_dir, "metab_white_eur_qced_no_outliers_log1p.csv"),
  row.names = FALSE
)

cat("\nDone. Outputs written to:", cfg$out_dir, "\n")
