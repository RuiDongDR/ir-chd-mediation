#!/usr/bin/env Rscript
# =============================================================================
# 13_serial_mediation.R
#
# Title:   Serial (two-mediator) mediation: PRS_IR -> protein -> metabolite -> CHD
#
# Purpose: Estimate serial mediation models (Hayes PROCESS Model 6) as
#          structural equation models with lavaan, for protein x metabolite
#          pairs where both the protein and the metabolite passed the
#          single-mediator screen. A two-pass scheme is used:
#            Prescreen  Retain proteins and metabolites whose single-mediator
#                       proportion-mediated p-value < 0.05 (from step 12).
#            Pass 1     Fit every prescreened protein x metabolite triad by
#                       maximum likelihood; screen the serial-indirect path
#                       (a1 * d21 * b2) using delta-method standard errors,
#                       keeping triads with p < 0.05.
#            Pass 2     Refit each surviving triad with nonparametric
#                       bootstrap (1000 resamples) for percentile CIs and
#                       bootstrap p-values on all indirect paths.
#          Before modelling, each mediator and the exposure are rank-inverse-
#          normal transformed and residualized on the covariates.
#
# Compute note: Compute-intensive; originally run as parallel array jobs on a
#          cluster (one array task per protein, PSOCK-parallel over metabolites).
#          Here proteins are looped within R, or a single protein can be
#          supplied as an argument; the metabolite loop remains PSOCK-parallel.
#
# Inputs:  cfg$prs_ir_scores    - PRS_IR per participant (IID -> eid, PRS_IR)
#          cfg$sample_qc_ids     - QC-passed participant ids to retain (eid)
#          cfg$chd_cases         - CHD case ids (eid)
#          cfg$covariate_data    - analysis covariates (eid + covariates)
#          cfg$protein_data      - analysis-ready Olink NPX (eid + proteins)
#          cfg$metabolite_data   - analysis-ready NMR metabolites (eid + metabolites)
#          <out_dir>/single_mediator_protein_mediation_results.csv     (step 12)
#          <out_dir>/single_mediator_metabolite_mediation_results.csv  (step 12)
#
# Outputs: <out_dir>/serial_pass1_delta_screen.csv    - all fitted triads (Pass 1)
#          <out_dir>/serial_pass2_bootstrap.csv         - bootstrap-confirmed triads (Pass 2)
#
# Requires: R (>= 4.1)
#           data.table (>= 1.14)
#           dplyr      (>= 1.1)
#           lavaan     (>= 0.6-17)
#           parallel   (base R)
#
# Run after: 12_single_mediator_mediation.R (its result CSVs drive the prescreen).
#
# Usage:   Rscript 13_serial_mediation.R                 # all prescreened proteins
#          Rscript 13_serial_mediation.R <protein_name>   # single protein only
#
# CONFIG KEYS USED:
#   cfg$prs_ir_scores
#   cfg$sample_qc_ids
#   cfg$chd_cases
#   cfg$covariate_data
#   cfg$protein_data
#   cfg$metabolite_data
#   cfg$out_dir
# =============================================================================

suppressMessages({
    library(data.table)
    library(dplyr)
    library(lavaan)
    library(parallel)
})

source("config.R")

set.seed(42)
options(digits = 6, scipen = 0)

# -----------------------------------------------------------------------------
# Analysis constants (kept identical to the original cluster scripts)
# -----------------------------------------------------------------------------
PRESCREEN_PVAL        <- 0.05    # single-mediator prop-mediated p to enter prescreen
STEP1_PVAL_THRESHOLD  <- 0.05    # Pass 1 delta-method p to carry to bootstrap
N_BOOTSTRAP           <- 1000    # Pass 2 nonparametric bootstrap resamples

PROT_SINGLE_RESULTS  <- file.path(cfg$out_dir, "single_mediator_protein_mediation_results.csv")
METAB_SINGLE_RESULTS <- file.path(cfg$out_dir, "single_mediator_metabolite_mediation_results.csv")
OUT_PASS1            <- file.path(cfg$out_dir, "serial_pass1_delta_screen.csv")
OUT_PASS2            <- file.path(cfg$out_dir, "serial_pass2_bootstrap.csv")

# Optional single-protein target
args           <- commandArgs(trailingOnly = TRUE)
target_protein <- if (length(args) >= 1) args[1] else NA_character_

# Number of PSOCK workers for the metabolite loop (conservative local default;
# raise it to match your machine, or parallelize across proteins as you prefer).
N_CORES <- max(1L, min(4L, parallel::detectCores() - 1L))

cat("R version:", R.version.string, "\n")
cat("lavaan version:", as.character(packageVersion("lavaan")), "\n")
cat(sprintf("Prescreen p-value:        %.2f\n", PRESCREEN_PVAL))
cat(sprintf("Pass 1 delta p-value:     %.2f\n", STEP1_PVAL_THRESHOLD))
cat(sprintf("Pass 2 bootstrap samples: %d\n", N_BOOTSTRAP))
cat(sprintf("PSOCK workers:            %d\n\n", N_CORES))

# =============================================================================
# Prescreen — proteins and metabolites significant in single-mediator step
# =============================================================================
cat("=== Prescreen: loading single-mediator results (step 12) ===\n")

if (!file.exists(PROT_SINGLE_RESULTS) || !file.exists(METAB_SINGLE_RESULTS)) {
    stop("Single-mediator result CSVs not found. Run 12_single_mediator_mediation.R first.")
}

prot_med  <- fread(PROT_SINGLE_RESULTS)
metab_med <- fread(METAB_SINGLE_RESULTS)
prot_med$prop_mediated_pvalue  <- as.numeric(prot_med$prop_mediated_pvalue)
metab_med$prop_mediated_pvalue <- as.numeric(metab_med$prop_mediated_pvalue)

sig_proteins    <- prot_med[!is.na(prop_mediated_pvalue) &
                            prop_mediated_pvalue < PRESCREEN_PVAL, protein]
sig_metabolites <- metab_med[!is.na(prop_mediated_pvalue) &
                            prop_mediated_pvalue < PRESCREEN_PVAL, metabolite]

cat(sprintf("  Proteins passing prescreen:    %d\n", length(sig_proteins)))
cat(sprintf("  Metabolites passing prescreen: %d\n", length(sig_metabolites)))
cat(sprintf("  Candidate triads:              %d\n\n",
            length(sig_proteins) * length(sig_metabolites)))

# =============================================================================
# Load data and build the modelling matrices (INRT + residualization)
# =============================================================================
cat("=== Loading data ===\n")
t_load <- Sys.time()

prs_data        <- fread(cfg$prs_ir_scores)
prs_data$eid    <- prs_data$IID
sample_qc       <- fread(cfg$sample_qc_ids)
chd_cases       <- fread(cfg$chd_cases)
covariates      <- fread(cfg$covariate_data)

proteomics      <- fread(cfg$protein_data)
proteomics      <- proteomics %>% dplyr::select(!ends_with("_dup"))
metabolomics    <- fread(cfg$metabolite_data)

cat(sprintf("  Data loaded in %.0f seconds\n",
            as.numeric(difftime(Sys.time(), t_load, units = "secs"))))

# Restrict to prescreened mediators that exist in the assay tables
protein_cols_all <- setdiff(names(proteomics), "eid")
protein_cols_all <- protein_cols_all[!endsWith(protein_cols_all, "_dup")]
metab_cols_all   <- setdiff(names(metabolomics), "eid")

proteins_keep <- intersect(sig_proteins,    protein_cols_all)
metab_cols    <- intersect(sig_metabolites, metab_cols_all)

if (!is.na(target_protein)) {
    if (!(target_protein %in% proteins_keep)) {
        stop(sprintf("Target protein '%s' is not a prescreened protein present in the data.",
                     target_protein))
    }
    proteins_keep <- target_protein
}
cat(sprintf("  Proteins to model:    %d\n", length(proteins_keep)))
cat(sprintf("  Metabolites to model: %d\n", length(metab_cols)))

# Mean-impute mediator columns (proteins + all prescreened metabolites)
proteomics_sub   <- proteomics[, c("eid", proteins_keep), with = FALSE]
metabolomics_sub <- metabolomics[, c("eid", metab_cols), with = FALSE]
for (col in proteins_keep) {
    proteomics_sub[is.na(get(col)), (col) := mean(proteomics_sub[[col]], na.rm = TRUE)]
}
for (col in metab_cols) {
    metabolomics_sub[is.na(get(col)), (col) := mean(metabolomics_sub[[col]], na.rm = TRUE)]
}

# Assemble analysis dataset
analysis_data <- as.data.frame(prs_data[, .(eid, PRS_IR)]) %>%
    inner_join(as.data.frame(sample_qc[, .(eid)]), by = "eid") %>%
    inner_join(as.data.frame(covariates),          by = "eid") %>%
    mutate(CHD = as.integer(eid %in% chd_cases$eid)) %>%
    inner_join(as.data.frame(proteomics_sub),      by = "eid") %>%
    inner_join(as.data.frame(metabolomics_sub),    by = "eid") %>%
    as.data.table()

cat(sprintf("  Analysis dataset: N = %d\n", nrow(analysis_data)))
cat(sprintf("  CHD cases: %d (%.1f%%)\n",
            sum(analysis_data$CHD), 100 * mean(analysis_data$CHD)))

# --- Factor coding and covariate imputation ---------------------------------
analysis_data$sex               <- as.factor(analysis_data$sex)
analysis_data$highest_education <- as.factor(analysis_data$highest_education)
analysis_data$employment        <- as.factor(analysis_data$employment)
analysis_data$smoking_status    <- as.factor(analysis_data$smoking_status)
analysis_data$sleeplessness     <- as.factor(analysis_data$sleeplessness)

get_mode <- function(x) {
    x <- x[!is.na(x)]
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
}
continuous_covars  <- c("age_at_recruitment", "townsend_deprivation_index", paste0("pc", 1:10))
categorical_covars <- c("sex", "highest_education", "employment",
                        "physical_activity", "smoking_status", "sleeplessness")
for (v in continuous_covars) {
    if (v %in% names(analysis_data) && any(is.na(analysis_data[[v]])))
        analysis_data[[v]][is.na(analysis_data[[v]])] <- mean(analysis_data[[v]], na.rm = TRUE)
}
for (v in categorical_covars) {
    if (v %in% names(analysis_data) && any(is.na(analysis_data[[v]])))
        analysis_data[[v]][is.na(analysis_data[[v]])] <- get_mode(analysis_data[[v]])
}
covariate_vars <- c("age_at_recruitment", "sex", paste0("pc", 1:10),
                    "townsend_deprivation_index", "highest_education", "employment",
                    "physical_activity", "smoking_status", "sleeplessness")
covariate_vars <- intersect(covariate_vars, names(analysis_data))

# --- INRT normalization and covariate residualization -----------------------
cat("\n=== INRT normalization and residualization ===\n")
t0 <- Sys.time()

inrt <- function(x) {
    qnorm((rank(x, na.last = "keep", ties.method = "average") - 0.5) / sum(!is.na(x)))
}

ad <- as.data.frame(analysis_data)

proteins_mat <- matrix(
    sapply(proteins_keep, function(p) inrt(ad[[p]])),
    nrow = nrow(ad), ncol = length(proteins_keep),
    dimnames = list(ad$eid, proteins_keep)
)
metabs_mat <- sapply(metab_cols, function(m) inrt(ad[[m]]))
if (!is.matrix(metabs_mat)) {
    metabs_mat <- matrix(metabs_mat, nrow = nrow(ad), ncol = length(metab_cols),
                         dimnames = list(ad$eid, metab_cols))
} else {
    rownames(metabs_mat) <- ad$eid
}

cov_formula <- as.formula(paste("~ ", paste(covariate_vars, collapse = " + ")))
cov_mat     <- model.matrix(cov_formula, data = ad)[, -1]

resid_matrix_safe <- function(Y_mat, X_mat) {
    Xd <- cbind(1, X_mat)
    apply(Y_mat, 2, function(y) {
        ok  <- !is.na(y)
        out <- rep(NA_real_, length(y))
        if (sum(ok) > ncol(Xd) + 1) {
            fit <- lm.fit(Xd[ok, , drop = FALSE], y[ok])
            out[ok] <- fit$residuals
        }
        out
    })
}

resid_prs   <- residuals(lm.fit(cbind(1, cov_mat), ad$PRS_IR))
resid_prot  <- matrix(resid_matrix_safe(proteins_mat, cov_mat),
                      nrow = nrow(proteins_mat), ncol = ncol(proteins_mat),
                      dimnames = dimnames(proteins_mat))
resid_metab <- resid_matrix_safe(metabs_mat, cov_mat)
chd_outcome <- ad$CHD

cat(sprintf("  Residualization completed in %.0f seconds\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# =============================================================================
# Serial SEM model template (Hayes PROCESS Model 6)
# =============================================================================
serial_model_template <- '
  protein_k    ~ a1 * PRS_IR
  metabolite_l ~ a2 * PRS_IR + d21 * protein_k
  CHD          ~ c_prime * PRS_IR + b1 * protein_k + b2 * metabolite_l

  serial_indirect := a1 * d21 * b2
  indirect_prot   := a1 * b1
  indirect_metab  := a2 * b2
  total_indirect  := serial_indirect + indirect_prot + indirect_metab
  total_effect    := total_indirect + c_prime
  prop_serial     := serial_indirect / total_effect
'

# =============================================================================
# PASS 1 — Delta-method screening (ML fit) over all prescreened triads
# =============================================================================
# One triad: exposure -> protein k -> metabolite l -> CHD, fitted by ML with
# analytic (delta-method) SEs; the serial_indirect path is screened at p < 0.05.
fit_triad_delta <- function(prot_col, metab_col,
                            prs_vec, prot_vec, metab_vec, chd_vec) {
    model_spec <- gsub("protein_k",    prot_col,  serial_model_template)
    model_spec <- gsub("metabolite_l", metab_col, model_spec)

    dat <- data.frame(PRS_IR = prs_vec, CHD = chd_vec)
    dat[[prot_col]]  <- prot_vec
    dat[[metab_col]] <- metab_vec

    fit <- tryCatch(
        suppressWarnings(lavaan::sem(model_spec, data = dat,
                                     se = "standard", estimator = "ML")),
        error = function(e) NULL
    )
    if (is.null(fit)) {
        return(data.frame(protein = prot_col, metabolite = metab_col,
                          est = NA_real_, se = NA_real_, pvalue = NA_real_,
                          converged = FALSE))
    }
    pe <- tryCatch(lavaan::parameterEstimates(fit), error = function(e) NULL)
    if (is.null(pe)) {
        return(data.frame(protein = prot_col, metabolite = metab_col,
                          est = NA_real_, se = NA_real_, pvalue = NA_real_,
                          converged = FALSE))
    }
    row <- pe[pe$label == "serial_indirect", ]
    data.frame(
        protein    = prot_col,
        metabolite = metab_col,
        est        = row$est,
        se         = row$se,
        pvalue     = row$pvalue,
        converged  = lavaan::lavInspect(fit, "converged")
    )
}

cat("\n=== Pass 1: delta-method screening ===\n")
cat(sprintf("  Screening %d proteins x %d metabolites\n",
            length(proteins_keep), length(metab_cols)))

pass1_list <- list()
for (pi in seq_along(proteins_keep)) {
    prot_col <- proteins_keep[pi]
    cat(sprintf("  [%d/%d] protein %s\n", pi, length(proteins_keep), prot_col))

    cl <- makeCluster(N_CORES, type = "PSOCK")
    clusterExport(cl, varlist = c("fit_triad_delta", "serial_model_template"),
                  envir = environment())
    invisible(clusterEvalQ(cl, suppressMessages(library(lavaan))))

    batch <- clusterMap(
        cl, fun = fit_triad_delta,
        metab_col = metab_cols,
        metab_vec = lapply(metab_cols, function(m) resid_metab[, m]),
        MoreArgs = list(
            prot_col  = prot_col,
            prs_vec   = resid_prs,
            prot_vec  = resid_prot[, prot_col],
            chd_vec   = chd_outcome
        ),
        SIMPLIFY = FALSE
    )
    stopCluster(cl)

    pass1_list[[pi]] <- rbindlist(batch, fill = TRUE)
}

pass1_dt <- rbindlist(pass1_list, fill = TRUE)
fwrite(pass1_dt, OUT_PASS1)
cat(sprintf("  Pass 1 triads written: %d -> %s\n", nrow(pass1_dt), OUT_PASS1))

sig_triads <- pass1_dt[!is.na(pvalue) & pvalue < STEP1_PVAL_THRESHOLD]
cat(sprintf("  Significant triads (p < %.2f): %d\n",
            STEP1_PVAL_THRESHOLD, nrow(sig_triads)))

if (nrow(sig_triads) == 0) {
    cat("No triads survived Pass 1. Nothing to bootstrap. Done.\n")
    quit(save = "no", status = 0)
}

# =============================================================================
# PASS 2 — Bootstrap confirmation of surviving triads
# =============================================================================
# Refit each surviving triad with nonparametric bootstrap (1000 resamples) to
# obtain percentile CIs and bootstrap p-values for every indirect path.
fit_triad_bootstrap <- function(prot_col, metab_col,
                                prs_vec, prot_vec, metab_vec, chd_vec, n_boot) {
    model_spec <- gsub("protein_k",    prot_col,  serial_model_template)
    model_spec <- gsub("metabolite_l", metab_col, model_spec)

    dat <- data.frame(PRS_IR = prs_vec, CHD = chd_vec)
    dat[[prot_col]]  <- prot_vec
    dat[[metab_col]] <- metab_vec

    na_row <- data.frame(
        protein = prot_col, metabolite = metab_col,
        serial_indirect_est = NA_real_, serial_indirect_se = NA_real_,
        serial_indirect_pval = NA_real_,
        serial_indirect_ci_lower = NA_real_, serial_indirect_ci_upper = NA_real_,
        indirect_prot_est = NA_real_,  indirect_prot_pval = NA_real_,
        indirect_metab_est = NA_real_, indirect_metab_pval = NA_real_,
        total_indirect_est = NA_real_, total_indirect_pval = NA_real_,
        total_effect_est = NA_real_,   prop_serial_est = NA_real_,
        converged = FALSE, n_boot_used = NA_integer_
    )

    fit <- tryCatch(
        suppressWarnings(lavaan::sem(model_spec, data = dat,
                                     se = "bootstrap", bootstrap = n_boot,
                                     estimator = "ML")),
        error = function(e) NULL
    )
    if (is.null(fit)) return(na_row)

    pe <- tryCatch(
        lavaan::parameterEstimates(fit, ci = TRUE, level = 0.95,
                                   boot.ci.type = "perc"),
        error = function(e) NULL
    )
    if (is.null(pe)) return(na_row)

    get_row <- function(label) pe[pe$label == label, ]
    r_si <- get_row("serial_indirect"); r_ip <- get_row("indirect_prot")
    r_im <- get_row("indirect_metab");  r_ti <- get_row("total_indirect")
    r_te <- get_row("total_effect");    r_ps <- get_row("prop_serial")

    data.frame(
        protein = prot_col, metabolite = metab_col,
        serial_indirect_est      = r_si$est,
        serial_indirect_se       = r_si$se,
        serial_indirect_pval     = r_si$pvalue,
        serial_indirect_ci_lower = r_si$ci.lower,
        serial_indirect_ci_upper = r_si$ci.upper,
        indirect_prot_est   = r_ip$est,  indirect_prot_pval  = r_ip$pvalue,
        indirect_metab_est  = r_im$est,  indirect_metab_pval = r_im$pvalue,
        total_indirect_est  = r_ti$est,  total_indirect_pval = r_ti$pvalue,
        total_effect_est    = r_te$est,
        prop_serial_est     = r_ps$est,
        converged   = lavaan::lavInspect(fit, "converged"),
        n_boot_used = n_boot
    )
}

cat("\n=== Pass 2: bootstrap confirmation ===\n")
cat(sprintf("  Bootstrapping %d triads (%d resamples each)\n",
            nrow(sig_triads), N_BOOTSTRAP))

boot_proteins <- sort(unique(sig_triads$protein))
pass2_list <- list()
for (pi in seq_along(boot_proteins)) {
    prot_col <- boot_proteins[pi]
    metabs_this <- sig_triads[protein == prot_col, metabolite]
    cat(sprintf("  [%d/%d] protein %s (%d metabolites)\n",
                pi, length(boot_proteins), prot_col, length(metabs_this)))

    cl <- makeCluster(N_CORES, type = "PSOCK")
    clusterExport(cl, varlist = c("fit_triad_bootstrap", "serial_model_template"),
                  envir = environment())
    invisible(clusterEvalQ(cl, suppressMessages(library(lavaan))))

    batch <- clusterMap(
        cl, fun = fit_triad_bootstrap,
        metab_col = metabs_this,
        metab_vec = lapply(metabs_this, function(m) resid_metab[, m]),
        MoreArgs = list(
            prot_col = prot_col,
            prs_vec  = resid_prs,
            prot_vec = resid_prot[, prot_col],
            chd_vec  = chd_outcome,
            n_boot   = N_BOOTSTRAP
        ),
        SIMPLIFY = FALSE
    )
    stopCluster(cl)

    pass2_list[[pi]] <- rbindlist(batch, fill = TRUE)
}

pass2_dt <- rbindlist(pass2_list, fill = TRUE)
fwrite(pass2_dt, OUT_PASS2)

# =============================================================================
# Summary
# =============================================================================
cat("\n=== Summary ===\n")
cat(sprintf("  Pass 2 triads bootstrapped: %d -> %s\n", nrow(pass2_dt), OUT_PASS2))
cat(sprintf("  Converged: %d (%.1f%%)\n",
            sum(pass2_dt$converged, na.rm = TRUE),
            100 * mean(pass2_dt$converged, na.rm = TRUE)))
cat(sprintf("  Bootstrap serial-indirect p < 0.05: %d\n",
            sum(!is.na(pass2_dt$serial_indirect_pval) &
                pass2_dt$serial_indirect_pval < 0.05)))
cat(sprintf("  95%% CI excludes 0: %d\n",
            sum((!is.na(pass2_dt$serial_indirect_ci_lower) &
                 pass2_dt$serial_indirect_ci_lower > 0) |
                (!is.na(pass2_dt$serial_indirect_ci_upper) &
                 pass2_dt$serial_indirect_ci_upper < 0), na.rm = TRUE)))
cat("\nDone:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
