#!/usr/bin/env Rscript
# =============================================================================
# 12_single_mediator_mediation.R
#
# Title:   Single-mediator causal mediation: IR-PRS -> (protein | metabolite) -> CHD
#
# Purpose: For each candidate mediator (one Olink protein or one NMR
#          metabolite), estimate a single-mediator causal mediation model in
#          the potential-outcomes framework of Imai, Keele & Tingley (the
#          `mediation` package). The exposure is the polygenic score for
#          insulin resistance (PRS_IR), the mediator is the individual
#          protein/metabolite, and the outcome is coronary heart disease
#          (CHD, binary). We report the average causal mediation (indirect)
#          effect (ACME), the average direct effect (ADE), the total effect,
#          and the proportion mediated, each with a nonparametric bootstrap
#          p-value (1000 resamples). Proteins are retained only if their
#          complete-case sample size exceeds 30,000; metabolites only if it
#          exceeds 350,000.
#
# Compute note: Compute-intensive; originally run as parallel array jobs on a
#          cluster (one array task per mediator). Here the mediators are looped
#          within R, or a single mediator can be supplied as an argument.
#
# Inputs:  cfg$prs_ir_scores    - PRS_IR per participant (IID -> eid, PRS_IR)
#          cfg$sample_qc_ids     - QC-passed participant ids to retain (eid)
#          cfg$chd_cases         - CHD case ids (eid)
#          cfg$covariate_data    - analysis covariates (eid + covariates)
#          cfg$protein_data      - analysis-ready Olink NPX (eid + proteins)
#          cfg$metabolite_data   - analysis-ready NMR metabolites (eid + metabolites)
#
# Outputs: <out_dir>/single_mediator_protein_mediation_results.csv
#          <out_dir>/single_mediator_metabolite_mediation_results.csv
#          (one row per retained mediator; combined across mediators)
#
# Requires: R (>= 4.1)
#           data.table (>= 1.14)
#           dplyr      (>= 1.1)
#           mediation  (>= 4.5.0)
#
# Run after: analysis-ready protein / metabolite / PRS / covariate tables are
#            prepared (step 11).
#
# Usage:   Rscript 12_single_mediator_mediation.R                  # both types, all retained mediators
#          Rscript 12_single_mediator_mediation.R protein          # proteins only
#          Rscript 12_single_mediator_mediation.R metabolite        # metabolites only
#          Rscript 12_single_mediator_mediation.R protein <index>   # single protein by column index
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
    library(mediation)
})

source("config.R")

# -----------------------------------------------------------------------------
# Analysis constants (kept identical to the original cluster scripts)
# -----------------------------------------------------------------------------
N_SIMS          <- 1000      # nonparametric bootstrap resamples for mediate()
MIN_N_PROTEIN   <- 30000     # retain proteins with complete-case N > this
MIN_N_METAB     <- 350000    # retain metabolites with complete-case N > this
SEED_BASE       <- 42        # per-mediator seed = SEED_BASE + column index

# -----------------------------------------------------------------------------
# Command-line arguments: [mediator_type] [optional single-mediator index]
# -----------------------------------------------------------------------------
args          <- commandArgs(trailingOnly = TRUE)
which_types   <- if (length(args) >= 1) tolower(args[1]) else "both"
single_index  <- if (length(args) >= 2) as.integer(args[2]) else NA_integer_
if (!which_types %in% c("protein", "metabolite", "both")) {
    stop("First argument must be 'protein', 'metabolite', or omitted (both).")
}

cat("R version:", R.version.string, "\n")
cat("mediation package version:", as.character(packageVersion("mediation")), "\n")
cat("Bootstrap resamples:", N_SIMS, "\n\n")

# =============================================================================
# Load shared exposure / outcome / covariate tables
# =============================================================================
cat("Loading exposure, outcome and covariate data...\n")

prs_data        <- fread(cfg$prs_ir_scores)
prs_data$eid    <- prs_data$IID
sample_qc       <- fread(cfg$sample_qc_ids)
chd_cases       <- fread(cfg$chd_cases)
covariates      <- fread(cfg$covariate_data)

# =============================================================================
# Covariate definition, factor coding and imputation
# =============================================================================
get_mode <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return(NA)
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
}

prepare_covariates <- function(dat) {
    dat$sex               <- as.factor(dat$sex)
    dat$highest_education <- as.factor(dat$highest_education)
    dat$employment        <- as.factor(dat$employment)
    dat$smoking_status    <- as.factor(dat$smoking_status)
    dat$sleeplessness     <- as.factor(dat$sleeplessness)

    continuous_covars  <- c("age_at_recruitment", "townsend_deprivation_index", paste0("pc", 1:10))
    categorical_covars <- c("sex", "highest_education", "employment",
                            "physical_activity", "smoking_status", "sleeplessness")

    for (v in continuous_covars) {
        if (v %in% names(dat) && any(is.na(dat[[v]]))) {
            dat[[v]][is.na(dat[[v]])] <- mean(dat[[v]], na.rm = TRUE)
        }
    }
    for (v in categorical_covars) {
        if (v %in% names(dat) && any(is.na(dat[[v]]))) {
            dat[[v]][is.na(dat[[v]])] <- get_mode(dat[[v]])
        }
    }
    dat
}

pc_vars        <- paste0("pc", 1:10)
covariate_vars <- c(
    "age_at_recruitment", "sex", pc_vars,
    "townsend_deprivation_index", "highest_education", "employment",
    "physical_activity", "smoking_status", "sleeplessness"
)
covariate_formula_str <- paste(covariate_vars, collapse = " + ")

# =============================================================================
# Single-mediator mediation model (Imai et al. potential-outcomes framework)
# =============================================================================
# Mediator model:  lm(MEDIATOR  ~ PRS_IR + covariates)
# Outcome model:   glm(CHD ~ PRS_IR + MEDIATOR + covariates, binomial logit)
# mediate(): bootstrap ACME (d0), ADE (z0), total effect (tau), prop. mediated (n0)
run_mediation_single <- function(mediator_name, mediator_label, data,
                                 covariate_str, n_sims = 1000) {

    error_result <- function(msg, n = NA) {
        setNames(
            data.frame(
                mediator = mediator_name,
                n_complete = n,
                ACME_estimate = NA, ACME_pvalue = NA,
                ADE_estimate = NA,  ADE_pvalue = NA,
                total_effect = NA,
                prop_mediated = NA, prop_mediated_pvalue = NA,
                error = msg
            ),
            c(mediator_label, "n_complete", "ACME_estimate", "ACME_pvalue",
              "ADE_estimate", "ADE_pvalue", "total_effect",
              "prop_mediated", "prop_mediated_pvalue", "error")
        )
    }

    temp_data_name <- paste0("med_data_", gsub("[^a-zA-Z0-9]", "_", mediator_name))

    tryCatch({
        required_cols <- c("CHD", "PRS_IR", mediator_name,
                          strsplit(covariate_str, " \\+ ")[[1]])

        df <- as.data.frame(data)[, required_cols, drop = FALSE]
        df <- df[complete.cases(df[[mediator_name]]), ]
        df <- droplevels(df)

        n_complete <- nrow(df)
        if (n_complete < 100) {
            return(error_result("Insufficient sample size", n_complete))
        }

        for (col in names(df)) {
            if (is.factor(df[[col]]) && nlevels(df[[col]]) < 2) {
                return(error_result(paste0("Factor '", col, "' has <2 levels"), n_complete))
            }
        }

        df$MEDIATOR <- df[[mediator_name]]
        assign(temp_data_name, df, envir = .GlobalEnv)

        med_call <- paste0("lm(MEDIATOR ~ PRS_IR + ", covariate_str,
                           ", data = ", temp_data_name, ")")
        out_call <- paste0("glm(CHD ~ PRS_IR + MEDIATOR + ", covariate_str,
                           ", data = ", temp_data_name, ", family = binomial(link = 'logit'))")

        mod_med <- eval(parse(text = med_call))
        mod_out <- eval(parse(text = out_call))

        med_res <- mediate(
            model.m  = mod_med,
            model.y  = mod_out,
            treat    = "PRS_IR",
            mediator = "MEDIATOR",
            boot     = TRUE,
            sims     = n_sims
        )

        rm(list = temp_data_name, envir = .GlobalEnv)

        setNames(
            data.frame(
                mediator = mediator_name,
                n_complete = n_complete,
                ACME_estimate = med_res$d0,   ACME_pvalue = med_res$d0.p,
                ADE_estimate  = med_res$z0,    ADE_pvalue  = med_res$z0.p,
                total_effect  = med_res$tau.coef,
                prop_mediated = med_res$n0,    prop_mediated_pvalue = med_res$n0.p,
                error = ""
            ),
            c(mediator_label, "n_complete", "ACME_estimate", "ACME_pvalue",
              "ADE_estimate", "ADE_pvalue", "total_effect",
              "prop_mediated", "prop_mediated_pvalue", "error")
        )

    }, error = function(e) {
        try(rm(list = temp_data_name, envir = .GlobalEnv), silent = TRUE)
        error_result(e$message)
    })
}

# =============================================================================
# Driver: run all retained mediators of one type
# =============================================================================
run_mediator_type <- function(type) {

    if (type == "protein") {
        mediator_label <- "protein"
        min_n          <- MIN_N_PROTEIN
        mediators      <- fread(cfg$protein_data)
        # drop duplicate-suffixed assay columns if present
        mediators      <- mediators %>% dplyr::select(!ends_with("_dup"))
        out_file       <- file.path(cfg$out_dir,
                                    "single_mediator_protein_mediation_results.csv")
    } else {
        mediator_label <- "metabolite"
        min_n          <- MIN_N_METAB
        mediators      <- fread(cfg$metabolite_data)
        out_file       <- file.path(cfg$out_dir,
                                    "single_mediator_metabolite_mediation_results.csv")
    }

    cat("\n=====================================================\n")
    cat("Mediator type:", type, "\n")
    cat("=====================================================\n")

    # Assemble analysis dataset -------------------------------------------------
    analysis_data <- merge(prs_data[, .(eid, PRS_IR)], covariates, by = "eid")
    analysis_data <- analysis_data[eid %in% sample_qc$eid]
    analysis_data$CHD <- as.integer(analysis_data$eid %in% chd_cases$eid)
    analysis_data <- merge(analysis_data, mediators, by = "eid")
    analysis_data <- prepare_covariates(analysis_data)

    cat("Sample size:", nrow(analysis_data), "\n")
    cat("CHD cases:", sum(analysis_data$CHD), "\n")

    mediator_cols <- setdiff(names(mediators), "eid")
    mediator_cols <- intersect(mediator_cols, names(analysis_data))
    cat("Total", type, "mediators:", length(mediator_cols), "\n")

    # Retain mediators exceeding the complete-case N threshold ------------------
    complete_n <- vapply(mediator_cols,
                         function(m) sum(!is.na(analysis_data[[m]])),
                         integer(1))
    retained <- mediator_cols[complete_n > min_n]
    cat(sprintf("Retained (complete-case N > %d): %d of %d\n",
                min_n, length(retained), length(mediator_cols)))

    # Optional: restrict to a single mediator (by index in the full column list)
    idx_seq <- seq_along(mediator_cols)
    if (!is.na(single_index)) {
        if (single_index < 1 || single_index > length(mediator_cols)) {
            stop(sprintf("Index %d out of range 1-%d", single_index, length(mediator_cols)))
        }
        idx_seq <- single_index
    }

    results <- list()
    for (i in idx_seq) {
        mediator_name <- mediator_cols[i]
        if (!(mediator_name %in% retained)) next

        # Per-mediator seed matches the original array-job seeding
        # (SEED_BASE + array task id, task id = column index).
        set.seed(SEED_BASE + i)

        cat(sprintf("[%d/%d] %s ...\n", i, length(mediator_cols), mediator_name))
        res <- run_mediation_single(
            mediator_name  = mediator_name,
            mediator_label = mediator_label,
            data           = analysis_data,
            covariate_str  = covariate_formula_str,
            n_sims         = N_SIMS
        )
        results[[length(results) + 1L]] <- res
    }

    if (length(results) == 0) {
        cat("No mediators of type", type, "were run.\n")
        return(invisible(NULL))
    }

    results_dt <- rbindlist(results, fill = TRUE)
    fwrite(results_dt, out_file)
    cat("Results saved:", out_file, "(", nrow(results_dt), "rows )\n")
    invisible(results_dt)
}

# =============================================================================
# Execute
# =============================================================================
types_to_run <- if (which_types == "both") c("protein", "metabolite") else which_types
for (tp in types_to_run) run_mediator_type(tp)

cat("\nDone:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
