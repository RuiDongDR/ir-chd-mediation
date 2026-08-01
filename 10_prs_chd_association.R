## =====================================================================
## 10_prs_chd_association.R
##
## Title:   Association of the insulin-resistance (IR) PRS with CHD
##
## Purpose: Logistic regression of coronary heart disease (CHD) on the
##          IR polygenic risk score (PRS) in the QC-passing European
##          cohort. The PRS is standardised to unit SD so effects are
##          reported as the odds ratio per 1-SD increase in the PRS, with
##          95% confidence intervals and P values. Three nested models are
##          fitted:
##            Model 0  CHD ~ PRS_IR                     (unadjusted)
##            Model 1  + age at recruitment, sex, 10 PCs (primary model)
##            Model 2  + Townsend, education, employment, physical
##                       activity, smoking, sleeplessness (fully adjusted)
##
## Inputs:  cfg$prs_ir_scores          - per-individual IR PRS (columns IID, PRS_IR).
##          cfg$covariate_data  - cleaned covariate table (from step 8).
##          cfg$chd_cases   - CHD case EID list (from step 8).
##
## Outputs: cfg$out_dir/prs_ir_chd_association.txt - tidy OR-per-SD table
##          (per-SD OR, 95% CI, P) for the PRS term across the three models.
##          Full coefficient tables are also printed to the console.
##
## Requires: R (>= 4.1)
##           data.table (>= 1.14)
##           tidyverse  (>= 1.3.0)   [dplyr]
##           broom      (>= 0.7)
##
## Run after: 8_phenotypes_covariates.R
##
## ---------------------------------------------------------------------
## CONFIG KEYS USED:
##   out_dir         (existing)  final-output directory
##   prs_ir          (NEW)       per-individual IR PRS scores
##   covariate_file  (NEW)       cleaned covariate table (step 8 output)
##   chd_case_list   (NEW)       CHD case EID list (step 8 output)
##   -> add the NEW keys to config.R before running.
## =====================================================================

source("config.R")

suppressMessages({
  library(data.table)
  library(tidyverse)
  library(broom)
})

## ---- Load inputs -----------------------------------------------------
## The PRS file defines the analysis sample (one row per genotyped
## individual); rename the PLINK IID column to eid.
PRS_IR <- fread(cfg$prs_ir_scores) %>%
  rename(eid = IID) %>%
  select(eid, PRS_IR)

covariates <- fread(cfg$covariate_data)
chd_cases  <- fread(cfg$chd_cases)$eid

## ---- Assemble analysis frame -----------------------------------------
## Merge covariates onto the PRS sample and add the binary CHD outcome
## (1 = case, 0 = control) so every PRS-carrying individual is classified.
df <- PRS_IR %>%
  left_join(covariates, by = "eid") %>%
  mutate(chd = as.integer(eid %in% chd_cases))

## Standardise the PRS to unit SD -> regression coefficients are per-SD.
df$PRS_IR <- as.numeric(scale(df$PRS_IR))

## ---- Recode categorical covariates with explicit reference levels ----
df$sex <- relevel(factor(df$sex), ref = "Female")
df$highest_education <- relevel(factor(df$highest_education), ref = "CSEs or equivalent")
df$employment        <- relevel(factor(df$employment),       ref = "Unemployed")
df$smoking_status    <- relevel(factor(df$smoking_status),   ref = "Never")
df$sleeplessness     <- relevel(factor(df$sleeplessness),    ref = "Never/rarely")

## ---- Helper: fit model and extract OR-per-SD for the PRS term ---------
fit_and_summarise <- function(formula_str, model_label, data) {
  fit <- glm(as.formula(formula_str),
             family = binomial(link = "logit"),
             data = data)
  cat("\n================================================================\n")
  cat(model_label, "\n")
  cat("================================================================\n")
  print(summary(fit)$coefficients)
  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term == "PRS_IR") %>%
    transmute(model = model_label,
              or_per_sd = estimate,
              ci_low = conf.low,
              ci_high = conf.high,
              p_value = p.value)
}

pc_vars <- paste0("pc", 1:10, collapse = " + ")

## ---- Model 0: unadjusted ---------------------------------------------
res0 <- fit_and_summarise(
  "chd ~ PRS_IR",
  "Model 0: CHD ~ PRS-IR (unadjusted)", df)

## ---- Model 1: age + sex + 10 PCs (primary) ---------------------------
res1 <- fit_and_summarise(
  paste("chd ~ PRS_IR + age_at_recruitment + sex +", pc_vars),
  "Model 1: CHD ~ PRS-IR + age + sex + PC1-10 (primary)", df)

## ---- Model 2: fully adjusted -----------------------------------------
res2 <- fit_and_summarise(
  paste("chd ~ PRS_IR + age_at_recruitment + sex +", pc_vars,
        "+ townsend_deprivation_index + highest_education + employment",
        "+ physical_activity + smoking_status + sleeplessness"),
  "Model 2: CHD ~ PRS-IR + age + sex + PC1-10 + lifestyle (fully adjusted)", df)

## ---- Combine and write OR-per-SD summary -----------------------------
prs_summary <- bind_rows(res0, res1, res2)
cat("\n================================================================\n")
cat("IR PRS - CHD association: OR per 1-SD increase (95% CI, P)\n")
cat("================================================================\n")
print(prs_summary)

out_file <- file.path(cfg$out_dir, "prs_ir_chd_association.txt")
write.table(prs_summary, out_file,
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("\nOR-per-SD summary written to: %s\n", out_file))
