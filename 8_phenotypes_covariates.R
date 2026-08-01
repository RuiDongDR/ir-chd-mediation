## =====================================================================
## 8_phenotypes_covariates.R
##
## Title:   Coronary heart disease (CHD) outcome and covariate cleaning
##
## Purpose: Two tasks, run over the QC-passing European cohort:
##          (A) Define the CHD outcome (prevalent + incident combined) from:
##              - self-reported non-cancer illness (Field 20002, codes
##                1074 angina / 1075 heart attack),
##              - self-reported vascular/heart problems diagnosed by doctor
##                (Field 6150, codes 1 heart attack / 2 angina),
##              - hospital-inpatient ICD-10 I20-I25 (Fields 41270 / 41202 /
##                41204) and ICD-9 410-414 (Fields 41271 / 41203 / 41205),
##              - algorithmically defined first-occurrence dates for
##                I20-I25 (Fields 131296 / 131298 / 131300 / 131302 /
##                131304 / 131306).
##              A per-case dataframe of first-occurrence dates and timing
##              relative to assessment-centre attendance is also written.
##          (B) Clean the analysis covariate set: age at recruitment, sex,
##              10 genetic principal components, Townsend deprivation index,
##              highest education, employment, physical activity, smoking
##              status, and sleeplessness.
##
## Inputs:  cfg$sample_qc_ids  - QC-passing EID list (from step 7).
##          cfg$ukb_pheno    - extracted UK Biobank phenotype table.
##          disease_case_definition_functions.R  - shared date helpers
##                             (is_valid_date, find_earliest_valid_date,
##                             create_case_dataframe_with_dates,
##                             save_case_dataframe); ship alongside this
##                             script.
##
## Outputs: cfg$chd_cases   - one-column ("eid") CHD case list.
##          cfg$chd_case_dates  - per-case first-occurrence dates and
##                                timing vs attendance.
##          cfg$covariate_data  - cleaned covariate table (one row per
##                                QC-passing participant).
##
## Requires: R (>= 4.1)
##           data.table (>= 1.14)
##           tidyverse  (>= 1.3.0)   [dplyr, tidyr, stringr]
##           lubridate  (>= 1.8)
##
## Run after: 7_ukb_cohort_qc.R
##
## ---------------------------------------------------------------------
## CONFIG KEYS USED:
##   ukb_pheno       (existing)  extracted UKB phenotype table
##   out_dir         (existing)  final-output directory
##   qc_eid_list     (NEW)       QC-passing EID list (step 7 output)
##   chd_case_list   (NEW)       CHD case EID list (output)
##   chd_case_dates  (NEW)       CHD per-case date dataframe (output)
##   covariate_file  (NEW)       cleaned covariate table (output)
##   -> add the NEW keys to config.R before running.
##
## UK Biobank public fields used (kept for reproducibility; see Methods):
##   Outcome: 20002, 6150, 41270, 41202, 41204, 41271, 41203, 41205,
##            131296, 131298, 131300, 131302, 131304, 131306, 53 (attendance)
##   Covariates: 21003 age at recruitment, 31 sex, 22009 PCs,
##            22189 Townsend, 6138 education, 6142 employment,
##            22040 physical activity (MET-min/wk), 20116 smoking,
##            1200 sleeplessness
## =====================================================================

source("config.R")

suppressMessages({
  library(data.table)
  library(tidyverse)
  library(lubridate)
})

## Shared date-validation / date-extraction helpers.
source("disease_case_definition_functions.R")

## ---- Load QC-passing cohort ------------------------------------------
cat("Loading QC-passing EID list...\n")
qced_eids <- fread(cfg$sample_qc_ids)$eid
cat(sprintf("QC-passing individuals: %d\n", length(qced_eids)))

cat("Loading UK Biobank phenotype table...\n")
ukbb_raw  <- fread(cfg$ukb_pheno, na.strings = c("", "NA"))
ukbb_qced <- ukbb_raw %>% filter(eid %in% qced_eids)
cat(sprintf("Cohort after subsetting: %d individuals\n\n", nrow(ukbb_qced)))

## =====================================================================
## (A) CHD OUTCOME DEFINITION  (prevalent + incident combined)
## =====================================================================
## Cases accumulate across sources; setdiff() tracks newly added EIDs per
## criterion, and the running union is the CHD case set.
cases <- character(0)

## Helper: append EIDs matched by a regex across a set of columns.
append_code_cases <- function(df, prefix, regex, cases, label) {
  cols <- grep(paste0("^", prefix), colnames(df), value = TRUE)
  if (length(cols) == 0) {
    cat(sprintf("Warning: %s (%s) columns not found; skipped.\n", label, prefix))
    return(cases)
  }
  hits <- df %>%
    filter(if_any(all_of(cols),
                  ~ !is.na(.) & str_detect(as.character(.), regex))) %>%
    pull(eid) %>% unique()
  n_new <- length(setdiff(hits, cases))
  cases <- unique(c(cases, hits))
  cat(sprintf("%s: %d matched, %d new, %d cumulative\n",
              label, length(hits), n_new, length(cases)))
  cases
}

## --- Self-report: Field 20002 (codes 1074 angina, 1075 heart attack) ---
cols_20002 <- grep("^p20002(_|$)", colnames(ukbb_qced), value = TRUE)
hits_20002 <- ukbb_qced %>%
  filter(if_any(all_of(cols_20002), ~ !is.na(.) & (. == "1074" | . == "1075"))) %>%
  pull(eid) %>% unique()
cases <- unique(c(cases, hits_20002))
cat(sprintf("Field 20002 (self-report CHD 1074/1075): %d matched, %d cumulative\n",
            length(hits_20002), length(cases)))

## --- Self-report: Field 6150 vascular/heart problems diagnosed by doctor
## (code 1 = Heart attack, 2 = Angina). Values are pipe-delimited arrays.
cases <- append_code_cases(ukbb_qced, "p6150",
                           "(^|\\|)(1|2)(\\||$)", cases,
                           "Field 6150 (heart attack/angina, codes 1/2)")

## --- Hospital ICD-10 I20-I25 (Fields 41270, 41202 main, 41204 secondary)
icd10_re <- "(^|\\|)(I20.*|I21.*|I22.*|I23.*|I24.*|I25.*)(\\||$)"
cases <- append_code_cases(ukbb_qced, "p41270", icd10_re, cases, "Field 41270 (ICD-10 I20-I25)")
cases <- append_code_cases(ukbb_qced, "p41202", icd10_re, cases, "Field 41202 (ICD-10 I20-I25)")
cases <- append_code_cases(ukbb_qced, "p41204", icd10_re, cases, "Field 41204 (ICD-10 I20-I25)")

## --- Hospital ICD-9 410-414 (Fields 41271, 41203 main, 41205 secondary) -
icd9_re <- "(^|\\|)(410.*|411.*|412.*|413.*|414.*)(\\||$)"
cases <- append_code_cases(ukbb_qced, "p41271", icd9_re, cases, "Field 41271 (ICD-9 410-414)")
cases <- append_code_cases(ukbb_qced, "p41203", icd9_re, cases, "Field 41203 (ICD-9 410-414)")
cases <- append_code_cases(ukbb_qced, "p41205", icd9_re, cases, "Field 41205 (ICD-9 410-414)")

## --- Algorithmic first-occurrence dates for I20-I25 --------------------
## Fields 131296 (I20), 131298 (I21), 131300 (I22), 131302 (I23),
## 131304 (I24), 131306 (I25). A participant is a case if any field holds a
## valid (non-placeholder) date; is_valid_date() screens UKB date codes.
foc_fields <- c("131296", "131298", "131300", "131302", "131304", "131306")
foc_labels <- c("I20", "I21", "I22", "I23", "I24", "I25")
for (i in seq_along(foc_fields)) {
  cols <- grep(paste0("^p", foc_fields[i]), colnames(ukbb_qced), value = TRUE)
  if (length(cols) == 0) {
    cat(sprintf("Warning: Field %s (%s first occurrence) not found; skipped.\n",
                foc_fields[i], foc_labels[i]))
    next
  }
  hits <- ukbb_qced %>%
    filter(if_any(all_of(cols), is_valid_date)) %>%
    pull(eid) %>% unique()
  n_new <- length(setdiff(hits, cases))
  cases <- unique(c(cases, hits))
  cat(sprintf("Field %s (first-occurrence %s): %d matched, %d new, %d cumulative\n",
              foc_fields[i], foc_labels[i], length(hits), n_new, length(cases)))
}

## --- CHD case summary and output --------------------------------------
cat("\n==============================================================\n")
cat("CHD CASE SUMMARY\n")
cat("==============================================================\n")
cat(sprintf("QC-passing individuals: %d\n", nrow(ukbb_qced)))
cat(sprintf("CHD cases:              %d (%.2f%%)\n",
            length(cases), 100 * length(cases) / nrow(ukbb_qced)))
cat(sprintf("Controls:               %d\n\n", nrow(ukbb_qced) - length(cases)))

write.table(data.frame(eid = cases),
            cfg$chd_cases,
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
cat(sprintf("CHD case EID list written to: %s\n\n", cfg$chd_cases))

## --- Per-case first-occurrence dates and timing vs attendance ----------
case_df <- create_case_dataframe_with_dates(
  df               = ukbb_qced,
  case_eids        = cases,
  date_field_ids   = foc_fields,
  date_field_names = foc_labels,
  attendance_field = "53",       # Field 53: date of attending assessment centre
  attendance_instance = 0
)
write.table(case_df,
            cfg$chd_case_dates,
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
cat(sprintf("CHD per-case date dataframe written to: %s\n\n", cfg$chd_case_dates))

## =====================================================================
## (B) COVARIATE CLEANING
## =====================================================================
## One row per QC-passing participant, labelled by CHD status, with the
## analysis covariate set derived from raw UKB fields.
pheno <- ukbb_qced %>%
  mutate(
    chd_status = factor(as.integer(eid %in% cases),
                        levels = c(1, 0), labels = c("Case", "Control")),
    sex                = factor(p31, labels = c("Female", "Male")),
    age_at_recruitment = p21003_i0                # Field 21003: age at recruitment
  )

## --- 10 genetic principal components (Field 22009, arrays 1-10) --------
pheno <- pheno %>%
  rename_with(~ paste0("pc", 1:10), .cols = paste0("p22009_a", 1:10))

## --- Townsend deprivation index (Field 22189) -------------------------
pheno <- pheno %>% mutate(townsend_deprivation_index = p22189)

## --- Highest education (Field 6138) -----------------------------------
## Take the highest qualification across instances/arrays (lowest code in
## the ordered list = highest qualification).
edu_order  <- c("1", "2", "3", "4", "5", "6", "-7", "-3")
edu_labels <- c(
  "College/University degree",
  "A levels/AS levels or equivalent",
  "O levels/GCSEs or equivalent",
  "CSEs or equivalent",
  "NVQ/HND/HNC or equivalent",
  "Other professional (e.g., nursing, teaching)",
  "None of the above",
  "Prefer not to answer"
)
p6138_cols <- grep("^p6138", names(pheno), value = TRUE)
p6138_highest_df <- pheno %>%
  select(eid, all_of(p6138_cols)) %>%
  pivot_longer(cols = all_of(p6138_cols), values_to = "p6138") %>%
  filter(!is.na(p6138) & p6138 != "") %>%
  separate_rows(p6138, sep = "\\|") %>%
  filter(p6138 %in% edu_order) %>%
  mutate(p6138 = factor(p6138, levels = edu_order, ordered = TRUE)) %>%
  group_by(eid) %>%
  summarize(p6138_highest = as.character(min(p6138)), .groups = "drop") %>%
  mutate(highest_education = factor(p6138_highest, levels = edu_order,
                                    labels = edu_labels, ordered = TRUE))
pheno <- pheno %>%
  select(-matches("^highest_education$")) %>%
  left_join(p6138_highest_df, by = "eid")

## --- Employment (Field 6142, instance 0) ------------------------------
pheno <- pheno %>%
  mutate(
    p6142_raw = as.character(p6142_i0),
    employment = case_when(
      is.na(p6142_raw) | p6142_raw == ""       ~ NA_character_,
      str_detect(p6142_raw, "\\b1\\b")         ~ "Paid employment/self-employed",
      str_detect(p6142_raw, "\\b2\\b")         ~ "Retired",
      str_detect(p6142_raw, "\\b(3|4|5|6)\\b") ~ "Unemployed",
      str_detect(p6142_raw, "\\b7\\b")         ~ "Student",
      p6142_raw == "-7"                        ~ "None of the above",
      p6142_raw == "-3"                        ~ "Prefer not to answer",
      TRUE                                     ~ "Other/Unknown"
    )
  )

## --- Physical activity (Field 22040: summed MET-min/week) -------------
pheno <- pheno %>% rename(physical_activity = p22040_i0)

## --- Smoking status (Field 20116, instance 0) -------------------------
pheno <- pheno %>%
  mutate(smoking_status = case_when(
    p20116_i0 == "2"  ~ "Current",
    p20116_i0 == "1"  ~ "Previous",
    p20116_i0 == "0"  ~ "Never",
    p20116_i0 == "-3" ~ "Prefer not to answer",
    TRUE ~ NA_character_
  ))

## --- Sleeplessness / insomnia (Field 1200, instance 0) ----------------
pheno <- pheno %>%
  mutate(sleeplessness = case_when(
    p1200_i0 == "1"  ~ "Never/rarely",
    p1200_i0 == "2"  ~ "Sometimes",
    p1200_i0 == "3"  ~ "Usually",
    p1200_i0 == "-3" ~ "Prefer not to answer",
    TRUE ~ NA_character_
  ))

## --- Assemble and write cleaned covariate table -----------------------
covariate_cols <- c("eid", "chd_status", "sex", "age_at_recruitment",
                    paste0("pc", 1:10), "townsend_deprivation_index",
                    "highest_education", "employment", "physical_activity",
                    "smoking_status", "sleeplessness")
output <- pheno %>% select(all_of(covariate_cols))

write.table(output,
            cfg$covariate_data,
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("Cleaned covariate table written to: %s\n", cfg$covariate_data))
cat(sprintf("Rows: %d  Columns: %d\n", nrow(output), ncol(output)))
