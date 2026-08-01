## =====================================================================
## 3_factor_model.R
##
## Purpose:
##   Assess suitability of the three IR indicators for factor analysis and
##   fit the common-factor model on their genetic covariance structure:
##     (1) Kaiser-Meyer-Olkin (KMO) sampling-adequacy test on the
##         GenomicSEM-derived genetic-correlation matrix;
##     (2) exploratory factor analysis (single factor, promax rotation)
##         reporting the proportion of variance explained;
##     (3) the single common-factor SEM (GenomicSEM commonfactor),
##         reporting the standardized F1 loadings.
##
## Inputs:
##   - cfg$work_dir/LDSCoutput.RData   (GenomicSEM covariance object, step 2)
##
## Outputs:
##   - cfg$out_dir/kmo_result.txt          (overall + per-item MSA)
##   - cfg$out_dir/efa_loadings.txt        (EFA loadings + proportion variance)
##   - cfg$out_dir/commonfactor_F1.txt     (standardized F1 loadings, SEM)
##
## Requires:
##   - R packages: GenomicSEM, psych (KMO), Matrix (nearPD), stats (factanal)
##
## Run after: 2_ldsc_genetic_correlation
## =====================================================================
##
## CONFIG KEYS USED:
##   cfg$work_dir  intermediate-file directory (LDSCoutput.RData)
##   cfg$out_dir   final-output directory
## =====================================================================

source("config.R")

suppressMessages({
  library(GenomicSEM)
  library(psych)
  library(Matrix)
})

trait.names <- c("fasting_insulin", "insulin_sensitivity_index", "pro_insulin")

## ---- load genetic covariance structure -------------------------------
load(file.path(cfg$work_dir, "LDSCoutput.RData"))  # -> LDSCoutput

## Genetic-correlation matrix implied by the multivariable LDSC fit. This is
## the GenomicSEM-derived rg matrix used for the KMO test (values reported in
## the manuscript were approximately:
##   rg(FI, ISI) = -0.835,  rg(FI, ProI) = 0.345,  rg(ISI, ProI) = -0.230).
cor_mat <- cov2cor(LDSCoutput$S)
dimnames(cor_mat) <- list(trait.names, trait.names)

## ---- (1) KMO test -----------------------------------------------------
cat("psych package version:", as.character(packageVersion("psych")), "\n")
kmo <- KMO(cor_mat)
print(kmo)

kmo_tbl <- data.frame(
  metric = c("Overall MSA", paste0("MSA ", trait.names)),
  MSA    = round(c(kmo$MSA, kmo$MSAi[trait.names]), 4)
)
kmo_tbl$above_0.5 <- ifelse(kmo_tbl$MSA > 0.5, "yes", "NO")
write.table(kmo_tbl, file = file.path(cfg$out_dir, "kmo_result.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("\nConclusion:", ifelse(kmo$MSA > 0.5,
    "Overall MSA exceeds the 0.5 cut-off; the three IR indicators are suitable for factor analysis.",
    "Overall MSA falls below the 0.5 cut-off."), "\n")

## ---- (2) Exploratory factor analysis ---------------------------------
## Smooth the genetic covariance matrix to the nearest positive-definite
## matrix, then run a single-factor EFA with promax rotation. The loadings
## print includes the proportion of variance explained by the factor.
Ssmooth <- as.matrix(nearPD(LDSCoutput$S, corr = FALSE)$mat)
dimnames(Ssmooth) <- list(trait.names, trait.names)

EFA <- factanal(covmat = Ssmooth, factors = 1, rotation = "promax")
print(EFA$loadings)

efa_load <- data.frame(trait = trait.names,
                       loading = round(as.numeric(EFA$loadings[, 1]), 4))
prop_var <- sum(EFA$loadings[, 1]^2) / length(trait.names)  # proportion variance
write.table(efa_load, file = file.path(cfg$out_dir, "efa_loadings.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("\nEFA proportion of variance explained by the single factor: %.4f\n",
            prop_var))

## ---- (3) Single common-factor SEM ------------------------------------
## GenomicSEM common-factor model (F1 =~ FI + ISI + ProI). DWLS is the
## primary estimator; ML is reported for comparison. The standardized
## solution (STD_Genotype) gives the standardized F1 loadings.
CommonFactor_DWLS <- commonfactor(covstruc = LDSCoutput, estimation = "DWLS")
CommonFactor_ML   <- commonfactor(covstruc = LDSCoutput, estimation = "ML")

cat("\n--- Common factor (DWLS) ---\n"); print(CommonFactor_DWLS)
cat("\n--- Common factor (ML) ---\n");   print(CommonFactor_ML)

## Extract the standardized factor loadings (F1 =~ indicator rows).
res      <- CommonFactor_DWLS$results
loadings <- res[res$op == "=~", c("lhs", "rhs", "STD_Genotype", "STD_Genotype_SE")]
write.table(loadings, file = file.path(cfg$out_dir, "commonfactor_F1.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("\nStandardized F1 loadings (DWLS):\n"); print(loadings, row.names = FALSE)

message("Step 3 complete: KMO, EFA, and common-factor SEM results written.")
