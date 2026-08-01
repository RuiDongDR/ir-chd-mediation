## =====================================================================
## 4_multivariate_gwas.R
##
## Purpose:
##   Run the GenomicSEM common-factor multivariate GWAS of the latent IR
##   factor (F1) across the genome, with per-SNP heterogeneity (Q_SNP)
##   testing whether each variant acts through the common factor. The
##   insulin-sensitivity-index (ISI) effect sizes are sign-flipped before
##   modelling so that all three indicators load in the same direction and
##   higher F1 corresponds to MORE insulin resistance.
##
## Inputs:
##   - cfg$work_dir/rsid/<trait>_with_rsid.tsv   (from step 1)
##   - cfg$ldsc_hm3        HapMap3 SNP list (munge)
##   - cfg$ldsc_ld_ref     LD scores + weights directory (eur_w_ld_chr)
##   - cfg$genomicsem_ref  GenomicSEM 1000G reference (reference.1000G.maf.0.005.txt)
##
## Outputs:
##   - cfg$work_dir/insulin_sensitivity_index_flipped_with_rsid.tsv
##   - cfg$work_dir/gwas_LDSCoutput.rds        (LDSC covariance for the GWAS)
##   - cfg$work_dir/gwas_sumstats.rds          (merged, formatted per-SNP input)
##   - cfg$out_dir/F1_multivariate_gwas.txt    (F1~SNP effects + Q_SNP heterogeneity)
##   - cfg$out_dir/F1_prscs_sst.txt            (PRS-CS sst file: SNP A1 A2 BETA SE)
##
## Requires:
##   - R package: GenomicSEM (munge, ldsc, sumstats, userGWAS)
##   - SNP filters kept from the analysis: munge info.filter = 0.9 /
##     maf.filter = 0.005; sumstats info.filter = 0.6 / maf.filter = 0.005.
##
## Run after: 1_gwas_qc  (and typically 2_ldsc_genetic_correlation / 3_factor_model)
## =====================================================================
##
## CONFIG KEYS USED:
##   cfg$ldsc_hm3        HapMap3 SNP list for munging (w_hm3.snplist)
##   cfg$ldsc_ld_ref     LD scores / weights directory (eur_w_ld_chr)
##   cfg$genomicsem_ref  GenomicSEM 1000G reference panel for sumstats()   [NEW]
##   cfg$work_dir        intermediate-file directory
##   cfg$out_dir         final-output directory
## =====================================================================

source("config.R")

suppressMessages({
  library(GenomicSEM)
  library(data.table)
})

set.seed(42)  # reproducibility (matches the source analysis)

rsid_dir    <- file.path(cfg$work_dir, "rsid")
trait.names <- c("fasting_insulin", "insulin_sensitivity_index", "pro_insulin")

## ---- (0) Flip the ISI effect sizes -----------------------------------
## Higher ISI = greater insulin SENSITIVITY (opposite of IR). Flipping the
## sign aligns ISI with fasting insulin / pro-insulin so that a higher
## common factor F1 means MORE insulin resistance.
isi_in  <- file.path(rsid_dir, "insulin_sensitivity_index_with_rsid.tsv")
isi_out <- file.path(cfg$work_dir, "insulin_sensitivity_index_flipped_with_rsid.tsv")
isi     <- fread(isi_in)
isi$BETA <- -isi$BETA
fwrite(isi, file = isi_out, sep = "\t", quote = FALSE, na = "NA")

files <- c(file.path(rsid_dir, "fasting_insulin_with_rsid.tsv"),
           isi_out,
           file.path(rsid_dir, "pro_insulin_with_rsid.tsv"))

## Run in the working directory so munge outputs land there.
old_wd <- getwd(); on.exit(setwd(old_wd), add = TRUE)
setwd(cfg$work_dir)

## ---- (1) Munge to HapMap3 --------------------------------------------
munge(files       = files,
      hm3         = cfg$ldsc_hm3,
      trait.names = trait.names,
      info.filter = 0.9,
      maf.filter  = 0.005,
      overwrite   = FALSE)

## ---- (2) Multivariable LDSC ------------------------------------------
munged <- paste0(trait.names, ".sumstats.gz")
LDSCoutput <- ldsc(traits          = munged,
                   sample.prev     = c(NA, NA, NA),
                   population.prev = c(NA, NA, NA),
                   ld              = cfg$ldsc_ld_ref,
                   wld             = cfg$ldsc_ld_ref,
                   trait.names     = trait.names)
setwd(old_wd)
saveRDS(LDSCoutput, file = file.path(cfg$work_dir, "gwas_LDSCoutput.rds"))

## ---- (3) Prepare per-SNP summary statistics --------------------------
## Merge the three GWAS on a common reference for the per-SNP model.
## OLS = TRUE (continuous traits); se.logit = FALSE; INFO >= 0.6, MAF >= 0.005.
gwas_sumstats <- sumstats(files       = as.list(files),
                          ref         = cfg$genomicsem_ref,
                          trait.names = trait.names,
                          se.logit    = c(FALSE, FALSE, FALSE),
                          OLS         = c(TRUE, TRUE, TRUE),
                          linprob     = NULL,
                          N           = c(NA, NA, NA),
                          betas       = NULL,
                          info.filter = 0.6,
                          maf.filter  = 0.005,
                          keep.indel  = FALSE,
                          parallel    = TRUE,
                          cores       = NULL)
saveRDS(gwas_sumstats, file = file.path(cfg$work_dir, "gwas_sumstats.rds"))

## ---- (4) Common-factor multivariate GWAS with per-SNP Q --------------
## F1 loads on the three indicators; each SNP predicts F1. fix_measurement
## holds the measurement model at the common-factor solution; Q_SNP tests
## per-SNP heterogeneity (departure from a pure common-factor effect).
## NOTE: this is genome-wide and computationally intensive; set `cores`
## to match the available hardware.
model <- 'F1 =~ fasting_insulin + insulin_sensitivity_index + pro_insulin
          F1 ~ SNP'

insulin_gwas <- userGWAS(covstruc        = LDSCoutput,
                         SNPs            = gwas_sumstats,
                         estimation      = "DWLS",
                         model           = model,
                         printwarn       = TRUE,
                         sub             = c("F1~SNP"),
                         parallel        = TRUE,
                         cores           = 10,
                         smooth_check    = TRUE,
                         fix_measurement = TRUE,
                         Q_SNP           = TRUE)

## userGWAS returns a list keyed by `sub`; element [[1]] is the F1~SNP table
## (includes est, SE, Z, Pval_Estimate and Q_SNP / Q_SNP_df / Q_SNP_pval).
F1_gwas <- insulin_gwas[[1]]
fwrite(F1_gwas, file = file.path(cfg$out_dir, "F1_multivariate_gwas.txt"),
       sep = "\t", quote = FALSE, na = "NA")

## ---- (5) Format the F1 summary statistics for PRS-CS -----------------
## PRS-CS expects an sst file with columns: SNP A1 A2 BETA SE. Take the
## F1~SNP effect (est) and its SE, drop non-finite / zero SE, and remove
## duplicate SNPs (keep first). This file is PRS-CS's --sst_file in step 6.
prscs_sst <- F1_gwas[, c("SNP", "A1", "A2", "est", "SE")]
colnames(prscs_sst) <- c("SNP", "A1", "A2", "BETA", "SE")
prscs_sst <- prscs_sst[is.finite(prscs_sst$SE) & prscs_sst$SE > 0, ]
prscs_sst <- prscs_sst[!duplicated(prscs_sst$SNP), ]
fwrite(prscs_sst, file = file.path(cfg$out_dir, "F1_prscs_sst.txt"),
       sep = "\t", quote = FALSE, na = "NA")

message("Step 4 complete: F1 multivariate GWAS + PRS-CS sst file written.")
