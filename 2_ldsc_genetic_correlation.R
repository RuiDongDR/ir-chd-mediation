## =====================================================================
## 2_ldsc_genetic_correlation.R
##
## Purpose:
##   Munge the three rsID-annotated IR GWAS to HapMap3 SNPs and run
##   multivariable LD-score regression in GenomicSEM to obtain the genetic
##   covariance structure of the three IR indicators: SNP heritabilities,
##   pairwise genetic correlations (rg), and the cross-trait (bivariate)
##   LDSC intercepts. The resulting covariance object (LDSCoutput) is the
##   input to the factor model (step 3) and the multivariate GWAS (step 4).
##
## Inputs:
##   - cfg$work_dir/rsid/<trait>_with_rsid.tsv   (from step 1)
##   - cfg$ldsc_hm3      HapMap3 SNP list (munge --merge / hm3)
##   - cfg$ldsc_ld_ref   LD scores + weights directory (eur_w_ld_chr)
##
## Outputs:
##   - cfg$work_dir/<trait>.sumstats.gz          (munged sumstats)
##   - cfg$work_dir/LDSCoutput.RData             (GenomicSEM covariance object)
##   - cfg$out_dir/ldsc_h2.txt                   (SNP h2 per trait)
##   - cfg$out_dir/ldsc_rg_matrix.txt            (pairwise genetic correlations)
##   - cfg$out_dir/ldsc_intercepts.txt           (cross-trait LDSC intercepts)
##
## Requires:
##   - R package: GenomicSEM (munge, ldsc)
##   - SNP filters kept from the analysis: info.filter = 0.9, maf.filter = 0.005
##
## Run after: 1_gwas_qc
## =====================================================================
##
## CONFIG KEYS USED:
##   cfg$ldsc_hm3     HapMap3 SNP list for munging (w_hm3.snplist)
##   cfg$ldsc_ld_ref  LD scores / weights directory (eur_w_ld_chr)
##   cfg$work_dir     intermediate-file directory
##   cfg$out_dir      final-output directory
## =====================================================================

source("config.R")

suppressMessages(library(GenomicSEM))

## ---- inputs -----------------------------------------------------------
rsid_dir    <- file.path(cfg$work_dir, "rsid")
trait.names <- c("fasting_insulin", "insulin_sensitivity_index", "pro_insulin")
files       <- file.path(rsid_dir, paste0(trait.names, "_with_rsid.tsv"))

## GenomicSEM writes "<trait.name>.sumstats.gz" into the working directory;
## run munge from cfg$work_dir so the outputs land there.
old_wd <- getwd(); on.exit(setwd(old_wd), add = TRUE)
setwd(cfg$work_dir)

## ---- munge to HapMap3 -------------------------------------------------
## rsID-annotated files carry columns rsID / EA / NEA / N / P (a1/a2 inferred
## by GenomicSEM). Filters: INFO >= 0.9, MAF >= 0.005.
munge(files       = files,
      hm3         = cfg$ldsc_hm3,
      trait.names = trait.names,
      info.filter = 0.9,
      maf.filter  = 0.005,
      overwrite   = FALSE)

## ---- multivariable LDSC ----------------------------------------------
munged          <- paste0(trait.names, ".sumstats.gz")
sample.prev     <- c(NA, NA, NA)   # continuous traits
population.prev <- c(NA, NA, NA)   # continuous traits

LDSCoutput <- ldsc(traits          = munged,
                   sample.prev     = sample.prev,
                   population.prev = population.prev,
                   ld              = cfg$ldsc_ld_ref,
                   wld             = cfg$ldsc_ld_ref,
                   trait.names     = trait.names)

setwd(old_wd)
save(LDSCoutput, file = file.path(cfg$work_dir, "LDSCoutput.RData"))
print(LDSCoutput)

## ---- extract and report -----------------------------------------------
## S is the genetic covariance matrix (SNP h2 on the diagonal); cov2cor(S)
## gives the pairwise genetic correlations. I holds the LDSC intercepts:
## diagonal = univariate (mean chi-square) intercepts, off-diagonal =
## cross-trait intercepts (index of sample overlap).
S  <- LDSCoutput$S
rg <- cov2cor(S)
dimnames(S) <- dimnames(rg) <- list(trait.names, trait.names)

h2 <- data.frame(trait = trait.names, SNP_h2 = diag(S))
write.table(h2, file = file.path(cfg$out_dir, "ldsc_h2.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(round(rg, 4), file = file.path(cfg$out_dir, "ldsc_rg_matrix.txt"),
            sep = "\t", quote = FALSE, col.names = NA)

I <- LDSCoutput$I
dimnames(I) <- list(trait.names, trait.names)
write.table(round(I, 4), file = file.path(cfg$out_dir, "ldsc_intercepts.txt"),
            sep = "\t", quote = FALSE, col.names = NA)

cat("\nSNP heritabilities:\n");            print(h2, row.names = FALSE)
cat("\nGenetic correlation matrix:\n");    print(round(rg, 4))
cat("\nLDSC intercepts (diag=univariate, off-diag=cross-trait):\n"); print(round(I, 4))

message("Step 2 complete: LDSCoutput and genetic-correlation summaries written.")
