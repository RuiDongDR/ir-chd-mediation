## =====================================================================
## config.R  —  central configuration for the IR -> CHD analysis pipeline
##
## Edit the placeholder paths below to match YOUR environment, then every
## analysis script does `source("config.R")` and reads paths from `cfg`.
## No absolute paths are shipped; every path is a "PATH/TO/..." placeholder
## you must set. Entries are grouped as (A) external references, (B) raw
## inputs, and (C) files PRODUCED by earlier pipeline steps (the producing
## step is noted). See README.md for the full run order.
## =====================================================================

cfg <- list(

  ## ==== (A) external reference panels & annotations =====================
  ref_1000g_eur = "PATH/TO/1000G_phase3_EUR",     # PLINK bfile prefix (.bed/.bim/.fam), EUR
  easyqc_1kg_ref= "PATH/TO/1000GP_p3v5_legends",  # 1000G allele-frequency reference for EasyQC MERGE
  ldsc_ld_ref   = "PATH/TO/eur_w_ld_chr",         # LD scores for (Genomic)SEM LDSC
  ldsc_hm3      = "PATH/TO/w_hm3.snplist",         # HapMap3 SNP list used for munging
  genomicsem_ref= "PATH/TO/reference.1000G.maf.0.005.txt", # GenomicSEM sumstats() reference
  dbsnp_ref     = "PATH/TO/dbsnp151_GRCh37.tsv",  # dbSNP151 (hg19) table for rsID annotation
  dbsnp_vcf     = "PATH/TO/dbsnp151_GRCh37.vcf.gz", # dbSNP151 (hg19) VCF for gwaslab rsID mapping
  prscs_ld_ref  = "PATH/TO/ldblk_1kg_eur",        # PRS-CS 1000G EUR LD reference
  nmr_annotation= "PATH/TO/nightingale_biomarker_groups.csv", # cols: Biomarker, Group, Sub.Group

  ## ==== (B) raw inputs ==================================================
  ## Discovery GWAS summary statistics (public; see README for accessions:
  ## fasting insulin GCST90002238, Modified Stumvoll ISI GCST90267573,
  ## pro-insulin Broadaway et al. 2023).
  gwas_fasting_insulin = "PATH/TO/fasting_insulin_sumstats.gz",
  gwas_isi             = "PATH/TO/modified_stumvoll_isi_sumstats.gz",
  gwas_proinsulin      = "PATH/TO/proinsulin_sumstats.gz",

  ## UK Biobank raw data (require an approved UK Biobank application).
  ukb_geno_dir  = "PATH/TO/ukb_imputed",          # imputed genotypes (per-chromosome)
  ukb_pheno     = "PATH/TO/ukb_phenotypes.tab",   # extracted phenotype table
  ukb_nmr       = "PATH/TO/ukb_nightingale_nmr.csv", # raw Nightingale NMR biomarkers (input to step 11)

  ## ==== (C) files produced by earlier pipeline steps ====================
  ## --- genetics (steps 4-5) ---
  f1_gwas_dir   = "results/f1_multivariate_gwas",     # step 4: F1 common-factor GWAS results dir
  f1_loci       = "results/f1_independent_loci.csv",  # step 5: F1 LD-clumped loci
  qsnp_loci     = "results/qsnp_independent_loci.csv",# step 5: Q_SNP LD-clumped loci

  ## --- cohort, phenotypes, PRS (steps 7-10) ---
  sample_qc_ids = "results/ukb_sample_qc_ids.txt",    # step 7: QC-passed European-ancestry ids (col: eid)
  chd_cases     = "results/ukb_chd_cases.txt",        # step 8: CHD case ids (col: eid)
  chd_case_dates= "results/ukb_chd_case_dates.csv",   # step 8: CHD first-occurrence dates
  t2d_cases     = "results/ukb_t2d_cases.txt",        # step 8: type 2 diabetes case ids (Table 1)
  htn_cases     = "results/ukb_htn_cases.txt",        # step 8: hypertension case ids (Table 1)
  hlp_cases     = "results/ukb_hlp_cases.txt",        # step 8: hyperlipidemia case ids (Table 1)
  covariate_data= "results/ukb_cleaned_covariates.txt", # step 8: analysis covariates (age, sex, PCs, ...)
  prs_ir_scores = "results/ukb_prs_ir_scores.txt",    # step 9: per-participant IR PRS (cols: IID, PRS_IR)

  ## --- analysis-ready molecular data (step 11) ---
  protein_data    = "PATH/TO/ukb_olink_proteomics.csv",   # per-sample Olink NPX (cols: eid, <proteins>)
  metabolite_data = "PATH/TO/ukb_nmr_metabolites.csv",    # per-sample NMR metabolites, log1p (cols: eid, <metabolites>)

  ## --- mediation results (steps 12-15) ---
  protein_mediation        = "results/single_mediator_proteins.csv",   # step 12
  metabolite_mediation     = "results/single_mediator_metabolites.csv",# step 12
  serial_mediation_results = "results/serial_mediation_results.csv",   # step 13 (converged pairs)
  level1_pathways          = "results/level1_pathways.csv",            # step 14 (prioritized pathways)
  core_proteins            = "results/core_proteins.csv",              # step 14 (core mediating proteins)
  hdl_pathways             = "results/hdl_subclass_pathways.csv",      # step 15

  ## ==== working / output directories ====================================
  work_dir      = "work",                          # intermediate files
  out_dir       = "results"                         # final outputs (tables, figures)
)

## Create working / output directories if missing.
for (d in c(cfg$work_dir, cfg$out_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
