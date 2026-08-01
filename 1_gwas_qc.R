## =====================================================================
## 1_gwas_qc.R
##
## Purpose:
##   Two-pass EasyQC quality control of the three insulin-resistance (IR)
##   discovery GWAS (fasting insulin, modified Stumvoll ISI, pro-insulin),
##   followed by post-QC filtering and rsID annotation (gwaslab), producing
##   clean, rsID-tagged summary statistics for downstream LDSC / GenomicSEM.
##
##   Pass 1 (standard QC):   drop non-ATCG / monomorphic / MAF<0.5% /
##                           duplicate variants and invalid P/SE/BETA.
##   Pass 2 (alignment QC):  same filters PLUS merge against the 1000
##                           Genomes Phase 3 EUR reference and flag allele
##                           matches (accounting for strand flips).
##   The pass-2 (aligned) output is carried forward, then filtered in R on
##   EAF in (0.005, 0.995), non-missing P, and AlleleMatch == TRUE, and
##   finally annotated with dbSNP rsIDs via gwaslab.
##
## Inputs:
##   - cfg$gwas_fasting_insulin, cfg$gwas_isi, cfg$gwas_proinsulin
##       raw GRCh37 GWAS summary statistics (tab-separated). NOTE: build
##       harmonisation to GRCh37 (liftover from GRCh38 where needed) is a
##       prerequisite performed upstream with MungeSumstats::liftover().
##   - cfg$easyqc_1kg_ref   1000G Phase3 EUR allele-frequency reference (EasyQC MERGE)
##   - cfg$dbsnp_ref        gwaslab dbSNP151 hg19 rsID tsv (1kg_dbsnp151_hg19_auto)
##   - cfg$dbsnp_vcf        dbSNP151 GRCh37 VCF for rsID mapping (gwaslab)
##
## Outputs:
##   - cfg$work_dir/1_EasyQC/<trait>.<...>.tsv.txt                (aligned, EasyQC)
##   - cfg$work_dir/1_EasyQC/<trait>_filtered.tsv.txt             (post-QC filtered)
##   - cfg$work_dir/rsid/<trait>_with_rsid.tsv                    (rsID-annotated)
##
## Requires:
##   - R packages: EasyQC, data.table, dplyr
##   - Python 3 with gwaslab (rsID annotation step)
##   - EasyQC methodological filters: MAF < 0.005, INFO < 0.9 (if present),
##     GRCh37 build, 1000G Phase3 EUR reference.
##
## Run after: build harmonisation / liftover to GRCh37 (upstream)
## =====================================================================
##
## CONFIG KEYS USED:
##   cfg$gwas_fasting_insulin  raw GRCh37 fasting-insulin GWAS sumstats
##   cfg$gwas_isi              raw GRCh37 modified-Stumvoll-ISI GWAS sumstats
##   cfg$gwas_proinsulin       raw GRCh37 pro-insulin GWAS sumstats
##   cfg$easyqc_1kg_ref        1000G Phase3 EUR allele-freq reference for EasyQC MERGE  [NEW]
##   cfg$dbsnp_ref             gwaslab dbSNP151 hg19 rsID reference tsv (1kg_dbsnp151_hg19_auto)
##   cfg$dbsnp_vcf             dbSNP151 GRCh37 VCF for gwaslab rsID mapping             [NEW]
##   cfg$work_dir              intermediate-file directory
## =====================================================================

source("config.R")

suppressMessages({
  library(EasyQC)
  library(data.table)
  library(dplyr)
})

## ---- output locations -------------------------------------------------
qc_dir   <- file.path(cfg$work_dir, "1_EasyQC")
rsid_dir <- file.path(cfg$work_dir, "rsid")
dir.create(qc_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(rsid_dir, recursive = TRUE, showWarnings = FALSE)

## ---- trait metadata ---------------------------------------------------
## Each GWAS uses slightly different column names for the non-effect allele
## and the sample-size column; capture those here so the ECF is built
## faithfully per trait.
traits <- list(
  fasting_insulin = list(
    file        = cfg$gwas_fasting_insulin,
    other_col   = "other_allele",      # non-effect-allele column name
    n_col       = "sample_size"
  ),
  insulin_sensitivity_index = list(
    file        = cfg$gwas_isi,
    other_col   = "other_allele",
    n_col       = "n"
  ),
  pro_insulin = list(
    file        = cfg$gwas_proinsulin,
    other_col   = "reference_allele",  # pro-insulin uses 'reference_allele'
    n_col       = "sample_size"
  )
)

## ---------------------------------------------------------------------
## Build an EasyQC configuration (.ecf) for one trait / one pass.
##   align = FALSE -> Pass 1 (standard QC, no reference merge)
##   align = TRUE  -> Pass 2 (adds 1000G Phase3 EUR alignment + allele match)
## Methodological filters kept verbatim from the analysis: MAF < 0.005,
## monomorphic / non-ATCG / duplicate removal, INFO < 0.9 if an INFO column
## is present.
## ---------------------------------------------------------------------
make_ecf <- function(name, meta, out_path, align) {
  other <- meta$other_col
  ncol  <- meta$n_col
  lines <- c(
    sprintf("DEFINE\t--pathOut %s", out_path),
    sprintf("\t--acolIn chromosome;base_pair_location;effect_allele;%s;effect_allele_frequency;beta;standard_error;p_value;%s", other, ncol),
    "\t--acolInClasses character;numeric;character;character;numeric;numeric;numeric;numeric;numeric",
    "\t--strMissing NA",
    "\t--strSeparator TAB",
    sprintf("EASYIN --fileIn %s", meta$file),
    "START EASYQC",
    ## (1) Non-ATCG alleles
    "CLEAN --rcdClean !(effect_allele%in%c('A','C','G','T')) --strCleanName numDrop_NonATCG_Effect",
    sprintf("CLEAN --rcdClean !(%s%%in%%c('A','C','G','T')) --strCleanName numDrop_NonATCG_Other", other),
    ## (2) Monomorphic SNPs
    "CLEAN --rcdClean (effect_allele_frequency==0)|(effect_allele_frequency==1) --strCleanName numDrop_Monomorph",
    ## (3) MAF < 0.5%
    "ADDCOL --rcdAddCol pmin(effect_allele_frequency,1-effect_allele_frequency,na.rm=TRUE) --colOut MAF",
    "CLEAN --rcdClean MAF < 0.005 --strCleanName numDrop_MAF_below_0.5pct",
    ## SNP identifier (chr:pos)
    "ADDCOL --rcdAddCol paste(chromosome,base_pair_location,sep=\":\") --colOut cptid",
    ## (4) Duplicates / multiallelic
    "CLEANDUPLICATES --colInMarker cptid --strMode keepfirst"
    ## (5) INFO < 0.9 applied automatically when an INFO column is present:
    ##     CLEAN --rcdClean INFO < 0.9 --strCleanName numDrop_LowInfo
  )
  if (align) {
    ## (6) Merge against 1000 Genomes Phase 3 EUR + allele matching
    lines <- c(lines,
      sprintf("MERGE --colInMarker cptid --fileRef %s --colRefMarker cptid --acolIn cptid;ea;oa;eaf --acolInClasses character;character;character;numeric --strSeparator TAB --strRefSuffix .ref --blnInAll 1 --blnRefAll 0 --blnWriteNotInRef 1 --strTag 1kGP3", cfg$easyqc_1kg_ref),
      "CLEAN --rcdClean is.na(cptid) --strCleanName numDrop_AbsentIn1kGP3",
      sprintf("ADDCOL --rcdAddCol ((effect_allele==ea.ref&%s==oa.ref)|(effect_allele==oa.ref&%s==ea.ref)) --colOut AlleleMatch", other, other),
      "CLEAN --rcdClean !AlleleMatch --strCleanName numDrop_AlleleMismatch_1kGP3"
    )
  }
  lines <- c(lines,
    ## Standard invalid-value QC
    "CLEAN --rcdClean is.na(p_value) --strCleanName numDrop_Missing_P",
    "CLEAN --rcdClean p_value<0|p_value>1 --strCleanName numDrop_Invalid_P",
    "CLEAN --rcdClean standard_error<=0 --strCleanName numDrop_Invalid_SE",
    "CLEAN --rcdClean abs(beta)==Inf --strCleanName numDrop_Invalid_BETA",
    "CLEAN --rcdClean is.na(beta) --strCleanName numDrop_Missing_BETA",
    "CLEAN --rcdClean is.na(standard_error) --strCleanName numDrop_Missing_SE",
    sprintf("WRITE --strPrefix %s --strMissing NA", name),
    "STOP EASYQC"
  )
  ecf <- file.path(out_path, sprintf("%s_%s.ecf", name, if (align) "alignment" else "standard"))
  writeLines(lines, ecf)
  ecf
}

## ---------------------------------------------------------------------
## Run both QC passes for every trait, then post-QC filter the aligned
## output and annotate rsIDs.
## ---------------------------------------------------------------------
for (name in names(traits)) {
  meta <- traits[[name]]
  message("==== EasyQC: ", name, " ====")

  ## Pass 1: standard QC (no reference alignment)
  EasyQC(make_ecf(name, meta, qc_dir, align = FALSE))

  ## Pass 2: alignment QC (1000G Phase3 EUR merge + allele match) -> used downstream
  EasyQC(make_ecf(name, meta, qc_dir, align = TRUE))

  ## ---- post-QC filtering in R (on the aligned output) ----------------
  ## EasyQC writes "<prefix>.<inputbasename>.tsv.txt"; pick the aligned file.
  aligned <- list.files(qc_dir, pattern = paste0("^", name, "\\..*\\.txt$"),
                        full.names = TRUE)
  aligned <- aligned[!grepl("notinref|\\.ecf$|\\.rep$", aligned, ignore.case = TRUE)]
  aligned <- aligned[which.max(file.info(aligned)$mtime)]  # most recent = pass 2

  dat <- fread(aligned)
  dat_filtered <- dat %>%
    filter(!is.na(effect_allele_frequency) &
             effect_allele_frequency > 0.005 &
             effect_allele_frequency < 0.995) %>%   # EAF in (0.005, 0.995)
    filter(!is.na(p_value)) %>%                     # non-missing P
    filter(AlleleMatch == TRUE)                     # matched to 1000G EUR

  filt_path <- file.path(qc_dir, sprintf("%s_filtered.tsv.txt", name))
  write.table(dat_filtered, file = filt_path, sep = "\t",
              row.names = FALSE, quote = FALSE, na = "NA")
  message(sprintf("  %s: %d -> %d variants after post-QC filtering",
                  name, nrow(dat), nrow(dat_filtered)))

  ## ---- rsID annotation with gwaslab (Python) -------------------------
  ## gwaslab standardises columns, infers/basic-checks the build, and
  ## assigns rsIDs from a dbSNP151 (hg19) reference. Output columns
  ## (rsID, EA, NEA, BETA, SE, P, N) feed the LDSC munge step.
  rsid_out <- file.path(rsid_dir, sprintf("%s_with_rsid.tsv", name))
  py <- file.path(rsid_dir, "annotate_rsid.py")
  writeLines(c(
    "import sys, gwaslab as gl",
    "infile, outfile, ref_tsv, ref_vcf = sys.argv[1:5]",
    "ss = gl.Sumstats(infile, chrom='chromosome', pos='base_pair_location',",
    "                 ea='effect_allele', nea='other_allele',",
    "                 eaf='effect_allele_frequency', beta='beta',",
    "                 se='standard_error', p='p_value', n='sample_size', sep='\\t')",
    "ss.infer_build()",
    "ss.basic_check()",
    "ss.assign_rsid(ref_rsid_tsv=ref_tsv, ref_rsid_vcf=ref_vcf,",
    "               chr_dict=gl.get_number_to_NC(build='19'), n_cores=2)",
    "ss.data.rename(columns={'EA':'EA','NEA':'NEA'}).to_csv(outfile, sep='\\t', index=False)"
  ), py)
  ## NOTE: pro-insulin/ISI non-effect-allele and N columns differ; the
  ## gwaslab call above assumes standardised 'other_allele'/'sample_size'
  ## names, which the filtered EasyQC output already provides.
  system2("python3", c(py, shQuote(filt_path), shQuote(rsid_out),
                       shQuote(cfg$dbsnp_ref), shQuote(cfg$dbsnp_vcf)))
  message("  rsID-annotated file written: ", rsid_out)
}

message("Step 1 complete: QC + rsID annotation finished for all three IR GWAS.")
