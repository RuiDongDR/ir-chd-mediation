# Analysis code — genetic insulin resistance → CHD via proteomic–metabolomic mediation

Release code for the study *"Dissecting the Genetic Link Between Insulin Resistance
and Coronary Heart Disease: A Multi-Trait GWAS and Proteomic–Metabolomic Mediation
Analysis."* The scripts reproduce the analysis end-to-end given access to the
required data. **No absolute paths or dataset-specific locations are shipped** — all
paths are placeholders you set once in `config.R` (R steps) and `config.sh`
(shell steps).

## How to run

1. Install the required software (below).
2. Edit `config.R` and `config.sh`, replacing every `PATH/TO/...` placeholder with a
   path on your system. Files produced by earlier steps default to `results/…`.
3. Run the scripts **in the numbered order**. R: `Rscript N_name.R`; shell:
   `bash N_name.sh`; Python: `python N_name.py`.

## Pipeline

| # | Script | Does |
|---|--------|------|
| 1 | `1_gwas_qc.R` | Two-pass EasyQC of the three IR GWAS + rsID annotation (gwaslab) |
| 2 | `2_ldsc_genetic_correlation.R` | Multivariable LDSC (GenomicSEM): SNP h², pairwise rg, cross-trait intercepts |
| 3 | `3_factor_model.R` | KMO test, exploratory factor analysis, single common-factor SEM (F1 loadings) |
| 4 | `4_multivariate_gwas.R` | GenomicSEM common-factor GWAS of F1 + per-SNP Q; latent-factor effective N |
| 5 | `5_ld_clumping.sh` | PLINK LD clumping → independent F1 and Q_SNP loci |
| 6 | `6_prscs_weights.sh` | PRS-CS posterior SNP weights for the IR factor (per autosome) |
| 7 | `7_ukb_cohort_qc.R` | UK Biobank European-ancestry sample QC |
| 8 | `8_phenotypes_covariates.R` | CHD outcome definition + covariate cleaning |
| 9 | `9_compute_prs.sh` (+ `9_compute_prs.R`) | Apply PRS-CS weights to UKB genotypes → per-individual IR PRS |
| 10 | `10_prs_chd_association.R` | Logistic regression: IR PRS → CHD (OR per SD) |
| 11 | `11_nmr_qc.R` | Nightingale NMR technical-variation removal + outlier flagging (ukbnmr) |
| 12 | `12_single_mediator_mediation.R` | Single-mediator causal mediation (per protein / per metabolite) |
| 13 | `13_serial_mediation.R` | Serial mediation (IR PRS → protein → metabolite → CHD), two-pass |
| 14 | `14_pathway_prioritization.R` | FDR + effect-size filter; rank core mediating proteins |
| 15 | `15_hdl_subclass_pathways.R` | Localize core-protein signal to Large/Very-large HDL subclasses |
| 16 | `16_figures_tables.R` | All manuscript figures (PNG + 300-dpi PDF) and tables |

`config.R` / `config.sh` — central path configuration (edit before running).

## Required software

- **R** (≥ 4.4) with: `GenomicSEM`, `lavaan`, `mediation`, `ukbnmr`, `data.table`,
  `dplyr`, `readr`, `tidyr`, `broom`, `ggplot2`, `reshape2`, `png`, `EasyQC`, `gwaslab` (Python, via reticulate/CLI).
- **PLINK 1.9** (clumping) and **PLINK 2.0** (`--score`).
- **PRS-CS** (`PRScs.py`) with its 1000 Genomes EUR LD reference.

## Data

- **Discovery GWAS (public):** fasting insulin — GWAS Catalog `GCST90002238`;
  Modified Stumvoll ISI — `GCST90267573`; pro-insulin — Broadaway et al. 2023,
  *Am J Hum Genet* (MAGIC). All European-ancestry, excluding UK Biobank.
- **Reference panel:** 1000 Genomes Phase 3 European.
- **UK Biobank** (individual-level genotypes, phenotypes, Olink proteomics,
  Nightingale NMR metabolomics): available to approved researchers via the UK
  Biobank Access Management System (https://www.ukbiobank.ac.uk/). An approved
  application is required; these data cannot be redistributed.

## Compute notes

- **Steps 12–13** (single-mediator and serial mediation) are **compute-intensive**
  (thousands of models, each with 1,000 bootstrap resamples). The scripts run the same
  models everywhere — loop over targets, or pass a single target via a command-line
  argument — and can be parallelized across targets to suit your setup.

## Notes

- Set seeds are preserved from the original analysis for reproducibility.

## Authors

Rui Dong¹\*, Jiao Fu²\*, Bowen Zhang³, Qingling Zhang⁴, Yang Yan³, Guoliang Li⁵

\* Rui Dong and Jiao Fu contributed equally.
Corresponding authors: Guoliang Li and Yang Yan.

1. Hetao Institute of Mathematics and Interdisciplinary Sciences, Shenzhen, China
2. Department of Endocrinology and Metabolism, the First Affiliated Hospital of Xi'an Jiaotong University, Xi'an, China
3. Department of Cardiovascular Surgery, the First Affiliated Hospital of Xi'an Jiaotong University, Xi'an, China
4. Tianjin Key Laboratory of Ionic-Molecular Function of Cardiovascular Disease, Department of Cardiology, Tianjin Institute of Cardiology, Second Hospital of Tianjin Medical University, Tianjin, China
5. Department of Cardiovascular Medicine, the First Affiliated Hospital of Xi'an Jiaotong University, Xi'an, China

## License

Released under the **Creative Commons Attribution-NonCommercial 4.0 International
License (CC BY-NC 4.0)** — see [`LICENSE`](LICENSE). You may share and adapt this
material for **non-commercial** purposes with appropriate attribution.

© 2026 Rui Dong, Jiao Fu, Bowen Zhang, Qingling Zhang, Yang Yan, Guoliang Li.
