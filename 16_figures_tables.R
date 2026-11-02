## =====================================================================
## 16_figures_tables.R
##
## Title:   Manuscript figures and tables (IR -> CHD study)
##
## Purpose: Generate every publication figure and table from the finalized
##          upstream result files. Faithfully reproduces the plotting and
##          table-formatting logic used in the analysis, writing all outputs
##          to a single results directory.
##
## Inputs   (set these paths in config.R):
##   cfg$f1_gwas_dir            directory of the multivariate (common-factor F1)
##                              GWAS .rds results; largest .rds is the full set
##   cfg$f1_loci                F1 independent GWS loci (LD-clumped) CSV
##   cfg$qsnp_loci              Q_SNP independent loci (LD-clumped) CSV
##   cfg$prs_ir_scores                 IR polygenic risk score table (per-individual)
##   cfg$covariate_data         cleaned cohort covariate table
##   cfg$chd_cases              coronary heart disease case IDs
##   cfg$t2d_cases              type 2 diabetes case IDs
##   cfg$htn_cases              hypertension case IDs
##   cfg$hlp_cases              hyperlipidemia case IDs
##   cfg$protein_mediation      single-mediator (protein) results CSV
##   cfg$metabolite_mediation   single-mediator (metabolite) results CSV
##   cfg$serial_mediation_results serial-mediation input (proteins carried forward)
##   cfg$core_proteins          core mediating proteins summary CSV
##   cfg$level1_pathways        level-1 serial-mediation pathways CSV
##   cfg$hdl_pathways           core-protein x HDL-subclass pathways CSV
##
## Outputs  (written to cfg$out_dir):
##   Figures: Figure 1 (genetic-correlation matrix), Figure 2 (Manhattan),
##            Figure 3 (IR PRS distribution), Figure 4 (core-protein landscape),
##            Figure 5 (core-protein x HDL-subclass heatmap), Figure S1 (QQ).
##            Each figure is saved as BOTH PNG and PDF at 300 dpi.
##   Tables:  Table 1 (baseline characteristics), Table 2a/2b (top mediators),
##            Table 3 (core proteins); supplementary Tables S1-S6. CSV.
##
## Requires: ggplot2, data.table, png, reshape2
##
## Run after: all prior steps.
##
## CONFIG KEYS USED:
##   out_dir, f1_gwas_dir, f1_loci, qsnp_loci, ir_prs, ukb_covariates,
##   chd_cases, t2d_cases, htn_cases, hlp_cases, protein_mediation,
##   metabolite_mediation, serial_mediation_input, core_proteins,
##   level1_pathways, hdl_pathways
## =====================================================================

source("config.R")

suppressMessages({
  library(ggplot2)
  library(data.table)
  library(png)
  library(reshape2)
})

set.seed(1)   # reproducible Manhattan baseline thinning

## single output location for ALL figures and tables
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)


## =====================================================================
## Figure 1 — genetic-correlation matrix (GenomicSEM ldsc() output)
## =====================================================================
cor_mat <- matrix(c(1,       -0.8351,  0.3450,
                    -0.8351,  1,      -0.2301,
                    0.3450,  -0.2301,  1),
                  nrow = 3, byrow = TRUE)
trait_labels <- c("Fasting insulin", "Insulin sensitivity\nindex", "Pro-insulin")
colnames(cor_mat) <- rownames(cor_mat) <- trait_labels

# Long format for ggplot; keep factor order so the matrix reads top-left to bottom-right
df <- melt(cor_mat, varnames = c("Var1", "Var2"), value.name = "rg")
df$Var1 <- factor(df$Var1, levels = trait_labels)
df$Var2 <- factor(df$Var2, levels = rev(trait_labels))

p <- ggplot(df, aes(x = Var1, y = Var2, fill = rg)) +
  geom_tile(color = "white", linewidth = 1.2) +
  geom_text(aes(label = sprintf("%.3f", rg),
                color = abs(rg) > 0.6),
            size = 5, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = c(`TRUE` = "white", `FALSE` = "black")) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(-1, 1),
    name = expression(r[g]),
    breaks = c(-1, -0.5, 0, 0.5, 1)
  ) +
  coord_fixed() +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.grid = element_blank(),
    legend.key.height = unit(1.4, "cm")
  )

ggsave(file.path(cfg$out_dir, "20260613_Figure1_genetic_correlation_matrix.png"),
       p, width = 7, height = 6, dpi = 300)
ggsave(file.path(cfg$out_dir, "20260613_Figure1_genetic_correlation_matrix.pdf"),
       p, width = 7, height = 6, dpi = 300, device = cairo_pdf)


## =====================================================================
## Load the multivariate GWAS results (largest .rds in the results folder)
## =====================================================================
rds_files <- list.files(cfg$f1_gwas_dir, pattern = "\\.rds$", full.names = TRUE, ignore.case = TRUE)
stopifnot("No .rds found in cfg$f1_gwas_dir" = length(rds_files) >= 1)
gwas_path <- rds_files[which.max(file.info(rds_files)$size)]   # largest file = full results
cat("Loading:", gwas_path, "\n")

res  <- readRDS(gwas_path)
gwas <- if (is.list(res) && !is.data.frame(res)) res[[1]] else res
rm(res); invisible(gc())

cat("dim:", nrow(gwas), "x", ncol(gwas), "\n")


## ---- Helper functions ----

# Genomic inflation factor (lambda_GC)
lambda_gc <- function(p, df = 1) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  median(qchisq(p, df = df, lower.tail = FALSE)) / qchisq(0.5, df = df)
}

# Greedy distance-based clumping: keep lead (lowest-P) SNP per +/- window
clump_df <- function(sig, pcol, window = 5e5) {
  sig <- sig[is.finite(sig[[pcol]]) & sig[[pcol]] > 0, ]
  sig <- sig[order(sig[[pcol]]), ]
  lc <- integer(0); lb <- numeric(0); keep <- logical(nrow(sig))
  for (i in seq_len(nrow(sig))) {
    if (!any(lc == sig$CHR[i] & abs(lb - sig$BP[i]) <= window)) {
      keep[i] <- TRUE; lc <- c(lc, sig$CHR[i]); lb <- c(lb, sig$BP[i])
    }
  }
  out <- sig[keep, ]
  out[order(out$CHR, out$BP), ]
}

# Build a thinned Manhattan plotting frame (all suggestive SNPs + sampled background)
manhattan_df <- function(CHR, BP, P, keep_below = 1e-3, n_bg = 3e5, seed = 1) {
  ok <- is.finite(P) & P > 0 & P <= 1
  CHR <- CHR[ok]; BP <- BP[ok]; P <- P[ok]
  set.seed(seed)
  top <- which(P <  keep_below)
  bg  <- which(P >= keep_below)
  if (length(bg) > n_bg) bg <- sample(bg, n_bg)
  idx <- c(top, bg)
  chrs    <- sort(unique(CHR))
  chr_len <- tapply(BP, CHR, max)[as.character(chrs)]
  offset  <- c(0, cumsum(as.numeric(chr_len))[-length(chrs)])
  names(offset) <- as.character(chrs)
  d <- data.frame(CHR = CHR[idx], BP = BP[idx], P = P[idx])
  d$pos  <- d$BP + offset[as.character(d$CHR)]
  d$logp <- -log10(d$P)
  axis_df <- data.frame(CHR = chrs, center = as.numeric(offset) + as.numeric(chr_len) / 2)
  list(d = d, axis = axis_df)
}

# QQ plotting frame (keep the full tail + a sample of the bulk for speed)
qq_df <- function(p, n_top = 5000, n_bg = 95000, seed = 1) {
  p <- sort(p[is.finite(p) & p > 0 & p <= 1])   # ascending: most significant first
  n <- length(p)
  obs <- -log10(p)
  exp <- -log10(ppoints(n))
  set.seed(seed)
  idx <- if (n > (n_top + n_bg)) c(1:n_top, sort(sample((n_top + 1):n, n_bg))) else seq_len(n)
  data.frame(exp = exp[idx], obs = obs[idx])
}


## ---- Significance counts and genomic inflation ----
P   <- gwas$Pval_Estimate
QP  <- gwas$Q_SNP_pval
gws <- 5e-8

lam_f1   <- lambda_gc(P)
mean_chi <- mean(qchisq(P[is.finite(P) & P > 0 & P <= 1], 1, lower.tail = FALSE))

cat("Total SNPs:            ", nrow(gwas), "\n")
cat("F1    GWS SNPs (P<5e-8):", sum(P  < gws, na.rm = TRUE), "\n")
cat("Q_SNP GWS SNPs (P<5e-8):", sum(QP < gws, na.rm = TRUE), "\n")
cat(sprintf("F1 lambda_GC = %.3f | mean chi^2 = %.3f\n", lam_f1, mean_chi))


## =====================================================================
## Figure 2 — Manhattan plot of the F1 common-factor association
## =====================================================================
mh <- manhattan_df(gwas$CHR, gwas$BP, gwas$Pval_Estimate)

# --- thin the overplotted Manhattan baseline for speed & small file ---
# keep every SNP with -log10(P) > 2 (all peaks); subsample the dense baseline to 200k
set.seed(1)
.keep <- mh$d$logp > 2
.bulk <- which(!.keep)
.samp <- if (length(.bulk) > 200000) sample(.bulk, 200000) else .bulk
mh_thin <- mh$d[sort(c(which(.keep), .samp)), ]
cat(sprintf("Manhattan points: %d -> %d after thinning\n", nrow(mh$d), nrow(mh_thin)))

p_manhattan <- ggplot(mh_thin, aes(pos, logp, color = factor(CHR %% 2))) +
  geom_point(size = 0.35, alpha = 0.8) +
  scale_color_manual(values = c(`0` = "#08519C", `1` = "#6BAED6"), guide = "none") +
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "red", linewidth = 0.4) +
  scale_x_continuous(breaks = mh$axis$center, labels = mh$axis$CHR, expand = c(0.01, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Chromosome", y = expression(-log[10](italic(P)))) +
  theme_classic(base_size = 12) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(size = 7))

# Manhattan: raster PNG at 300 dpi; PDF embeds the same 300-dpi PNG
ggsave(file.path(cfg$out_dir, "20260613_Figure2_multivariate_GWAS_manhattan.png"),
       p_manhattan, width = 11, height = 4.5, dpi = 300)
pdf(file.path(cfg$out_dir, "20260613_Figure2_multivariate_GWAS_manhattan.pdf"), width = 11, height = 4.5)
grid::grid.raster(png::readPNG(file.path(cfg$out_dir, "20260613_Figure2_multivariate_GWAS_manhattan.png")),
                  width = grid::unit(1, "npc"), height = grid::unit(1, "npc"))
invisible(dev.off())


## =====================================================================
## Figure S1 — QQ plot with genomic inflation lambda
## =====================================================================
qd      <- qq_df(gwas$Pval_Estimate)
lam_lab <- sprintf("lambda[GC] == %.2f", lam_f1)

p_qq <- ggplot(qd, aes(exp, obs)) +
  geom_abline(slope = 1, intercept = 0, color = "grey60", linewidth = 0.4) +
  geom_point(size = 0.5, color = "#08519C", alpha = 0.7) +
  annotate("text", x = 0, y = max(qd$obs) * 0.95, hjust = 0,
           label = lam_lab, parse = TRUE, size = 4.5) +
  labs(x = expression(Expected~-log[10](italic(P))),
       y = expression(Observed~-log[10](italic(P)))) +
  theme_classic(base_size = 12) +
  theme(panel.grid = element_blank())

# QQ: raster PNG at 300 dpi; PDF embeds the same 300-dpi PNG
ggsave(file.path(cfg$out_dir, "20260613_FigureS1_QQ_plot.png"),
       p_qq, width = 5.5, height = 5.5, dpi = 300)
pdf(file.path(cfg$out_dir, "20260613_FigureS1_QQ_plot.pdf"), width = 5.5, height = 5.5)
grid::grid.raster(png::readPNG(file.path(cfg$out_dir, "20260613_FigureS1_QQ_plot.png")),
                  width = grid::unit(1, "npc"), height = grid::unit(1, "npc"))
invisible(dev.off())


## =====================================================================
## Table S1 — F1 independent genome-wide-significant loci (LD-clumped)
## Loci computed upstream by PLINK LD-clumping; here we format and save.
## =====================================================================
f1_loci <- fread(cfg$f1_loci)
fwrite(f1_loci, file.path(cfg$out_dir, "20260613_TableS1_F1_loci.csv"))
cat("Table S1:", nrow(f1_loci), "independent F1 loci (LD-clumped)\n")


## =====================================================================
## Table S2 — Q_SNP independent loci (LD-clumped)
## =====================================================================
q_loci <- fread(cfg$qsnp_loci)
fwrite(q_loci, file.path(cfg$out_dir, "20260613_TableS2_QSNP_loci.csv"))
cat("Table S2:", nrow(q_loci), "independent Q_SNP loci (LD-clumped)\n")


## =====================================================================
## Figure 3 — distribution of the IR polygenic risk score
## =====================================================================
prs <- fread(cfg$prs_ir_scores)
cat("N individuals:", nrow(prs), "\n")
prs_mean <- mean(prs$PRS_IR); prs_sd <- sd(prs$PRS_IR)
cat(sprintf("mean = %.4f | SD = %.4f | median = %.4f\n", prs_mean, prs_sd, median(prs$PRS_IR)))

p_prs <- ggplot(prs, aes(x = PRS_IR)) +
  geom_histogram(aes(y = after_stat(density)), bins = 80,
                 fill = "#6BAED6", color = "white", linewidth = 0.1) +
  stat_function(fun = dnorm, args = list(mean = prs_mean, sd = prs_sd),
                color = "#08519C", linewidth = 0.8) +
  geom_vline(xintercept = prs_mean, linetype = "dashed", color = "red", linewidth = 0.4) +
  labs(x = "IR polygenic risk score", y = "Density") +
  theme_classic(base_size = 12) +
  theme(panel.grid = element_blank())

ggsave(file.path(cfg$out_dir, "20260613_Figure3_IR_PRS_distribution.png"),
       p_prs, width = 6, height = 4.5, dpi = 300)
ggsave(file.path(cfg$out_dir, "20260613_Figure3_IR_PRS_distribution.pdf"),
       p_prs, width = 6, height = 4.5, dpi = 300, device = cairo_pdf)


## =====================================================================
## Table 1 — baseline characteristics by CHD status (full cohort)
## =====================================================================
cov <- fread(cfg$covariate_data)
prs <- fread(cfg$prs_ir_scores)[, .(eid = IID, PRS_IR)]
chd <- fread(cfg$chd_cases)
t2d <- fread(cfg$t2d_cases)
htn <- fread(cfg$htn_cases)
hlp <- fread(cfg$hlp_cases)

dt <- merge(cov, prs, by = "eid", all.x = TRUE)
dt[, physical_activity := suppressWarnings(as.numeric(physical_activity))]
dt[, CHD := factor(fifelse(eid %in% chd$eid, "CHD case", "Control"), levels = c("Control", "CHD case"))]
dt[, T2D            := as.integer(eid %in% t2d$eid)]
dt[, Hypertension   := as.integer(eid %in% htn$eid)]
dt[, Hyperlipidemia := as.integer(eid %in% hlp$eid)]
grp <- dt$CHD
cat("N =", nrow(dt), "| CHD cases =", sum(grp == "CHD case"),
    sprintf("(%.1f%%)\n", 100 * mean(grp == "CHD case")))

# ---- formatting helpers ----
fmt_cont <- function(x) sprintf("%.1f (%.1f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
fmt_cat  <- function(x, lvl) { d <- x[!is.na(x)]; n <- sum(d == lvl); sprintf("%d (%.1f%%)", n, 100 * n / length(d)) }
fmt_p    <- function(p) if (is.na(p)) "" else if (p < 1e-300) "<1e-300" else if (p < 0.001) sprintf("%.1e", p) else sprintf("%.3f", p)

rows <- list(); add <- function(...) rows[[length(rows) + 1]] <<- data.frame(..., stringsAsFactors = FALSE)

cont_vars <- c("Age at recruitment (years)" = "age_at_recruitment", "BMI (kg/m2)" = "BMI",
               "Townsend deprivation index" = "townsend_deprivation_index",
               "Physical activity (MET-min/wk)" = "physical_activity", "IR PRS" = "PRS_IR")
cat_vars  <- c("Sex" = "sex", "Smoking status" = "smoking_status",
               "Alcohol intake frequency" = "alcohol_intake_frequency", "Highest education" = "highest_education")
bin_vars  <- c("Type 2 diabetes" = "T2D", "Hypertension" = "Hypertension", "Hyperlipidemia" = "Hyperlipidemia")

for (nm in names(cont_vars)) {
  v <- dt[[cont_vars[nm]]]
  p <- tryCatch(t.test(v ~ grp)$p.value, error = function(e) NA)
  add(Variable = paste0(nm, ", mean (SD)"), Level = "", Overall = fmt_cont(v),
      Control = fmt_cont(v[grp == "Control"]), CHD_case = fmt_cont(v[grp == "CHD case"]), P = fmt_p(p))
}
for (nm in names(cat_vars)) {
  v <- dt[[cat_vars[nm]]]
  p <- tryCatch(chisq.test(table(v, grp))$p.value, error = function(e) NA)
  add(Variable = paste0(nm, ", n (%)"), Level = "", Overall = "", Control = "", CHD_case = "", P = fmt_p(p))
  for (lvl in sort(unique(v[!is.na(v)])))
    add(Variable = "", Level = lvl, Overall = fmt_cat(v, lvl),
        Control = fmt_cat(v[grp == "Control"], lvl), CHD_case = fmt_cat(v[grp == "CHD case"], lvl), P = "")
}
for (nm in names(bin_vars)) {
  v <- dt[[bin_vars[nm]]]
  p <- tryCatch(chisq.test(table(v, grp))$p.value, error = function(e) NA)
  add(Variable = paste0(nm, ", n (%)"), Level = "", Overall = fmt_cat(v, 1),
      Control = fmt_cat(v[grp == "Control"], 1), CHD_case = fmt_cat(v[grp == "CHD case"], 1), P = fmt_p(p))
}

table1 <- do.call(rbind, rows)
colnames(table1) <- c("Variable", "Level", sprintf("Overall (N=%d)", nrow(dt)),
                      sprintf("Control (n=%d)", sum(grp == "Control")),
                      sprintf("CHD case (n=%d)", sum(grp == "CHD case")), "P")

write.csv(table1, file.path(cfg$out_dir, "20260613_Table1_baseline_characteristics.csv"), row.names = FALSE)
cat("Wrote Table 1 (", nrow(table1), "rows )\n")


## =====================================================================
## Load single-mediator results
## =====================================================================
prot  <- fread(cfg$protein_mediation)
metab <- fread(cfg$metabolite_mediation)

## ---- Table 2a/2b — top mediators (bootstrap p < 0.05 AND |proportion mediated| > 0.30) ----
PROP_CUT <- 0.30
mk <- function(d, idcol, nmin) {
  d <- d[is.finite(prop_mediated) & prop_mediated_pvalue < 0.05 & n_complete > nmin & abs(prop_mediated) > PROP_CUT]
  d <- d[order(-abs(prop_mediated))]
  data.table(
    Mediator      = toupper(d[[idcol]]),
    N             = d$n_complete,
    ACME          = signif(d$ACME_estimate, 3),
    Total_effect  = signif(d$total_effect, 3),
    Prop_mediated = round(d$prop_mediated, 3),
    P_value       = fifelse(d$prop_mediated_pvalue == 0, "<1e-3",
                            formatC(d$prop_mediated_pvalue, format = "g", digits = 2))
  )
}
table2a <- mk(prot,  "protein",    30000)   # Table 2a — proteins
table2b <- mk(metab, "metabolite", 350000)  # Table 2b — metabolites

write.csv(table2a, file.path(cfg$out_dir, "20260613_Table2a_top_protein_mediators.csv"),    row.names = FALSE)
write.csv(table2b, file.path(cfg$out_dir, "20260613_Table2b_top_metabolite_mediators.csv"), row.names = FALSE)
cat("Table 2a:", nrow(table2a), "top protein mediators | Table 2b:", nrow(table2b), "top metabolite mediators\n")


## =====================================================================
## Supplementary Tables S3 (proteins) & S4 (metabolites): full single-mediator results
## =====================================================================
# proteins carried into serial mediation (distinct proteins in the serial-mediation input)
serial_prot <- unique(fread(cfg$serial_mediation_results, select = "protein")$protein)

tidyfull <- function(d, idcol, nmin, carried = NULL) {
  d <- copy(d)
  if ("error" %in% names(d)) d[, error := NULL]
  setnames(d, idcol, "Mediator")
  d[, Mediator         := toupper(Mediator)]
  d[, significant_p05  := prop_mediated_pvalue < 0.05]
  d[, passed_n_filter  := n_complete > nmin]
  if (!is.null(carried)) d[, carried_to_serial := tolower(Mediator) %in% carried]
  setorder(d, Mediator)
  d[]
}
s3 <- tidyfull(prot,  "protein",    30000, carried = serial_prot)
s4 <- tidyfull(metab, "metabolite", 350000)

write.csv(s3, file.path(cfg$out_dir, "20260613_TableS3_protein_single_mediator_full.csv"), row.names = FALSE)
write.csv(s4, file.path(cfg$out_dir, "20260613_TableS4_metabolite_single_mediator_full.csv"), row.names = FALSE)
cat("Table S3 (proteins):  ", nrow(s3), "rows |", sum(s3$significant_p05), "significant (p<0.05) |",
    sum(s3$carried_to_serial), "carried to serial\n")
cat("Table S4 (metabolites):", nrow(s4), "rows |", sum(s4$significant_p05), "significant (p<0.05)\n")


## =====================================================================
## Figure 4 — core-protein landscape (core proteins)
## =====================================================================
core <- fread(cfg$core_proteins)
core[, Protein   := toupper(protein)]
core[, Mechanism := factor(mechanism_type,
                           levels = c("Amplifying-dominant", "Compensatory-dominant", "Mixed"))]
setorder(core, -max_abs_prop)   # label priority: largest-effect proteins first (check_overlap)

p_core <- ggplot(core, aes(n_pathways, max_abs_prop, color = Mechanism)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_text(aes(label = Protein), size = 2.6, vjust = -0.8, color = "black",
            check_overlap = TRUE, show.legend = FALSE) +
  scale_color_manual(values = c("Amplifying-dominant"   = "#B2182B",
                                "Compensatory-dominant" = "#2166AC",
                                "Mixed"                 = "grey55")) +
  labs(x = "Number of independent pathways",
       y = "Max |proportion mediated|", color = "Mechanism") +
  theme_classic(base_size = 12) +
  theme(panel.grid = element_blank(), legend.position = "right")

ggsave(file.path(cfg$out_dir, "20260613_Figure4_core_protein_landscape.png"),
       p_core, width = 8, height = 5.5, dpi = 300)
ggsave(file.path(cfg$out_dir, "20260613_Figure4_core_protein_landscape.pdf"),
       p_core, width = 8, height = 5.5, dpi = 300, device = cairo_pdf)
cat("Core proteins:", nrow(core), "| by mechanism:\n"); print(table(core$Mechanism))


## =====================================================================
## Table 3 — core mediating proteins
## =====================================================================
core <- fread(cfg$core_proteins)
core[, mech := fcase(
  mechanism_type == "Amplifying-dominant",   "Amplifying",
  mechanism_type == "Compensatory-dominant", "Compensatory",
  default = "Mixed")]
core[, mech := factor(mech, levels = c("Amplifying", "Compensatory", "Mixed"))]
table3 <- core[order(mech, -max_abs_prop), .(
  Protein              = toupper(protein),
  Dominant_mechanism   = as.character(mech),
  N_pathways           = n_pathways,
  N_metabolite_classes = n_metab_class,
  Max_prop_mediated    = round(max_abs_prop, 3),
  Mean_prop_mediated   = round(mean_abs_prop, 3),
  Pct_amplifying       = round(100 * prop_amplifying, 1)
)]
setnames(table3, "Dominant_mechanism", "Dominant Mechanism")
write.csv(table3, file.path(cfg$out_dir, "20260613_Table3_core_proteins.csv"), row.names = FALSE)
cat("Table 3:", nrow(table3), "core proteins\n")


## =====================================================================
## Supplementary Table S5 — all level-1 serial-mediation pathways
## =====================================================================
s5 <- fread(cfg$level1_pathways)
s5_out <- s5[order(protein, metabolite), .(
  Protein           = toupper(protein),
  Metabolite        = metabolite,
  Metabolite_class  = metabolite_class,
  Serial_indirect   = signif(serial_indirect_est, 3),
  CI_lower          = signif(serial_indirect_ci_lower, 3),
  CI_upper          = signif(serial_indirect_ci_upper, 3),
  Prop_mediated     = round(prop_serial_est, 3),
  FDR_serial        = signif(fdr_serial, 3),
  FDR_protein       = signif(fdr_prot, 3),
  FDR_metabolite    = signif(fdr_metab, 3),
  Direction         = direction
)]
write.csv(s5_out, file.path(cfg$out_dir, "20260613_TableS5_level1_pathways.csv"), row.names = FALSE)
cat("Table S5:", nrow(s5_out), "level-1 pathways |",
    sum(s5_out$Direction == "Amplifying"), "amplifying,",
    sum(s5_out$Direction == "Compensatory"), "compensatory\n")


## =====================================================================
## Figure 5 — core protein x HDL-subclass mediation heatmap
## =====================================================================
hdl <- fread(cfg$hdl_pathways)

# strongest pathway per protein x metabolite cell
hm <- hdl[, .SD[which.max(abs(prop_serial_est))], by = .(protein, metabolite)]
hm[, Protein := toupper(protein)]

# row order: proteins by max |proportion| (strongest at top); column order: Large then Very-large HDL
prot_order  <- hm[, .(m = max(abs(prop_serial_est))), by = Protein][order(m), Protein]
metab_order <- c("L_HDL_P","L_HDL_L","L_HDL_PL","L_HDL_C","L_HDL_CE","L_HDL_FC","L_HDL_TG",
                 "XL_HDL_P","XL_HDL_L","XL_HDL_PL","XL_HDL_C","XL_HDL_CE","XL_HDL_FC","XL_HDL_TG")
metab_order <- metab_order[metab_order %in% hm$metabolite]
hm[, Protein    := factor(Protein, levels = prot_order)]
hm[, metabolite := factor(metabolite, levels = metab_order)]

p_hdl <- ggplot(hm, aes(metabolite, Protein, fill = prop_serial_est)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.0f", prop_serial_est * 100)), size = 2) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       name = "Proportion\nmediated (%)", labels = function(x) x * 100) +
  labs(x = "HDL-subclass metabolite", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), panel.grid = element_blank())

ggsave(file.path(cfg$out_dir, "20260613_Figure5_HDL_heatmap.png"),
       p_hdl, width = 10, height = 8, dpi = 300)
ggsave(file.path(cfg$out_dir, "20260613_Figure5_HDL_heatmap.pdf"),
       p_hdl, width = 10, height = 8, dpi = 300, device = cairo_pdf)
cat("Figure 5: heatmap of", nrow(hm), "protein x HDL-metabolite cells (", uniqueN(hm$Protein), "proteins )\n")


## =====================================================================
## Supplementary Table S6 — core-protein x HDL-subclass serial-mediation pathways
## =====================================================================
s6 <- hdl[order(protein, metabolite), .(
  Protein         = toupper(protein),
  Metabolite      = metabolite,
  HDL_subclass    = metabolite_class,
  Serial_indirect = signif(serial_indirect_est, 3),
  CI_lower        = signif(serial_indirect_ci_lower, 3),
  CI_upper        = signif(serial_indirect_ci_upper, 3),
  Prop_mediated   = round(prop_serial_est, 3),
  FDR_serial      = signif(fdr_serial, 3),
  Direction       = direction
)]
write.csv(s6, file.path(cfg$out_dir, "20260613_TableS6_HDL_pathways.csv"), row.names = FALSE)
cat("Table S6:", nrow(s6), "HDL pathways |",
    sum(s6$Direction == "Amplifying"), "amplifying,",
    sum(s6$Direction == "Compensatory"), "compensatory\n")

cat("\nAll manuscript figures and tables written to:", cfg$out_dir, "\n")
