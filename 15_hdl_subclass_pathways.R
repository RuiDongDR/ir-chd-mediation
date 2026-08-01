## =====================================================================
## 15_hdl_subclass_pathways.R
##
## Purpose : Localize the core-protein mediation signal to the Large and
##           Very-large HDL particle subclasses.
## Inputs  : cfg$serial_mediation_results  (converged serial-mediation results, step 13)
##           cfg$nmr_annotation            (Nightingale biomarker group table)
##           cfg$core_proteins             (core proteins from step 14)
## Outputs : cfg$hdl_pathways   (core-protein x HDL-subclass pathways, ranked)
## Requires: R (readr, dplyr)
## Run after: 14_pathway_prioritization.R
## CONFIG KEYS USED:
##   serial_mediation_results - converged serial-mediation result table
##   nmr_annotation           - Nightingale biomarker -> class annotation
##   core_proteins            - core mediating proteins (from step 14)
##   hdl_pathways             - output: core-protein x HDL-subclass pathways
## =====================================================================

suppressMessages({
  library(readr)
  library(dplyr)
})
source("config.R")

## ---- inputs ----------------------------------------------------------
raw_data      <- read_csv(cfg$serial_mediation_results, show_col_types = FALSE)
nmr_annot     <- read_csv(cfg$nmr_annotation,           show_col_types = FALSE)
core_proteins <- read_csv(cfg$core_proteins,            show_col_types = FALSE)$protein

nmr_lookup <- nmr_annot %>%
  mutate(group_info = ifelse(Group == "Lipoprotein subclasses", Sub.Group, Group)) %>%
  select(Biomarker, group_info) %>%
  rename(metabolite = Biomarker, metabolite_class = group_info)

## The two largest HDL subclasses (14 lipoprotein measures in total).
hdl_groups <- c(
  "Large HDL (average diameter 12.1 nm)",
  "Very large HDL (average diameter 14.3 nm)"
)
hdl_metabs <- nmr_lookup %>% filter(metabolite_class %in% hdl_groups)

cat(sprintf("Core proteins: %d | HDL groups: %d | HDL metabolites: %d\n",
            length(core_proteins), length(hdl_groups), nrow(hdl_metabs)))

annotated <- raw_data %>% left_join(nmr_lookup, by = "metabolite")

## ---- level-1 filter (identical to step 14) ---------------------------
level1 <- annotated %>%
  filter(
    fdr_serial < 0.05,
    abs(prop_serial_est) > 0.05,
    fdr_prot  < 0.05,
    fdr_metab < 0.05
  ) %>%
  mutate(direction = ifelse(prop_serial_est > 0, "Amplifying", "Compensatory"))

## ---- restrict to core proteins x the two HDL subclasses --------------
hdl_pathways <- level1 %>%
  filter(protein %in% core_proteins, metabolite_class %in% hdl_groups) %>%
  arrange(desc(abs(prop_serial_est)))

write_csv(hdl_pathways, cfg$hdl_pathways)

cat(sprintf("\nCore-protein x HDL pathways: %d\n", nrow(hdl_pathways)))
cat(sprintf("  Unique proteins: %d of %d | Unique metabolites: %d of %d\n",
            n_distinct(hdl_pathways$protein), length(core_proteins),
            n_distinct(hdl_pathways$metabolite), nrow(hdl_metabs)))
cat("\nDirection breakdown:\n");        print(table(hdl_pathways$direction))
cat("\nPathways per HDL subclass:\n");   print(table(hdl_pathways$metabolite_class))

## ---- lipoprotein-class representation across all level-1 pathways -----
cat("\nMetabolite-class representation (all level-1 pathways):\n")
print(sort(table(level1$metabolite_class), decreasing = TRUE))
