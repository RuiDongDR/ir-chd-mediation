## =====================================================================
## 14_pathway_prioritization.R
##
## Purpose : Prioritize serial-mediation pathways and identify the
##           broadly-acting "core" mediating proteins.
## Inputs  : cfg$serial_mediation_results  (converged serial-mediation
##             results from step 13; one row per protein x metabolite pair)
##           cfg$nmr_annotation            (Nightingale biomarker group table
##             with columns Biomarker, Group, Sub.Group)
## Outputs : cfg$level1_pathways   (prioritized "level-1" pathways)
##           cfg$core_proteins     (core proteins with breadth/effect metrics)
## Requires: R (readr, dplyr)
## Run after: 13_serial_mediation.R
## CONFIG KEYS USED:
##   serial_mediation_results - converged serial-mediation result table
##   nmr_annotation           - Nightingale biomarker -> class annotation
##   level1_pathways          - output: prioritized pathways
##   core_proteins            - output: core mediating proteins
## =====================================================================

suppressMessages({
  library(readr)
  library(dplyr)
})
source("config.R")

## ---- inputs ----------------------------------------------------------
raw_data  <- read_csv(cfg$serial_mediation_results, show_col_types = FALSE)
nmr_annot <- read_csv(cfg$nmr_annotation,           show_col_types = FALSE)

## metabolite -> class lookup (lipoprotein subclasses keep their sub-group label)
nmr_lookup <- nmr_annot %>%
  mutate(group_info = ifelse(Group == "Lipoprotein subclasses", Sub.Group, Group)) %>%
  select(Biomarker, group_info) %>%
  rename(metabolite = Biomarker, metabolite_class = group_info)

annotated_data <- raw_data %>% left_join(nmr_lookup, by = "metabolite")

## ---- level-1 prioritization filter -----------------------------------
## Keep FDR-significant serial pathways that explain > 5% of the total
## IR-PRS -> CHD effect, and whose protein-only and metabolite-only
## indirect effects are BOTH FDR-significant.
level1 <- annotated_data %>%
  filter(
    fdr_serial < 0.05,
    abs(prop_serial_est) > 0.05,
    fdr_prot  < 0.05,
    fdr_metab < 0.05
  ) %>%
  mutate(direction = ifelse(prop_serial_est > 0, "Amplifying", "Compensatory")) %>%
  arrange(desc(abs(prop_serial_est)))

write_csv(level1, cfg$level1_pathways)
cat(sprintf("Level-1 pathways: %d | proteins: %d | metabolites: %d\n",
            nrow(level1), n_distinct(level1$protein), n_distinct(level1$metabolite)))
cat(sprintf("  Amplifying: %d | Compensatory: %d\n",
            sum(level1$direction == "Amplifying"),
            sum(level1$direction == "Compensatory")))

## ---- protein-level breadth / effect-size summary ---------------------
protein_summary <- level1 %>%
  group_by(protein) %>%
  summarise(
    n_pathways        = n(),
    n_metab_class     = n_distinct(metabolite_class),
    max_abs_prop      = max(abs(prop_serial_est),  na.rm = TRUE),
    mean_abs_prop     = mean(abs(prop_serial_est), na.rm = TRUE),
    prop_amplifying   = mean(prop_serial_est > 0,  na.rm = TRUE),
    prop_compensatory = mean(prop_serial_est < 0,  na.rm = TRUE),
    .groups = "drop"
  )

## ---- select "core" proteins ------------------------------------------
## Retain proteins exceeding the 75th percentile on ALL FOUR breadth /
## effect-size metrics, then classify each protein's dominant mechanism.
metrics <- c("n_pathways", "n_metab_class", "max_abs_prop", "mean_abs_prop")
thresh  <- sapply(metrics, function(m) quantile(protein_summary[[m]], 0.75, na.rm = TRUE))

core_proteins <- protein_summary %>%
  filter(
    n_pathways    >= thresh["n_pathways"],
    n_metab_class >= thresh["n_metab_class"],
    max_abs_prop  >= thresh["max_abs_prop"],
    mean_abs_prop >= thresh["mean_abs_prop"]
  ) %>%
  mutate(
    mechanism_type = case_when(
      prop_amplifying   > 0.7 & max_abs_prop > 0.1 ~ "Amplifying-dominant",
      prop_compensatory > 0.7 & max_abs_prop > 0.1 ~ "Compensatory-dominant",
      TRUE                                         ~ "Mixed"
    )
  ) %>%
  arrange(desc(n_pathways))

write_csv(core_proteins, cfg$core_proteins)
cat(sprintf("\nCore proteins: %d\n", nrow(core_proteins)))
print(table(core_proteins$mechanism_type))
