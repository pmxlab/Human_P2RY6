# ============================================================================#
# S87_Fibroblast_Abundance.R
# Fibroblast abundance per participant by condition, to check whether the
# higher fibroblast PDGFRB expression can be read as proliferation.
# Reviewer comment R2.16 (Life, life-4655118)
#
# Per participant (KPMP annotation, same participants and conditions as the
# compartment analyses of S50/S60):
#   scRNA-seq: fibroblasts (subclass.level1 == "FIB") and the stromal
#              compartment (class == "stroma cells") as % of all cells
#   snRNA-seq: fibroblasts (subclass.l1 == "FIB") as % of all nuclei
# All CKD and DKD-focused (CKD restricted to adjudicated DKD, as in S52/S62).
# Each disease group vs Healthy: two-sided Wilcoxon rank-sum test,
# BH-adjusted over the two comparisons of each measure and set (as in S80).
#
# Only the cell metadata of the h5Seurat files are read (hdf5r).
#
# Outputs: 04-Outputs/S87_Fibroblast_Abundance/
#   FIB_abundance_per_participant.csv
#   Table_FIB_abundance.csv          median (IQR) by condition + BH p
# ============================================================================#

rm(list = ls())

library(hdf5r)
library(dplyr)
library(tidyr)

# ============================================================================#
# SETTINGS ----
# ============================================================================#

sc_file <- "../../KPMP_scRNA/KPMP_PREMIERE_SC_version2_ForExplorer_RemovedBatchEffect_Final2025.h5seurat"
sn_file <- "../../KPMP_snRNA/SingleNucleus_KPMP_Explorer_05182025.h5Seurat"
sc_clin <- "../../KPMP_scRNA/OpenAccessClinicalData/20251205_OpenAccessClinicalData.csv"
sn_clin <- "../../KPMP_snRNA/6d3047c2-588f-4ad4-8eba-f4c46a1e5c14_20251205_OpenAccessClinicalData/20251205_OpenAccessClinicalData.csv"
out_dir <- "../04-Outputs/S87_Fibroblast_Abundance"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

groups <- c("Healthy", "AKI", "CKD")
sc_lab <- c("Healthy Reference Tissue" = "Healthy", "Acute Kidney Injury" = "AKI",
            "Chronic Kidney Disease" = "CKD")
sn_lab <- c("Healthy Reference" = "Healthy", "AKI" = "AKI", "CKD" = "CKD")

# ============================================================================#
# HELPER FUNCTIONS ----
# ============================================================================#

# metadata columns of an h5Seurat file (factors are stored as levels/values)
read_meta <- function(file, cols) {
  h <- H5File$new(file, mode = "r")
  on.exit(h$close_all())
  m <- h[["meta.data"]]
  as_tibble(lapply(setNames(cols, cols), function(col) {
    x <- m[[col]]
    if (inherits(x, "H5Group")) x[["levels"]][][x[["values"]][]] else x[]
  }))
}

med_iqr <- function(v) {
  q <- quantile(v, c(0.5, 0.25, 0.75), names = FALSE)
  sprintf("%.2f (%.2f–%.2f)", q[1], q[2], q[3])
}

# each disease group vs Healthy, BH over the comparisons of one measure
vs_healthy_bh <- function(v, g) {
  p <- vapply(c("AKI", "CKD"), function(x) wilcox.test(v[g == x], v[g == "Healthy"])$p.value,
              numeric(1))
  p.adjust(p, method = "BH")
}

# ============================================================================#
# 1. PER-PARTICIPANT ABUNDANCE ----
# ============================================================================#

clin_sc <- read.csv(sc_clin)
sc <- read_meta(sc_file, c("SampleID", "EnrollmentCategory", "class", "subclass.level1")) %>%
  transmute(id = as.character(SampleID), group = unname(sc_lab[EnrollmentCategory]),
            fib = subclass.level1 == "FIB", stroma = class == "stroma cells") %>%
  filter(!is.na(group)) %>%
  group_by(id, group) %>%
  summarise(n_all = n(), pct_fib = 100 * mean(fib), pct_stroma = 100 * mean(stroma), .groups = "drop") %>%
  mutate(dataset = "scRNA-seq",
         adjudicated = clin_sc$Primary.Adjudicated.Category[match(id, clin_sc$Participant.ID)])

clin_sn <- read.csv(sn_clin)
sn <- read_meta(sn_file, c("patient", "subclass.l1")) %>%
  transmute(id = patient,
            group = unname(sn_lab[clin_sn$Enrollment.Category[match(patient, clin_sn$Participant.ID)]]),
            fib = subclass.l1 == "FIB") %>%
  filter(!is.na(group)) %>%
  group_by(id, group) %>%
  summarise(n_all = n(), pct_fib = 100 * mean(fib), .groups = "drop") %>%
  mutate(dataset = "snRNA-seq",
         adjudicated = clin_sn$Primary.Adjudicated.Category[match(id, clin_sn$Participant.ID)])

fib_all <- bind_rows(sc, sn) %>% mutate(group = factor(group, levels = groups))
write.csv(fib_all, file.path(out_dir, "FIB_abundance_per_participant.csv"), row.names = FALSE)

# ============================================================================#
# 2. SUMMARY TABLE ----
# ============================================================================#

fib_sets <- bind_rows(
  fib_all %>% mutate(set = "all CKD"),
  fib_all %>% filter(!(group == "CKD" & adjudicated != "Diabetic Kidney Disease")) %>%
    mutate(set = "DKD-focused"))

tab <- fib_sets %>%
  pivot_longer(c(pct_fib, pct_stroma), names_to = "measure", values_to = "value") %>%
  filter(!is.na(value)) %>%
  group_by(dataset, set, measure) %>%
  group_modify(function(d, k) {
    p <- vs_healthy_bh(d$value, d$group)
    tibble(n = paste(as.integer(table(d$group)), collapse = " / "),
           Healthy = med_iqr(d$value[d$group == "Healthy"]),
           AKI     = med_iqr(d$value[d$group == "AKI"]),
           CKD     = med_iqr(d$value[d$group == "CKD"]),
           p_AKI_vs_Healthy = p[["AKI"]], p_CKD_vs_Healthy = p[["CKD"]])
  }) %>%
  ungroup()
print(as.data.frame(tab), digits = 3)
write.csv(tab, file.path(out_dir, "Table_FIB_abundance.csv"), row.names = FALSE)
