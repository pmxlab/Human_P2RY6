# ============================================================================#
# S90_P2RY6_All_CellTypes.R
# P2RY6 expression across all annotated cell types of the KPMP snRNA-seq atlas
# by condition, to judge which cell types carry the P2RY6 increase.
# Reviewer comment R2.15 (Life, life-4655118)
#
# HEAVY: loads the snRNA-seq h5Seurat object once (raw counts only, as in
# S83). Run on the Windows host, not in the background on WSL.
#
# Cell types: KPMP subclass.l1 (all 18 subclasses). Conditions: Healthy, AKI,
# CKD (all causes) as in S60; CKD (DKD) = CKD with adjudicated DKD.
# P2RY6 is log-normalized per nucleus as in S60/S83 (counts / total counts of
# the nucleus x 10,000, natural log of 1 + x).
#
# 1. Per cell type and condition (all nuclei pooled):
#    nuclei, % of all nuclei, % P2RY6-positive nuclei, mean log-normalized
#    P2RY6, % of all P2RY6 UMIs of the condition
# 2. Per participant and cell type:
#    pseudobulk mean P2RY6 (participants with >= 20 nuclei of the cell type)
#    tissue contribution = (nuclei of the cell type / all nuclei) x mean
#      P2RY6 of the cell type; summed over cell types it gives the mean P2RY6
#      over all nuclei of the participant (a pseudo-bulk P2RY6 level), so the
#      contributions show which cell types carry a change of total P2RY6
# 3. Group comparisons per cell type: AKI and CKD vs Healthy, two-sided
#    Wilcoxon rank-sum test, BH-adjusted over the two comparisons (as in S80)
# Check: the fibroblast pseudobulk reproduces Figure 2d (S60 pb_fib).
#
# Inputs : KPMP snRNA-seq h5Seurat object, clinical CSV,
#          04-Outputs/S60_KPMP_snRNA/S60_pseudobulk.rds
# Outputs: 04-Outputs/S90_P2RY6_All_CellTypes/
#          S90_per_participant.rds
#          Table_P2RY6_celltypes_pooled.csv       (1)
#          Table_P2RY6_celltypes_pseudobulk.csv   (2, 3; medians and BH p)
#          Fig_P2RY6_celltypes.png/.pdf           dot plot cell type x condition
# ============================================================================#

rm(list = ls())

library(Seurat)
library(SeuratDisk)
library(Matrix)
library(dplyr)
library(tidyr)
library(ggplot2)

# ============================================================================#
# SETTINGS ----
# ============================================================================#

sn_file <- "../../KPMP_snRNA/SingleNucleus_KPMP_Explorer_05182025.h5Seurat"
sn_clin <- "../../KPMP_snRNA/6d3047c2-588f-4ad4-8eba-f4c46a1e5c14_20251205_OpenAccessClinicalData/20251205_OpenAccessClinicalData.csv"
s60_rds <- "../04-Outputs/S60_KPMP_snRNA/S60_pseudobulk.rds"
out_dir <- "../04-Outputs/S90_P2RY6_All_CellTypes"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

min_cells <- 20
groups    <- c("Healthy", "AKI", "CKD")
sn_lab    <- c("Healthy Reference" = "Healthy", "AKI" = "AKI", "CKD" = "CKD")

# ============================================================================#
# HELPER FUNCTIONS ----
# ============================================================================#

med_iqr <- function(v, digits = 3) {
  v <- v[!is.na(v)]
  if (!length(v)) return("–")
  q <- quantile(v, c(0.5, 0.25, 0.75), names = FALSE)
  sprintf(paste0("%.", digits, "f (%.", digits, "f–%.", digits, "f)"), q[1], q[2], q[3])
}

# AKI and CKD vs Healthy, BH over the two comparisons
vs_healthy_bh <- function(v, g) {
  p <- vapply(c("AKI", "CKD"), function(x) {
    a <- v[g == x]; b <- v[g == "Healthy"]
    if (sum(!is.na(a)) < 3 || sum(!is.na(b)) < 3) return(NA_real_)
    suppressWarnings(wilcox.test(a, b)$p.value)
  }, numeric(1))
  p.adjust(p, method = "BH")
}

# ============================================================================#
# 1. LOAD snRNA-seq (RAW COUNTS ONLY) AND LOG-NORMALIZE P2RY6 ----
# ============================================================================#

cat("Loading snRNA object...\n")
obj <- LoadH5Seurat(sn_file, assays = list(RNA = "counts"), reductions = FALSE,
                    graphs = FALSE, neighbors = FALSE, images = FALSE, misc = FALSE, tools = FALSE)
clin <- read.csv(sn_clin)
counts <- LayerData(obj, assay = "RNA", layer = "counts")
stopifnot("P2RY6" %in% rownames(counts))
umi <- as.numeric(counts["P2RY6", ])
tot <- Matrix::colSums(counts)
rm(counts); gc()

nuc <- tibble(
  id        = obj$patient,
  group     = unname(sn_lab[clin$Enrollment.Category[match(obj$patient, clin$Participant.ID)]]),
  adjud     = clin$Primary.Adjudicated.Category[match(obj$patient, clin$Participant.ID)],
  cell_type = as.character(obj$subclass.l1),
  umi       = umi,
  p2ry6     = log1p(umi / tot * 1e4)) %>%
  filter(!is.na(group))
rm(obj); gc()

# ============================================================================#
# 2. POOLED TABLE PER CELL TYPE AND CONDITION ----
# ============================================================================#

pooled <- nuc %>%
  group_by(group, cell_type) %>%
  summarise(nuclei = n(), pct_pos = 100 * mean(umi > 0), mean_p2ry6 = mean(p2ry6),
            umi = sum(umi), .groups = "drop") %>%
  group_by(group) %>%
  mutate(pct_of_nuclei = 100 * nuclei / sum(nuclei), pct_of_P2RY6_umi = 100 * umi / sum(umi)) %>%
  ungroup() %>%
  mutate(group = factor(group, levels = groups)) %>%
  arrange(cell_type, group)
print(as.data.frame(pooled), digits = 3)
write.csv(pooled, file.path(out_dir, "Table_P2RY6_celltypes_pooled.csv"), row.names = FALSE)

# ============================================================================#
# 3. PER PARTICIPANT: PSEUDOBULK AND TISSUE CONTRIBUTION ----
# ============================================================================#

n_all <- nuc %>% count(id, name = "n_all")
pp <- nuc %>%
  group_by(id, group, adjud, cell_type) %>%
  summarise(n = n(), mean_p2ry6 = mean(p2ry6), .groups = "drop") %>%
  complete(nesting(id, group, adjud), cell_type, fill = list(n = 0L)) %>%
  left_join(n_all, by = "id") %>%
  mutate(contribution = ifelse(n > 0, n / n_all * mean_p2ry6, 0),
         pseudobulk   = ifelse(n >= min_cells, mean_p2ry6, NA_real_),
         group        = factor(group, levels = groups))

# check: fibroblast pseudobulk = Figure 2d (S60)
s60 <- readRDS(s60_rds)
fib <- pp %>% filter(cell_type == "FIB", !is.na(pseudobulk))
m <- match(s60$pb_fib$patient, fib$id)
stopifnot(!anyNA(m), max(abs(fib$pseudobulk[m] - s60$pb_fib$mean_P2RY6)) < 1e-6)

# whole-tissue (all nuclei) P2RY6 per participant = sum of contributions
tissue <- pp %>% group_by(id, group, adjud) %>%
  summarise(tissue_p2ry6 = sum(contribution), .groups = "drop")
saveRDS(list(per_participant = pp, tissue = tissue), file.path(out_dir, "S90_per_participant.rds"))

# ============================================================================#
# 4. GROUP COMPARISONS PER CELL TYPE ----
# ============================================================================#

summ <- function(d, value, measure, set) {
  d %>% group_by(cell_type) %>%
    group_modify(function(x, k) {
      p <- vs_healthy_bh(x[[value]], x$group)
      tibble(set = set, measure = measure,
             n = paste(vapply(groups, function(g) sum(!is.na(x[[value]][x$group == g])), 1L), collapse = " / "),
             Healthy = med_iqr(x[[value]][x$group == "Healthy"]),
             AKI     = med_iqr(x[[value]][x$group == "AKI"]),
             CKD     = med_iqr(x[[value]][x$group == "CKD"]),
             p_AKI_vs_Healthy = p[["AKI"]], p_CKD_vs_Healthy = p[["CKD"]])
    }) %>% ungroup()
}
dkd <- pp %>% filter(!(group == "CKD" & !adjud %in% "Diabetic Kidney Disease"))
tab <- bind_rows(
  summ(pp,  "pseudobulk",   "Mean P2RY6 per cell type (pseudobulk)",       "all CKD"),
  summ(pp,  "contribution", "Contribution to whole-tissue P2RY6",           "all CKD"),
  summ(dkd, "pseudobulk",   "Mean P2RY6 per cell type (pseudobulk)",       "DKD-focused"),
  summ(dkd, "contribution", "Contribution to whole-tissue P2RY6",           "DKD-focused"))
tissue_tab <- bind_rows(
  summ(tissue %>% mutate(cell_type = "All nuclei"), "tissue_p2ry6", "Whole-tissue P2RY6 (all nuclei)", "all CKD"),
  summ(tissue %>% filter(!(group == "CKD" & !adjud %in% "Diabetic Kidney Disease")) %>%
         mutate(cell_type = "All nuclei"), "tissue_p2ry6", "Whole-tissue P2RY6 (all nuclei)", "DKD-focused"))
tab <- bind_rows(tissue_tab, tab)
print(as.data.frame(tab), digits = 3)
write.csv(tab, file.path(out_dir, "Table_P2RY6_celltypes_pseudobulk.csv"), row.names = FALSE)

# ============================================================================#
# 5. DOT PLOT ----
# ============================================================================#

ct_order <- pooled %>% group_by(cell_type) %>% summarise(m = max(mean_p2ry6)) %>% arrange(m) %>% pull(cell_type)
p <- ggplot(pooled, aes(x = group, y = factor(cell_type, levels = ct_order))) +
  geom_point(aes(size = pct_pos, colour = mean_p2ry6)) +
  scale_colour_gradient(low = "grey85", high = "#08306b", name = "Mean P2RY6\n(log-normalized)") +
  scale_size_area(max_size = 8, name = "% P2RY6+\nnuclei") +
  labs(x = NULL, y = NULL, title = "KPMP snRNA-seq: P2RY6 by cell type and condition") +
  theme_bw(base_size = 12)
ggsave(file.path(out_dir, "Fig_P2RY6_celltypes.png"), p, width = 6.5, height = 7, dpi = 300)
ggsave(file.path(out_dir, "Fig_P2RY6_celltypes.pdf"), p, width = 6.5, height = 7, device = cairo_pdf)
cat("S90 done.\n")
