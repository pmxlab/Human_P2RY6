# ============================================================================#
# S83_P2RY6_Stromal_CellTypes.R
# Consistency check of P2RY6 between KPMP scRNA-seq and snRNA-seq for the same
# stromal cell types (reviewer comment R1.04).
#
# The scRNA-seq "stroma cells" class (Figure 2c) is ~82% VSMC/P and ~18%
# fibroblasts (subtypes MYOF and pvFIB only); the snRNA-seq fibroblast subclass
# (Figure 2e) contains all fibroblast subtypes. Cell sets (figure rows:
# scRNA-seq, snRNA-seq; columns: stroma, VSMC/P, fibroblasts):
#   (1) scRNA-seq stroma as in Figure 2c (checked against S50)
#   (2) scRNA-seq VSMC/P only (subclass.level1 == "VSMC/P")
#   (3) scRNA-seq fibroblasts only (subclass.level1 == "FIB")
#   (4) snRNA-seq stroma with scRNA-seq composition: per participant,
#       weighted mean of the FIB and VSM/P pseudobulk values with the
#       scRNA-seq stroma weights (FIB ~18%, VSMC/P ~82%)
#   (5) snRNA-seq VSM/P only  (subclass.l1 == "VSM/P")
#   (6) snRNA-seq fibroblasts as in Figure 2e (checked against S60)
# Per participant (>= 20 cells/nuclei): mean log-normalized P2RY6 and % of
# P2RY6-positive cells/nuclei; the table also gives the number of cells/nuclei
# per participant before the filter. All conditions; Wilcoxon rank-sum test vs
# Healthy, BH-adjusted within each cell set (as in the manuscript figures).
#
# Inputs : KPMP scRNA-seq and snRNA-seq h5Seurat objects, clinical CSVs,
#          04-Outputs/S50_KPMP_scRNA/S50_pseudobulk.rds
#          04-Outputs/S60_KPMP_snRNA/S60_pseudobulk.rds
# Output : 04-Outputs/S83_P2RY6_Stromal_CellTypes/
#          S83_pseudobulk.rds, Table_P2RY6_cell_types.csv/.md,
#          Table_P2RY6_by_class.csv (share of cells and P2RY6 UMIs per class),
#          Table_S2_cell_composition.csv,
#          Fig_P2RY6_cell_types.png/.pdf
# ============================================================================#

rm(list = ls())

library(Seurat)
library(SeuratDisk)
library(dplyr)
library(tidyr)
library(ggplot2)

# ============================================================================#
# SETTINGS ----
# ============================================================================#

sc_file <- "../../KPMP_scRNA/KPMP_PREMIERE_SC_version2_ForExplorer_RemovedBatchEffect_Final2025.h5seurat"
sn_file <- "../../KPMP_snRNA/SingleNucleus_KPMP_Explorer_05182025.h5Seurat"
sn_clin <- "../../KPMP_snRNA/6d3047c2-588f-4ad4-8eba-f4c46a1e5c14_20251205_OpenAccessClinicalData/20251205_OpenAccessClinicalData.csv"
s50_rds <- "../04-Outputs/S50_KPMP_scRNA/S50_pseudobulk.rds"
s60_rds <- "../04-Outputs/S60_KPMP_snRNA/S60_pseudobulk.rds"
out_dir <- "../04-Outputs/S83_P2RY6_Stromal_CellTypes"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

min_cells <- 20
groups    <- c("Healthy", "AKI", "CKD")
sc_lab    <- c("Healthy Reference Tissue" = "Healthy", "Acute Kidney Injury" = "AKI",
               "Chronic Kidney Disease" = "CKD")
sn_lab    <- c("Healthy Reference" = "Healthy", "AKI" = "AKI", "CKD" = "CKD")
group_colors <- c(Healthy = "#2166ac", AKI = "#e07b39", CKD = "#756bb1")

set_levels <- c(
  "scRNA-seq stroma (Figure 2c)",
  "scRNA-seq VSMC/P",
  "scRNA-seq fibroblasts",
  "snRNA-seq stroma, scRNA-seq composition",
  "snRNA-seq VSM/P",
  "snRNA-seq fibroblasts (Figure 2e)"
)

# Only raw counts and cell metadata are loaded (lower memory); the subsets are
# log-normalized here, which gives the same values as in S50/S60.

# ============================================================================#
# HELPER FUNCTIONS ----
# ============================================================================#

# per-participant mean log-normalized P2RY6, % positive cells and cell count
# (all participants with >= 1 cell; the min_cells filter is applied later)
p2ry6_pseudobulk <- function(obj, cells, id_col, group_col) {
  sub <- subset(obj, cells = cells)
  sub <- NormalizeData(sub, normalization.method = "LogNormalize",
                       scale.factor = 10000, verbose = FALSE)
  df  <- FetchData(sub, vars = c("P2RY6", id_col, group_col))
  names(df) <- c("P2RY6", "id", "group")
  df %>%
    group_by(id, group) %>%
    summarise(mean_P2RY6    = mean(P2RY6),
              pct_P2RY6_pos = mean(P2RY6 > 0) * 100,
              n_cells       = n(),
              .groups = "drop")
}

# share of all cells/nuclei and of all P2RY6 UMIs per cell class (whole tissue)
p2ry6_by_class <- function(obj, keep, cls, dataset) {
  umi <- FetchData(obj, vars = "P2RY6", layer = "counts")[keep, 1]   # one gene, no matrix copy
  tibble(dataset = dataset, class = cls[keep], umi = as.numeric(umi)) %>%
    group_by(dataset, class) %>%
    summarise(cells = n(), pos = sum(umi > 0), umi = sum(umi), .groups = "drop") %>%
    mutate(pct_of_cells = 100 * cells / sum(cells), pct_of_P2RY6_umi = 100 * umi / sum(umi),
           pct_P2RY6_pos = 100 * pos / cells) %>%
    arrange(desc(pct_of_cells))
}

med_iqr <- function(v, digits) {
  q <- quantile(v, c(0.5, 0.25, 0.75), names = FALSE)
  sprintf(paste0("%.", digits, "f (%.", digits, "f–%.", digits, "f)"), q[1], q[2], q[3])
}

# ============================================================================#
# 1. scRNA-seq: STROMA, VSMC/P AND FIBROBLASTS ----
# ============================================================================#

cat("Loading scRNA object...\n")
obj <- LoadH5Seurat(sc_file, assays = list(RNA = "counts"), reductions = FALSE,
                   graphs = FALSE, neighbors = FALSE, images = FALSE, misc = FALSE, tools = FALSE)
obj$group <- unname(sc_lab[as.character(obj$EnrollmentCategory)])
is_stroma <- !is.na(obj$group) & obj$class == "stroma cells"
stroma_comp <- table(as.character(obj$subclass.level1[is_stroma]))
w_fib <- unname(stroma_comp["FIB"] / sum(stroma_comp))
cat(sprintf("scRNA stroma composition: FIB %.1f%%, VSMC/P %.1f%%\n",
            100 * w_fib, 100 * (1 - w_fib)))
stopifnot(setequal(names(stroma_comp), c("FIB", "VSMC/P")))

all_sc <- list(
  stroma = p2ry6_pseudobulk(obj, colnames(obj)[is_stroma], "SampleID", "group"),
  vsmc   = p2ry6_pseudobulk(obj, colnames(obj)[is_stroma & obj$subclass.level1 == "VSMC/P"],
                            "SampleID", "group"),
  fib    = p2ry6_pseudobulk(obj, colnames(obj)[is_stroma & obj$subclass.level1 == "FIB"],
                            "SampleID", "group"))
sc_cls <- ifelse(obj$cell_type == "epithelial cell of proximal tubule", "PT epithelial cells",
                 ifelse(obj$class == "epithelial cells", "other epithelial cells", as.character(obj$class)))
by_class_sc <- p2ry6_by_class(obj, !is.na(obj$group), sc_cls, "scRNA-seq")
rm(obj); gc()

# ============================================================================#
# 2. snRNA-seq: VSM/P AND FIBROBLASTS ----
# ============================================================================#

cat("Loading snRNA object...\n")
obj <- LoadH5Seurat(sn_file, assays = list(RNA = "counts"), reductions = FALSE,
                   graphs = FALSE, neighbors = FALSE, images = FALSE, misc = FALSE, tools = FALSE)
clin <- read.csv(sn_clin)
obj$group <- unname(sn_lab[clin$Enrollment.Category[match(obj$patient, clin$Participant.ID)]])

all_sn <- list(
  vsmp = p2ry6_pseudobulk(obj, colnames(obj)[!is.na(obj$group) & obj$subclass.l1 == "VSM/P"],
                          "patient", "group"),
  fib  = p2ry6_pseudobulk(obj, colnames(obj)[!is.na(obj$group) & obj$subclass.l1 == "FIB"],
                          "patient", "group"))
sn_cls <- ifelse(obj$subclass.l1 == "PT", "PT epithelial cells",
                 ifelse(obj$class == "epithelial cells", "other epithelial cells", as.character(obj$class)))
by_class_sn <- p2ry6_by_class(obj, !is.na(obj$group), sn_cls, "snRNA-seq")
rm(obj); gc()

# ============================================================================#
# 3. REFERENCE CHECK AND COMPOSITION-MATCHED snRNA-seq STROMA ----
# ============================================================================#

# the stroma (scRNA-seq) and fibroblast (snRNA-seq) sets must reproduce Figure 2c/2e
s50 <- readRDS(s50_rds)
s60 <- readRDS(s60_rds)
chk <- function(pb, ref, id_col) {
  m <- inner_join(filter(pb, n_cells >= min_cells), ref, by = c(id = id_col))
  stopifnot(nrow(m) == nrow(ref), isTRUE(all.equal(m$mean_P2RY6.x, m$mean_P2RY6.y, tolerance = 1e-6)))
}
chk(all_sc$stroma, s50$pb_fib %>% select(SampleID, mean_P2RY6), "SampleID")
chk(all_sn$fib,    s60$pb_fib %>% select(patient,  mean_P2RY6), "patient")

# snRNA-seq fibroblasts and VSM/P combined per participant with the scRNA-seq
# stroma weights (FIB ~18%, VSMC/P ~82%); participants need both components
all_sn$matched <- inner_join(all_sn$fib, all_sn$vsmp, by = c("id", "group"), suffix = c("_fib", "_vsmp")) %>%
  transmute(id, group,
            mean_P2RY6    = w_fib * mean_P2RY6_fib    + (1 - w_fib) * mean_P2RY6_vsmp,
            pct_P2RY6_pos = w_fib * pct_P2RY6_pos_fib + (1 - w_fib) * pct_P2RY6_pos_vsmp,
            n_cells       = n_cells_fib + n_cells_vsmp,
            n_min         = pmin(n_cells_fib, n_cells_vsmp))

# all participants with >= 1 cell/nucleus of the set (used for the coverage columns)
pb_all <- bind_rows(
  mutate(all_sc$stroma,  set = set_levels[1]),
  mutate(all_sc$vsmc,    set = set_levels[2]),
  mutate(all_sc$fib,     set = set_levels[3]),
  mutate(all_sn$matched, set = set_levels[4]),
  mutate(all_sn$vsmp,    set = set_levels[5]),
  mutate(all_sn$fib,     set = set_levels[6])
) %>%
  mutate(set = factor(set, levels = set_levels), group = factor(group, levels = groups))
stopifnot(!anyNA(pb_all$group))
pb <- filter(pb_all, n_cells >= min_cells, is.na(n_min) | n_min >= min_cells)
by_class <- bind_rows(by_class_sc, by_class_sn)
print(as.data.frame(by_class), digits = 3)
write.csv(by_class, file.path(out_dir, "Table_P2RY6_by_class.csv"), row.names = FALSE)
saveRDS(list(pb = pb, pb_all = pb_all, w_fib = w_fib, by_class = by_class),
        file.path(out_dir, "S83_pseudobulk.rds"))

# ============================================================================#
# 4. TESTS AND TABLE ----
# ============================================================================#

tests <- pb %>%
  group_by(set) %>%
  group_modify(function(d, k) {
    # no test if a group has fewer than 3 participants
    p <- sapply(c("AKI", "CKD"), function(g) {
      x <- d$mean_P2RY6[d$group == g]; y <- d$mean_P2RY6[d$group == "Healthy"]
      if (length(x) < 3 || length(y) < 3) NA_real_ else wilcox.test(x, y)$p.value
    })
    tibble(group = c("AKI", "CKD"), p = p, p_BH = p.adjust(p, "BH"))
  }) %>%
  ungroup() %>%
  mutate(stars = cut(p_BH, c(-Inf, 1e-4, 1e-3, 1e-2, 0.05, Inf),
                     labels = c("****", "***", "**", "*", "ns")))

coverage <- pb_all %>%
  group_by(set, group) %>%
  summarise(participants_any = n(),
            cells_all_med_iqr = med_iqr(n_cells, 0),
            .groups = "drop")

tab <- pb %>%
  group_by(set, group) %>%
  summarise(participants     = n(),
            cells_med_iqr    = med_iqr(n_cells, 0),
            mean_P2RY6       = med_iqr(mean_P2RY6, 4),
            pct_pos_med_iqr  = med_iqr(pct_P2RY6_pos, 2),
            .groups = "drop") %>%
  right_join(coverage, by = c("set", "group")) %>%
  arrange(set, group) %>%
  left_join(tests %>% mutate(group = factor(group, levels = groups)), by = c("set", "group")) %>%
  mutate(participants = coalesce(participants, 0L),
         p_BH_vs_Healthy = case_when(group == "Healthy" ~ "",
                                     is.na(p_BH) ~ "not tested",
                                     TRUE ~ sprintf("%.3g (%s)", p_BH, stars))) %>%
  select(set, group, participants_any, cells_all_med_iqr, participants, cells_med_iqr,
         mean_P2RY6, pct_pos_med_iqr, p_BH_vs_Healthy)

print(as.data.frame(tab))
write.csv(tab, file.path(out_dir, "Table_P2RY6_cell_types.csv"), row.names = FALSE)
hdr <- c("Cell set", "Group", "Participants with >= 1 cell/nucleus",
         "Cells/nuclei per participant (all), median (IQR)",
         "Participants with >= 20 cells/nuclei (analyzed)",
         "Cells/nuclei per analyzed participant, median (IQR)",
         "Mean P2RY6 (log-normalized), median (IQR)",
         "P2RY6+ cells/nuclei per participant, %, median (IQR)", "BH-adjusted p vs Healthy")
writeLines(c(paste("|", paste(hdr, collapse = " | "), "|"),
             paste("|", paste(rep("---", length(hdr)), collapse = " | "), "|"),
             apply(tab, 1, function(r) paste("|", paste(r, collapse = " | "), "|"))),
           file.path(out_dir, "Table_P2RY6_cell_types.md"))

# ============================================================================#
# 5. FIGURE ----
# ============================================================================#

ymax <- pb %>% group_by(set) %>% summarise(top = max(mean_P2RY6), rng = diff(range(mean_P2RY6)))
lab  <- tests %>% left_join(ymax, by = "set") %>%
  mutate(y = top + rng * ifelse(group == "AKI", 0.10, 0.22),
         x = ifelse(group == "AKI", 1.5, 2.0), label = paste0(group, ": ", ifelse(is.na(stars), "not tested", as.character(stars))))
n_lab <- pb %>% count(set, group) %>%
  group_by(set) %>% summarise(n = paste(sprintf("%s=%d", group, n), collapse = ", "))
pb_plot <- pb %>% left_join(n_lab, by = "set") %>% mutate(facet = paste0(set, "\n", n))
lab <- lab %>% left_join(n_lab, by = "set") %>% mutate(facet = paste0(set, "\n", n))
facet_levels <- unique(pb_plot$facet[order(pb_plot$set)])

p <- ggplot(pb_plot %>% mutate(facet = factor(facet, levels = facet_levels)),
            aes(x = group, y = mean_P2RY6, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.45, width = 0.5) +
  geom_jitter(width = 0.12, size = 1.4, alpha = 0.8) +
  geom_text(data = lab %>% mutate(facet = factor(facet, levels = facet_levels)),
            aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 3.2) +
  facet_wrap(~facet, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = group_colors, guide = "none") +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(x = NULL, y = "P2RY6 mean expression (pseudobulk, log-normalized)") +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(face = "bold", size = 9))

ggsave(file.path(out_dir, "Fig_P2RY6_cell_types.png"), p, width = 11, height = 7.5, dpi = 300)
ggsave(file.path(out_dir, "Fig_P2RY6_cell_types.pdf"), p, width = 11, height = 7.5)
cat("Written to", out_dir, "\n")

# ============================================================================#
# 6. CELL COMPOSITION AND P2RY6 TRANSCRIPTS (Table S2) ----
# ============================================================================#

# participants of the Healthy, AKI and CKD groups; fibroblast counts from the
# fibroblast sets of all participants with >= 1 fibroblast; P2RY6 transcripts
# are the raw P2RY6 counts (UMIs) summed over the cells/nuclei of a class
by_class <- readRDS(file.path(out_dir, "S83_pseudobulk.rds"))$by_class
pb_all   <- readRDS(file.path(out_dir, "S83_pseudobulk.rds"))$pb_all
totals <- by_class %>% group_by(dataset) %>% summarise(cells = sum(cells), umi = sum(umi))
share  <- function(sel, what) {
  by_class %>% filter(sel(class)) %>% group_by(dataset) %>% summarise(v = sum(.data[[what]])) %>%
    left_join(totals, by = "dataset") %>% mutate(pct = 100 * v / .data[[what]]) %>% pull(pct)
}
n_fib <- pb_all %>% filter(set %in% c("scRNA-seq fibroblasts", "snRNA-seq fibroblasts (Figure 2e)")) %>%
  group_by(set) %>% summarise(n = sum(n_cells)) %>% pull(n)
is_epi <- function(cl) grepl("epithelial", cl)
is_str <- function(cl) cl == "stroma cells"
vals <- rbind(share(is_epi, "cells"), share(is_epi, "umi"), share(is_str, "cells"), share(is_str, "umi"),
              100 * n_fib / totals$cells)
comp <- tibble(measure = c("Epithelial cells, % of all cells or nuclei", "Epithelial cells, % of all P2RY6 transcripts",
                           "Stroma, % of all cells or nuclei", "Stroma, % of all P2RY6 transcripts",
                           "Fibroblasts, % of all cells or nuclei"),
               `scRNA-seq` = sprintf("%.1f%%", vals[, 1]), `snRNA-seq` = sprintf("%.1f%%", vals[, 2]))
print(as.data.frame(comp))
print(totals)
write.csv(comp, file.path(out_dir, "Table_S2_cell_composition.csv"), row.names = FALSE)
