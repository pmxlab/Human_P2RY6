# ============================================================================#
# S50_KPMP_scRNA_loadOnce.R
# Tier 2 — KPMP single-cell RNA-seq (integrated, batch-corrected)
#
# Dataset : KPMP_PREMIERE_SC_version2_ForExplorer_RemovedBatchEffect_Final2025
#           This integrated object includes Lake 2023 scRNA cells.
# Conditions : Healthy Reference Tissue + CKD (+ AKI if present after table check)
# PT analysis : cell_type == "epithelial cell of proximal tubule"
#               pyrimidine module score, pseudobulk by SampleID
# FIB analysis: class == "stroma cells" — pooled pseudobulk + per-subtype
# Outputs    : 04-Outputs/S50_KPMP_scRNA/
#              S50_pseudobulk.rds  (used by S51 and S52)
# ============================================================================#

rm(list = ls())

library(Seurat)
library(SeuratDisk)
library(hdf5r)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

# ============================================================================#
# SETTINGS ----
# ============================================================================#

sc_file   <- "../../KPMP_scRNA/KPMP_PREMIERE_SC_version2_ForExplorer_RemovedBatchEffect_Final2025.h5seurat"
clin_file <- "../../KPMP_scRNA/OpenAccessClinicalData/20251205_OpenAccessClinicalData.csv"
out_dir   <- "../04-Outputs/S50_KPMP_scRNA"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

pyrim_genes   <- c("CAD", "DHODH", "UMPS", "UCK1", "CDA", "DCTPP1")
calcium_genes <- c("GNAQ", "ITPR2", "ITPR3", "ORAI1")
ecm_genes     <- c("COL1A1", "COL3A1", "FN1", "LUM", "DCN")

min_cells_per_sample <- 20

# Expected EnrollmentCategory values in this dataset
# (AKI will be included only if actually present — checked below)
conditions_keep  <- c("Healthy Reference Tissue",
                       "Acute Kidney Injury",
                       "Chronic Kidney Disease")
condition_levels <- c("Healthy Reference Tissue",
                       "Acute Kidney Injury",
                       "Chronic Kidney Disease")
condition_short  <- c("Healthy", "AKI", "CKD")

condition_colors <- c(
  "Healthy Reference Tissue" = "#2166ac",
  "Acute Kidney Injury"      = "#e07b39",
  "Chronic Kidney Disease"   = "#756bb1"
)

comparisons_main <- list(
  c("Healthy Reference Tissue", "Acute Kidney Injury"),
  c("Healthy Reference Tissue", "Chronic Kidney Disease")
)

dpi <- 300

# ============================================================================#
# HELPER FUNCTIONS ----
# ============================================================================#

spearman_ci <- function(x, y) {
  keep <- !is.na(x) & !is.na(y)
  x <- x[keep]; y <- y[keep]
  if (length(x) < 4) {
    return(list(rho = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
                p = NA_real_, n = length(x),
                label = "rho = NA [NA, NA]\np = NA"))
  }
  ct  <- cor.test(x, y, method = "spearman", exact = FALSE)
  rho <- unname(ct$estimate)
  n   <- length(x)
  if (is.na(rho) || abs(rho) >= 1) {
    return(list(rho = rho, ci_lo = NA_real_, ci_hi = NA_real_,
                p = ct$p.value, n = n,
                label = sprintf("rho = %.2f [NA, NA]\np = %s", rho, format_p(ct$p.value))))
  }
  z  <- atanh(rho); se <- 1 / sqrt(n - 3)
  ci <- tanh(c(z - 1.96 * se, z + 1.96 * se))
  list(rho = rho, ci_lo = ci[1], ci_hi = ci[2], p = ct$p.value, n = n,
       label = sprintf("rho = %.2f [%.2f, %.2f]\np = %s",
                       rho, ci[1], ci[2], format_p(ct$p.value)))
}

format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("< 0.001")
  sprintf("%.3f", p)
}

module_mean <- function(df, genes) {
  found <- intersect(genes, colnames(df))
  if (length(found) == 0) return(rep(NA_real_, nrow(df)))
  rowMeans(df[, found, drop = FALSE], na.rm = TRUE)
}

theme_box <- function() {
  theme_bw(base_size = 13) +
    theme(axis.text.x  = element_text(angle = 20, hjust = 1),
          plot.title   = element_text(face = "bold"),
          strip.text   = element_text(face = "bold"))
}

make_subtitle <- function(pb, cond_col = "EnrollmentCategory") {
  cond_vals <- pb[[cond_col]]
  active_levels <- condition_levels[condition_levels %in% unique(cond_vals)]
  active_short  <- condition_short [match(active_levels, condition_levels)]
  sprintf("n samples: %s",
    paste(mapply(function(lv, sh) {
      n <- sum(as.character(cond_vals) == lv, na.rm = TRUE)
      sprintf("%s=%d", sh, n)
    }, active_levels, active_short), collapse = ", "))
}

# ============================================================================#
# 1. LOAD OBJECT ONCE + ADD CLINICAL METADATA ----
# ============================================================================#

cat("Loading scRNA object from h5Seurat...\n")
obj <- LoadH5Seurat(sc_file)

cat("Total cells:", ncol(obj), "| Metadata cols:", ncol(obj@meta.data), "\n")
cat("Metadata columns:\n"); print(colnames(obj@meta.data))

# Inspect available enrollment categories
cat("\nEnrollmentCategory distribution (all cells):\n")
print(table(obj@meta.data$EnrollmentCategory, useNA = "ifany"))

# Restrict conditions_keep to those actually present
available_conds <- intersect(conditions_keep, unique(obj@meta.data$EnrollmentCategory))
cat("\nConditions available for analysis:", paste(available_conds, collapse = ", "), "\n")

# Adjust comparisons to only include present conditions
comparisons_main <- Filter(
  function(pair) all(pair %in% available_conds),
  comparisons_main
)
cat("Comparisons used:", length(comparisons_main), "\n")

# Add adjudicated category from clinical CSV (needed for S52 DKD filter)
clin <- read.csv(clin_file)
names(clin)[names(clin) == "Participant.ID"]               <- "participant"
names(clin)[names(clin) == "Primary.Adjudicated.Category"] <- "adjudicated"

# KPMPID in the scRNA object maps to Participant.ID in clinical CSV
if ("KPMPID" %in% colnames(obj@meta.data)) {
  obj$adjudicated <- clin$adjudicated[match(obj$KPMPID, clin$participant)]
  cat("\nAdjudicated categories added via KPMPID.\n")
} else {
  cat("WARNING: KPMPID column not found. Adjudicated will be NA.\n")
  obj$adjudicated <- NA_character_
}

cat("Adjudicated distribution (CKD samples only):\n")
ckd_mask <- obj$EnrollmentCategory == "Chronic Kidney Disease"
print(table(unique(obj@meta.data[ckd_mask, c("SampleID", "adjudicated")])$adjudicated,
            useNA = "ifany"))

# Cell barcodes for each compartment (filtered to target conditions)
in_target  <- obj$EnrollmentCategory %in% available_conds
cells_pt   <- colnames(obj)[in_target & obj$cell_type == "epithelial cell of proximal tubule"]
cells_fib  <- colnames(obj)[in_target & obj$class     == "stroma cells"]

cat("\nTarget-condition cells — PT:", length(cells_pt), "| Stroma:", length(cells_fib), "\n")

# Patient counts per condition
sample_counts <- table(unique(
  obj@meta.data[in_target, c("SampleID", "EnrollmentCategory")])$EnrollmentCategory)
cat("Samples per condition:", paste(names(sample_counts), sample_counts,
                                    sep = "=", collapse = ", "), "\n")

# ============================================================================#
# 2. PT ANALYSIS — Pyrimidine module score ----
# ============================================================================#

cat("\n--- PT analysis ---\n")
pt <- subset(obj, cells = cells_pt)
pt <- NormalizeData(pt, normalization.method = "LogNormalize",
                    scale.factor = 10000, verbose = FALSE)

cat("PT cells:", ncol(pt), "| Samples:", length(unique(pt$SampleID)), "\n")

pyrim_found <- intersect(pyrim_genes, rownames(pt))
cat("Pyrimidine genes found:", paste(pyrim_found, collapse = ", "), "\n")

df_pt <- FetchData(pt, vars = c(pyrim_found, "SampleID", "EnrollmentCategory"))
df_pt$pyrim_score <- module_mean(df_pt, pyrim_found)

pb_pt <- df_pt %>%
  group_by(SampleID, EnrollmentCategory) %>%
  summarise(
    pyrim_score = mean(pyrim_score, na.rm = TRUE),
    across(all_of(pyrim_found), mean),
    n_cells = n(),
    .groups = "drop"
  ) %>%
  filter(n_cells >= min_cells_per_sample) %>%
  mutate(EnrollmentCategory = factor(EnrollmentCategory, levels = condition_levels))

cat("PT samples in pseudobulk:", nrow(pb_pt), "\n")
pw_pt <- pairwise.wilcox.test(pb_pt$pyrim_score, pb_pt$EnrollmentCategory,
                               p.adjust.method = "BH")
cat("PT pyrimidine pairwise Wilcoxon (BH-FDR):\n"); print(pw_pt$p.value)

subtitle_pt <- make_subtitle(pb_pt)

# . Pyrimidine score boxplot ----
p_pt_score <- ggplot(pb_pt, aes(x = EnrollmentCategory, y = pyrim_score,
                                 fill = EnrollmentCategory)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.45, width = 0.5) +
  geom_jitter(width = 0.12, size = 2.2, alpha = 0.85) +
  stat_compare_means(comparisons = comparisons_main, method = "wilcox.test",
                     label = "p.signif", hide.ns = FALSE) +
  scale_fill_manual(values = condition_colors, guide = "none") +
  scale_x_discrete(labels = setNames(condition_short, condition_levels), drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.30))) +
  labs(x = NULL, y = "Pyrimidine module score (mean expression)",
       title = "KPMP scRNA — PT pyrimidine score (pseudobulk by sample)",
       subtitle = subtitle_pt) +
  theme_box()

ggsave(file.path(out_dir, "Fig_PT_pyrimidine_score.png"), p_pt_score,
       width = 6, height = 5, dpi = dpi)

# . Individual pyrimidine genes ----
df_pt_long <- pb_pt %>%
  select(SampleID, EnrollmentCategory, n_cells, all_of(pyrim_found)) %>%
  pivot_longer(cols = all_of(pyrim_found), names_to = "gene", values_to = "expr")

p_pt_genes <- ggplot(df_pt_long, aes(x = EnrollmentCategory, y = expr,
                                      fill = EnrollmentCategory)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.45, width = 0.5) +
  geom_jitter(width = 0.12, size = 1.4, alpha = 0.85) +
  stat_compare_means(comparisons = comparisons_main, method = "wilcox.test",
                     label = "p.signif", hide.ns = FALSE) +
  scale_fill_manual(values = condition_colors, guide = "none") +
  scale_x_discrete(labels = setNames(condition_short, condition_levels), drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.55))) +
  facet_wrap(~gene, scales = "free_y", ncol = 3) +
  labs(x = NULL, y = "Mean expression (pseudobulk by sample)",
       title = "KPMP scRNA — Individual pyrimidine genes in PT",
       subtitle = subtitle_pt) +
  theme_box()

ggsave(file.path(out_dir, "Fig_PT_pyrimidine_genes.png"), p_pt_genes,
       width = 11, height = 7, dpi = dpi)

# ============================================================================#
# 3. FIBROBLAST / STROMA ANALYSIS — class == "stroma cells" ----
# ============================================================================#

cat("\n--- Stroma / fibroblast analysis ---\n")
fib <- subset(obj, cells = cells_fib)
fib <- NormalizeData(fib, normalization.method = "LogNormalize",
                     scale.factor = 10000, verbose = FALSE)

cat("Stroma cells:", ncol(fib), "| Samples:", length(unique(fib$SampleID)), "\n")
cat("Stroma subtype counts (SubclassLevel2_FullName):\n")
print(sort(table(fib$SubclassLevel2_FullName), decreasing = TRUE))

all_fib_genes <- c("P2RY6", "PDGFRB", calcium_genes, ecm_genes)
fib_found     <- intersect(all_fib_genes, rownames(fib))
cat("Stroma genes found:", paste(fib_found, collapse = ", "), "\n")

df_fib <- FetchData(fib, vars = c(fib_found, "SubclassLevel2_FullName",
                                   "SampleID", "EnrollmentCategory", "adjudicated"))

df_fib$calcium_score <- module_mean(df_fib, calcium_genes)
df_fib$ecm_score     <- module_mean(df_fib, ecm_genes)
if (!"PDGFRB" %in% fib_found) df_fib$PDGFRB <- NA_real_
if (!"P2RY6"  %in% fib_found) df_fib$P2RY6  <- NA_real_

# . 3a. Combined stroma pseudobulk (all subtypes pooled) ----
pb_fib <- df_fib %>%
  group_by(SampleID, EnrollmentCategory, adjudicated) %>%
  summarise(
    mean_P2RY6         = mean(.data$P2RY6,         na.rm = TRUE),
    pct_P2RY6_pos      = mean(.data$P2RY6 > 0,     na.rm = TRUE) * 100,
    mean_PDGFRB        = mean(.data$PDGFRB,         na.rm = TRUE),
    mean_calcium_score = mean(.data$calcium_score,  na.rm = TRUE),
    mean_ecm_score     = mean(.data$ecm_score,      na.rm = TRUE),
    n_cells            = n(),
    .groups            = "drop"
  ) %>%
  filter(n_cells >= min_cells_per_sample) %>%
  mutate(EnrollmentCategory = factor(EnrollmentCategory, levels = condition_levels))

cat("Stroma samples in pseudobulk:", nrow(pb_fib), "\n")

pw_fib_p2ry6 <- pairwise.wilcox.test(pb_fib$mean_P2RY6, pb_fib$EnrollmentCategory,
                                       p.adjust.method = "BH")
cat("Stroma P2RY6 pairwise Wilcoxon (BH-FDR):\n"); print(pw_fib_p2ry6$p.value)

subtitle_fib <- make_subtitle(pb_fib)

plot_fib_box <- function(df, y_col, y_label, title_str) {
  ymax <- max(df[[y_col]], na.rm = TRUE)
  ymin <- min(df[[y_col]], na.rm = TRUE)
  rng  <- ymax - ymin
  if (!is.finite(rng) || rng == 0) rng <- max(abs(c(ymax, ymin)), 1)
  ggplot(df, aes(x = EnrollmentCategory, y = .data[[y_col]],
                 fill = EnrollmentCategory)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.45, width = 0.5) +
    geom_jitter(width = 0.12, size = 2.2, alpha = 0.85) +
    stat_compare_means(comparisons = comparisons_main, method = "wilcox.test",
                       label = "p.signif", hide.ns = FALSE) +
    scale_fill_manual(values = condition_colors, guide = "none") +
    scale_x_discrete(labels = setNames(condition_short, condition_levels), drop = FALSE) +
    scale_y_continuous(limits = c(ymin - 0.05 * rng, ymax + 0.26 * rng)) +
    labs(x = NULL, y = y_label, title = title_str, subtitle = subtitle_fib) +
    theme_box()
}

p_p2ry6   <- plot_fib_box(pb_fib, "mean_P2RY6",         "P2RY6 mean expression (pseudobulk)",        "KPMP scRNA — P2RY6 in stroma cells")
p_pct     <- plot_fib_box(pb_fib, "pct_P2RY6_pos",      "% P2RY6+ cells per sample",                 "KPMP scRNA — P2RY6 detection in stroma")
p_calcium <- plot_fib_box(pb_fib, "mean_calcium_score",  "Calcium signaling score (mean expression)",  "KPMP scRNA — Calcium score in stroma")
p_ecm     <- plot_fib_box(pb_fib, "mean_ecm_score",      "ECM score (mean expression)",                "KPMP scRNA — ECM score in stroma")
p_pdgfrb  <- plot_fib_box(pb_fib, "mean_PDGFRB",         "PDGFRB mean expression (pseudobulk)",        "KPMP scRNA — PDGFRB in stroma cells")

ggsave(file.path(out_dir, "Fig_FB_P2RY6.png"),       p_p2ry6,   width = 6, height = 5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_FB_P2RY6_pct.png"),   p_pct,     width = 6, height = 5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_FB_calcium.png"),      p_calcium, width = 6, height = 5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_FB_ECM.png"),          p_ecm,     width = 6, height = 5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_FB_PDGFRB.png"),       p_pdgfrb,  width = 6, height = 5, dpi = dpi)

# . 3b. Per-subtype pseudobulk — P2RY6 ----
# Identify subtypes with sufficient representation
subtype_counts <- sort(table(fib$SubclassLevel2_FullName), decreasing = TRUE)
cat("\nSubtype cell counts:\n"); print(subtype_counts)

# Keep subtypes with at least 3 samples meeting the min_cells threshold
pb_fib_sub <- df_fib %>%
  group_by(SampleID, EnrollmentCategory, SubclassLevel2_FullName) %>%
  summarise(
    mean_P2RY6 = mean(.data$P2RY6, na.rm = TRUE),
    n_cells    = n(),
    .groups    = "drop"
  ) %>%
  filter(n_cells >= min_cells_per_sample) %>%
  mutate(EnrollmentCategory    = factor(EnrollmentCategory, levels = condition_levels),
         SubclassLevel2_FullName = factor(SubclassLevel2_FullName))

# Only keep subtypes with data in ≥2 conditions
subtype_cond_count <- pb_fib_sub %>%
  group_by(SubclassLevel2_FullName) %>%
  summarise(n_conds = n_distinct(EnrollmentCategory), .groups = "drop") %>%
  filter(n_conds >= 2)
pb_fib_sub <- pb_fib_sub %>%
  filter(SubclassLevel2_FullName %in% subtype_cond_count$SubclassLevel2_FullName)

cat("Subtypes retained for per-subtype plot:",
    paste(levels(droplevels(pb_fib_sub$SubclassLevel2_FullName)), collapse = ", "), "\n")

if (nrow(pb_fib_sub) > 0) {
  n_sub <- n_distinct(pb_fib_sub$SubclassLevel2_FullName)
  p_sub_p2ry6 <- ggplot(pb_fib_sub, aes(x = EnrollmentCategory, y = mean_P2RY6,
                                          fill = EnrollmentCategory)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.45, width = 0.5) +
    geom_jitter(width = 0.12, size = 1.4, alpha = 0.85) +
    stat_compare_means(comparisons = comparisons_main, method = "wilcox.test",
                       label = "p.signif", hide.ns = FALSE) +
    scale_fill_manual(values = condition_colors, guide = "none") +
    scale_x_discrete(labels = setNames(condition_short, condition_levels), drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.55))) +
    facet_wrap(~SubclassLevel2_FullName, scales = "free_y",
               ncol = min(4, ceiling(sqrt(n_sub)))) +
    labs(x = NULL, y = "P2RY6 mean expression (pseudobulk)",
         title = "KPMP scRNA — P2RY6 per stroma subtype",
         subtitle = "Min 20 cells per sample per subtype") +
    theme_box()

  ggsave(file.path(out_dir, "Fig_FB_P2RY6_by_subtype.png"), p_sub_p2ry6,
         width = 14, height = 8, dpi = dpi)
}

# ============================================================================#
# 4. CROSS-COMPARTMENT CORRELATION ----
# ============================================================================#

cat("\n--- Cross-compartment correlation ---\n")

pb_cross <- inner_join(
  pb_pt  %>% select(SampleID, EnrollmentCategory, pyrim_score),
  pb_fib %>% select(SampleID, mean_P2RY6, pct_P2RY6_pos, mean_ecm_score),
  by = "SampleID"
) %>%
  mutate(EnrollmentCategory = factor(as.character(EnrollmentCategory),
                                     levels = condition_levels))

cat("Matched samples for correlation:", nrow(pb_cross), "\n")

ct_p2ry6 <- spearman_ci(pb_cross$pyrim_score, pb_cross$mean_P2RY6)
ct_ecm   <- spearman_ci(pb_cross$pyrim_score, pb_cross$mean_ecm_score)
cat("PT pyrimidine vs Stroma P2RY6:", ct_p2ry6$label, "\n")
cat("PT pyrimidine vs Stroma ECM:  ", ct_ecm$label, "\n")

plot_cor <- function(df, x_col, y_col, x_lab, y_lab, title_str, ct) {
  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]],
                 color = EnrollmentCategory)) +
    geom_point(size = 2.8, alpha = 0.85) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8,
                inherit.aes = FALSE, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    annotate("text", x = Inf, y = Inf,
             label = ct$label, hjust = 1.05, vjust = 1.2, size = 3.8) +
    scale_color_manual(values = condition_colors, name = NULL) +
    labs(x = x_lab, y = y_lab, title = title_str,
         subtitle = sprintf("n = %d matched samples", nrow(df))) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold"), legend.position = "right")
}

p_cor_p2ry6 <- plot_cor(pb_cross, "pyrim_score", "mean_P2RY6",
  "PT pyrimidine score (pseudobulk)", "Stroma P2RY6 (pseudobulk)",
  "KPMP scRNA — PT pyrimidine vs Stroma P2RY6", ct_p2ry6)
p_cor_ecm <- plot_cor(pb_cross, "pyrim_score", "mean_ecm_score",
  "PT pyrimidine score (pseudobulk)", "Stroma ECM score (pseudobulk)",
  "KPMP scRNA — PT pyrimidine vs Stroma ECM", ct_ecm)

ggsave(file.path(out_dir, "Fig_Cor_PT_Pyrimidine_vs_FB_P2RY6.png"),
       p_cor_p2ry6, width = 7, height = 4.5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_Cor_PT_Pyrimidine_vs_FB_ECM.png"),
       p_cor_ecm, width = 7, height = 4.5, dpi = dpi)

ct_p2ry6_calcium <- spearman_ci(x = pb_fib$mean_P2RY6, y = pb_fib$mean_calcium_score)
p_cor_p2ry6_calcium <- plot_cor(
  df        = pb_fib,
  x_col     = "mean_P2RY6",
  y_col     = "mean_calcium_score",
  x_lab     = "Stroma P2RY6 (pseudobulk)",
  y_lab     = "Stroma calcium score (pseudobulk)",
  title_str = "KPMP scRNA — Stroma P2RY6 vs Calcium",
  ct        = ct_p2ry6_calcium
) + labs(subtitle = sprintf("n = %d stroma samples", ct_p2ry6_calcium$n))
ggsave(filename = file.path(out_dir, "Fig_Cor_FB_P2RY6_vs_Calcium.png"),
       plot = p_cor_p2ry6_calcium, width = 7, height = 4.5, dpi = dpi)

ct_p2ry6_ecm <- spearman_ci(x = pb_fib$mean_P2RY6, y = pb_fib$mean_ecm_score)
p_cor_p2ry6_ecm <- plot_cor(
  df        = pb_fib,
  x_col     = "mean_P2RY6",
  y_col     = "mean_ecm_score",
  x_lab     = "Stroma P2RY6 (pseudobulk)",
  y_lab     = "Stroma ECM score (pseudobulk)",
  title_str = "KPMP scRNA — Stroma P2RY6 vs ECM",
  ct        = ct_p2ry6_ecm
) + labs(subtitle = sprintf("n = %d stroma samples", ct_p2ry6_ecm$n))
ggsave(filename = file.path(out_dir, "Fig_Cor_FB_P2RY6_vs_ECM.png"),
       plot = p_cor_p2ry6_ecm, width = 7, height = 4.5, dpi = dpi)

ct_pdgfrb_p2ry6 <- spearman_ci(x = pb_fib$mean_PDGFRB, y = pb_fib$mean_P2RY6)
p_cor_pdgfrb_p2ry6 <- plot_cor(
  df        = pb_fib,
  x_col     = "mean_PDGFRB",
  y_col     = "mean_P2RY6",
  x_lab     = "Stroma PDGFRB (pseudobulk)",
  y_lab     = "Stroma P2RY6 (pseudobulk)",
  title_str = "KPMP scRNA — Stroma PDGFRB vs P2RY6",
  ct        = ct_pdgfrb_p2ry6
) + labs(subtitle = sprintf("n = %d stroma samples", ct_pdgfrb_p2ry6$n))
ggsave(filename = file.path(out_dir, "Fig_Cor_FB_PDGFRB_vs_P2RY6.png"),
       plot = p_cor_pdgfrb_p2ry6, width = 7, height = 4.5, dpi = dpi)

# ============================================================================#
# 5. SAVE PSEUDOBULK DATA FOR S51 / S52 ----
# ============================================================================#

saveRDS(
  list(
    pb_pt      = pb_pt,
    pb_fib     = pb_fib,
    pb_fib_sub = pb_fib_sub,
    pb_cross   = pb_cross
  ),
  file.path(out_dir, "S50_pseudobulk.rds")
)

cat("\nS50 complete. Outputs saved to:", out_dir, "\n")
cat("PT samples:", nrow(pb_pt),
    "| Stroma samples:", nrow(pb_fib),
    "| Matched:", nrow(pb_cross), "\n")
