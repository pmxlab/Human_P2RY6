# ============================================================================#
# S53_scRNA_Fibroblast.R
# Tier 2 — KPMP scRNA, fibroblast-specific stroma analysis
#
# Problem with S50: using class == "stroma cells" (all subtypes) dilutes the
# PDGFRB+ fibroblast signal. P2RY6 and PDGFRB were ns because VSMC, mesangial,
# and pericyte cells wash out the fibroblast-specific expression.
#
# Solution: restrict stroma to "Perivascular Fibroblast" + "Myofibroblast"
# — the PDGFRB+ fibrogenic equivalents in this scRNA annotation, analogous to
# subclass.l1 == "FIB" in the KPMP snRNA.
#
# PT analysis: SKIPPED — reuses pb_pt from S50_pseudobulk.rds (no re-loading needed)
# FIB analysis: Perivascular Fibroblast + Myofibroblast pseudobulk by SampleID
# Outputs    : 04-Outputs/S53_scRNA_Fibroblast/
#              S53_pseudobulk.rds  (used by S54 DKD-focused variant)
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
s50_rds   <- "../04-Outputs/S50_KPMP_scRNA/S50_pseudobulk.rds"
out_dir   <- "../04-Outputs/S53_scRNA_Fibroblast"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Fibroblast subtypes to keep (PDGFRB+ fibrogenic cells)
fib_subtypes <- c("Perivascular Fibroblast", "Myofibroblast")

calcium_genes <- c("GNAQ", "ITPR2", "ITPR3", "ORAI1")
ecm_genes     <- c("COL1A1", "COL3A1", "FN1", "LUM", "DCN")

min_cells_per_sample <- 20

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

make_subtitle <- function(pb) {
  active_levels <- condition_levels[condition_levels %in% unique(as.character(pb$EnrollmentCategory))]
  active_short  <- condition_short [match(active_levels, condition_levels)]
  sprintf("n samples: %s",
    paste(mapply(function(lv, sh) {
      sprintf("%s=%d", sh, sum(as.character(pb$EnrollmentCategory) == lv, na.rm = TRUE))
    }, active_levels, active_short), collapse = ", "))
}

get_comparisons <- function(df) {
  present <- unique(as.character(df$EnrollmentCategory))
  Filter(function(p) all(p %in% present), comparisons_main)
}

# ============================================================================#
# 1. LOAD PT PSEUDOBULK FROM S50 (skip re-running PT analysis) ----
# ============================================================================#

cat("Loading PT pseudobulk from S50...\n")
s50   <- readRDS(s50_rds)
pb_pt <- s50$pb_pt
cat("PT donors from S50:", nrow(pb_pt), "\n")
rm(s50); gc()

# ============================================================================#
# 2. LOAD SCRNA OBJECT + IDENTIFY FIBROBLAST CELLS ----
# ============================================================================#

cat("\nLoading scRNA object from h5Seurat...\n")
obj <- LoadH5Seurat(sc_file)
cat("Total cells:", ncol(obj), "\n")

# Add adjudicated category from clinical CSV
clin <- read.csv(clin_file)
names(clin)[names(clin) == "Participant.ID"]               <- "participant"
names(clin)[names(clin) == "Primary.Adjudicated.Category"] <- "adjudicated"

if ("KPMPID" %in% colnames(obj@meta.data)) {
  obj$adjudicated <- clin$adjudicated[match(obj$KPMPID, clin$participant)]
} else {
  cat("WARNING: KPMPID not found — adjudicated will be NA\n")
  obj$adjudicated <- NA_character_
}

# Restrict to target conditions (reuse same logic as S50)
available_conds <- intersect(condition_levels, unique(obj$EnrollmentCategory))
cat("Conditions available:", paste(available_conds, collapse = ", "), "\n")

comparisons_main <- Filter(
  function(pair) all(pair %in% available_conds),
  comparisons_main
)

in_target  <- obj$EnrollmentCategory %in% available_conds
cells_fib  <- colnames(obj)[
  in_target &
  obj$class == "stroma cells" &
  obj$SubclassLevel2_FullName %in% fib_subtypes
]

cat("\nFibroblast cells (", paste(fib_subtypes, collapse=" + "), "):", length(cells_fib), "\n")

# Condition and subtype breakdown
meta_fib <- obj@meta.data[cells_fib, ]
cat("Fibroblast cells per condition:\n")
print(table(meta_fib$EnrollmentCategory))
cat("Fibroblast cells per subtype:\n")
print(table(meta_fib$SubclassLevel2_FullName))

# ============================================================================#
# 3. FIBROBLAST ANALYSIS ----
# ============================================================================#

cat("\n--- Fibroblast analysis (", paste(fib_subtypes, collapse=" + "), ") ---\n")
fib <- subset(obj, cells = cells_fib)
fib <- NormalizeData(fib, normalization.method = "LogNormalize",
                     scale.factor = 10000, verbose = FALSE)

cat("FIB cells:", ncol(fib), "| Samples:", length(unique(fib$SampleID)), "\n")

all_fib_genes <- c("P2RY6", "PDGFRB", calcium_genes, ecm_genes)
fib_found     <- intersect(all_fib_genes, rownames(fib))
cat("Fibroblast genes found:", paste(fib_found, collapse = ", "), "\n")

df_fib <- FetchData(fib, vars = c(fib_found, "SubclassLevel2_FullName",
                                   "SampleID", "EnrollmentCategory", "adjudicated"))

df_fib$calcium_score <- module_mean(df_fib, calcium_genes)
df_fib$ecm_score     <- module_mean(df_fib, ecm_genes)
if (!"PDGFRB" %in% fib_found) df_fib$PDGFRB <- NA_real_
if (!"P2RY6"  %in% fib_found) df_fib$P2RY6  <- NA_real_

# . 3a. Combined pseudobulk (both subtypes pooled) ----
pb_fib <- df_fib %>%
  group_by(SampleID, EnrollmentCategory, adjudicated) %>%
  summarise(
    mean_P2RY6         = mean(.data$P2RY6,        na.rm = TRUE),
    pct_P2RY6_pos      = mean(.data$P2RY6 > 0,    na.rm = TRUE) * 100,
    mean_PDGFRB        = mean(.data$PDGFRB,        na.rm = TRUE),
    mean_calcium_score = mean(.data$calcium_score, na.rm = TRUE),
    mean_ecm_score     = mean(.data$ecm_score,     na.rm = TRUE),
    n_cells            = n(),
    .groups            = "drop"
  ) %>%
  filter(n_cells >= min_cells_per_sample) %>%
  mutate(EnrollmentCategory = factor(EnrollmentCategory, levels = condition_levels))

cat("FIB donors in pseudobulk:", nrow(pb_fib), "\n")

pw_p2ry6 <- pairwise.wilcox.test(pb_fib$mean_P2RY6, pb_fib$EnrollmentCategory,
                                   p.adjust.method = "BH")
cat("FIB P2RY6 pairwise Wilcoxon (BH-FDR):\n"); print(pw_p2ry6$p.value)

pw_pdgfrb <- pairwise.wilcox.test(pb_fib$mean_PDGFRB, pb_fib$EnrollmentCategory,
                                   p.adjust.method = "BH")
cat("FIB PDGFRB pairwise Wilcoxon (BH-FDR):\n"); print(pw_pdgfrb$p.value)

subtitle_fib <- make_subtitle(pb_fib)

plot_fib_box <- function(df, y_col, y_label, title_str) {
  cmps <- get_comparisons(df)
  ymax <- max(df[[y_col]], na.rm = TRUE)
  ymin <- min(df[[y_col]], na.rm = TRUE)
  rng  <- ymax - ymin
  if (!is.finite(rng) || rng == 0) rng <- max(abs(c(ymax, ymin)), 1)
  ggplot(df, aes(x = EnrollmentCategory, y = .data[[y_col]],
                 fill = EnrollmentCategory)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.45, width = 0.5) +
    geom_jitter(width = 0.12, size = 2.2, alpha = 0.85) +
    {if (length(cmps) > 0)
       stat_compare_means(comparisons = cmps, method = "wilcox.test",
                          label = "p.signif", hide.ns = FALSE)} +
    scale_fill_manual(values = condition_colors, guide = "none") +
    scale_x_discrete(labels = setNames(condition_short, condition_levels), drop = FALSE) +
    scale_y_continuous(limits = c(ymin - 0.05 * rng, ymax + 0.26 * rng)) +
    labs(x = NULL, y = y_label, title = title_str, subtitle = subtitle_fib) +
    theme_box()
}

p_p2ry6   <- plot_fib_box(pb_fib, "mean_P2RY6",         "P2RY6 mean expression (pseudobulk)",        "KPMP scRNA — P2RY6 in fibroblasts (pvFIB + MYOF)")
p_pct     <- plot_fib_box(pb_fib, "pct_P2RY6_pos",      "% P2RY6+ cells per sample",                 "KPMP scRNA — P2RY6 detection in fibroblasts")
p_calcium <- plot_fib_box(pb_fib, "mean_calcium_score",  "Calcium signaling score (mean expression)",  "KPMP scRNA — Calcium score in fibroblasts")
p_ecm     <- plot_fib_box(pb_fib, "mean_ecm_score",      "ECM score (mean expression)",                "KPMP scRNA — ECM score in fibroblasts")
p_pdgfrb  <- plot_fib_box(pb_fib, "mean_PDGFRB",         "PDGFRB mean expression (pseudobulk)",        "KPMP scRNA — PDGFRB in fibroblasts (pvFIB + MYOF)")

ggsave(file.path(out_dir, "Fig_FB_P2RY6.png"),       p_p2ry6,   width = 6, height = 5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_FB_P2RY6_pct.png"),   p_pct,     width = 6, height = 5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_FB_calcium.png"),      p_calcium, width = 6, height = 5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_FB_ECM.png"),          p_ecm,     width = 6, height = 5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_FB_PDGFRB.png"),       p_pdgfrb,  width = 6, height = 5, dpi = dpi)

# . 3b. Per-subtype pseudobulk — P2RY6 ----
pb_fib_sub <- df_fib %>%
  group_by(SampleID, EnrollmentCategory, SubclassLevel2_FullName) %>%
  summarise(
    mean_P2RY6 = mean(.data$P2RY6, na.rm = TRUE),
    n_cells    = n(),
    .groups    = "drop"
  ) %>%
  filter(n_cells >= min_cells_per_sample) %>%
  mutate(EnrollmentCategory    = factor(EnrollmentCategory, levels = condition_levels),
         SubclassLevel2_FullName = factor(SubclassLevel2_FullName, levels = fib_subtypes))

if (nrow(pb_fib_sub) > 0 && n_distinct(pb_fib_sub$SubclassLevel2_FullName) > 1) {
  p_sub <- ggplot(pb_fib_sub, aes(x = EnrollmentCategory, y = mean_P2RY6,
                                   fill = EnrollmentCategory)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.45, width = 0.5) +
    geom_jitter(width = 0.12, size = 1.4, alpha = 0.85) +
    {cmps_sub <- get_comparisons(pb_fib_sub)
     if (length(cmps_sub) > 0)
       stat_compare_means(comparisons = cmps_sub, method = "wilcox.test",
                          label = "p.signif", hide.ns = FALSE)} +
    scale_fill_manual(values = condition_colors, guide = "none") +
    scale_x_discrete(labels = setNames(condition_short, condition_levels), drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.55))) +
    facet_wrap(~SubclassLevel2_FullName, scales = "free_y", ncol = 2) +
    labs(x = NULL, y = "P2RY6 mean expression (pseudobulk)",
         title = "KPMP scRNA — P2RY6: Perivascular Fibroblast vs Myofibroblast",
         subtitle = "Min 20 cells per sample per subtype") +
    theme_box()

  ggsave(file.path(out_dir, "Fig_FB_P2RY6_by_subtype.png"), p_sub,
         width = 9, height = 5, dpi = dpi)
}

# ============================================================================#
# 4. CROSS-COMPARTMENT CORRELATION ----
# ============================================================================#

cat("\n--- Cross-compartment correlation ---\n")

pb_cross <- inner_join(
  pb_pt  %>% select(SampleID, EnrollmentCategory, pyrim_score),
  pb_fib %>% select(SampleID, adjudicated, mean_P2RY6, pct_P2RY6_pos, mean_ecm_score),
  by = "SampleID"
) %>%
  mutate(EnrollmentCategory = factor(as.character(EnrollmentCategory),
                                     levels = condition_levels))

cat("Matched samples for correlation:", nrow(pb_cross), "\n")

ct_p2ry6 <- spearman_ci(pb_cross$pyrim_score, pb_cross$mean_P2RY6)
ct_ecm   <- spearman_ci(pb_cross$pyrim_score, pb_cross$mean_ecm_score)
cat("PT pyrimidine vs FIB P2RY6:", ct_p2ry6$label, "\n")
cat("PT pyrimidine vs FIB ECM:  ", ct_ecm$label, "\n")

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
  "PT pyrimidine score (pseudobulk)", "FIB P2RY6 (pseudobulk)",
  "KPMP scRNA — PT pyrimidine vs FIB P2RY6 (pvFIB + MYOF)", ct_p2ry6)
p_cor_ecm <- plot_cor(pb_cross, "pyrim_score", "mean_ecm_score",
  "PT pyrimidine score (pseudobulk)", "FIB ECM score (pseudobulk)",
  "KPMP scRNA — PT pyrimidine vs FIB ECM (pvFIB + MYOF)", ct_ecm)

ggsave(file.path(out_dir, "Fig_Cor_PT_Pyrimidine_vs_FB_P2RY6.png"),
       p_cor_p2ry6, width = 7, height = 4.5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_Cor_PT_Pyrimidine_vs_FB_ECM.png"),
       p_cor_ecm, width = 7, height = 4.5, dpi = dpi)

ct_p2ry6_calcium <- spearman_ci(x = pb_fib$mean_P2RY6, y = pb_fib$mean_calcium_score)
p_cor_p2ry6_calcium <- plot_cor(
  df        = pb_fib,
  x_col     = "mean_P2RY6",
  y_col     = "mean_calcium_score",
  x_lab     = "FB P2RY6 (pseudobulk)",
  y_lab     = "FB calcium score (pseudobulk)",
  title_str = "KPMP scRNA — FB P2RY6 vs Calcium (pvFIB + MYOF)",
  ct        = ct_p2ry6_calcium
) + labs(subtitle = sprintf("n = %d fibroblast samples", ct_p2ry6_calcium$n))
ggsave(filename = file.path(out_dir, "Fig_Cor_FB_P2RY6_vs_Calcium.png"),
       plot = p_cor_p2ry6_calcium, width = 7, height = 4.5, dpi = dpi)

ct_p2ry6_ecm <- spearman_ci(x = pb_fib$mean_P2RY6, y = pb_fib$mean_ecm_score)
p_cor_p2ry6_ecm <- plot_cor(
  df        = pb_fib,
  x_col     = "mean_P2RY6",
  y_col     = "mean_ecm_score",
  x_lab     = "FB P2RY6 (pseudobulk)",
  y_lab     = "FB ECM score (pseudobulk)",
  title_str = "KPMP scRNA — FB P2RY6 vs ECM (pvFIB + MYOF)",
  ct        = ct_p2ry6_ecm
) + labs(subtitle = sprintf("n = %d fibroblast samples", ct_p2ry6_ecm$n))
ggsave(filename = file.path(out_dir, "Fig_Cor_FB_P2RY6_vs_ECM.png"),
       plot = p_cor_p2ry6_ecm, width = 7, height = 4.5, dpi = dpi)

ct_pdgfrb_p2ry6 <- spearman_ci(x = pb_fib$mean_PDGFRB, y = pb_fib$mean_P2RY6)
p_cor_pdgfrb_p2ry6 <- plot_cor(
  df        = pb_fib,
  x_col     = "mean_PDGFRB",
  y_col     = "mean_P2RY6",
  x_lab     = "FB PDGFRB (pseudobulk)",
  y_lab     = "FB P2RY6 (pseudobulk)",
  title_str = "KPMP scRNA — FB PDGFRB vs P2RY6 (pvFIB + MYOF)",
  ct        = ct_pdgfrb_p2ry6
) + labs(subtitle = sprintf("n = %d fibroblast samples", ct_pdgfrb_p2ry6$n))
ggsave(filename = file.path(out_dir, "Fig_Cor_FB_PDGFRB_vs_P2RY6.png"),
       plot = p_cor_pdgfrb_p2ry6, width = 7, height = 4.5, dpi = dpi)

# ============================================================================#
# 5. SAVE PSEUDOBULK FOR S54 ----
# ============================================================================#

saveRDS(
  list(
    pb_pt      = pb_pt,
    pb_fib     = pb_fib,
    pb_fib_sub = pb_fib_sub,
    pb_cross   = pb_cross
  ),
  file.path(out_dir, "S53_pseudobulk.rds")
)

cat("\nS53 complete. Outputs saved to:", out_dir, "\n")
cat("PT samples:", nrow(pb_pt),
    "| FIB samples:", nrow(pb_fib),
    "| Matched:", nrow(pb_cross), "\n")
