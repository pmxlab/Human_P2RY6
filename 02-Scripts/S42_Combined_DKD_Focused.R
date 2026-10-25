# =============================================================================
# S42_Combined_DKD_Focused.R
# Variant of S41: KPMP CKD patients restricted to Diabetic Kidney Disease only.
# KPMP AKI patients retained as-is (all adjudicated categories).
# Fan samples retained as-is (Normal / Early_DKD / Advanced_DKD).
#
# Loads: S40_corrected.rds  (batch-corrected matrix + unified meta)
#        S30_normalised.rds  (adds primary_adjudicated_category for KPMP samples)
# Output: 04-Outputs/S42_Combined_DKD_Focused/
# =============================================================================

rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

# =============================================================================
# SETTINGS ----
# =============================================================================

s40_rds <- "../04-Outputs/S40_Combined_BulkRNA/S40_corrected.rds"
s30_rds <- "../04-Outputs/S30_KPMP_BulkRNA/S30_normalised.rds"
out_dir <- "../04-Outputs/S42_Combined_DKD_Focused"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

pyrim_genes   <- c("CAD", "DHODH", "UMPS", "UCK1", "CDA", "DCTPP1")
calcium_genes <- c("GNAQ", "ITPR2", "ITPR3", "ORAI1")
ecm_genes     <- c("COL1A1", "COL3A1", "FN1", "LUM", "DCN")

# Fan: coloured by condition
fan_colors <- c(
  "Normal"       = "#2166ac",
  "Early_DKD"    = "#f4a582",
  "Advanced_DKD" = "#d6604d"
)

# KPMP: coloured by primary adjudicated category
# (only categories remaining after CKD filter: ATI, AIN, DKD)
adj_colors <- c(
  "Acute Tubular Injury"         = "#f39c12",
  "Acute Interstitial Nephritis" = "#e74c3c",
  "Diabetic Kidney Disease"      = "#2980b9",
  "Hypertensive Kidney Disease"  = "#8e44ad",
  "Other"                        = "#27ae60",
  "Cannot be determined"         = "#95a5a6"
)

all_colors <- c(fan_colors, adj_colors)

cohort_shapes <- c("Fan" = 16, "KPMP" = 17)

condition_levels <- c("Normal", "Early_DKD", "Advanced_DKD", "AKI", "CKD")
condition_labels <- c("Normal", "Early DKD", "Adv. DKD", "AKI", "CKD\n(DKD only)")

condition_colors <- c(
  "Normal"       = "#2166ac",
  "Early_DKD"    = "#f4a582",
  "Advanced_DKD" = "#d6604d",
  "AKI"          = "#e07b39",
  "CKD"          = "#756bb1"
)

comparisons_main <- list(
  c("Normal", "Early_DKD"),
  c("Normal", "Advanced_DKD"),
  c("Normal", "AKI"),
  c("Normal", "CKD")
)

dpi <- 300

# =============================================================================
# HELPER FUNCTIONS ----
# =============================================================================

spearman_ci <- function(x, y) {
  ct  <- cor.test(x, y, method = "spearman", exact = FALSE)
  rho <- unname(ct$estimate)
  n   <- sum(!is.na(x) & !is.na(y))
  z   <- atanh(rho); se <- 1 / sqrt(n - 3)
  ci  <- tanh(c(z - 1.96 * se, z + 1.96 * se))
  list(rho = rho, ci_lo = ci[1], ci_hi = ci[2], p = ct$p.value, n = n,
       label = sprintf("rho = %.2f [%.2f, %.2f]\np = %s",
                       rho, ci[1], ci[2], format_p(ct$p.value)))
}

format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("< 0.001")
  sprintf("%.3f", p)
}

module_score <- function(mat, genes) {
  found <- intersect(genes, rownames(mat))
  if (length(found) == 0) return(rep(NA_real_, ncol(mat)))
  colMeans(mat[found, , drop = FALSE], na.rm = TRUE)
}

theme_box <- function() {
  theme_bw(base_size = 13) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1),
          plot.title  = element_text(face = "bold"),
          strip.text  = element_text(face = "bold"))
}

# =============================================================================
# 1. LOAD DATA AND APPLY DKD FILTER
# =============================================================================

s40  <- readRDS(s40_rds)
mat  <- s40$corrected_mat
meta <- s40$meta   # sample_id, condition, cohort

# Add primary_adjudicated_category for KPMP samples
s30      <- readRDS(s30_rds)
adj_info <- s30$meta %>% select(sample_id, primary_adjudicated_category)
meta     <- left_join(meta, adj_info, by = "sample_id")

# Filter: keep only KPMP CKD patients with Diabetic Kidney Disease
# (Fan samples and KPMP AKI samples are unaffected)
excluded_ckd <- meta %>%
  filter(cohort == "KPMP", condition == "CKD",
         primary_adjudicated_category != "Diabetic Kidney Disease")
if (nrow(excluded_ckd) > 0) {
  cat("CKD exclusion (non-DKD adjudicated category):\n")
  print(excluded_ckd[, c("sample_id", "condition", "primary_adjudicated_category")])
}
meta <- meta %>%
  filter(!(cohort == "KPMP" & condition == "CKD" &
             primary_adjudicated_category != "Diabetic Kidney Disease"))

mat <- mat[, meta$sample_id]

cat("\nSamples after DKD filter:", nrow(meta), "\n")
cat("Conditions:", paste(names(table(meta$condition)), table(meta$condition),
                         sep = "=", collapse = ", "), "\n")

# Build color_group: Fan -> condition, KPMP -> adjudicated category
meta$color_group <- ifelse(
  meta$cohort == "Fan",
  as.character(meta$condition),
  meta$primary_adjudicated_category
)
meta$color_group <- factor(meta$color_group, levels = names(all_colors))

# =============================================================================
# 2. MODULE SCORES
# =============================================================================

meta$pyrim_score   <- module_score(mat, pyrim_genes)
meta$calcium_score <- module_score(mat, calcium_genes)
meta$ecm_score     <- module_score(mat, ecm_genes)
meta$P2RY6  <- if ("P2RY6"  %in% rownames(mat)) as.numeric(mat["P2RY6",  ]) else NA_real_
meta$PDGFRB <- if ("PDGFRB" %in% rownames(mat)) as.numeric(mat["PDGFRB", ]) else NA_real_

subtitle_n <- {
  ns  <- table(meta$condition)
  lbs <- setNames(condition_labels, condition_levels)
  paste(sapply(condition_levels[condition_levels %in% names(ns)],
               function(lv) sprintf("%s (n=%d)", lbs[lv], ns[lv])),
        collapse = ", ")
}

# =============================================================================
# 3. BOXPLOTS
# =============================================================================

plot_score_box <- function(df, y_col, y_label, title_str,
                           cmp = comparisons_main) {
  ymax        <- max(df[[y_col]], na.rm = TRUE)
  ymin        <- min(df[[y_col]], na.rm = TRUE)
  rng         <- ymax - ymin
  upper_extra <- (0.08 * length(cmp) + 0.10) * rng

  ggplot(df, aes(x = condition, y = .data[[y_col]])) +
    geom_boxplot(aes(fill = condition), outlier.shape = NA,
                 alpha = 0.45, width = 0.55) +
    geom_jitter(aes(color = color_group, shape = cohort),
                width = 0.12, size = 2.4, alpha = 0.95) +
    stat_compare_means(comparisons = cmp, method = "wilcox.test",
                       label = "p.signif", hide.ns = FALSE) +
    scale_fill_manual(values = condition_colors, guide = "none") +
    scale_color_manual(values = all_colors, name = NULL,
                       guide = guide_legend(ncol = 1)) +
    scale_shape_manual(values = cohort_shapes, name = "Cohort") +
    scale_x_discrete(labels = setNames(condition_labels, condition_levels)) +
    scale_y_continuous(limits = c(ymin - 0.05 * rng, ymax + upper_extra)) +
    labs(x = NULL, y = y_label, title = title_str, subtitle = subtitle_n) +
    theme_box() +
    theme(legend.position = "right",
          legend.text  = element_text(size = 9),
          legend.title = element_text(size = 10))
}

p_pyrim   <- plot_score_box(meta, "pyrim_score",   "Pyrimidine module score (log2)",  "Combined (DKD-focused) — Pyrimidine score")
p_p2ry6   <- plot_score_box(meta, "P2RY6",         "P2RY6 expression (log2)",         "Combined (DKD-focused) — P2RY6")
p_calcium <- plot_score_box(meta, "calcium_score", "Calcium signaling score (log2)",  "Combined (DKD-focused) — Calcium score")
p_pdgfrb  <- plot_score_box(meta, "PDGFRB",        "PDGFRB expression (log2)",        "Combined (DKD-focused) — PDGFRB")
p_ecm     <- plot_score_box(meta, "ecm_score",     "ECM / fibrotic score (log2)",     "Combined (DKD-focused) — ECM score")

ggsave(file.path(out_dir, "Fig_pyrimidine_score.png"),   p_pyrim,   width = 7, height = 5.5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_P2RY6.png"),              p_p2ry6,   width = 7, height = 5.5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_calcium_score.png"),      p_calcium, width = 7, height = 5.5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_PDGFRB.png"),             p_pdgfrb,  width = 7, height = 5.5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_ECM_score.png"),          p_ecm,     width = 7, height = 5.5, dpi = dpi)

# --- Individual pyrimidine genes ---
pyrim_found   <- intersect(pyrim_genes, rownames(mat))
df_pyrim_long <- meta %>%
  select(sample_id, condition, cohort, color_group) %>%
  bind_cols(as.data.frame(t(mat[pyrim_found, ]))) %>%
  pivot_longer(cols = all_of(pyrim_found), names_to = "gene", values_to = "expr")

p_pyrim_genes <- ggplot(df_pyrim_long, aes(x = condition, y = expr)) +
  geom_boxplot(aes(fill = condition), outlier.shape = NA,
               alpha = 0.45, width = 0.55) +
  geom_jitter(aes(color = color_group, shape = cohort),
              width = 0.12, size = 1.4, alpha = 0.9) +
  stat_compare_means(comparisons = comparisons_main, method = "wilcox.test",
                     label = "p.signif", hide.ns = FALSE) +
  scale_fill_manual(values = condition_colors, guide = "none") +
  scale_color_manual(values = all_colors, name = NULL,
                     guide = guide_legend(ncol = 1)) +
  scale_shape_manual(values = cohort_shapes, name = "Cohort") +
  scale_x_discrete(labels = setNames(condition_labels, condition_levels)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.55))) +
  facet_wrap(~gene, scales = "free_y", ncol = 3) +
  labs(x = NULL, y = "Expression (log2)",
       title = "Combined (DKD-focused) — Individual pyrimidine genes",
       subtitle = subtitle_n) +
  theme_box() +
  theme(legend.position = "right",
        legend.text  = element_text(size = 9),
        legend.title = element_text(size = 10))

ggsave(file.path(out_dir, "Fig_pyrimidine_genes.png"), p_pyrim_genes,
       width = 13, height = 7, dpi = dpi)

# =============================================================================
# 4. CORRELATION SCATTER PLOTS
# =============================================================================

plot_scatter <- function(df, x_col, y_col, x_label, y_label, title_str) {
  ct <- spearman_ci(df[[x_col]], df[[y_col]])
  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]],
                 color = color_group, shape = cohort)) +
    geom_point(size = 2.8, alpha = 0.95) +
    geom_smooth(method = "lm", se = FALSE, color = "black",
                linewidth = 0.8, inherit.aes = FALSE,
                aes(x = .data[[x_col]], y = .data[[y_col]])) +
    scale_color_manual(values = all_colors, name = NULL,
                       guide = guide_legend(ncol = 1)) +
    scale_shape_manual(values = cohort_shapes, name = "Cohort") +
    annotate("text", x = Inf, y = Inf,
             label = ct$label, hjust = 1.05, vjust = 1.2, size = 3.8) +
    labs(x = x_label, y = y_label, title = title_str, subtitle = subtitle_n) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "right",
          legend.text  = element_text(size = 9),
          legend.title = element_text(size = 10))
}

p_cor_pyrim_p2ry6   <- plot_scatter(meta, "pyrim_score", "P2RY6",
  "Pyrimidine score", "P2RY6 (log2)", "Combined (DKD-focused) — Pyrimidine vs P2RY6")
p_cor_pyrim_pdgfrb  <- plot_scatter(meta, "pyrim_score", "PDGFRB",
  "Pyrimidine score", "PDGFRB (log2)", "Combined (DKD-focused) — Pyrimidine vs PDGFRB")
p_cor_pyrim_ecm     <- plot_scatter(meta, "pyrim_score", "ecm_score",
  "Pyrimidine score", "ECM score", "Combined (DKD-focused) — Pyrimidine vs ECM")
p_cor_p2ry6_ecm     <- plot_scatter(meta, "P2RY6", "ecm_score",
  "P2RY6 (log2)", "ECM score", "Combined (DKD-focused) — P2RY6 vs ECM")
p_cor_p2ry6_calcium <- plot_scatter(meta, "P2RY6", "calcium_score",
  "P2RY6 (log2)", "Calcium score", "Combined (DKD-focused) — P2RY6 vs Calcium")
p_cor_pdgfrb_p2ry6  <- plot_scatter(meta, "PDGFRB", "P2RY6",
  "PDGFRB (log2)", "P2RY6 (log2)", "Combined (DKD-focused) — PDGFRB vs P2RY6")

ggsave(file.path(out_dir, "Fig_cor_Pyrimidine_P2RY6.png"),   p_cor_pyrim_p2ry6,   width = 7, height = 4.5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_cor_Pyrimidine_PDGFRB.png"),  p_cor_pyrim_pdgfrb,  width = 7, height = 4.5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_cor_Pyrimidine_ECM.png"),     p_cor_pyrim_ecm,     width = 7, height = 4.5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_cor_P2RY6_ECM.png"),          p_cor_p2ry6_ecm,     width = 7, height = 4.5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_cor_P2RY6_Calcium.png"),      p_cor_p2ry6_calcium, width = 7, height = 4.5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_cor_PDGFRB_P2RY6.png"),       p_cor_pdgfrb_p2ry6,  width = 7, height = 4.5, dpi = dpi)

cat("\nS42 complete. Outputs saved to:", out_dir, "\n")
cat("Final N =", nrow(meta), "| Conditions:",
    paste(names(table(meta$condition)), table(meta$condition),
          sep = "=", collapse = ", "), "\n")
