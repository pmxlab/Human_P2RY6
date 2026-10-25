# =========================================================================== #
# S40_Combined_BulkRNA.R
# Combined bulk RNA cohort: Fan DKD (GSE142025) + KPMP BulkRNA
# Tier 1, Pooled: ComBat batch correction -> joint analysis
#
# Fan  (Cohort A): Normal n=9, Early_DKD n=6, Advanced_DKD n=21  [N=36]
# KPMP (Cohort B): AKI n=12, CKD n=15                            [N=27]
# Total: N=63 across 5 disease states
#
# Batch correction: ComBat (log-domain; sva package)
#   batch     = cohort (Fan / KPMP)
#   covariate = Normal vs Diseased (binary, avoids confounding since
#               Normal samples exist only in Fan)
#
# NORMALIZATION LIMITATION:
#   Fan (S20) values are log2(quantile-normalized counts), pre-processed by the
#   authors before GEO deposit; raw counts are not publicly available and TMM
#   normalization cannot be applied retrospectively.
#   KPMP (S30) values are TMM log2 CPM, computed from raw counts in-script.
#   The two cohorts therefore enter ComBat on different normalization scales.
#   ComBat's per-gene location-and-scale correction absorbs this difference as
#   part of the batch effect. This is the best feasible approach given data
#   availability and should be stated as a limitation in the Methods section.
# =========================================================================== #

rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(sva)

# =========================================================================== #
# SETTINGS ----
# =========================================================================== #

s20_rds <- "../04-Outputs/S20_FanDKD_BulkRNA/S20_normalised.rds"
s30_rds <- "../04-Outputs/S30_KPMP_BulkRNA/S30_normalised.rds"
out_dir <- "../04-Outputs/S40_Combined_BulkRNA"
dir.create(file.path(out_dir, "QC"), showWarnings = FALSE, recursive = TRUE)

pyrim_genes   <- c("CAD", "DHODH", "UMPS", "UCK1", "CDA", "DCTPP1")
calcium_genes <- c("GNAQ", "ITPR2", "ITPR3", "ORAI1")
ecm_genes     <- c("COL1A1", "COL3A1", "FN1", "LUM", "DCN")

condition_levels <- c("Normal", "Early_DKD", "Advanced_DKD", "AKI", "CKD")
condition_labels <- c("Normal", "Early DKD", "Adv. DKD", "AKI", "CKD")
condition_colors <- c(
  "Normal"       = "#2166ac",
  "Early_DKD"    = "#f4a582",
  "Advanced_DKD" = "#d6604d",
  "AKI"          = "#e07b39",
  "CKD"          = "#756bb1"
)
cohort_shapes <- c("Fan" = 16, "KPMP" = 17)   # circle vs triangle

# Comparisons: Normal as reference vs each disease state
comparisons_main <- list(
  c("Normal", "Early_DKD"),
  c("Normal", "Advanced_DKD"),
  c("Normal", "AKI"),
  c("Normal", "CKD")
)

dpi <- 300

# =========================================================================== #
# HELPER FUNCTIONS ----
# =========================================================================== #

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

wilcox_fdr <- function(df, value_col, group_col, comparisons) {
  pvals <- sapply(comparisons, function(cmp) {
    a <- df[[value_col]][df[[group_col]] == cmp[1]]
    b <- df[[value_col]][df[[group_col]] == cmp[2]]
    if (length(a) < 2 || length(b) < 2) return(NA_real_)
    wilcox.test(a, b, exact = FALSE)$p.value
  })
  padj <- p.adjust(pvals, method = "BH")
  data.frame(
    group1  = sapply(comparisons, `[`, 1),
    group2  = sapply(comparisons, `[`, 2),
    p_value = pvals, p_adj = padj,
    signif  = ifelse(is.na(padj), "n/a",
               ifelse(padj < 0.001, "***",
                ifelse(padj < 0.01, "**",
                 ifelse(padj < 0.05, "*", "ns"))))
  )
}

module_score <- function(mat, genes) {
  found <- intersect(genes, rownames(mat))
  if (length(found) == 0) return(rep(NA_real_, ncol(mat)))
  colMeans(mat[found, , drop = FALSE], na.rm = TRUE)
}

theme_box <- function() {
  theme_bw(base_size = 13) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 20, hjust = 1),
          plot.title  = element_text(face = "bold"),
          strip.text  = element_text(face = "bold"))
}

n_subtitle <- function(meta) {
  ns  <- table(meta$condition)
  lbs <- setNames(condition_labels, condition_levels)
  paste(sapply(condition_levels[condition_levels %in% names(ns)],
               function(lv) sprintf("%s (n=%d)", lbs[lv], ns[lv])),
        collapse = ", ")
}

# =========================================================================== #
# 1. LOAD NORMALISED MATRICES FROM S20 AND S30 ----
# =========================================================================== #

s20 <- readRDS(s20_rds)
s30 <- readRDS(s30_rds)

fan_mat  <- s20$expr_mat
fan_meta <- s20$meta %>%
  rename(sample_id = sample) %>%
  mutate(cohort = "Fan",
         condition = as.character(condition))

kpmp_mat  <- s30$expr_mat
kpmp_meta <- s30$meta %>%
  mutate(cohort    = "KPMP",
         condition = as.character(enrollment_category)) %>%
  select(sample_id, condition, cohort)

# Common genes
common_genes <- intersect(rownames(fan_mat), rownames(kpmp_mat))
cat("Common genes:", length(common_genes), "\n")

fan_mat  <- fan_mat[common_genes, ]
kpmp_mat <- kpmp_mat[common_genes, ]

# =========================================================================== #
# 2. MERGE AND BUILD UNIFIED METADATA ----
# =========================================================================== #

combined_mat <- cbind(fan_mat, kpmp_mat)

meta <- bind_rows(
  fan_meta  %>% select(sample_id, condition, cohort),
  kpmp_meta %>% select(sample_id, condition, cohort)
)
meta$condition <- factor(meta$condition, levels = condition_levels)
meta$cohort    <- factor(meta$cohort, levels = c("Fan", "KPMP"))

stopifnot(identical(colnames(combined_mat), meta$sample_id))

cat("Total samples:", nrow(meta), "\n")
cat("Conditions:", paste(names(table(meta$condition)), table(meta$condition),
                         sep = "=", collapse = ", "), "\n")

# =========================================================================== #
# 3. COMBAT BATCH CORRECTION (log-domain) ----
# =========================================================================== #
# Biological covariate: binary Normal vs Diseased (Normal only in Fan —
# using the full condition would create a rank-deficient design matrix).

is_normal <- as.integer(meta$condition == "Normal")
mod_combat <- model.matrix(~ is_normal)
batch      <- as.character(meta$cohort)

cat("Running ComBat...\n")
corrected_mat <- ComBat(
  dat   = combined_mat,
  batch = batch,
  mod   = mod_combat,
  par.prior = TRUE,
  prior.plots = FALSE
)
cat("ComBat complete.\n")

# =========================================================================== #
# 4. PCA — PRE vs POST CORRECTION (QC figure) ----
# =========================================================================== #

top_var_genes <- function(mat, n = 2000) {
  rv <- apply(mat, 1, var)
  head(order(rv, decreasing = TRUE), n)
}

pca_plot <- function(mat, meta, title_str, color_by = "condition",
                     shape_by = "cohort") {
  idx  <- top_var_genes(mat, 2000)
  pca  <- prcomp(t(mat[idx, ]), scale. = TRUE)
  pve  <- round(summary(pca)$importance[2, 1:2] * 100, 1)
  df   <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                     condition = meta$condition, cohort = meta$cohort)
  ggplot(df, aes(x = PC1, y = PC2,
                 color = .data[[color_by]],
                 shape = .data[[shape_by]])) +
    geom_point(size = 3, alpha = 0.9) +
    scale_color_manual(values = condition_colors) +
    scale_shape_manual(values = cohort_shapes) +
    labs(x = sprintf("PC1 (%.1f%%)", pve[1]),
         y = sprintf("PC2 (%.1f%%)", pve[2]),
         title = title_str, color = NULL, shape = NULL) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))
}

p_pca_pre  <- pca_plot(combined_mat, meta,
                       "PCA — pre-ComBat (top 2000 variable genes)")
p_pca_post <- pca_plot(corrected_mat, meta,
                       "PCA — post-ComBat (top 2000 variable genes)")

ggsave(file.path(out_dir, "QC", "PCA_pre_ComBat.png"),  p_pca_pre,
       width = 6, height = 5, dpi = dpi)
ggsave(file.path(out_dir, "QC", "PCA_post_ComBat.png"), p_pca_post,
       width = 6, height = 5, dpi = dpi)

# =========================================================================== #
# 5. MODULE SCORES ON BATCH-CORRECTED MATRIX ----
# =========================================================================== #

meta$pyrim_score   <- module_score(corrected_mat, pyrim_genes)
meta$calcium_score <- module_score(corrected_mat, calcium_genes)
meta$ecm_score     <- module_score(corrected_mat, ecm_genes)
meta$P2RY6  <- if ("P2RY6"  %in% rownames(corrected_mat)) as.numeric(corrected_mat["P2RY6",  ]) else NA_real_
meta$PDGFRB <- if ("PDGFRB" %in% rownames(corrected_mat)) as.numeric(corrected_mat["PDGFRB", ]) else NA_real_

# =========================================================================== #
# 6. BOXPLOTS — ALL 5 CONDITIONS, NORMAL AS REFERENCE ----
# =========================================================================== #

subtitle_n <- n_subtitle(meta)

plot_score_box <- function(df, y_col, y_label, title_str,
                           cmp = comparisons_main) {
  ymax        <- max(df[[y_col]], na.rm = TRUE)
  ymin        <- min(df[[y_col]], na.rm = TRUE)
  rng         <- ymax - ymin
  upper_extra <- (0.08 * length(cmp) + 0.10) * rng   # ~0.10 per bracket

  ggplot(df, aes(x = condition, y = .data[[y_col]], fill = condition)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.55, width = 0.55) +
    geom_jitter(aes(shape = cohort), width = 0.12, size = 2.2, alpha = 0.9) +
    stat_compare_means(comparisons = cmp, method = "wilcox.test",
                       label = "p.signif", hide.ns = FALSE) +
    scale_fill_manual(values = condition_colors) +
    scale_shape_manual(values = cohort_shapes) +
    scale_x_discrete(labels = setNames(condition_labels, condition_levels)) +
    scale_y_continuous(limits = c(ymin - 0.05 * rng, ymax + upper_extra)) +
    labs(x = NULL, y = y_label, title = title_str, subtitle = subtitle_n,
         shape = "Cohort") +
    theme_box() +
    theme(legend.position = "right",
          legend.title    = element_text(size = 10))
}

p_pyrim   <- plot_score_box(meta, "pyrim_score",   "Pyrimidine module score (log2)",  "Combined — Pyrimidine score")
p_p2ry6   <- plot_score_box(meta, "P2RY6",         "P2RY6 expression (log2)",         "Combined — P2RY6")
p_calcium <- plot_score_box(meta, "calcium_score", "Calcium signaling score (log2)",  "Combined — Calcium score")
p_pdgfrb  <- plot_score_box(meta, "PDGFRB",        "PDGFRB expression (log2)",        "Combined — PDGFRB")
p_ecm     <- plot_score_box(meta, "ecm_score",     "ECM / fibrotic score (log2)",     "Combined — ECM score")

ggsave(file.path(out_dir, "Fig_pyrimidine_score.png"),   p_pyrim,   width = 6, height = 5.5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_P2RY6.png"),              p_p2ry6,   width = 6, height = 5.5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_calcium_score.png"),      p_calcium, width = 6, height = 5.5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_PDGFRB.png"),             p_pdgfrb,  width = 6, height = 5.5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_ECM_score.png"),          p_ecm,     width = 6, height = 5.5, dpi = dpi)

# --- Individual pyrimidine genes ---
pyrim_found <- intersect(pyrim_genes, rownames(corrected_mat))
df_pyrim_long <- meta %>%
  select(sample_id, condition, cohort) %>%
  bind_cols(as.data.frame(t(corrected_mat[pyrim_found, ]))) %>%
  pivot_longer(cols = all_of(pyrim_found), names_to = "gene", values_to = "expr")

p_pyrim_genes <- ggplot(df_pyrim_long,
                        aes(x = condition, y = expr, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.55, width = 0.55) +
  geom_jitter(aes(shape = cohort), width = 0.12, size = 1.2, alpha = 0.85) +
  stat_compare_means(comparisons = comparisons_main, method = "wilcox.test",
                     label = "p.signif", hide.ns = FALSE) +
  scale_fill_manual(values = condition_colors) +
  scale_shape_manual(values = cohort_shapes) +
  scale_x_discrete(labels = setNames(condition_labels, condition_levels)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.55))) +
  facet_wrap(~gene, scales = "free_y", ncol = 3) +
  labs(x = NULL, y = "Expression (log2)", shape = "Cohort",
       title = "Combined — Individual pyrimidine genes", subtitle = subtitle_n) +
  theme_box() + theme(legend.position = "right")

ggsave(file.path(out_dir, "Fig_pyrimidine_genes.png"), p_pyrim_genes,
       width = 12, height = 7, dpi = dpi)

# =========================================================================== #
# 7. STATISTICAL TABLES — BH-FDR CORRECTED ----
# =========================================================================== #

make_stats_table <- function(df, y_col, label) {
  tbl <- wilcox_fdr(df, y_col, "condition", comparisons_main)
  ns  <- table(df$condition)
  tbl$module   <- label
  tbl$n_group1 <- ns[tbl$group1]
  tbl$n_group2 <- ns[tbl$group2]
  tbl[, c("module", "group1", "group2", "n_group1", "n_group2",
          "p_value", "p_adj", "signif")]
}

stats_list <- list(
  make_stats_table(meta, "pyrim_score",   "Pyrimidine score"),
  make_stats_table(meta, "P2RY6",         "P2RY6"),
  make_stats_table(meta, "calcium_score", "Calcium score"),
  make_stats_table(meta, "PDGFRB",        "PDGFRB"),
  make_stats_table(meta, "ecm_score",     "ECM score")
)

stats_df <- do.call(rbind, stats_list)
write.csv(stats_df, file.path(out_dir, "Table_group_comparisons.csv"), row.names = FALSE)
cat("\n--- Group comparisons (vs Normal) ---\n")
print(stats_df)

# =========================================================================== #
# 8. CORRELATION PANELS — POOLED (ALL 63 SAMPLES) ----
# =========================================================================== #

plot_scatter <- function(df, x_col, y_col, x_label, y_label, title_str) {
  ct <- spearman_ci(df[[x_col]], df[[y_col]])
  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]],
                 color = condition, shape = cohort)) +
    geom_point(size = 2.5, alpha = 0.9) +
    geom_smooth(method = "lm", se = FALSE, color = "black",
                linewidth = 0.8, inherit.aes = FALSE,
                aes(x = .data[[x_col]], y = .data[[y_col]])) +
    scale_color_manual(values = condition_colors) +
    scale_shape_manual(values = cohort_shapes) +
    annotate("text", x = Inf, y = Inf,
             label = ct$label, hjust = 1.05, vjust = 1.2, size = 3.8) +
    labs(x = x_label, y = y_label, title = title_str, subtitle = subtitle_n,
         color = NULL, shape = "Cohort") +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold"), legend.position = "right")
}

p_cor_pyrim_p2ry6   <- plot_scatter(meta, "pyrim_score", "P2RY6",
  "Pyrimidine score", "P2RY6 (log2)", "Combined — Pyrimidine vs P2RY6")
p_cor_pyrim_pdgfrb  <- plot_scatter(meta, "pyrim_score", "PDGFRB",
  "Pyrimidine score", "PDGFRB (log2)", "Combined — Pyrimidine vs PDGFRB")
p_cor_pyrim_ecm     <- plot_scatter(meta, "pyrim_score", "ecm_score",
  "Pyrimidine score", "ECM score", "Combined — Pyrimidine vs ECM")
p_cor_p2ry6_ecm     <- plot_scatter(meta, "P2RY6", "ecm_score",
  "P2RY6 (log2)", "ECM score", "Combined — P2RY6 vs ECM")
p_cor_p2ry6_calcium <- plot_scatter(meta, "P2RY6", "calcium_score",
  "P2RY6 (log2)", "Calcium score", "Combined — P2RY6 vs Calcium")
p_cor_pdgfrb_p2ry6  <- plot_scatter(meta, "PDGFRB", "P2RY6",
  "PDGFRB (log2)", "P2RY6 (log2)", "Combined — PDGFRB vs P2RY6")

ggsave(file.path(out_dir, "Fig_cor_Pyrimidine_P2RY6.png"),   p_cor_pyrim_p2ry6,   width = 6, height = 4.5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_cor_Pyrimidine_PDGFRB.png"),  p_cor_pyrim_pdgfrb,  width = 6, height = 4.5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_cor_Pyrimidine_ECM.png"),     p_cor_pyrim_ecm,     width = 6, height = 4.5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_cor_P2RY6_ECM.png"),          p_cor_p2ry6_ecm,     width = 6, height = 4.5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_cor_P2RY6_Calcium.png"),      p_cor_p2ry6_calcium, width = 6, height = 4.5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_cor_PDGFRB_P2RY6.png"),       p_cor_pdgfrb_p2ry6,  width = 6, height = 4.5, dpi = dpi)

# =========================================================================== #
# 9. WITHIN-CONDITION CORRELATIONS (supplementary) ----
# =========================================================================== #

cor_pairs <- list(
  list(x = "pyrim_score", y = "P2RY6"),
  list(x = "pyrim_score", y = "ecm_score"),
  list(x = "P2RY6",       y = "ecm_score"),
  list(x = "P2RY6",       y = "calcium_score")
)

within_cor <- lapply(cor_pairs, function(pair) {
  meta %>%
    group_by(condition) %>%
    group_modify(~{
      if (nrow(.x) < 4) {
        return(data.frame(x_var = pair$x, y_var = pair$y,
                          rho = NA, ci_lo = NA, ci_hi = NA, p = NA, n = nrow(.x)))
      }
      ct <- spearman_ci(.x[[pair$x]], .x[[pair$y]])
      data.frame(x_var = pair$x, y_var = pair$y,
                 rho = ct$rho, ci_lo = ct$ci_lo, ci_hi = ct$ci_hi,
                 p = ct$p, n = ct$n)
    }) %>% ungroup()
}) |> do.call(what = rbind)

within_cor$p_adj <- p.adjust(within_cor$p, method = "BH")
write.csv(within_cor, file.path(out_dir, "Table_within_condition_correlations.csv"),
          row.names = FALSE)
cat("\n--- Within-condition correlations ---\n")
print(within_cor)

# =========================================================================== #
# 10. POOLED CORRELATION TABLE ----
# =========================================================================== #

all_pairs <- expand.grid(
  x = c("pyrim_score", "P2RY6", "calcium_score", "PDGFRB", "ecm_score"),
  y = c("pyrim_score", "P2RY6", "calcium_score", "PDGFRB", "ecm_score"),
  stringsAsFactors = FALSE
) |> filter(x < y)

pooled_cor <- lapply(seq_len(nrow(all_pairs)), function(i) {
  ct <- spearman_ci(meta[[all_pairs$x[i]]], meta[[all_pairs$y[i]]])
  data.frame(x_var = all_pairs$x[i], y_var = all_pairs$y[i],
             rho = ct$rho, ci_lo = ct$ci_lo, ci_hi = ct$ci_hi,
             p = ct$p, n = ct$n)
}) |> do.call(what = rbind)

pooled_cor$p_adj <- p.adjust(pooled_cor$p, method = "BH")
write.csv(pooled_cor, file.path(out_dir, "Table_pooled_correlations.csv"),
          row.names = FALSE)
cat("\n--- Pooled correlations ---\n")
print(pooled_cor)

# =========================================================================== #
# 11. SAVE BATCH-CORRECTED MATRIX ----
# =========================================================================== #

saveRDS(
  list(corrected_mat = corrected_mat, meta = meta),
  file = file.path(out_dir, "S40_corrected.rds")
)

cat("\nS40 complete. Outputs saved to:", out_dir, "\n")
cat("N =", nrow(meta), "| Conditions:",
    paste(names(table(meta$condition)), table(meta$condition), sep = "=", collapse = ", "), "\n")
