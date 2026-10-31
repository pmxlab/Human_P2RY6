# =========================================================================== #
# S20_FanDKD_BulkRNA.R
# Fan et al. 2020 (GSE142025) — Bulk RNA-seq DKD analysis
# Tier 1, Cohort A: Normal / Early DKD / Advanced DKD
#
# DATA — Source: GSE142025_RAW.tar (NCBI GEO, acc. GSE142025)
#
# Normalization (from !Sample_data_processing, GSE142025_series_matrix.txt):
#   "The reads ... were log2 transformed and quantile-normalized at an equal
#    global median value in order to compare transcription level across samples."
#
# Values in the uploaded files are therefore log2(quantile-normalized counts).
# They are NOT raw counts, log2(counts+1), or log2 CPM.
#
# Verified by inspection of all 36 sample files:
#   - All samples contain exactly 17,184 genes, reflecting a global low-count
#     filter (transcripts with < 100 reads across all samples) applied before QN.
#   - The maximum value (19.2853815603874) is identical across all 36 samples:
#     it is the rank-mean for rank 17,184 and is shared regardless of which gene
#     occupies that rank (e.g. GPX3 in Normal/Early DN, TTN in Advanced DN).
#     This is definitive evidence of quantile normalization.
#
# No further normalization is applied in this script.
# =========================================================================== #

rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(edgeR)

# =========================================================================== #
# SETTINGS ----
# =========================================================================== #

raw_dir    <- "../../Fan_DKD_2020/GSE142025_RAW"
out_dir    <- "../04-Outputs/S20_FanDKD_BulkRNA"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

pyrim_genes   <- c("CAD", "DHODH", "UMPS", "UCK1", "CDA", "DCTPP1")
calcium_genes <- c("GNAQ", "ITPR2", "ITPR3", "ORAI1")
ecm_genes     <- c("COL1A1", "COL3A1", "FN1", "LUM", "DCN")

condition_levels  <- c("Normal", "Early_DKD", "Advanced_DKD")
condition_labels  <- c("Normal", "Early DKD", "Advanced DKD")
condition_colors  <- c("Normal" = "#4393c3", "Early_DKD" = "#f4a582", "Advanced_DKD" = "#d6604d")

comparisons_all   <- list(c("Normal", "Early_DKD"), c("Normal", "Advanced_DKD"))
comparisons_full  <- list(c("Normal", "Early_DKD"), c("Normal", "Advanced_DKD"),
                          c("Early_DKD", "Advanced_DKD"))

dpi <- 300

# =========================================================================== #
# HELPER FUNCTIONS ----
# =========================================================================== #

# Spearman correlation with 95% CI (Fisher z-transform) and formatted label
spearman_ci <- function(x, y) {
  ct  <- cor.test(x, y, method = "spearman", exact = FALSE)
  rho <- unname(ct$estimate)
  n   <- sum(!is.na(x) & !is.na(y))
  # Fisher z CI
  z    <- atanh(rho)
  se   <- 1 / sqrt(n - 3)
  ci   <- tanh(c(z - 1.96 * se, z + 1.96 * se))
  list(
    rho = rho, ci_lo = ci[1], ci_hi = ci[2],
    p   = ct$p.value, n = n,
    label = sprintf("rho = %.2f [%.2f, %.2f]\np = %s",
                    rho, ci[1], ci[2], format_p(ct$p.value))
  )
}

# Format p-value for display
format_p <- function(p) {
  if (p < 0.001) return("< 0.001")
  sprintf("%.3f", p)
}

# Pairwise Wilcoxon with BH-FDR correction; returns tidy data.frame
wilcox_fdr <- function(df, value_col, group_col, comparisons) {
  pvals <- sapply(comparisons, function(cmp) {
    a <- df[[value_col]][df[[group_col]] == cmp[1]]
    b <- df[[value_col]][df[[group_col]] == cmp[2]]
    wilcox.test(a, b, exact = FALSE)$p.value
  })
  padj <- p.adjust(pvals, method = "BH")
  data.frame(
    group1  = sapply(comparisons, `[`, 1),
    group2  = sapply(comparisons, `[`, 2),
    p_value = pvals,
    p_adj   = padj,
    signif  = ifelse(padj < 0.001, "***",
               ifelse(padj < 0.01,  "**",
                ifelse(padj < 0.05,  "*", "ns")))
  )
}

# Module score: column means over found genes
module_score <- function(mat, genes) {
  found <- intersect(genes, rownames(mat))
  if (length(found) == 0) return(rep(NA_real_, ncol(mat)))
  colMeans(mat[found, , drop = FALSE], na.rm = TRUE)
}

# Standard boxplot theme
theme_box <- function() {
  theme_bw(base_size = 13) +
    theme(
      legend.position   = "none",
      axis.text.x       = element_text(angle = 20, hjust = 1),
      plot.title        = element_text(face = "bold"),
      strip.text        = element_text(face = "bold")
    )
}

# n label for subtitle: "Normal (n=12), Early DKD (n=10), Advanced DKD (n=14)"
n_subtitle <- function(meta) {
  ns <- table(meta$condition)
  lbs <- condition_labels
  names(lbs) <- condition_levels
  paste(
    mapply(function(lvl, lbl) sprintf("%s (n=%d)", lbl, ns[lvl]),
           condition_levels, lbs[condition_levels]),
    collapse = ", "
  )
}

# =========================================================================== #
# 1. DATA LOADING ----
# =========================================================================== #
# NOTE: GSE142025 files contain log2(quantile-normalized counts) — see header.
# No further normalization is applied; values are used directly.

all_files <- list.files(raw_dir, pattern = "\\.txt\\.gz$", full.names = TRUE)
stopifnot(length(all_files) > 0)

expr_list <- lapply(all_files, function(f) {
  x <- read.delim(gzfile(f), header = TRUE, sep = "\t", check.names = FALSE)
  colnames(x)[2] <- sub("\\.txt\\.gz$", "", basename(f))
  x
})

expr_merged <- Reduce(function(a, b) merge(a, b, by = "Symbol", all = TRUE), expr_list)
expr_mat <- as.matrix(expr_merged[, -1])
rownames(expr_mat) <- expr_merged$Symbol
mode(expr_mat) <- "numeric"

# Metadata
sample_names <- colnames(expr_mat)
meta <- data.frame(
  sample   = sample_names,
  short_id = sub("^GSM[0-9]+_", "", sample_names),
  stringsAsFactors = FALSE
)
meta$group <- ifelse(grepl("^N", meta$short_id), "Normal",
               ifelse(grepl("^A", meta$short_id), "DKD_A",
                ifelse(grepl("^B", meta$short_id), "DKD_B", NA)))
meta$condition <- factor(
  ifelse(meta$group == "Normal", "Normal",
   ifelse(meta$group == "DKD_A", "Advanced_DKD",
    ifelse(meta$group == "DKD_B", "Early_DKD", NA))),
  levels = condition_levels
)

# Align matrix and remove genes with any missing values
expr_mat <- expr_mat[, meta$sample]
expr_mat <- expr_mat[complete.cases(expr_mat), ]

cat("Samples:", nrow(meta), "\n")
cat("Conditions:", paste(table(meta$condition), collapse = " / "), "\n")
cat("Genes (complete cases):", nrow(expr_mat), "\n")
cat(n_subtitle(meta), "\n")
key_genes   <- unique(c(pyrim_genes, calcium_genes, ecm_genes, "P2RY6", "PDGFRB"))
missing_keys <- key_genes[!key_genes %in% rownames(expr_mat)]
if (length(missing_keys)) cat("WARNING - key genes absent from matrix:", paste(missing_keys, collapse = ", "), "\n")

# =========================================================================== #
# 2. MODULE SCORES ----
# =========================================================================== #

meta$pyrim_score   <- module_score(expr_mat, pyrim_genes)
meta$calcium_score <- module_score(expr_mat, calcium_genes)
meta$ecm_score     <- module_score(expr_mat, ecm_genes)
meta$P2RY6         <- if ("P2RY6"  %in% rownames(expr_mat)) as.numeric(expr_mat["P2RY6",  ]) else NA_real_
meta$PDGFRB        <- if ("PDGFRB" %in% rownames(expr_mat)) as.numeric(expr_mat["PDGFRB", ]) else NA_real_

# =========================================================================== #
# 3. BOXPLOTS — GROUP COMPARISONS ----
# =========================================================================== #

subtitle_n <- n_subtitle(meta)

# --- helper: generic score boxplot ---
plot_score_box <- function(df, y_col, y_label, title_str, cmp = comparisons_all) {
  pvals <- wilcox_fdr(df, y_col, "condition", cmp)
  ymax  <- max(df[[y_col]], na.rm = TRUE)
  ymin  <- min(df[[y_col]], na.rm = TRUE)
  rng   <- ymax - ymin

  ggplot(df, aes(x = condition, y = .data[[y_col]], fill = condition)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.55, width = 0.55) +
    geom_jitter(width = 0.12, size = 2.2, alpha = 0.9) +
    stat_compare_means(
      comparisons = cmp,
      method      = "wilcox.test",
      label       = "p.signif",
      hide.ns     = FALSE
    ) +
    scale_fill_manual(values = condition_colors) +
    scale_x_discrete(labels = setNames(condition_labels, condition_levels)) +
    scale_y_continuous(limits = c(ymin - 0.05 * rng, ymax + 0.25 * rng)) +
    labs(x = NULL, y = y_label, title = title_str, subtitle = subtitle_n) +
    theme_box()
}

p_pyrim   <- plot_score_box(meta, "pyrim_score",   "Pyrimidine module score (log2, QN)",
                            "Fan DKD — Pyrimidine score")
p_p2ry6   <- plot_score_box(meta, "P2RY6",         "P2RY6 expression (log2, QN)",
                            "Fan DKD — P2RY6")
p_calcium  <- plot_score_box(meta, "calcium_score", "Calcium signaling score (log2, QN)",
                             "Fan DKD — Calcium score")
p_pdgfrb   <- plot_score_box(meta, "PDGFRB",        "PDGFRB expression (log2, QN)",
                             "Fan DKD — PDGFRB")
p_ecm      <- plot_score_box(meta, "ecm_score",     "ECM / fibrotic score (log2, QN)",
                             "Fan DKD — ECM score")

ggsave(file.path(out_dir, "Fig_pyrimidine_score.png"),   p_pyrim,   width = 4, height = 5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_P2RY6.png"),              p_p2ry6,   width = 4, height = 5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_calcium_score.png"),      p_calcium,  width = 4, height = 5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_PDGFRB.png"),             p_pdgfrb,   width = 4, height = 5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_ECM_score.png"),          p_ecm,      width = 4, height = 5, dpi = dpi)

# --- Individual pyrimidine genes ---
df_pyrim_long <- meta %>%
  select(sample, condition) %>%
  bind_cols(as.data.frame(t(expr_mat[intersect(pyrim_genes, rownames(expr_mat)), ]))) %>%
  pivot_longer(cols = all_of(intersect(pyrim_genes, rownames(expr_mat))),
               names_to = "gene", values_to = "expr")

p_pyrim_genes <- ggplot(df_pyrim_long, aes(x = condition, y = expr, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.55, width = 0.55) +
  geom_jitter(width = 0.12, size = 1.2, alpha = 0.85) +
  stat_compare_means(
    comparisons = comparisons_all,
    method = "wilcox.test", label = "p.signif", hide.ns = FALSE
  ) +
  scale_fill_manual(values = condition_colors) +
  scale_x_discrete(labels = setNames(condition_labels, condition_levels)) +
  facet_wrap(~gene, scales = "free_y", ncol = 3) +
  labs(x = NULL, y = "Expression (log2, QN)",
       title = "Fan DKD — Individual pyrimidine genes", subtitle = subtitle_n) +
  theme_box()

ggsave(file.path(out_dir, "Fig_pyrimidine_genes.png"), p_pyrim_genes,
       width = 10, height = 7, dpi = dpi)

# =========================================================================== #
# 4. STATISTICAL TABLES — ALL PAIRWISE COMPARISONS WITH BH-FDR ----
# =========================================================================== #

make_stats_table <- function(df, y_col, label) {
  tbl <- wilcox_fdr(df, y_col, "condition", comparisons_full)
  ns  <- table(df$condition)
  tbl$module <- label
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
print(stats_df)

# =========================================================================== #
# 5. CORRELATION PANELS — POOLED AND WITHIN-CONDITION ----
# =========================================================================== #

# --- helper: scatter with CI annotation ---
plot_scatter <- function(df, x_col, y_col, x_label, y_label, title_str, color_by = "condition") {
  ct  <- spearman_ci(df[[x_col]], df[[y_col]])
  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]], color = .data[[color_by]])) +
    geom_point(size = 2.8, alpha = 0.9) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8) +
    scale_color_manual(values = condition_colors) +
    annotate("text", x = Inf, y = Inf,
             label = ct$label, hjust = 1.05, vjust = 1.2, size = 3.8) +
    labs(x = x_label, y = y_label, title = title_str, subtitle = subtitle_n,
         color = NULL) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "right")
}

p_cor_pyrim_p2ry6   <- plot_scatter(meta, "pyrim_score", "P2RY6",
  "Pyrimidine score", "P2RY6 (log2, QN)", "Fan DKD — Pyrimidine vs P2RY6")
p_cor_pyrim_pdgfrb  <- plot_scatter(meta, "pyrim_score", "PDGFRB",
  "Pyrimidine score", "PDGFRB (log2, QN)", "Fan DKD — Pyrimidine vs PDGFRB")
p_cor_pyrim_ecm     <- plot_scatter(meta, "pyrim_score", "ecm_score",
  "Pyrimidine score", "ECM score", "Fan DKD — Pyrimidine vs ECM")
p_cor_p2ry6_ecm     <- plot_scatter(meta, "P2RY6", "ecm_score",
  "P2RY6 (log2, QN)", "ECM score", "Fan DKD — P2RY6 vs ECM")
p_cor_p2ry6_calcium <- plot_scatter(meta, "P2RY6", "calcium_score",
  "P2RY6 (log2, QN)", "Calcium score", "Fan DKD — P2RY6 vs Calcium")
p_cor_pdgfrb_p2ry6  <- plot_scatter(meta, "PDGFRB", "P2RY6",
  "PDGFRB (log2, QN)", "P2RY6 (log2, QN)", "Fan DKD — PDGFRB vs P2RY6")

ggsave(file.path(out_dir, "Fig_cor_Pyrimidine_P2RY6.png"),   p_cor_pyrim_p2ry6,   width = 5.5, height = 4.5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_cor_Pyrimidine_PDGFRB.png"),  p_cor_pyrim_pdgfrb,  width = 5.5, height = 4.5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_cor_Pyrimidine_ECM.png"),     p_cor_pyrim_ecm,     width = 5.5, height = 4.5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_cor_P2RY6_ECM.png"),          p_cor_p2ry6_ecm,     width = 5.5, height = 4.5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_cor_P2RY6_Calcium.png"),      p_cor_p2ry6_calcium, width = 5.5, height = 4.5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_cor_PDGFRB_P2RY6.png"),       p_cor_pdgfrb_p2ry6,  width = 5.5, height = 4.5, dpi = dpi)

# --- Within-condition correlations (supplementary) ---
cor_pairs <- list(
  list(x = "pyrim_score", y = "P2RY6",         label_x = "Pyrimidine score", label_y = "P2RY6"),
  list(x = "pyrim_score", y = "ecm_score",      label_x = "Pyrimidine score", label_y = "ECM score"),
  list(x = "P2RY6",       y = "ecm_score",      label_x = "P2RY6",            label_y = "ECM score"),
  list(x = "P2RY6",       y = "calcium_score",  label_x = "P2RY6",            label_y = "Calcium score")
)

within_cor <- lapply(cor_pairs, function(pair) {
  meta %>%
    group_by(condition) %>%
    group_modify(~{
      ct <- spearman_ci(.x[[pair$x]], .x[[pair$y]])
      data.frame(
        x_var  = pair$x,   y_var = pair$y,
        rho    = ct$rho,   ci_lo = ct$ci_lo, ci_hi = ct$ci_hi,
        p      = ct$p,     n     = ct$n
      )
    }) %>%
    ungroup()
}) |> do.call(what = rbind)

within_cor$p_adj <- p.adjust(within_cor$p, method = "BH")
write.csv(within_cor, file.path(out_dir, "Table_within_condition_correlations.csv"), row.names = FALSE)
print(within_cor)

# =========================================================================== #
# 6. CORRELATION TABLE — ALL PAIRS, POOLED ----
# =========================================================================== #

all_pairs <- expand.grid(
  x = c("pyrim_score", "P2RY6", "calcium_score", "PDGFRB", "ecm_score"),
  y = c("pyrim_score", "P2RY6", "calcium_score", "PDGFRB", "ecm_score"),
  stringsAsFactors = FALSE
) |>
  filter(x < y)  # upper triangle only

pooled_cor <- lapply(seq_len(nrow(all_pairs)), function(i) {
  ct <- spearman_ci(meta[[all_pairs$x[i]]], meta[[all_pairs$y[i]]])
  data.frame(x_var = all_pairs$x[i], y_var = all_pairs$y[i],
             rho = ct$rho, ci_lo = ct$ci_lo, ci_hi = ct$ci_hi,
             p = ct$p, n = ct$n)
}) |> do.call(what = rbind)

pooled_cor$p_adj <- p.adjust(pooled_cor$p, method = "BH")
write.csv(pooled_cor, file.path(out_dir, "Table_pooled_correlations.csv"), row.names = FALSE)
print(pooled_cor)

# =========================================================================== #
# 7. SAVE NORMALISED EXPRESSION (used by S70, S75, S80 and later scripts) ----
# =========================================================================== #

saveRDS(
  list(expr_mat = expr_mat, meta = meta),
  file = file.path(out_dir, "S20_normalised.rds")
)

cat("\nS20 complete. Outputs saved to:", out_dir, "\n")
