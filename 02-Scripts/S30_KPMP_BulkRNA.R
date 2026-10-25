# =========================================================================== #
# S30_KPMP_BulkRNA.R
# KPMP BulkRNA — Bulk RNA-seq analysis
# Tier 1, Cohort B: AKI / CKD (no healthy controls)
#
# Condition assignment: derived from 0_ATLAS_PATIENT_FILE_METADATA.csv
#   enrollment_category          -> AKI / CKD
#   primary_adjudicated_category -> DKD, Acute Tubular Injury,
#                                   Hypertensive Kidney Disease,
#                                   Acute Interstitial Nephritis,
#                                   Other, Cannot be determined
#
# Deduplication: for participants with multiple count files, the file with the
#   highest total library size (sum of all mapped read counts) is retained.
#   This criterion should be stated in the Methods section.
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

raw_dir  <- "../../KPMP_BulkRNA"
meta_csv <- file.path(raw_dir, "0_ATLAS_PATIENT_FILE_METADATA.csv")
out_dir  <- "../04-Outputs/S30_KPMP_BulkRNA"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

pyrim_genes   <- c("CAD", "DHODH", "UMPS", "UCK1", "CDA", "DCTPP1")
calcium_genes <- c("GNAQ", "ITPR2", "ITPR3", "ORAI1")
ecm_genes     <- c("COL1A1", "COL3A1", "FN1", "LUM", "DCN")

# Primary comparison: enrollment category
enroll_levels <- c("AKI", "CKD")
enroll_colors <- c("AKI" = "#e67e22", "CKD" = "#3498db")

# Adjudicated sub-categories (for within-group figures)
adj_colors <- c(
  "Acute Tubular Injury"         = "#f39c12",
  "Acute Interstitial Nephritis" = "#e74c3c",
  "Diabetic Kidney Disease"      = "#2980b9",
  "Hypertensive Kidney Disease"  = "#8e44ad",
  "Other"                        = "#27ae60",
  "Cannot be determined"         = "#95a5a6"
)

comparisons_enroll <- list(c("AKI", "CKD"))

dpi <- 300

# =========================================================================== #
# HELPER FUNCTIONS ----
# =========================================================================== #

spearman_ci <- function(x, y) {
  ct  <- cor.test(x, y, method = "spearman", exact = FALSE)
  rho <- unname(ct$estimate)
  n   <- sum(!is.na(x) & !is.na(y))
  z   <- atanh(rho)
  se  <- 1 / sqrt(n - 3)
  ci  <- tanh(c(z - 1.96 * se, z + 1.96 * se))
  list(
    rho = rho, ci_lo = ci[1], ci_hi = ci[2],
    p   = ct$p.value, n = n,
    label = sprintf("rho = %.2f [%.2f, %.2f]\np = %s",
                    rho, ci[1], ci[2], format_p(ct$p.value))
  )
}

format_p <- function(p) {
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
    p_value = pvals,
    p_adj   = padj,
    signif  = ifelse(is.na(padj), "n/a",
               ifelse(padj < 0.001, "***",
                ifelse(padj < 0.01,  "**",
                 ifelse(padj < 0.05,  "*", "ns"))))
  )
}

module_score <- function(mat, genes) {
  found <- intersect(genes, rownames(mat))
  if (length(found) == 0) return(rep(NA_real_, ncol(mat)))
  colMeans(mat[found, , drop = FALSE], na.rm = TRUE)
}

theme_box <- function() {
  theme_bw(base_size = 13) +
    theme(
      legend.position = "none",
      axis.text.x     = element_text(angle = 20, hjust = 1),
      plot.title      = element_text(face = "bold"),
      strip.text      = element_text(face = "bold")
    )
}

# Compact n label, e.g. "AKI (n=12), CKD (n=15)"
n_subtitle_enroll <- function(meta) {
  ns <- table(meta$enrollment_category)
  paste(sapply(enroll_levels, function(lv) sprintf("%s (n=%d)", lv, ns[lv])),
        collapse = ", ")
}

# =========================================================================== #
# 1. BUILD FILE MANIFEST FROM METADATA CSV ----
# =========================================================================== #

atlas_meta <- read.csv(meta_csv, stringsAsFactors = FALSE, check.names = FALSE)
names(atlas_meta) <- trimws(names(atlas_meta))

# Keep only count files (not metadata XLSX, not raw-reads)
counts_meta <- atlas_meta %>%
  filter(grepl("counts\\.csv$|Counts\\.csv$", `File Name`, ignore.case = TRUE),
         !grepl("raw-reads", `File Name`, ignore.case = TRUE))

cat("Count file entries in metadata:", nrow(counts_meta), "\n")

# Resolve full file path
counts_meta$file_path <- file.path(raw_dir, counts_meta$`File Name`)

# Verify files exist
missing_files <- counts_meta$file_path[!file.exists(counts_meta$file_path)]
if (length(missing_files) > 0) {
  cat("WARNING - files listed in metadata but not found on disk:\n")
  cat(paste(" ", missing_files, collapse = "\n"), "\n")
}
counts_meta <- counts_meta[file.exists(counts_meta$file_path), ]

# =========================================================================== #
# 2. DEDUPLICATE: one file per participant (highest total library size) ----
# =========================================================================== #

read_counts_file <- function(f) {
  d <- read.csv(f, header = TRUE, check.names = FALSE)
  if (ncol(d) < 2) {
    # Try tab-separated fallback
    d <- read.table(f, header = TRUE, sep = "\t", check.names = FALSE)
  }
  if (ncol(d) < 2) stop("Cannot parse count file (< 2 columns): ", basename(f))
  d
}

counts_meta$total_counts <- sapply(counts_meta$file_path, function(f) {
  d <- read_counts_file(f)
  sum(d[, 2], na.rm = TRUE)
})

# Keep the file with the highest total counts per participant
counts_meta <- counts_meta %>%
  group_by(`Participant ID`) %>%
  slice_max(order_by = total_counts, n = 1, with_ties = FALSE) %>%
  ungroup()

# QC filter 1: exclude samples with total library size < 1M reads (failed runs)
min_library_size <- 1e6
excluded_qc <- counts_meta %>% filter(total_counts < min_library_size)
if (nrow(excluded_qc) > 0) {
  cat("QC exclusion (library size <", min_library_size, "):\n")
  print(excluded_qc[, c("Participant ID", "Enrollment Category",
                         "Primary Adjudicated Category", "total_counts")])
}
counts_meta <- counts_meta %>% filter(total_counts >= min_library_size)

# QC filter 2: exclude samples without a confirmed primary adjudicated category
excluded_adj <- counts_meta %>%
  filter(`Primary Adjudicated Category` == "Cannot be determined")
if (nrow(excluded_adj) > 0) {
  cat("Adjudicated category exclusion ('Cannot be determined'):\n")
  print(excluded_adj[, c("Participant ID", "Enrollment Category",
                          "Primary Adjudicated Category", "total_counts")])
}
counts_meta <- counts_meta %>%
  filter(`Primary Adjudicated Category` != "Cannot be determined")

cat("Unique participants after deduplication:", nrow(counts_meta), "\n")
cat("Enrollment categories:\n")
print(table(counts_meta$`Enrollment Category`))
cat("Primary adjudicated categories:\n")
print(table(counts_meta$`Primary Adjudicated Category`))

# =========================================================================== #
# 3. LOAD AND MERGE COUNT MATRIX ----
# =========================================================================== #

# Build unified gene list (intersect to be safe)
gene_sets <- lapply(counts_meta$file_path, function(f) {
  d <- read_counts_file(f)
  d[[1]]
})
common_genes <- Reduce(intersect, gene_sets)
cat("Genes in common across all samples:", length(common_genes), "\n")

# Assemble matrix: genes x samples
expr_raw <- do.call(cbind, lapply(seq_len(nrow(counts_meta)), function(i) {
  d    <- read_counts_file(counts_meta$file_path[i])
  idx  <- match(common_genes, d[[1]])
  as.integer(d[[2]][idx])
}))
rownames(expr_raw) <- common_genes
colnames(expr_raw) <- counts_meta$`Participant ID`

# =========================================================================== #
# 4. METADATA TABLE ----
# =========================================================================== #

meta <- data.frame(
  sample_id                    = counts_meta$`Participant ID`,
  enrollment_category          = counts_meta$`Enrollment Category`,
  primary_adjudicated_category = counts_meta$`Primary Adjudicated Category`,
  total_counts                 = counts_meta$total_counts,
  stringsAsFactors             = FALSE
)
meta$enrollment_category <- factor(meta$enrollment_category, levels = enroll_levels)

# Align matrix column order to meta row order
expr_raw <- expr_raw[, meta$sample_id]

# =========================================================================== #
# 5. TMM NORMALISATION ----
# =========================================================================== #

key_genes <- unique(c(pyrim_genes, calcium_genes, ecm_genes, "P2RY6", "PDGFRB"))
dge  <- DGEList(counts = expr_raw)
keep <- filterByExpr(dge, group = meta$enrollment_category) | (rownames(dge) %in% key_genes)
dge  <- dge[keep, ]
dge  <- calcNormFactors(dge, method = "TMM")
expr_mat <- cpm(dge, log = TRUE, prior.count = 1)   # log2 CPM

cat("\nSamples:", nrow(meta), "\n")
cat("Genes after filtering:", nrow(expr_mat), "\n")
cat(n_subtitle_enroll(meta), "\n")
missing_keys <- key_genes[!key_genes %in% rownames(expr_mat)]
if (length(missing_keys)) {
  cat("WARNING - key genes absent from matrix:", paste(missing_keys, collapse = ", "), "\n")
}

# =========================================================================== #
# 6. MODULE SCORES ----
# =========================================================================== #

meta$pyrim_score   <- module_score(expr_mat, pyrim_genes)
meta$calcium_score <- module_score(expr_mat, calcium_genes)
meta$ecm_score     <- module_score(expr_mat, ecm_genes)
meta$P2RY6         <- if ("P2RY6"  %in% rownames(expr_mat)) as.numeric(expr_mat["P2RY6",  ]) else NA_real_
meta$PDGFRB        <- if ("PDGFRB" %in% rownames(expr_mat)) as.numeric(expr_mat["PDGFRB", ]) else NA_real_

# =========================================================================== #
# 7. BOXPLOTS - AKI vs CKD (enrollment category) ----
# =========================================================================== #

subtitle_n <- n_subtitle_enroll(meta)

plot_score_box <- function(df, y_col, y_label, title_str,
                           group_col = "enrollment_category",
                           colors    = enroll_colors,
                           cmp       = comparisons_enroll) {
  pvals <- wilcox_fdr(df, y_col, group_col, cmp)
  ymax  <- max(df[[y_col]], na.rm = TRUE)
  ymin  <- min(df[[y_col]], na.rm = TRUE)
  rng   <- ymax - ymin

  ggplot(df, aes(x = .data[[group_col]], y = .data[[y_col]], fill = .data[[group_col]])) +
    geom_boxplot(outlier.shape = NA, alpha = 0.55, width = 0.55) +
    geom_jitter(width = 0.12, size = 2.2, alpha = 0.9) +
    stat_compare_means(
      comparisons = cmp,
      method      = "wilcox.test",
      label       = "p.signif",
      hide.ns     = FALSE
    ) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(limits = c(ymin - 0.05 * rng, ymax + 0.25 * rng)) +
    labs(x = NULL, y = y_label, title = title_str, subtitle = subtitle_n) +
    theme_box()
}

p_pyrim   <- plot_score_box(meta, "pyrim_score",   "Pyrimidine module score (log2 CPM)", "KPMP - Pyrimidine score")
p_p2ry6   <- plot_score_box(meta, "P2RY6",         "P2RY6 expression (log2 CPM)",        "KPMP - P2RY6")
p_calcium <- plot_score_box(meta, "calcium_score", "Calcium signaling score (log2 CPM)", "KPMP - Calcium score")
p_pdgfrb  <- plot_score_box(meta, "PDGFRB",        "PDGFRB expression (log2 CPM)",       "KPMP - PDGFRB")
p_ecm     <- plot_score_box(meta, "ecm_score",     "ECM / fibrotic score (log2 CPM)",    "KPMP - ECM score")

ggsave(file.path(out_dir, "Fig_pyrimidine_score.png"),   p_pyrim,   width = 4, height = 5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_P2RY6.png"),              p_p2ry6,   width = 4, height = 5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_calcium_score.png"),      p_calcium, width = 4, height = 5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_PDGFRB.png"),             p_pdgfrb,  width = 4, height = 5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_ECM_score.png"),          p_ecm,     width = 4, height = 5, dpi = dpi)

# --- Individual pyrimidine genes ---
df_pyrim_long <- meta %>%
  select(sample_id, enrollment_category) %>%
  bind_cols(as.data.frame(t(expr_mat[intersect(pyrim_genes, rownames(expr_mat)), ]))) %>%
  pivot_longer(cols = all_of(intersect(pyrim_genes, rownames(expr_mat))),
               names_to = "gene", values_to = "expr")

p_pyrim_genes <- ggplot(df_pyrim_long,
                        aes(x = enrollment_category, y = expr, fill = enrollment_category)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.55, width = 0.55) +
  geom_jitter(width = 0.12, size = 1.2, alpha = 0.85) +
  stat_compare_means(
    comparisons = comparisons_enroll,
    method = "wilcox.test", label = "p.signif", hide.ns = FALSE
  ) +
  scale_fill_manual(values = enroll_colors) +
  facet_wrap(~gene, scales = "free_y", ncol = 3) +
  labs(x = NULL, y = "Expression (log2 CPM)",
       title = "KPMP - Individual pyrimidine genes", subtitle = subtitle_n) +
  theme_box()

ggsave(file.path(out_dir, "Fig_pyrimidine_genes.png"), p_pyrim_genes,
       width = 10, height = 7, dpi = dpi)

# =========================================================================== #
# 8. BOXPLOTS - BY PRIMARY ADJUDICATED CATEGORY (exploratory) ----
# =========================================================================== #

# Only run if at least 2 adjudicated categories have n >= 3

adj_counts <- table(meta$primary_adjudicated_category)
adj_keep   <- names(adj_counts)[adj_counts >= 3]

if (length(adj_keep) >= 2) {
  meta_adj <- meta %>%
    filter(primary_adjudicated_category %in% adj_keep) %>%
    mutate(adj_cat = factor(primary_adjudicated_category, levels = adj_keep))

  adj_cmp <- combn(adj_keep, 2, simplify = FALSE)
  adj_n   <- paste(sapply(adj_keep, function(a) sprintf("%s (n=%d)", a, adj_counts[a])),
                   collapse = ", ")
  adj_col_used <- setNames(
    sapply(adj_keep, function(a) if (a %in% names(adj_colors)) adj_colors[a] else "#7f8c8d"),
    adj_keep
  )

  plot_adj <- function(y_col, y_label, title_str) {
    ymax <- max(meta_adj[[y_col]], na.rm = TRUE)
    ymin <- min(meta_adj[[y_col]], na.rm = TRUE)
    rng  <- ymax - ymin
    ggplot(meta_adj, aes(x = adj_cat, y = .data[[y_col]], fill = adj_cat)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.55, width = 0.55) +
      geom_jitter(width = 0.12, size = 2.2, alpha = 0.9) +
      stat_compare_means(
        comparisons = adj_cmp, method = "wilcox.test",
        label = "p.signif", hide.ns = FALSE
      ) +
      scale_fill_manual(values = adj_col_used) +
      scale_y_continuous(limits = c(ymin - 0.05 * rng, ymax + 0.35 * rng)) +
      labs(x = NULL, y = y_label, title = title_str, subtitle = adj_n) +
      theme_box() +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
  }

  w <- max(3 + length(adj_keep) * 0.8, 6)
  ggsave(file.path(out_dir, "Fig_adj_pyrimidine_score.png"),
         plot_adj("pyrim_score",   "Pyrimidine module score (log2 CPM)",
                  "KPMP - Pyrimidine by adjudicated category"),
         width = w, height = 5, dpi = dpi)
  ggsave(file.path(out_dir, "Fig_adj_P2RY6.png"),
         plot_adj("P2RY6",         "P2RY6 expression (log2 CPM)",
                  "KPMP - P2RY6 by adjudicated category"),
         width = w, height = 5, dpi = dpi)
  ggsave(file.path(out_dir, "Fig_adj_calcium_score.png"),
         plot_adj("calcium_score", "Calcium signaling score (log2 CPM)",
                  "KPMP - Calcium by adjudicated category"),
         width = w, height = 5, dpi = dpi)
  ggsave(file.path(out_dir, "Fig_adj_ecm_score.png"),
         plot_adj("ecm_score",     "ECM / fibrotic score (log2 CPM)",
                  "KPMP - ECM by adjudicated category"),
         width = w, height = 5, dpi = dpi)
} else {
  cat("NOTE: fewer than 2 adjudicated categories have n >= 3; skipping adjudicated-category plots.\n")
}

# =========================================================================== #
# 9. STATISTICAL TABLES ----
# =========================================================================== #

make_stats_table <- function(df, y_col, label, group_col, cmps) {
  tbl <- wilcox_fdr(df, y_col, group_col, cmps)
  ns  <- table(df[[group_col]])
  tbl$module   <- label
  tbl$n_group1 <- ns[tbl$group1]
  tbl$n_group2 <- ns[tbl$group2]
  tbl[, c("module", "group1", "group2", "n_group1", "n_group2",
          "p_value", "p_adj", "signif")]
}

stats_list <- list(
  make_stats_table(meta, "pyrim_score",   "Pyrimidine score", "enrollment_category", comparisons_enroll),
  make_stats_table(meta, "P2RY6",         "P2RY6",            "enrollment_category", comparisons_enroll),
  make_stats_table(meta, "calcium_score", "Calcium score",    "enrollment_category", comparisons_enroll),
  make_stats_table(meta, "PDGFRB",        "PDGFRB",           "enrollment_category", comparisons_enroll),
  make_stats_table(meta, "ecm_score",     "ECM score",        "enrollment_category", comparisons_enroll)
)

stats_df <- do.call(rbind, stats_list)
write.csv(stats_df, file.path(out_dir, "Table_group_comparisons.csv"), row.names = FALSE)
cat("\n--- Group comparisons (AKI vs CKD) ---\n")
print(stats_df)

# =========================================================================== #
# 10. CORRELATION PANELS - POOLED ----
# =========================================================================== #

plot_scatter <- function(df, x_col, y_col, x_label, y_label, title_str,
                         color_by = "enrollment_category",
                         colors   = enroll_colors) {
  ct <- spearman_ci(df[[x_col]], df[[y_col]])
  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]], color = .data[[color_by]])) +
    geom_point(size = 2.8, alpha = 0.9) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8) +
    scale_color_manual(values = colors) +
    annotate("text", x = Inf, y = Inf,
             label = ct$label, hjust = 1.05, vjust = 1.2, size = 3.8) +
    labs(x = x_label, y = y_label, title = title_str, subtitle = subtitle_n,
         color = NULL) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold"), legend.position = "right")
}

p_cor_pyrim_p2ry6   <- plot_scatter(meta, "pyrim_score", "P2RY6",
  "Pyrimidine score", "P2RY6 (log2 CPM)", "KPMP - Pyrimidine vs P2RY6")
p_cor_pyrim_pdgfrb  <- plot_scatter(meta, "pyrim_score", "PDGFRB",
  "Pyrimidine score", "PDGFRB (log2 CPM)", "KPMP - Pyrimidine vs PDGFRB")
p_cor_pyrim_ecm     <- plot_scatter(meta, "pyrim_score", "ecm_score",
  "Pyrimidine score", "ECM score", "KPMP - Pyrimidine vs ECM")
p_cor_p2ry6_ecm     <- plot_scatter(meta, "P2RY6", "ecm_score",
  "P2RY6 (log2 CPM)", "ECM score", "KPMP - P2RY6 vs ECM")
p_cor_p2ry6_calcium <- plot_scatter(meta, "P2RY6", "calcium_score",
  "P2RY6 (log2 CPM)", "Calcium score", "KPMP - P2RY6 vs Calcium")
p_cor_pdgfrb_p2ry6  <- plot_scatter(meta, "PDGFRB", "P2RY6",
  "PDGFRB (log2 CPM)", "P2RY6 (log2 CPM)", "KPMP - PDGFRB vs P2RY6")

ggsave(file.path(out_dir, "Fig_cor_Pyrimidine_P2RY6.png"),   p_cor_pyrim_p2ry6,   width = 5.5, height = 4.5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_cor_Pyrimidine_PDGFRB.png"),  p_cor_pyrim_pdgfrb,  width = 5.5, height = 4.5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_cor_Pyrimidine_ECM.png"),     p_cor_pyrim_ecm,     width = 5.5, height = 4.5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_cor_P2RY6_ECM.png"),          p_cor_p2ry6_ecm,     width = 5.5, height = 4.5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_cor_P2RY6_Calcium.png"),      p_cor_p2ry6_calcium, width = 5.5, height = 4.5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_cor_PDGFRB_P2RY6.png"),       p_cor_pdgfrb_p2ry6,  width = 5.5, height = 4.5, dpi = dpi)

# =========================================================================== #
# 11. WITHIN-ENROLLMENT-CATEGORY CORRELATIONS (supplementary) ----
# =========================================================================== #

cor_pairs <- list(
  list(x = "pyrim_score", y = "P2RY6",        label_x = "Pyrimidine score", label_y = "P2RY6"),
  list(x = "pyrim_score", y = "ecm_score",     label_x = "Pyrimidine score", label_y = "ECM score"),
  list(x = "P2RY6",       y = "ecm_score",     label_x = "P2RY6",            label_y = "ECM score"),
  list(x = "P2RY6",       y = "calcium_score", label_x = "P2RY6",            label_y = "Calcium score")
)

within_cor <- lapply(cor_pairs, function(pair) {
  meta %>%
    group_by(enrollment_category) %>%
    group_modify(~{
      if (nrow(.x) < 4) {
        return(data.frame(x_var = pair$x, y_var = pair$y,
                          rho = NA, ci_lo = NA, ci_hi = NA, p = NA, n = nrow(.x)))
      }
      ct <- spearman_ci(.x[[pair$x]], .x[[pair$y]])
      data.frame(x_var = pair$x, y_var = pair$y,
                 rho = ct$rho, ci_lo = ct$ci_lo, ci_hi = ct$ci_hi, p = ct$p, n = ct$n)
    }) %>%
    ungroup()
}) |> do.call(what = rbind)

within_cor$p_adj <- p.adjust(within_cor$p, method = "BH")
write.csv(within_cor, file.path(out_dir, "Table_within_category_correlations.csv"), row.names = FALSE)
cat("\n--- Within-enrollment-category correlations ---\n")
print(within_cor)

# =========================================================================== #
# 12. POOLED CORRELATION TABLE ----
# =========================================================================== #

all_pairs <- expand.grid(
  x = c("pyrim_score", "P2RY6", "calcium_score", "PDGFRB", "ecm_score"),
  y = c("pyrim_score", "P2RY6", "calcium_score", "PDGFRB", "ecm_score"),
  stringsAsFactors = FALSE
) |>
  filter(x < y)

pooled_cor <- lapply(seq_len(nrow(all_pairs)), function(i) {
  ct <- spearman_ci(meta[[all_pairs$x[i]]], meta[[all_pairs$y[i]]])
  data.frame(x_var = all_pairs$x[i], y_var = all_pairs$y[i],
             rho = ct$rho, ci_lo = ct$ci_lo, ci_hi = ct$ci_hi, p = ct$p, n = ct$n)
}) |> do.call(what = rbind)

pooled_cor$p_adj <- p.adjust(pooled_cor$p, method = "BH")
write.csv(pooled_cor, file.path(out_dir, "Table_pooled_correlations.csv"), row.names = FALSE)
cat("\n--- Pooled correlations ---\n")
print(pooled_cor)

# =========================================================================== #
# 13. SAVE NORMALISED EXPRESSION FOR S40 (combined cohort) ----
# =========================================================================== #

saveRDS(
  list(expr_mat = expr_mat, meta = meta),
  file = file.path(out_dir, "S30_normalised.rds")
)

cat("\nS30 complete. Outputs saved to:", out_dir, "\n")
cat("Samples saved for S40: AKI =",
    sum(meta$enrollment_category == "AKI"), "| CKD =",
    sum(meta$enrollment_category == "CKD"), "\n")
