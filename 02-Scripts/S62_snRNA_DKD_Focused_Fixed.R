# =============================================================================
# S62_snRNA_DKD_Focused_Fixed.R
# Variant of S60: CKD patients restricted to Diabetic Kidney Disease only.
# Healthy Reference and AKI patients retained as-is.
#
# Loads: S60_pseudobulk.rds  (pseudobulk data from S60)
# Output: 04-Outputs/S62_snRNA_DKD_Focused/
# =============================================================================

rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

# =============================================================================
# SETTINGS ----
# =============================================================================

s60_rds <- "../04-Outputs/S60_KPMP_snRNA/S60_pseudobulk.rds"
out_dir <- "../04-Outputs/S62_snRNA_DKD_Focused_Fixed"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

condition_levels <- c("Healthy Reference", "AKI", "CKD")
condition_short  <- c("Healthy", "AKI", "CKD\n(DKD only)")

condition_colors <- c(
  "Healthy Reference" = "#2166ac",
  "AKI"               = "#e07b39",
  "CKD"               = "#756bb1"
)

comparisons_main <- list(
  c("Healthy Reference", "AKI"),
  c("Healthy Reference", "CKD")
)

dpi <- 300

# =============================================================================
# HELPER FUNCTIONS ----
# =============================================================================

format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("< 0.001")
  sprintf("%.3f", p)
}

spearman_ci <- function(x, y) {
  keep <- !is.na(x) & !is.na(y)
  x <- x[keep]
  y <- y[keep]
  
  if (length(x) < 4) {
    return(list(
      rho = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
      p = NA_real_, n = length(x),
      label = "rho = NA [NA, NA]\np = NA"
    ))
  }
  
  ct  <- cor.test(x, y, method = "spearman", exact = FALSE)
  rho <- unname(ct$estimate)
  n   <- length(x)
  
  if (is.na(rho) || abs(rho) >= 1) {
    return(list(
      rho = rho, ci_lo = NA_real_, ci_hi = NA_real_,
      p = ct$p.value, n = n,
      label = sprintf("rho = %.2f [NA, NA]\np = %s", rho, format_p(ct$p.value))
    ))
  }
  
  z  <- atanh(rho)
  se <- 1 / sqrt(n - 3)
  ci <- tanh(c(z - 1.96 * se, z + 1.96 * se))
  
  list(
    rho = rho,
    ci_lo = ci[1],
    ci_hi = ci[2],
    p = ct$p.value,
    n = n,
    label = sprintf("rho = %.2f [%.2f, %.2f]\np = %s",
                    rho, ci[1], ci[2], format_p(ct$p.value))
  )
}

theme_box <- function() {
  theme_bw(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 20, hjust = 1),
      plot.title  = element_text(face = "bold"),
      strip.text  = element_text(face = "bold")
    )
}

make_subtitle <- function(pb) {
  sprintf(
    "n patients: %s",
    paste(
      mapply(function(lv, sh) {
        n <- sum(as.character(pb$enrollment) == lv, na.rm = TRUE)
        sprintf("%s=%d", gsub("\n", " ", sh), n)
      }, condition_levels, condition_short),
      collapse = ", "
    )
  )
}

plot_box <- function(df, y_col, y_label, title_str, subtitle_str) {
  ymax <- max(df[[y_col]], na.rm = TRUE)
  ymin <- min(df[[y_col]], na.rm = TRUE)
  rng  <- ymax - ymin
  
  if (!is.finite(rng) || rng == 0) {
    rng <- max(abs(c(ymax, ymin)), na.rm = TRUE)
    if (!is.finite(rng) || rng == 0) rng <- 1
  }
  
  ggplot(df, aes(x = enrollment, y = .data[[y_col]], fill = enrollment)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.45, width = 0.5) +
    geom_jitter(width = 0.12, size = 2.2, alpha = 0.85) +
    stat_compare_means(
      comparisons = comparisons_main,
      method = "wilcox.test",
      label = "p.signif",
      hide.ns = FALSE
    ) +
    scale_fill_manual(values = condition_colors, guide = "none") +
    scale_x_discrete(labels = setNames(condition_short, condition_levels), drop = FALSE) +
    scale_y_continuous(limits = c(ymin - 0.05 * rng, ymax + 0.26 * rng)) +
    labs(x = NULL, y = y_label, title = title_str, subtitle = subtitle_str) +
    theme_box()
}

plot_cor <- function(df, x_col, y_col, x_lab, y_lab, title_str, ct) {
  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]], color = enrollment)) +
    geom_point(size = 2.8, alpha = 0.85) +
    geom_smooth(
      method = "lm", se = FALSE, color = "black", linewidth = 0.8,
      inherit.aes = FALSE,
      aes(x = .data[[x_col]], y = .data[[y_col]])
    ) +
    annotate(
      "text", x = Inf, y = Inf,
      label = ct$label, hjust = 1.05, vjust = 1.2, size = 3.8
    ) +
    scale_color_manual(values = condition_colors, name = NULL, drop = FALSE) +
    labs(
      x = x_lab, y = y_lab, title = title_str,
      subtitle = sprintf("n = %d matched patients", nrow(df))
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
}

# =============================================================================
# 1. LOAD AND APPLY DKD FILTER ----
# =============================================================================

s60 <- readRDS(s60_rds)

cat("Loaded S60 pseudobulk object.\n")
cat("Objects available:", paste(names(s60), collapse = ", "), "\n\n")

cat("Columns in pb_pt:\n")
print(colnames(s60$pb_pt))
cat("\nColumns in pb_fib:\n")
print(colnames(s60$pb_fib))
cat("\nColumns in pb_fib_sub:\n")
print(colnames(s60$pb_fib_sub))

# CKD adjudication exists only in pb_fib
pb_ckd_all <- s60$pb_fib %>% filter(enrollment == "CKD")
cat("\nCKD patients before DKD filter:", nrow(pb_ckd_all), "\n")
cat("Primary adjudicated categories:\n")
print(table(pb_ckd_all$adjudicated, useNA = "ifany"))

# 1) Filter fibroblast table directly by adjudicated
pb_fib <- s60$pb_fib %>%
  filter(!(enrollment == "CKD" & adjudicated != "Diabetic Kidney Disease")) %>%
  mutate(enrollment = factor(as.character(enrollment), levels = condition_levels))

# 2) Use retained fibroblast patients as DKD-only CKD set
keep_patients <- unique(pb_fib$patient)

cat("\nPatients retained after DKD filter:", length(keep_patients), "\n")

# 3) Apply same patient-level restriction to PT
pb_pt <- s60$pb_pt %>%
  filter(patient %in% keep_patients) %>%
  mutate(enrollment = factor(as.character(enrollment), levels = condition_levels))

# 4) Apply same patient-level restriction to subtype table
pb_fib_sub <- s60$pb_fib_sub %>%
  filter(patient %in% keep_patients) %>%
  mutate(enrollment = factor(as.character(enrollment), levels = condition_levels))

cat("CKD patients after DKD filter:", sum(as.character(pb_fib$enrollment) == "CKD", na.rm = TRUE), "\n")
cat("Final N — PT:", nrow(pb_pt), "| FB:", nrow(pb_fib), "| FB subtype rows:", nrow(pb_fib_sub), "\n")

# Cross-compartment: re-join after filter
pb_cross <- inner_join(
  pb_pt  %>% select(patient, enrollment, pyrim_score),
  pb_fib %>% select(patient, mean_P2RY6, pct_P2RY6_pos, mean_ecm_score),
  by = "patient"
) %>%
  mutate(enrollment = factor(as.character(enrollment), levels = condition_levels))

cat("Matched patients for correlation:", nrow(pb_cross), "\n\n")

# =============================================================================
# 2. BOXPLOTS — PT ----
# =============================================================================

subtitle_pt  <- make_subtitle(pb_pt)
subtitle_fib <- make_subtitle(pb_fib)

p_pt_score <- plot_box(
  pb_pt, "pyrim_score",
  "Pyrimidine module score (mean expression)",
  "KPMP snRNA (DKD-focused) — PT pyrimidine score",
  subtitle_pt
)

ggsave(
  file.path(out_dir, "Fig_PT_pyrimidine_score.png"),
  p_pt_score, width = 6, height = 5, dpi = dpi
)

pyrim_found <- intersect(
  c("CAD", "DHODH", "UMPS", "UCK1", "CDA", "DCTPP1"),
  colnames(pb_pt)
)

if (length(pyrim_found) > 0) {
  df_pt_long <- pb_pt %>%
    select(patient, enrollment, all_of(pyrim_found)) %>%
    pivot_longer(
      cols = all_of(pyrim_found),
      names_to = "gene",
      values_to = "expr"
    ) %>%
    mutate(enrollment = factor(as.character(enrollment), levels = condition_levels))
  
  p_pt_genes <- ggplot(df_pt_long, aes(x = enrollment, y = expr, fill = enrollment)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.45, width = 0.5) +
    geom_jitter(width = 0.12, size = 1.4, alpha = 0.85) +
    stat_compare_means(
      comparisons = comparisons_main,
      method = "wilcox.test",
      label = "p.signif",
      hide.ns = FALSE
    ) +
    scale_fill_manual(values = condition_colors, guide = "none") +
    scale_x_discrete(labels = setNames(condition_short, condition_levels), drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.55))) +
    facet_wrap(~gene, scales = "free_y", ncol = 3) +
    labs(
      x = NULL,
      y = "Mean expression (pseudobulk)",
      subtitle = subtitle_pt,
      title = "KPMP snRNA (DKD-focused) — Individual pyrimidine genes in PT"
    ) +
    theme_box()
  
  ggsave(
    file.path(out_dir, "Fig_PT_pyrimidine_genes.png"),
    p_pt_genes, width = 11, height = 7, dpi = dpi
  )
}

# =============================================================================
# 3. BOXPLOTS — FIBROBLASTS ----
# =============================================================================

p_p2ry6 <- plot_box(
  pb_fib, "mean_P2RY6",
  "P2RY6 mean expression (pseudobulk)",
  "KPMP snRNA (DKD-focused) — P2RY6 in fibroblasts",
  subtitle_fib
)

p_pct <- plot_box(
  pb_fib, "pct_P2RY6_pos",
  "% P2RY6+ cells per patient",
  "KPMP snRNA (DKD-focused) — P2RY6 detection",
  subtitle_fib
)

p_calcium <- plot_box(
  pb_fib, "mean_calcium_score",
  "Calcium signaling score (mean expression)",
  "KPMP snRNA (DKD-focused) — Calcium score in FBs",
  subtitle_fib
)

p_ecm <- plot_box(
  pb_fib, "mean_ecm_score",
  "ECM score (mean expression)",
  "KPMP snRNA (DKD-focused) — ECM score in FBs",
  subtitle_fib
)

p_pdgfrb <- plot_box(
  pb_fib, "mean_PDGFRB",
  "PDGFRB mean expression (pseudobulk)",
  "KPMP snRNA (DKD-focused) — PDGFRB in fibroblasts",
  subtitle_fib
)

ggsave(file.path(out_dir, "Fig_FB_P2RY6.png"),     p_p2ry6,   width = 6, height = 5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_FB_P2RY6_pct.png"), p_pct,     width = 6, height = 5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_FB_calcium.png"),   p_calcium, width = 6, height = 5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_FB_ECM.png"),       p_ecm,     width = 6, height = 5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_FB_PDGFRB.png"),    p_pdgfrb,  width = 6, height = 5, dpi = dpi)

if (nrow(pb_fib_sub) > 0) {
  p_sub <- ggplot(pb_fib_sub, aes(x = enrollment, y = mean_P2RY6, fill = enrollment)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.45, width = 0.5) +
    geom_jitter(width = 0.12, size = 1.4, alpha = 0.85) +
    stat_compare_means(
      comparisons = comparisons_main,
      method = "wilcox.test",
      label = "p.signif",
      hide.ns = FALSE
    ) +
    scale_fill_manual(values = condition_colors, guide = "none") +
    scale_x_discrete(labels = setNames(condition_short, condition_levels), drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.55))) +
    facet_wrap(~subclass.l2, scales = "free_y", ncol = 4) +
    labs(
      x = NULL,
      y = "P2RY6 mean expression (pseudobulk)",
      title = "KPMP snRNA (DKD-focused) — P2RY6 per fibroblast subtype",
      subtitle = "Min 20 cells per patient per subtype"
    ) +
    theme_box()
  
  ggsave(
    file.path(out_dir, "Fig_FB_P2RY6_by_subtype.png"),
    p_sub, width = 14, height = 8, dpi = dpi
  )
}

# =============================================================================
# 4. CROSS-COMPARTMENT CORRELATION ----
# =============================================================================

ct_p2ry6 <- spearman_ci(pb_cross$pyrim_score, pb_cross$mean_P2RY6)
ct_ecm   <- spearman_ci(pb_cross$pyrim_score, pb_cross$mean_ecm_score)

cat("PT pyrimidine vs FB P2RY6:", ct_p2ry6$label, "\n")
cat("PT pyrimidine vs FB ECM:  ", ct_ecm$label, "\n")

p_cor_p2ry6 <- plot_cor(
  pb_cross, "pyrim_score", "mean_P2RY6",
  "PT pyrimidine score (pseudobulk)",
  "FB P2RY6 (pseudobulk)",
  "KPMP snRNA (DKD-focused) — PT pyrimidine vs FB P2RY6",
  ct_p2ry6
)

p_cor_ecm <- plot_cor(
  pb_cross, "pyrim_score", "mean_ecm_score",
  "PT pyrimidine score (pseudobulk)",
  "FB ECM score (pseudobulk)",
  "KPMP snRNA (DKD-focused) — PT pyrimidine vs FB ECM",
  ct_ecm
)

ggsave(
  file.path(out_dir, "Fig_Cor_PT_Pyrimidine_vs_FB_P2RY6.png"),
  p_cor_p2ry6, width = 7, height = 4.5, dpi = dpi
)

ggsave(
  file.path(out_dir, "Fig_Cor_PT_Pyrimidine_vs_FB_ECM.png"),
  p_cor_ecm, width = 7, height = 4.5, dpi = dpi
)

ct_p2ry6_calcium <- spearman_ci(x = pb_fib$mean_P2RY6, y = pb_fib$mean_calcium_score)
p_cor_p2ry6_calcium <- plot_cor(
  df        = pb_fib,
  x_col     = "mean_P2RY6",
  y_col     = "mean_calcium_score",
  x_lab     = "FB P2RY6 (pseudobulk)",
  y_lab     = "FB calcium score (pseudobulk)",
  title_str = "KPMP snRNA (DKD-focused) — FB P2RY6 vs Calcium",
  ct        = ct_p2ry6_calcium
) + labs(subtitle = sprintf("n = %d fibroblast patients", ct_p2ry6_calcium$n))
ggsave(
  filename = file.path(out_dir, "Fig_Cor_FB_P2RY6_vs_Calcium.png"),
  plot = p_cor_p2ry6_calcium, width = 7, height = 4.5, dpi = dpi
)

ct_p2ry6_ecm <- spearman_ci(x = pb_fib$mean_P2RY6, y = pb_fib$mean_ecm_score)
p_cor_p2ry6_ecm <- plot_cor(
  df        = pb_fib,
  x_col     = "mean_P2RY6",
  y_col     = "mean_ecm_score",
  x_lab     = "FB P2RY6 (pseudobulk)",
  y_lab     = "FB ECM score (pseudobulk)",
  title_str = "KPMP snRNA (DKD-focused) — FB P2RY6 vs ECM",
  ct        = ct_p2ry6_ecm
) + labs(subtitle = sprintf("n = %d fibroblast patients", ct_p2ry6_ecm$n))
ggsave(
  filename = file.path(out_dir, "Fig_Cor_FB_P2RY6_vs_ECM.png"),
  plot = p_cor_p2ry6_ecm, width = 7, height = 4.5, dpi = dpi
)

ct_pdgfrb_p2ry6 <- spearman_ci(x = pb_fib$mean_PDGFRB, y = pb_fib$mean_P2RY6)
p_cor_pdgfrb_p2ry6 <- plot_cor(
  df        = pb_fib,
  x_col     = "mean_PDGFRB",
  y_col     = "mean_P2RY6",
  x_lab     = "FB PDGFRB (pseudobulk)",
  y_lab     = "FB P2RY6 (pseudobulk)",
  title_str = "KPMP snRNA (DKD-focused) — FB PDGFRB vs P2RY6",
  ct        = ct_pdgfrb_p2ry6
) + labs(subtitle = sprintf("n = %d fibroblast patients", ct_pdgfrb_p2ry6$n))
ggsave(
  filename = file.path(out_dir, "Fig_Cor_FB_PDGFRB_vs_P2RY6.png"),
  plot = p_cor_pdgfrb_p2ry6, width = 7, height = 4.5, dpi = dpi
)

# =============================================================================
# 5. SAVE FILTERED OBJECTS ----
# =============================================================================

saveRDS(
  list(
    pb_pt      = pb_pt,
    pb_fib     = pb_fib,
    pb_fib_sub = pb_fib_sub,
    pb_cross   = pb_cross
  ),
  file.path(out_dir, "S62_pseudobulk_DKD_only.rds")
)

cat("\nS62 complete. Outputs saved to:", out_dir, "\n")
cat("Final N — PT:", nrow(pb_pt), "| FB:", nrow(pb_fib), "| Matched:", nrow(pb_cross), "\n")
