# ============================================================================#
# S51_scRNA_AdjudicatedColors.R
# Variant of S50: same pseudobulk data, CKD points coloured by primary
# adjudicated category (Diabetic KD, Hypertensive KD, etc.) to reveal which
# CKD sub-diagnosis drives each signal.
#
# Loads : S50_pseudobulk.rds
# Output: 04-Outputs/S51_scRNA_AdjudicatedColors/
# ============================================================================#

rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

# ============================================================================#
# SETTINGS ----
# ============================================================================#

s50_rds <- "../04-Outputs/S50_KPMP_scRNA/S50_pseudobulk.rds"
out_dir <- "../04-Outputs/S51_scRNA_AdjudicatedColors"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

condition_levels <- c("Healthy Reference Tissue",
                       "Acute Kidney Injury",
                       "Chronic Kidney Disease")
condition_short  <- c("Healthy", "AKI", "CKD")

# Box fill: by enrollment condition
condition_colors <- c(
  "Healthy Reference Tissue" = "#2166ac",
  "Acute Kidney Injury"      = "#e07b39",
  "Chronic Kidney Disease"   = "#756bb1"
)

# Point color: by adjudicated category (CKD subtype)
adj_colors <- c(
  "Diabetic Kidney Disease"       = "#1f78b4",
  "Hypertensive Kidney Disease"   = "#8e44ad",
  "Focal Segmental Glomerulosclerosis" = "#27ae60",
  "Lupus Nephritis"               = "#e74c3c",
  "Acute Tubular Injury"          = "#f39c12",
  "Acute Interstitial Nephritis"  = "#e67e22",
  "IgA Nephropathy"               = "#16a085",
  "Cannot be determined"          = "#95a5a6",
  "Other"                         = "#7f8c8d"
)

# For non-CKD samples: inherit condition color
non_ckd_colors <- c(
  "Healthy Reference Tissue" = "#2166ac",
  "Acute Kidney Injury"      = "#e07b39"
)

comparisons_main <- list(
  c("Healthy Reference Tissue", "Acute Kidney Injury"),
  c("Healthy Reference Tissue", "Chronic Kidney Disease")
)

dpi <- 300

# ============================================================================#
# HELPER FUNCTIONS ----
# ============================================================================#

format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("< 0.001")
  sprintf("%.3f", p)
}

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

theme_box <- function() {
  theme_bw(base_size = 13) +
    theme(axis.text.x  = element_text(angle = 20, hjust = 1),
          plot.title   = element_text(face = "bold"),
          strip.text   = element_text(face = "bold"),
          legend.position  = "right",
          legend.text  = element_text(size = 9),
          legend.title = element_text(size = 10, face = "bold"))
}

make_subtitle <- function(pb) {
  active_levels <- condition_levels[condition_levels %in% unique(as.character(pb$EnrollmentCategory))]
  active_short  <- condition_short [match(active_levels, condition_levels)]
  sprintf("n samples: %s",
    paste(mapply(function(lv, sh) {
      sprintf("%s=%d", sh, sum(as.character(pb$EnrollmentCategory) == lv, na.rm = TRUE))
    }, active_levels, active_short), collapse = ", "))
}

# ============================================================================#
# 1. LOAD PSEUDOBULK AND BUILD COLOR MAPPING ----
# ============================================================================#

s50 <- readRDS(s50_rds)
pb_pt    <- s50$pb_pt
pb_fib   <- s50$pb_fib
pb_cross <- s50$pb_cross

cat("Loaded S50 pseudobulk.\n")
cat("Adjudicated categories present:\n")
print(table(pb_fib$adjudicated, useNA = "ifany"))

# pb_cross and pb_pt were saved without adjudicated — join from pb_fib by SampleID
adj_lookup <- pb_fib %>% select(SampleID, adjudicated)
if (!"adjudicated" %in% colnames(pb_cross)) {
  pb_cross <- left_join(pb_cross, adj_lookup, by = "SampleID")
}
if (!"adjudicated" %in% colnames(pb_pt)) {
  pb_pt <- left_join(pb_pt, adj_lookup, by = "SampleID")
}

# Build point_color column:
# - Non-CKD: use condition color
# - CKD: use adjudicated category color (unknown → "Other")
build_point_color <- function(df) {
  df$point_color <- ifelse(
    as.character(df$EnrollmentCategory) == "Chronic Kidney Disease",
    ifelse(
      !is.na(df$adjudicated) & df$adjudicated %in% names(adj_colors),
      df$adjudicated,
      "Other"
    ),
    as.character(df$EnrollmentCategory)
  )
  df$point_color <- factor(
    df$point_color,
    levels = c(names(adj_colors), names(non_ckd_colors))
  )
  df
}

all_point_colors <- c(adj_colors, non_ckd_colors)

pb_fib   <- build_point_color(pb_fib)
pb_cross <- build_point_color(pb_cross)

# Add point_color to PT data using enrollment only (no adjudicated for PT table)
pb_pt <- pb_pt %>%
  mutate(point_color = factor(
    ifelse(as.character(EnrollmentCategory) %in% names(non_ckd_colors),
           as.character(EnrollmentCategory), "Other"),
    levels = c(names(adj_colors), names(non_ckd_colors))
  ))

# ============================================================================#
# 2. BOXPLOTS WITH ADJUDICATED POINT COLORS ----
# ============================================================================#

# Build dynamic comparisons (only include conditions present in each df)
get_comparisons <- function(df) {
  present <- unique(as.character(df$EnrollmentCategory))
  Filter(function(p) all(p %in% present), comparisons_main)
}

plot_adj_box <- function(df, y_col, y_label, title_str) {
  subtitle_str <- make_subtitle(df)
  cmps  <- get_comparisons(df)
  ymax  <- max(df[[y_col]], na.rm = TRUE)
  ymin  <- min(df[[y_col]], na.rm = TRUE)
  rng   <- ymax - ymin
  if (!is.finite(rng) || rng == 0) rng <- max(abs(c(ymax, ymin)), 1)

  ggplot(df, aes(x = EnrollmentCategory, y = .data[[y_col]])) +
    geom_boxplot(aes(fill = EnrollmentCategory),
                 outlier.shape = NA, alpha = 0.45, width = 0.5) +
    geom_jitter(aes(color = point_color),
                width = 0.12, size = 2.4, alpha = 0.95) +
    {if (length(cmps) > 0)
       stat_compare_means(comparisons = cmps, method = "wilcox.test",
                          label = "p.signif", hide.ns = FALSE)} +
    scale_fill_manual(values = condition_colors, guide = "none") +
    scale_color_manual(values = all_point_colors, name = "Adjudicated",
                       drop = TRUE,
                       guide = guide_legend(ncol = 1, override.aes = list(size = 3))) +
    scale_x_discrete(labels = setNames(condition_short, condition_levels), drop = FALSE) +
    scale_y_continuous(limits = c(ymin - 0.05 * rng, ymax + 0.26 * rng)) +
    labs(x = NULL, y = y_label, title = title_str, subtitle = subtitle_str) +
    theme_box()
}

subtitle_fib <- make_subtitle(pb_fib)

p_p2ry6   <- plot_adj_box(pb_fib, "mean_P2RY6",         "P2RY6 mean expression (pseudobulk)",       "KPMP scRNA — P2RY6 in stroma (adjudicated)")
p_pct     <- plot_adj_box(pb_fib, "pct_P2RY6_pos",      "% P2RY6+ cells per sample",                "KPMP scRNA — P2RY6 detection (adjudicated)")
p_calcium <- plot_adj_box(pb_fib, "mean_calcium_score",  "Calcium signaling score (mean expression)", "KPMP scRNA — Calcium score (adjudicated)")
p_ecm     <- plot_adj_box(pb_fib, "mean_ecm_score",      "ECM score (mean expression)",               "KPMP scRNA — ECM score (adjudicated)")
p_pdgfrb  <- plot_adj_box(pb_fib, "mean_PDGFRB",         "PDGFRB mean expression (pseudobulk)",       "KPMP scRNA — PDGFRB (adjudicated)")
p_pt      <- plot_adj_box(pb_pt,  "pyrim_score",         "Pyrimidine module score (mean expression)", "KPMP scRNA — PT pyrimidine (adjudicated)")

ggsave(file.path(out_dir, "Fig_FB_P2RY6.png"),          p_p2ry6,   width = 7, height = 5.5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_FB_P2RY6_pct.png"),      p_pct,     width = 7, height = 5.5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_FB_calcium.png"),         p_calcium, width = 7, height = 5.5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_FB_ECM.png"),             p_ecm,     width = 7, height = 5.5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_FB_PDGFRB.png"),          p_pdgfrb,  width = 7, height = 5.5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_PT_pyrimidine_score.png"), p_pt,     width = 7, height = 5.5, dpi = dpi)

# ============================================================================#
# 3. CORRELATION SCATTER — POINTS COLOURED BY ADJUDICATED ----
# ============================================================================#

plot_cor_adj <- function(df, x_col, y_col, x_lab, y_lab, title_str, ct) {
  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]],
                 color = point_color)) +
    geom_point(size = 2.8, alpha = 0.9) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8,
                inherit.aes = FALSE, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    annotate("text", x = Inf, y = Inf,
             label = ct$label, hjust = 1.05, vjust = 1.2, size = 3.8) +
    scale_color_manual(values = all_point_colors, name = "Adjudicated",
                       drop = TRUE,
                       guide = guide_legend(ncol = 1, override.aes = list(size = 3))) +
    labs(x = x_lab, y = y_lab, title = title_str,
         subtitle = sprintf("n = %d matched samples", nrow(df))) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "right",
          legend.text  = element_text(size = 9),
          legend.title = element_text(size = 10, face = "bold"))
}

ct_p2ry6 <- spearman_ci(pb_cross$pyrim_score, pb_cross$mean_P2RY6)
ct_ecm   <- spearman_ci(pb_cross$pyrim_score, pb_cross$mean_ecm_score)

p_cor_p2ry6 <- plot_cor_adj(pb_cross, "pyrim_score", "mean_P2RY6",
  "PT pyrimidine score (pseudobulk)", "Stroma P2RY6 (pseudobulk)",
  "KPMP scRNA — PT pyrimidine vs Stroma P2RY6 (adjudicated)", ct_p2ry6)
p_cor_ecm <- plot_cor_adj(pb_cross, "pyrim_score", "mean_ecm_score",
  "PT pyrimidine score (pseudobulk)", "Stroma ECM score (pseudobulk)",
  "KPMP scRNA — PT pyrimidine vs Stroma ECM (adjudicated)", ct_ecm)

ggsave(file.path(out_dir, "Fig_Cor_PT_Pyrimidine_vs_FB_P2RY6.png"),
       p_cor_p2ry6, width = 7.5, height = 4.5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_Cor_PT_Pyrimidine_vs_FB_ECM.png"),
       p_cor_ecm, width = 7.5, height = 4.5, dpi = dpi)

ct_p2ry6_calcium <- spearman_ci(x = pb_fib$mean_P2RY6, y = pb_fib$mean_calcium_score)
p_cor_p2ry6_calcium <- plot_cor_adj(
  df        = pb_fib,
  x_col     = "mean_P2RY6",
  y_col     = "mean_calcium_score",
  x_lab     = "Stroma P2RY6 (pseudobulk)",
  y_lab     = "Stroma calcium score (pseudobulk)",
  title_str = "KPMP scRNA — Stroma P2RY6 vs Calcium (adjudicated)",
  ct        = ct_p2ry6_calcium
) + labs(subtitle = sprintf("n = %d stroma samples", ct_p2ry6_calcium$n))
ggsave(filename = file.path(out_dir, "Fig_Cor_FB_P2RY6_vs_Calcium.png"),
       plot = p_cor_p2ry6_calcium, width = 7.5, height = 4.5, dpi = dpi)

ct_p2ry6_ecm <- spearman_ci(x = pb_fib$mean_P2RY6, y = pb_fib$mean_ecm_score)
p_cor_p2ry6_ecm <- plot_cor_adj(
  df        = pb_fib,
  x_col     = "mean_P2RY6",
  y_col     = "mean_ecm_score",
  x_lab     = "Stroma P2RY6 (pseudobulk)",
  y_lab     = "Stroma ECM score (pseudobulk)",
  title_str = "KPMP scRNA — Stroma P2RY6 vs ECM (adjudicated)",
  ct        = ct_p2ry6_ecm
) + labs(subtitle = sprintf("n = %d stroma samples", ct_p2ry6_ecm$n))
ggsave(filename = file.path(out_dir, "Fig_Cor_FB_P2RY6_vs_ECM.png"),
       plot = p_cor_p2ry6_ecm, width = 7.5, height = 4.5, dpi = dpi)

ct_pdgfrb_p2ry6 <- spearman_ci(x = pb_fib$mean_PDGFRB, y = pb_fib$mean_P2RY6)
p_cor_pdgfrb_p2ry6 <- plot_cor_adj(
  df        = pb_fib,
  x_col     = "mean_PDGFRB",
  y_col     = "mean_P2RY6",
  x_lab     = "Stroma PDGFRB (pseudobulk)",
  y_lab     = "Stroma P2RY6 (pseudobulk)",
  title_str = "KPMP scRNA — Stroma PDGFRB vs P2RY6 (adjudicated)",
  ct        = ct_pdgfrb_p2ry6
) + labs(subtitle = sprintf("n = %d stroma samples", ct_pdgfrb_p2ry6$n))
ggsave(filename = file.path(out_dir, "Fig_Cor_FB_PDGFRB_vs_P2RY6.png"),
       plot = p_cor_pdgfrb_p2ry6, width = 7.5, height = 4.5, dpi = dpi)

cat("\nS51 complete. Outputs saved to:", out_dir, "\n")
