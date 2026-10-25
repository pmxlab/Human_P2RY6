# ============================================================================#
# S80_Manuscript_Figures.R
# Assembles the multi-panel figures of the revised manuscript (MDPI Life,
# Submission 2) from the panel code of the analysis scripts S40-S62.
#
# Replaces the manually collated slides in
#   07-Manuscript/Submission 2 - MDPI-LIFE/Figures.pptx
#
# Panel -> source mapping (verified against the images embedded in the pptx):
#
#   Figure     (a)        (b)        (c)        (d)        (e)        (f)
#   ---------  ---------  ---------  ---------  ---------  ---------  ---------
#   Figure 1   S40        S42        S50        S52        S60        S62      Pyrimidine score
#   Figure 2   S40        S42        S50        S52        S60        S62      P2RY6
#   Figure 3   S40        S42        S50        S52        S60        S62      Calcium score
#   Figure 4   S40        S42        S50        S52        S60        S62      PDGFRB
#   Figure 5   S40        S42        S50        S52        S60        S62      ECM score
#   Figure 6   S41        S42        S50        S52        S60        S62      Pyrimidine vs ECM
#   Figure S1  S40 (pre-ComBat PCA)  S40 (post-ComBat PCA)                    Batch correction QC
#   Figure S2  S40        S42        S51        S52        S60        S62      P2RY6 vs Calcium
#   Figure S3  S40        S42        S51        S52        S60        S62      P2RY6 vs PDGFRB (**)
#   Figure S4  S40        S42        S51        S52        S60        S62      P2RY6 vs ECM (*)
#
#   Figure 2c shows mean P2RY6 expression like the other P2RY6 panels
#       (Figures.pptx showed % P2RY6+ stroma cells; changed in the v5 revision).
#   (*) In Figures.pptx the slide "Suppl. Figure 3" is an exact copy of
#       "Suppl. Figure 2". The manuscript text describes this figure as the
#       P2RY6 vs ECM correlations, which are assembled here following the
#       layout of the other correlation figures.
#   (**) Axes swapped relative to S40-S62 (P2RY6 on x, PDGFRB on y) so that
#       P2RY6 is on the x-axis in all supplementary correlation figures.
#   Supplementary numbering follows the order of citation in the manuscript:
#       S1 (PCA) is cited in the Methods, S2-S4 in the Results.
#
# Inputs (no heavy re-computation; Seurat objects and ComBat are not needed):
#   S20_normalised.rds, S30_normalised.rds, S40_corrected.rds,
#   S50_pseudobulk.rds, S60_pseudobulk.rds
#
# Differences to the panels saved by the source scripts:
#   - significance stars show Benjamini-Hochberg-adjusted p-values (two-sided
#     Wilcoxon rank-sum test, adjusted over the comparisons of each panel);
#     the source scripts showed unadjusted p-values
#   - jitter uses a fixed seed (identical points in PNG and PDF, reproducible)
#   - long titles/subtitles are wrapped instead of being cut at the panel edge
#   - extra top margin per panel to make room for the panel label
#
# Outputs: 04-Outputs/S80_Manuscript_Figures/Figure_<n>.png / .pdf
# ============================================================================#

rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(cowplot)
library(stringr)

# cowplot converts panels to grobs on a "null" device; use Cairo so that its
# font metrics match the Cairo-based PNG/PDF devices (avoids clipped legends)
set_null_device("cairo")

# ============================================================================#
# SETTINGS ----
# ============================================================================#

s20_rds <- "../04-Outputs/S20_FanDKD_BulkRNA/S20_normalised.rds"
s40_rds <- "../04-Outputs/S40_Combined_BulkRNA/S40_corrected.rds"
s30_rds <- "../04-Outputs/S30_KPMP_BulkRNA/S30_normalised.rds"
s50_rds <- "../04-Outputs/S50_KPMP_scRNA/S50_pseudobulk.rds"
s60_rds <- "../04-Outputs/S60_KPMP_snRNA/S60_pseudobulk.rds"
out_dir <- "../04-Outputs/S80_Manuscript_Figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dpi         <- 600   # MDPI: figures preferably >= 600 dpi
jitter_seed <- 2026

# Jitter with fixed seed (replaces geom_jitter(width = w) of the sources)
jitter_pos <- function(width = 0.12) {
  position_jitter(width = width, seed = jitter_seed)
}

# ============================================================================#
# SHARED HELPER FUNCTIONS ----
# ============================================================================#

format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("< 0.001")
  sprintf("%.3f", p)
}

# NA-safe version (S50/S51/S52/S62); gives identical results to the
# S40/S41/S42/S60 version on complete data.
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
    theme(axis.text.x = element_text(angle = 20, hjust = 1),
          plot.title  = element_text(face = "bold"),
          strip.text  = element_text(face = "bold"))
}

# Significance brackets with Benjamini-Hochberg-adjusted p-values ----
# Two-sided Wilcoxon rank-sum test with R's defaults (exact for small groups
# without ties, normal approximation otherwise), i.e. the same test that
# ggpubr::stat_compare_means() runs; the p-values of all comparisons drawn in
# one panel are then BH-adjusted together. Brackets are drawn with
# ggsignif::geom_signif() using ggpubr's default layout.
p_stars <- function(p) {
  as.character(cut(p, c(-Inf, 1e-4, 1e-3, 1e-2, 0.05, Inf),
                   c("****", "***", "**", "*", "ns")))
}
bh_test <- function(df, x_col, y_col, comparisons) {
  x <- as.character(df[[x_col]]); y <- df[[y_col]]
  p <- vapply(comparisons, function(cc)
    suppressWarnings(wilcox.test(y[x == cc[1]], y[x == cc[2]])$p.value), numeric(1))
  data.frame(group1 = vapply(comparisons, `[`, "", 1),
             group2 = vapply(comparisons, `[`, "", 2),
             p = p, p_adj = p.adjust(p, method = "BH"),
             label = p_stars(p.adjust(p, method = "BH")))
}
bh_signif <- function(df, x_col, y_col, comparisons) {
  st <- bh_test(df, x_col, y_col, comparisons)
  ggsignif::geom_signif(comparisons = comparisons, annotations = st$label,
                        step_increase = 0.12, tip_length = 0.03, size = 0.3,
                        textsize = 3.88, vjust = 0)
}

# Layout polish applied to every panel before assembly:
# break the title at " — ", wrap the subtitle, align both to the panel edge
# and leave room at the top for the panel label.
tidy_panel <- function(p, subtitle_width = 52) {
  if (!is.null(p$labels$title)) {
    p$labels$title <- sub(" — ", "\n", p$labels$title, fixed = TRUE)
  }
  if (!is.null(p$labels$subtitle)) {
    p$labels$subtitle <- str_wrap(gsub("\n", " ", p$labels$subtitle),
                                  width = subtitle_width)
  }
  p + theme(plot.title.position = "plot",
            plot.margin = margin(t = 26, r = 8, b = 6, l = 8))
}

# Assemble a 3 x 2 figure, label panels (a)-(f) and save PNG + PDF
save_figure <- function(panels, name, width, heights, ncol = 2) {
  panels <- lapply(panels, tidy_panel)
  fig <- plot_grid(plotlist = panels, ncol = ncol,
                   labels = sprintf("(%s)", letters[seq_along(panels)]),
                   label_size = 18, label_fontface = "bold",
                   rel_heights = heights)
  fig <- ggdraw(fig) + theme(plot.background = element_rect(fill = "white", colour = NA))
  h <- sum(heights)
  ggsave(file.path(out_dir, paste0(name, ".png")), fig,
         width = width, height = h, dpi = dpi, bg = "white")
  # Render the PDF to a temp file first: on Windows/WSL a PDF that is open in a
  # viewer cannot be overwritten; warn instead of aborting the whole script.
  tmp_pdf <- tempfile(fileext = ".pdf")
  ggsave(tmp_pdf, fig, width = width, height = h, device = cairo_pdf)
  pdf_ok <- suppressWarnings(file.copy(tmp_pdf, file.path(out_dir, paste0(name, ".pdf")),
                                       overwrite = TRUE))
  if (!pdf_ok) warning(name, ".pdf could not be written (file open in a viewer?)",
                       call. = FALSE)
  cat("Saved", name, if (!pdf_ok) "(PNG only, PDF locked)", "\n")
  invisible(fig)
}

# ============================================================================#
# DATA ----
# ============================================================================#

# . Bulk: S40 (combined, all conditions) ----
s40      <- readRDS(s40_rds)
meta_s40 <- s40$meta   # sample_id, condition, cohort + module scores (S40 §5)

# . Bulk: S41 (KPMP points coloured by adjudicated category) ----
s30      <- readRDS(s30_rds)
adj_info <- s30$meta %>% select(sample_id, primary_adjudicated_category)
meta_s41 <- left_join(meta_s40, adj_info, by = "sample_id")

# . Bulk: S42 (S41 + KPMP CKD restricted to Diabetic Kidney Disease) ----
meta_s42 <- meta_s41 %>%
  filter(!(cohort == "KPMP" & condition == "CKD" &
             primary_adjudicated_category != "Diabetic Kidney Disease"))

# Module scores are per-sample means of the batch-corrected matrix, so the
# S40 values are identical to those recomputed in S41/S42.

# . scRNA: S50 (all CKD) ----
s50         <- readRDS(s50_rds)
pb_pt_s50   <- s50$pb_pt
pb_fib_s50  <- s50$pb_fib
pb_cross_s50 <- s50$pb_cross

# . scRNA: S52 (CKD restricted to DKD; S52 §1) ----
s52_levels <- c("Healthy Reference Tissue", "Acute Kidney Injury", "Chronic Kidney Disease")
pb_fib_s52 <- s50$pb_fib %>%
  filter(!(as.character(EnrollmentCategory) == "Chronic Kidney Disease" &
           adjudicated != "Diabetic Kidney Disease")) %>%
  mutate(EnrollmentCategory = factor(as.character(EnrollmentCategory),
                                     levels = s52_levels))
keep_samples_s52 <- unique(pb_fib_s52$SampleID)
pb_pt_s52 <- s50$pb_pt %>%
  filter(SampleID %in% keep_samples_s52) %>%
  mutate(EnrollmentCategory = factor(as.character(EnrollmentCategory),
                                     levels = s52_levels))
pb_cross_s52 <- inner_join(
  pb_pt_s52  %>% select(SampleID, EnrollmentCategory, pyrim_score),
  pb_fib_s52 %>% select(SampleID, mean_P2RY6, pct_P2RY6_pos, mean_ecm_score),
  by = "SampleID"
) %>%
  mutate(EnrollmentCategory = factor(as.character(EnrollmentCategory),
                                     levels = s52_levels))

# . snRNA: S60 (all CKD) ----
s60          <- readRDS(s60_rds)
pb_pt_s60    <- s60$pb_pt
pb_fib_s60   <- s60$pb_fib
pb_cross_s60 <- s60$pb_cross

# . snRNA: S62 (CKD restricted to DKD; S62 §1) ----
s62_levels <- c("Healthy Reference", "AKI", "CKD")
pb_fib_s62 <- s60$pb_fib %>%
  filter(!(enrollment == "CKD" & adjudicated != "Diabetic Kidney Disease")) %>%
  mutate(enrollment = factor(as.character(enrollment), levels = s62_levels))
keep_patients_s62 <- unique(pb_fib_s62$patient)
pb_pt_s62 <- s60$pb_pt %>%
  filter(patient %in% keep_patients_s62) %>%
  mutate(enrollment = factor(as.character(enrollment), levels = s62_levels))
pb_cross_s62 <- inner_join(
  pb_pt_s62  %>% select(patient, enrollment, pyrim_score),
  pb_fib_s62 %>% select(patient, mean_P2RY6, pct_P2RY6_pos, mean_ecm_score),
  by = "patient"
) %>%
  mutate(enrollment = factor(as.character(enrollment), levels = s62_levels))

# ============================================================================#
# PANEL CODE BY SOURCE SCRIPT ----
# ============================================================================#
# Plot functions are copied from the source scripts. Script-level globals
# (colours, labels, subtitles) are kept in a settings list per source so the
# functions of different scripts do not clash.

# . S40_Combined_BulkRNA.R ----
s40_set <- list(
  condition_levels = c("Normal", "Early_DKD", "Advanced_DKD", "AKI", "CKD"),
  condition_labels = c("Normal", "Early DKD", "Adv. DKD", "AKI", "CKD"),
  condition_colors = c("Normal" = "#2166ac", "Early_DKD" = "#f4a582",
                       "Advanced_DKD" = "#d6604d", "AKI" = "#e07b39",
                       "CKD" = "#756bb1"),
  cohort_shapes    = c("Fan" = 16, "KPMP" = 17),
  comparisons_main = list(c("Normal", "Early_DKD"), c("Normal", "Advanced_DKD"),
                          c("Normal", "AKI"), c("Normal", "CKD"))
)

n_subtitle <- function(meta, levels, labels) {
  ns  <- table(meta$condition)
  lbs <- setNames(labels, levels)
  paste(sapply(levels[levels %in% names(ns)],
               function(lv) sprintf("%s (n=%d)", lbs[lv], ns[lv])),
        collapse = ", ")
}
s40_set$subtitle_n <- n_subtitle(meta_s40, s40_set$condition_levels,
                                 s40_set$condition_labels)

s40_box <- function(y_col, y_label, title_str, df = meta_s40, S = s40_set) {
  cmp         <- S$comparisons_main
  ymax        <- max(df[[y_col]], na.rm = TRUE)
  ymin        <- min(df[[y_col]], na.rm = TRUE)
  rng         <- ymax - ymin
  upper_extra <- (0.08 * length(cmp) + 0.10) * rng
  ggplot(df, aes(x = condition, y = .data[[y_col]], fill = condition)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.55, width = 0.55) +
    geom_jitter(aes(shape = cohort), position = jitter_pos(0.12),
                size = 2.2, alpha = 0.9) +
    bh_signif(df, "condition", y_col, cmp) +
    scale_fill_manual(values = S$condition_colors) +
    scale_shape_manual(values = S$cohort_shapes) +
    scale_x_discrete(labels = setNames(S$condition_labels, S$condition_levels)) +
    scale_y_continuous(limits = c(ymin - 0.05 * rng, ymax + upper_extra)) +
    labs(x = NULL, y = y_label, title = title_str, subtitle = S$subtitle_n,
         shape = "Cohort") +
    theme_box() +
    theme(legend.position = "right",
          legend.title    = element_text(size = 10))
}

s40_scatter <- function(x_col, y_col, x_label, y_label, title_str,
                        df = meta_s40, S = s40_set) {
  ct <- spearman_ci(df[[x_col]], df[[y_col]])
  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]],
                 color = condition, shape = cohort)) +
    geom_point(size = 2.5, alpha = 0.9) +
    geom_smooth(method = "lm", se = FALSE, color = "black",
                linewidth = 0.8, inherit.aes = FALSE,
                aes(x = .data[[x_col]], y = .data[[y_col]])) +
    scale_color_manual(values = S$condition_colors) +
    scale_shape_manual(values = S$cohort_shapes) +
    annotate("text", x = Inf, y = Inf,
             label = ct$label, hjust = 1.05, vjust = 1.2, size = 3.8) +
    labs(x = x_label, y = y_label, title = title_str, subtitle = S$subtitle_n,
         color = NULL, shape = "Cohort") +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold"), legend.position = "right")
}

# . S41_Combined_AdjudicatedColors.R / S42_Combined_DKD_Focused.R ----
# Both scripts share the same plot functions; they differ only in the data
# (S42: KPMP CKD = DKD only) and in the CKD axis label.
adj_bulk_colors <- c(
  # Fan: coloured by condition
  "Normal"       = "#2166ac",
  "Early_DKD"    = "#f4a582",
  "Advanced_DKD" = "#d6604d",
  # KPMP: coloured by primary adjudicated category
  "Acute Tubular Injury"         = "#f39c12",
  "Acute Interstitial Nephritis" = "#e74c3c",
  "Diabetic Kidney Disease"      = "#2980b9",
  "Hypertensive Kidney Disease"  = "#8e44ad",
  "Other"                        = "#27ae60",
  "Cannot be determined"         = "#95a5a6"
)

add_color_group <- function(meta) {
  meta$color_group <- ifelse(meta$cohort == "Fan",
                             as.character(meta$condition),
                             meta$primary_adjudicated_category)
  meta$color_group <- factor(meta$color_group, levels = names(adj_bulk_colors))
  meta
}
meta_s41 <- add_color_group(meta_s41)
meta_s42 <- add_color_group(meta_s42)

s41_set <- s40_set
s41_set$subtitle_n <- n_subtitle(meta_s41, s41_set$condition_levels,
                                 s41_set$condition_labels)

s42_set <- s40_set
s42_set$condition_labels <- c("Normal", "Early DKD", "Adv. DKD", "AKI", "CKD\n(DKD only)")
s42_set$subtitle_n <- n_subtitle(meta_s42, s42_set$condition_levels,
                                 s42_set$condition_labels)

adj_box <- function(y_col, y_label, title_str, df, S) {
  cmp         <- S$comparisons_main
  ymax        <- max(df[[y_col]], na.rm = TRUE)
  ymin        <- min(df[[y_col]], na.rm = TRUE)
  rng         <- ymax - ymin
  upper_extra <- (0.08 * length(cmp) + 0.10) * rng
  ggplot(df, aes(x = condition, y = .data[[y_col]])) +
    geom_boxplot(aes(fill = condition), outlier.shape = NA,
                 alpha = 0.45, width = 0.55) +
    geom_jitter(aes(color = color_group, shape = cohort),
                position = jitter_pos(0.12), size = 2.4, alpha = 0.95) +
    bh_signif(df, "condition", y_col, cmp) +
    scale_fill_manual(values = S$condition_colors, guide = "none") +
    scale_color_manual(values = adj_bulk_colors, name = NULL,
                       guide = guide_legend(ncol = 1)) +
    scale_shape_manual(values = S$cohort_shapes, name = "Cohort") +
    scale_x_discrete(labels = setNames(S$condition_labels, S$condition_levels)) +
    scale_y_continuous(limits = c(ymin - 0.05 * rng, ymax + upper_extra)) +
    labs(x = NULL, y = y_label, title = title_str, subtitle = S$subtitle_n) +
    theme_box() +
    theme(legend.position = "right",
          legend.text  = element_text(size = 9),
          legend.title = element_text(size = 10))
}

adj_scatter <- function(x_col, y_col, x_label, y_label, title_str, df, S) {
  ct <- spearman_ci(df[[x_col]], df[[y_col]])
  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]],
                 color = color_group, shape = cohort)) +
    geom_point(size = 2.8, alpha = 0.95) +
    geom_smooth(method = "lm", se = FALSE, color = "black",
                linewidth = 0.8, inherit.aes = FALSE,
                aes(x = .data[[x_col]], y = .data[[y_col]])) +
    scale_color_manual(values = adj_bulk_colors, name = NULL,
                       guide = guide_legend(ncol = 1)) +
    scale_shape_manual(values = S$cohort_shapes, name = "Cohort") +
    annotate("text", x = Inf, y = Inf,
             label = ct$label, hjust = 1.05, vjust = 1.2, size = 3.8) +
    labs(x = x_label, y = y_label, title = title_str, subtitle = S$subtitle_n) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "right",
          legend.text  = element_text(size = 9),
          legend.title = element_text(size = 10))
}

s42_box     <- function(...) adj_box(..., df = meta_s42, S = s42_set)
s41_scatter <- function(...) adj_scatter(..., df = meta_s41, S = s41_set)
s42_scatter <- function(...) adj_scatter(..., df = meta_s42, S = s42_set)

# . S50_KPMP_scRNA_loadOnce.R ----
s50_set <- list(
  condition_levels = c("Healthy Reference Tissue", "Acute Kidney Injury",
                       "Chronic Kidney Disease"),
  condition_short  = c("Healthy", "AKI", "CKD"),
  condition_colors = c("Healthy Reference Tissue" = "#2166ac",
                       "Acute Kidney Injury"      = "#e07b39",
                       "Chronic Kidney Disease"   = "#756bb1"),
  comparisons_main = list(c("Healthy Reference Tissue", "Acute Kidney Injury"),
                          c("Healthy Reference Tissue", "Chronic Kidney Disease"))
)

make_subtitle_sc <- function(pb, S, unit = "samples", cond_col = "EnrollmentCategory") {
  cond_vals     <- as.character(pb[[cond_col]])
  active_levels <- S$condition_levels[S$condition_levels %in% unique(cond_vals)]
  active_short  <- S$condition_short[match(active_levels, S$condition_levels)]
  sprintf("n %s: %s", unit,
    paste(mapply(function(lv, sh) {
      sprintf("%s=%d", gsub("\n", " ", sh), sum(cond_vals == lv, na.rm = TRUE))
    }, active_levels, active_short), collapse = ", "))
}

s50_pt_box <- function(S = s50_set) {
  ggplot(pb_pt_s50, aes(x = EnrollmentCategory, y = pyrim_score,
                        fill = EnrollmentCategory)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.45, width = 0.5) +
    geom_jitter(position = jitter_pos(0.12), size = 2.2, alpha = 0.85) +
    bh_signif(pb_pt_s50, "EnrollmentCategory", "pyrim_score", S$comparisons_main) +
    scale_fill_manual(values = S$condition_colors, guide = "none") +
    scale_x_discrete(labels = setNames(S$condition_short, S$condition_levels), drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.30))) +
    labs(x = NULL, y = "Pyrimidine module score (mean expression)",
         title = "KPMP scRNA — PT pyrimidine score (pseudobulk by sample)",
         subtitle = make_subtitle_sc(pb_pt_s50, S)) +
    theme_box()
}

s50_fib_box <- function(y_col, y_label, title_str, df = pb_fib_s50, S = s50_set) {
  ymax <- max(df[[y_col]], na.rm = TRUE)
  ymin <- min(df[[y_col]], na.rm = TRUE)
  rng  <- ymax - ymin
  if (!is.finite(rng) || rng == 0) rng <- max(abs(c(ymax, ymin)), 1)
  ggplot(df, aes(x = EnrollmentCategory, y = .data[[y_col]],
                 fill = EnrollmentCategory)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.45, width = 0.5) +
    geom_jitter(position = jitter_pos(0.12), size = 2.2, alpha = 0.85) +
    bh_signif(df, "EnrollmentCategory", y_col, S$comparisons_main) +
    scale_fill_manual(values = S$condition_colors, guide = "none") +
    scale_x_discrete(labels = setNames(S$condition_short, S$condition_levels), drop = FALSE) +
    scale_y_continuous(limits = c(ymin - 0.05 * rng, ymax + 0.26 * rng)) +
    labs(x = NULL, y = y_label, title = title_str,
         subtitle = make_subtitle_sc(df, S)) +
    theme_box()
}

# Used by S50, S52 (EnrollmentCategory) and S60, S62 (enrollment)
sc_cor <- function(df, x_col, y_col, x_lab, y_lab, title_str, colors,
                   color_col = "EnrollmentCategory",
                   subtitle = NULL, drop = TRUE) {
  ct <- spearman_ci(df[[x_col]], df[[y_col]])
  if (is.null(subtitle)) subtitle <- sprintf("n = %d matched samples", nrow(df))
  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]],
                 color = .data[[color_col]])) +
    geom_point(size = 2.8, alpha = 0.85) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8,
                inherit.aes = FALSE, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    annotate("text", x = Inf, y = Inf,
             label = ct$label, hjust = 1.05, vjust = 1.2, size = 3.8) +
    scale_color_manual(values = colors, name = NULL, drop = drop) +
    labs(x = x_lab, y = y_lab, title = title_str, subtitle = subtitle) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold"), legend.position = "right")
}

# . S51_scRNA_AdjudicatedColors.R ----
s51_adj_colors <- c(
  "Diabetic Kidney Disease"            = "#1f78b4",
  "Hypertensive Kidney Disease"        = "#8e44ad",
  "Focal Segmental Glomerulosclerosis" = "#27ae60",
  "Lupus Nephritis"                    = "#e74c3c",
  "Acute Tubular Injury"               = "#f39c12",
  "Acute Interstitial Nephritis"       = "#e67e22",
  "IgA Nephropathy"                    = "#16a085",
  "Cannot be determined"               = "#95a5a6",
  "Other"                              = "#7f8c8d"
)
s51_non_ckd_colors <- c("Healthy Reference Tissue" = "#2166ac",
                        "Acute Kidney Injury"      = "#e07b39")
s51_all_colors <- c(s51_adj_colors, s51_non_ckd_colors)

pb_fib_s51 <- pb_fib_s50
pb_fib_s51$point_color <- ifelse(
  as.character(pb_fib_s51$EnrollmentCategory) == "Chronic Kidney Disease",
  ifelse(!is.na(pb_fib_s51$adjudicated) & pb_fib_s51$adjudicated %in% names(s51_adj_colors),
         pb_fib_s51$adjudicated, "Other"),
  as.character(pb_fib_s51$EnrollmentCategory)
)
pb_fib_s51$point_color <- factor(pb_fib_s51$point_color,
                                 levels = names(s51_all_colors))

s51_cor_adj <- function(x_col, y_col, x_lab, y_lab, title_str, df = pb_fib_s51) {
  ct <- spearman_ci(df[[x_col]], df[[y_col]])
  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]], color = point_color)) +
    geom_point(size = 2.8, alpha = 0.9) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8,
                inherit.aes = FALSE, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    annotate("text", x = Inf, y = Inf,
             label = ct$label, hjust = 1.05, vjust = 1.2, size = 3.8) +
    scale_color_manual(values = s51_all_colors, name = "Adjudicated", drop = TRUE,
                       guide = guide_legend(ncol = 1, override.aes = list(size = 3))) +
    labs(x = x_lab, y = y_lab, title = title_str,
         subtitle = sprintf("n = %d stroma samples", ct$n)) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "right",
          legend.text  = element_text(size = 9),
          legend.title = element_text(size = 10, face = "bold"))
}

# . S52_scRNA_DKD_Focused.R ----
s52_set <- s50_set
s52_set$condition_short <- c("Healthy", "AKI", "CKD\n(DKD only)")

# S52 plot_box (also used by S62 with the enrollment column)
dkd_box <- function(df, y_col, y_label, title_str, subtitle_str, S,
                    cond_col = "EnrollmentCategory") {
  ymax <- max(df[[y_col]], na.rm = TRUE)
  ymin <- min(df[[y_col]], na.rm = TRUE)
  rng  <- ymax - ymin
  if (!is.finite(rng) || rng == 0) {
    rng <- max(abs(c(ymax, ymin)), na.rm = TRUE)
    if (!is.finite(rng) || rng == 0) rng <- 1
  }
  present <- unique(as.character(df[[cond_col]]))
  cmps    <- Filter(function(p) all(p %in% present), S$comparisons_main)
  ggplot(df, aes(x = .data[[cond_col]], y = .data[[y_col]],
                 fill = .data[[cond_col]])) +
    geom_boxplot(outlier.shape = NA, alpha = 0.45, width = 0.5) +
    geom_jitter(position = jitter_pos(0.12), size = 2.2, alpha = 0.85) +
    {if (length(cmps) > 0)
      bh_signif(df, cond_col, y_col, cmps)} +
    scale_fill_manual(values = S$condition_colors, guide = "none") +
    scale_x_discrete(labels = setNames(S$condition_short, S$condition_levels), drop = FALSE) +
    scale_y_continuous(limits = c(ymin - 0.05 * rng, ymax + 0.26 * rng)) +
    labs(x = NULL, y = y_label, title = title_str, subtitle = subtitle_str) +
    theme_box()
}

s52_fib_box <- function(y_col, y_label, title_str) {
  dkd_box(pb_fib_s52, y_col, y_label, title_str,
          make_subtitle_sc(pb_fib_s52, s52_set), s52_set)
}

# . S60_KPMP_snRNA_loadOnce.R ----
s60_set <- list(
  condition_levels = c("Healthy Reference", "AKI", "CKD"),
  condition_short  = c("Healthy", "AKI", "CKD"),
  condition_colors = c("Healthy Reference" = "#2166ac", "AKI" = "#e07b39",
                       "CKD" = "#756bb1"),
  comparisons_main = list(c("Healthy Reference", "AKI"),
                          c("Healthy Reference", "CKD"))
)

s60_pt_box <- function(S = s60_set) {
  ggplot(pb_pt_s60, aes(x = enrollment, y = pyrim_score, fill = enrollment)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.45, width = 0.5) +
    geom_jitter(position = jitter_pos(0.12), size = 2.2, alpha = 0.85) +
    bh_signif(pb_pt_s60, "enrollment", "pyrim_score", S$comparisons_main) +
    scale_fill_manual(values = S$condition_colors, guide = "none") +
    scale_x_discrete(labels = setNames(S$condition_short, S$condition_levels)) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.30))) +
    labs(x = NULL, y = "Pyrimidine module score (mean expression)",
         title = "KPMP snRNA — PT pyrimidine score (pseudobulk by patient)",
         subtitle = make_subtitle_sc(pb_pt_s60, S, "patients", "enrollment")) +
    theme_box()
}

s60_fib_box <- function(y_col, y_label, title_str, df = pb_fib_s60, S = s60_set) {
  ymax <- max(df[[y_col]], na.rm = TRUE)
  ymin <- min(df[[y_col]], na.rm = TRUE)
  rng  <- ymax - ymin
  ggplot(df, aes(x = enrollment, y = .data[[y_col]], fill = enrollment)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.45, width = 0.5) +
    geom_jitter(position = jitter_pos(0.12), size = 2.2, alpha = 0.85) +
    bh_signif(df, "enrollment", y_col, S$comparisons_main) +
    scale_fill_manual(values = S$condition_colors, guide = "none") +
    scale_x_discrete(labels = setNames(S$condition_short, S$condition_levels)) +
    scale_y_continuous(limits = c(ymin - 0.05 * rng, ymax + 0.26 * rng)) +
    labs(x = NULL, y = y_label, title = title_str,
         subtitle = make_subtitle_sc(df, S, "patients", "enrollment")) +
    theme_box()
}

# . S62_snRNA_DKD_Focused_Fixed.R ----
s62_set <- s60_set
s62_set$condition_short <- c("Healthy", "AKI", "CKD\n(DKD only)")

s62_fib_box <- function(y_col, y_label, title_str) {
  dkd_box(pb_fib_s62, y_col, y_label, title_str,
          make_subtitle_sc(pb_fib_s62, s62_set, "patients", "enrollment"),
          s62_set, cond_col = "enrollment")
}

# ============================================================================#
# FIGURE 1 — Pyrimidine score ----
# ============================================================================#

# (a) S40 Fig_pyrimidine_score.png
f1a <- s40_box("pyrim_score", "Pyrimidine module score (log2)",
               "Combined — Pyrimidine score")
# (b) S42 Fig_pyrimidine_score.png
f1b <- s42_box("pyrim_score", "Pyrimidine module score (log2)",
               "Combined (DKD-focused) — Pyrimidine score")
# (c) S50 Fig_PT_pyrimidine_score.png
f1c <- s50_pt_box()
# (d) S52 Fig_PT_pyrimidine_score.png
f1d <- dkd_box(pb_pt_s52, "pyrim_score",
               "Pyrimidine module score (mean expression)",
               "KPMP scRNA (DKD-focused) — PT pyrimidine score",
               make_subtitle_sc(pb_pt_s52, s52_set), s52_set)
# (e) S60 Fig_PT_pyrimidine_score.png
f1e <- s60_pt_box()
# (f) S62 Fig_PT_pyrimidine_score.png
f1f <- dkd_box(pb_pt_s62, "pyrim_score",
               "Pyrimidine module score (mean expression)",
               "KPMP snRNA (DKD-focused) — PT pyrimidine score",
               make_subtitle_sc(pb_pt_s62, s62_set, "patients", "enrollment"),
               s62_set, cond_col = "enrollment")

save_figure(list(f1a, f1b, f1c, f1d, f1e, f1f), "Figure_1",
            width = 15, heights = c(6.2, 5.6, 5.6))

# ============================================================================#
# FIGURE 2 — P2RY6 ----
# ============================================================================#

# (a) S40 Fig_P2RY6.png
f2a <- s40_box("P2RY6", "P2RY6 expression (log2)", "Combined — P2RY6")
# (b) S42 Fig_P2RY6.png
f2b <- s42_box("P2RY6", "P2RY6 expression (log2)", "Combined (DKD-focused) — P2RY6")
# (c) S50 Fig_FB_P2RY6.png
#     (Figures.pptx used Fig_FB_P2RY6_pct.png, % P2RY6+ cells; replaced by
#      mean expression for consistency with panels (d)-(f), v5 revision)
f2c <- s50_fib_box("mean_P2RY6", "P2RY6 mean expression (pseudobulk)",
                   "KPMP scRNA — P2RY6 in stroma cells")
# (d) S52 Fig_FB_P2RY6.png
f2d <- s52_fib_box("mean_P2RY6", "P2RY6 mean expression (pseudobulk)",
                   "KPMP scRNA (DKD-focused) — P2RY6 in stroma")
# (e) S60 Fig_FB_P2RY6.png
f2e <- s60_fib_box("mean_P2RY6", "P2RY6 mean expression (pseudobulk)",
                   "KPMP snRNA — P2RY6 in fibroblasts")
# (f) S62 Fig_FB_P2RY6.png
f2f <- s62_fib_box("mean_P2RY6", "P2RY6 mean expression (pseudobulk)",
                   "KPMP snRNA (DKD-focused) — P2RY6 in fibroblasts")

save_figure(list(f2a, f2b, f2c, f2d, f2e, f2f), "Figure_2",
            width = 15, heights = c(6.2, 5.6, 5.6))

# ============================================================================#
# FIGURE 3 — Calcium signaling score ----
# ============================================================================#

# (a) S40 Fig_calcium_score.png
f3a <- s40_box("calcium_score", "Calcium signaling score (log2)",
               "Combined — Calcium score")
# (b) S42 Fig_calcium_score.png
f3b <- s42_box("calcium_score", "Calcium signaling score (log2)",
               "Combined (DKD-focused) — Calcium score")
# (c) S50 Fig_FB_calcium.png
f3c <- s50_fib_box("mean_calcium_score", "Calcium signaling score (mean expression)",
                   "KPMP scRNA — Calcium score in stroma")
# (d) S52 Fig_FB_calcium.png
f3d <- s52_fib_box("mean_calcium_score", "Calcium signaling score (mean expression)",
                   "KPMP scRNA (DKD-focused) — Calcium score in stroma")
# (e) S60 Fig_FB_calcium.png
f3e <- s60_fib_box("mean_calcium_score", "Calcium signaling score (mean expression)",
                   "KPMP snRNA — Calcium score in fibroblasts")
# (f) S62 Fig_FB_calcium.png
f3f <- s62_fib_box("mean_calcium_score", "Calcium signaling score (mean expression)",
                   "KPMP snRNA (DKD-focused) — Calcium score in FBs")

save_figure(list(f3a, f3b, f3c, f3d, f3e, f3f), "Figure_3",
            width = 15, heights = c(6.2, 5.6, 5.6))

# ============================================================================#
# FIGURE 4 — PDGFRB ----
# ============================================================================#

# (a) S40 Fig_PDGFRB.png
f4a <- s40_box("PDGFRB", "PDGFRB expression (log2)", "Combined — PDGFRB")
# (b) S42 Fig_PDGFRB.png
f4b <- s42_box("PDGFRB", "PDGFRB expression (log2)", "Combined (DKD-focused) — PDGFRB")
# (c) S50 Fig_FB_PDGFRB.png
f4c <- s50_fib_box("mean_PDGFRB", "PDGFRB mean expression (pseudobulk)",
                   "KPMP scRNA — PDGFRB in stroma cells")
# (d) S52 Fig_FB_PDGFRB.png
f4d <- s52_fib_box("mean_PDGFRB", "PDGFRB mean expression (pseudobulk)",
                   "KPMP scRNA (DKD-focused) — PDGFRB in stroma")
# (e) S60 Fig_FB_PDGFRB.png
f4e <- s60_fib_box("mean_PDGFRB", "PDGFRB mean expression (pseudobulk)",
                   "KPMP snRNA — PDGFRB in fibroblasts")
# (f) S62 Fig_FB_PDGFRB.png
f4f <- s62_fib_box("mean_PDGFRB", "PDGFRB mean expression (pseudobulk)",
                   "KPMP snRNA (DKD-focused) — PDGFRB in fibroblasts")

save_figure(list(f4a, f4b, f4c, f4d, f4e, f4f), "Figure_4",
            width = 15, heights = c(6.2, 5.6, 5.6))

# ============================================================================#
# FIGURE 5 — ECM score ----
# ============================================================================#

# (a) S40 Fig_ECM_score.png
f5a <- s40_box("ecm_score", "ECM / fibrotic score (log2)", "Combined — ECM score")
# (b) S42 Fig_ECM_score.png
f5b <- s42_box("ecm_score", "ECM / fibrotic score (log2)",
               "Combined (DKD-focused) — ECM score")
# (c) S50 Fig_FB_ECM.png
f5c <- s50_fib_box("mean_ecm_score", "ECM score (mean expression)",
                   "KPMP scRNA — ECM score in stroma")
# (d) S52 Fig_FB_ECM.png
f5d <- s52_fib_box("mean_ecm_score", "ECM score (mean expression)",
                   "KPMP scRNA (DKD-focused) — ECM score in stroma")
# (e) S60 Fig_FB_ECM.png
f5e <- s60_fib_box("mean_ecm_score", "ECM score (mean expression)",
                   "KPMP snRNA — ECM score in fibroblasts")
# (f) S62 Fig_FB_ECM.png
f5f <- s62_fib_box("mean_ecm_score", "ECM score (mean expression)",
                   "KPMP snRNA (DKD-focused) — ECM score in FBs")

save_figure(list(f5a, f5b, f5c, f5d, f5e, f5f), "Figure_5",
            width = 15, heights = c(6.2, 5.6, 5.6))

# ============================================================================#
# FIGURE 6 — PT pyrimidine score vs ECM score ----
# ============================================================================#

# (a) S41 Fig_cor_Pyrimidine_ECM.png
f6a <- s41_scatter("pyrim_score", "ecm_score", "Pyrimidine score", "ECM score",
                   "Combined — Pyrimidine vs ECM")
# (b) S42 Fig_cor_Pyrimidine_ECM.png
f6b <- s42_scatter("pyrim_score", "ecm_score", "Pyrimidine score", "ECM score",
                   "Combined (DKD-focused) — Pyrimidine vs ECM")
# (c) S50 Fig_Cor_PT_Pyrimidine_vs_FB_ECM.png
f6c <- sc_cor(pb_cross_s50, "pyrim_score", "mean_ecm_score",
              "PT pyrimidine score (pseudobulk)", "Stroma ECM score (pseudobulk)",
              "KPMP scRNA — PT pyrimidine vs Stroma ECM",
              colors = s50_set$condition_colors)
# (d) S52 Fig_Cor_PT_Pyrimidine_vs_FB_ECM.png
f6d <- sc_cor(pb_cross_s52, "pyrim_score", "mean_ecm_score",
              "PT pyrimidine score (pseudobulk)", "Stroma ECM score (pseudobulk)",
              "KPMP scRNA (DKD-focused) — PT pyrimidine vs Stroma ECM",
              colors = s52_set$condition_colors, drop = FALSE)
# (e) S60 Fig_Cor_PT_Pyrimidine_vs_FB_ECM.png
f6e <- sc_cor(pb_cross_s60, "pyrim_score", "mean_ecm_score",
              "PT pyrimidine score (pseudobulk)", "FB ECM score (pseudobulk)",
              "KPMP snRNA — PT pyrimidine vs FB ECM",
              colors = s60_set$condition_colors, color_col = "enrollment",
              subtitle = sprintf("n = %d matched patients", nrow(pb_cross_s60)))
# (f) S62 Fig_Cor_PT_Pyrimidine_vs_FB_ECM.png
f6f <- sc_cor(pb_cross_s62, "pyrim_score", "mean_ecm_score",
              "PT pyrimidine score (pseudobulk)", "FB ECM score (pseudobulk)",
              "KPMP snRNA (DKD-focused) — PT pyrimidine vs FB ECM",
              colors = s62_set$condition_colors, color_col = "enrollment",
              subtitle = sprintf("n = %d matched patients", nrow(pb_cross_s62)),
              drop = FALSE)

save_figure(list(f6a, f6b, f6c, f6d, f6e, f6f), "Figure_6",
            width = 16, heights = c(5.6, 5.2, 5.2))

# ============================================================================#
# SUPPLEMENTARY FIGURE S1 — PCA before / after ComBat ----
# ============================================================================#
# Code from S40_Combined_BulkRNA.R §1-§4 (merge of S20/S30 matrices and
# PCA QC plots); the ComBat-corrected matrix is taken from S40_corrected.rds.

s20 <- readRDS(s20_rds)
fan_mat  <- s20$expr_mat
kpmp_mat <- s30$expr_mat
common_genes <- intersect(rownames(fan_mat), rownames(kpmp_mat))
combined_mat <- cbind(fan_mat[common_genes, ], kpmp_mat[common_genes, ])
corrected_mat <- s40$corrected_mat
stopifnot(identical(colnames(combined_mat), meta_s40$sample_id),
          identical(colnames(corrected_mat), meta_s40$sample_id))

top_var_genes <- function(mat, n = 2000) {
  rv <- apply(mat, 1, var)
  head(order(rv, decreasing = TRUE), n)
}

pca_plot <- function(mat, meta, title_str, S = s40_set) {
  idx <- top_var_genes(mat, 2000)
  pca <- prcomp(t(mat[idx, ]), scale. = TRUE)
  pve <- round(summary(pca)$importance[2, 1:2] * 100, 1)
  df  <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                    condition = meta$condition, cohort = meta$cohort)
  ggplot(df, aes(x = PC1, y = PC2, color = condition, shape = cohort)) +
    geom_point(size = 3, alpha = 0.9) +
    scale_color_manual(values = S$condition_colors) +
    scale_shape_manual(values = S$cohort_shapes) +
    labs(x = sprintf("PC1 (%.1f%%)", pve[1]),
         y = sprintf("PC2 (%.1f%%)", pve[2]),
         title = title_str, color = NULL, shape = NULL) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))
}

# (a) S40 QC/PCA_pre_ComBat.png
fs1a <- pca_plot(combined_mat,  meta_s40, "PCA — pre-ComBat (top 2000 variable genes)")
# (b) S40 QC/PCA_post_ComBat.png
fs1b <- pca_plot(corrected_mat, meta_s40, "PCA — post-ComBat (top 2000 variable genes)")

save_figure(list(fs1a, fs1b), "Figure_S1", width = 13, heights = 5.8)

# ============================================================================#
# SUPPLEMENTARY FIGURE S2 — P2RY6 vs Calcium score ----
# ============================================================================#

# (a) S40 Fig_cor_P2RY6_Calcium.png
fs2a <- s40_scatter("P2RY6", "calcium_score", "P2RY6 (log2)", "Calcium score",
                    "Combined — P2RY6 vs Calcium")
# (b) S42 Fig_cor_P2RY6_Calcium.png
fs2b <- s42_scatter("P2RY6", "calcium_score", "P2RY6 (log2)", "Calcium score",
                    "Combined (DKD-focused) — P2RY6 vs Calcium")
# (c) S51 Fig_Cor_FB_P2RY6_vs_Calcium.png
fs2c <- s51_cor_adj("mean_P2RY6", "mean_calcium_score",
                    "Stroma P2RY6 (pseudobulk)", "Stroma calcium score (pseudobulk)",
                    "KPMP scRNA — Stroma P2RY6 vs Calcium (adjudicated)")
# (d) S52 Fig_Cor_FB_P2RY6_vs_Calcium.png
fs2d <- sc_cor(pb_fib_s52, "mean_P2RY6", "mean_calcium_score",
               "Stroma P2RY6 (pseudobulk)", "Stroma calcium score (pseudobulk)",
               "KPMP scRNA (DKD-focused) — Stroma P2RY6 vs Calcium",
               colors = s52_set$condition_colors, drop = FALSE,
               subtitle = sprintf("n = %d stroma samples", nrow(pb_fib_s52)))
# (e) S60 Fig_Cor_FB_P2RY6_vs_Calcium.png
fs2e <- sc_cor(pb_fib_s60, "mean_P2RY6", "mean_calcium_score",
               "FB P2RY6 (pseudobulk)", "FB calcium score (pseudobulk)",
               "KPMP snRNA — FB P2RY6 vs Calcium",
               colors = s60_set$condition_colors, color_col = "enrollment",
               subtitle = sprintf("n = %d fibroblast patients", nrow(pb_fib_s60)))
# (f) S62 Fig_Cor_FB_P2RY6_vs_Calcium.png
fs2f <- sc_cor(pb_fib_s62, "mean_P2RY6", "mean_calcium_score",
               "FB P2RY6 (pseudobulk)", "FB calcium score (pseudobulk)",
               "KPMP snRNA (DKD-focused) — FB P2RY6 vs Calcium",
               colors = s62_set$condition_colors, color_col = "enrollment",
               subtitle = sprintf("n = %d fibroblast patients", nrow(pb_fib_s62)),
               drop = FALSE)

save_figure(list(fs2a, fs2b, fs2c, fs2d, fs2e, fs2f), "Figure_S2",
            width = 16, heights = c(5.6, 5.2, 5.2))

# ============================================================================#
# SUPPLEMENTARY FIGURE S3 — P2RY6 vs PDGFRB ----
# ============================================================================#
# Axes swapped relative to the source panels (P2RY6 on x, PDGFRB on y), as in
# Figures S2 and S4. Spearman rho is symmetric; only the linear fit line changes.

# (a) S40 Fig_cor_PDGFRB_P2RY6.png
fs3a <- s40_scatter("P2RY6", "PDGFRB", "P2RY6 (log2)", "PDGFRB (log2)",
                    "Combined — P2RY6 vs PDGFRB")
# (b) S42 Fig_cor_PDGFRB_P2RY6.png
fs3b <- s42_scatter("P2RY6", "PDGFRB", "P2RY6 (log2)", "PDGFRB (log2)",
                    "Combined (DKD-focused) — P2RY6 vs PDGFRB")
# (c) S51 Fig_Cor_FB_PDGFRB_vs_P2RY6.png
fs3c <- s51_cor_adj("mean_P2RY6", "mean_PDGFRB",
                    "Stroma P2RY6 (pseudobulk)", "Stroma PDGFRB (pseudobulk)",
                    "KPMP scRNA — Stroma P2RY6 vs PDGFRB (adjudicated)")
# (d) S52 Fig_Cor_FB_PDGFRB_vs_P2RY6.png
fs3d <- sc_cor(pb_fib_s52, "mean_P2RY6", "mean_PDGFRB",
               "Stroma P2RY6 (pseudobulk)", "Stroma PDGFRB (pseudobulk)",
               "KPMP scRNA (DKD-focused) — Stroma P2RY6 vs PDGFRB",
               colors = s52_set$condition_colors, drop = FALSE,
               subtitle = sprintf("n = %d stroma samples", nrow(pb_fib_s52)))
# (e) S60 Fig_Cor_FB_PDGFRB_vs_P2RY6.png
fs3e <- sc_cor(pb_fib_s60, "mean_P2RY6", "mean_PDGFRB",
               "FB P2RY6 (pseudobulk)", "FB PDGFRB (pseudobulk)",
               "KPMP snRNA — FB P2RY6 vs PDGFRB",
               colors = s60_set$condition_colors, color_col = "enrollment",
               subtitle = sprintf("n = %d fibroblast patients", nrow(pb_fib_s60)))
# (f) S62 Fig_Cor_FB_PDGFRB_vs_P2RY6.png
fs3f <- sc_cor(pb_fib_s62, "mean_P2RY6", "mean_PDGFRB",
               "FB P2RY6 (pseudobulk)", "FB PDGFRB (pseudobulk)",
               "KPMP snRNA (DKD-focused) — FB P2RY6 vs PDGFRB",
               colors = s62_set$condition_colors, color_col = "enrollment",
               subtitle = sprintf("n = %d fibroblast patients", nrow(pb_fib_s62)),
               drop = FALSE)

save_figure(list(fs3a, fs3b, fs3c, fs3d, fs3e, fs3f), "Figure_S3",
            width = 16, heights = c(5.6, 5.2, 5.2))

# ============================================================================#
# SUPPLEMENTARY FIGURE S4 — P2RY6 vs ECM score ----
# ============================================================================#
# Not present in Figures.pptx (slide duplicated Suppl. Figure 2); assembled
# from the P2RY6 vs ECM outputs following the layout of Figures S2-S3.

# (a) S40 Fig_cor_P2RY6_ECM.png
fs4a <- s40_scatter("P2RY6", "ecm_score", "P2RY6 (log2)", "ECM score",
                    "Combined — P2RY6 vs ECM")
# (b) S42 Fig_cor_P2RY6_ECM.png
fs4b <- s42_scatter("P2RY6", "ecm_score", "P2RY6 (log2)", "ECM score",
                    "Combined (DKD-focused) — P2RY6 vs ECM")
# (c) S51 Fig_Cor_FB_P2RY6_vs_ECM.png
fs4c <- s51_cor_adj("mean_P2RY6", "mean_ecm_score",
                    "Stroma P2RY6 (pseudobulk)", "Stroma ECM score (pseudobulk)",
                    "KPMP scRNA — Stroma P2RY6 vs ECM (adjudicated)")
# (d) S52 Fig_Cor_FB_P2RY6_vs_ECM.png
fs4d <- sc_cor(pb_fib_s52, "mean_P2RY6", "mean_ecm_score",
               "Stroma P2RY6 (pseudobulk)", "Stroma ECM score (pseudobulk)",
               "KPMP scRNA (DKD-focused) — Stroma P2RY6 vs ECM",
               colors = s52_set$condition_colors, drop = FALSE,
               subtitle = sprintf("n = %d stroma samples", nrow(pb_fib_s52)))
# (e) S60 Fig_Cor_FB_P2RY6_vs_ECM.png
fs4e <- sc_cor(pb_fib_s60, "mean_P2RY6", "mean_ecm_score",
               "FB P2RY6 (pseudobulk)", "FB ECM score (pseudobulk)",
               "KPMP snRNA — FB P2RY6 vs ECM",
               colors = s60_set$condition_colors, color_col = "enrollment",
               subtitle = sprintf("n = %d fibroblast patients", nrow(pb_fib_s60)))
# (f) S62 Fig_Cor_FB_P2RY6_vs_ECM.png
fs4f <- sc_cor(pb_fib_s62, "mean_P2RY6", "mean_ecm_score",
               "FB P2RY6 (pseudobulk)", "FB ECM score (pseudobulk)",
               "KPMP snRNA (DKD-focused) — FB P2RY6 vs ECM",
               colors = s62_set$condition_colors, color_col = "enrollment",
               subtitle = sprintf("n = %d fibroblast patients", nrow(pb_fib_s62)),
               drop = FALSE)

save_figure(list(fs4a, fs4b, fs4c, fs4d, fs4e, fs4f), "Figure_S4",
            width = 16, heights = c(5.6, 5.2, 5.2))

cat("\nS80 complete. Figures saved to:", out_dir, "\n")
