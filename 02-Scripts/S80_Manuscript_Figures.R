# ============================================================================#
# S80_Manuscript_Figures.R
# Assembles the multi-panel figures of the revised manuscript (MDPI Life,
# Submission 2) from the panel code of the analysis scripts S20-S62.
#
# Replaces the manually collated slides in
#   07-Manuscript/Submission 2 - MDPI-LIFE/Figures.pptx
#
# Panel -> source mapping (Revision 1):
#
#   Figure     top left  (a)        (b)        (c)        (d)        (e)
#   ---------  --------  ---------  ---------  ---------  ---------  ---------
#   Figure 1   S20 (*1)  S20        S50        S52        S60        S62      Pyrimidine score
#   Figure 2   empty     S20        S50        S52        S60        S62      P2RY6
#   Figure 3   empty     S20        S50        S52        S60        S62      Calcium score
#   Figure 4   empty     S20        S50        S52        S60        S62      PDGFRB
#   Figure 5   empty     S20        S50        S52        S60        S62      ECM score
#   Figure 6   empty     S20        S50        S52        S60        S62      Pyrimidine vs ECM
#   Figure S1  empty     S20        S51        S52        S60        S62      P2RY6 vs Calcium
#   Figure S2  empty     S20        S51        S52        S60        S62      P2RY6 vs PDGFRB (**)
#   Figure S3  empty     S20        S51        S52        S60        S62      P2RY6 vs ECM (*)
#
#   (a) Fan DKD bulk RNA-seq only (Normal, Early DKD, Advanced DKD; values
#       of S20 without batch correction), centred in the top row (A. Figurek,
#       after ms_rev1_v32; previously top right with the top left empty).
#       Revision 1 (R1.10, R2.02, R2.03): the KPMP
#       bulk data and the ComBat-corrected combined Fan + KPMP panels
#       (S40-S42) and the PCA figure (former Figure S1) were removed, because
#       the AKI/CKD profiles of KPMP cannot be compared reliably with the
#       Normal samples of another cohort (sensitivity analysis: S43). S30 and
#       S40-S43 were removed from the repository (see git history).
#   Figure 2b shows mean P2RY6 expression like the other P2RY6 panels
#       (Figures.pptx showed % P2RY6+ stroma cells; changed in the v5 revision).
#   (*) In Figures.pptx the slide "Suppl. Figure 3" is an exact copy of
#       "Suppl. Figure 2". The manuscript text describes this figure as the
#       P2RY6 vs ECM correlations, which are assembled here following the
#       layout of the other correlation figures.
#   (*1) Revision 1 (comment of A. Figurek on ms_rev1_v31_AF): Figure 1 has
#       no empty position; the bulk PT marker score (former Supplementary
#       Figure S5a, S86) is panel (a) and the panels below are (b)-(f).
#   Supplementary figures were renumbered in order of their citation in the
#   Results (same revision): Figures S1, S2, S3 of this script are Figures
#   S2, S3, S4 of the manuscript; the output file names are unchanged.
#   Revision 1 (comments of A. Figurek on ms_rev1_v42_AF): Figures 2-5 keep
#   only the snRNA-seq fibroblast panels (d, e -> a, b); their bulk and
#   scRNA-seq panels (a-c) form new Supplementary Figures S2, S3, S5 and S7.
#   The former Figure 6 is Supplementary Figure S10, and the new Figure 6 shows
#   the within-AKI and within-CKD correlations of its snRNA-seq panel (S10d).
#   With the new figures, the supplementary figures are renumbered in order of
#   citation, and the output files are named by the manuscript numbers:
#     Figure_S2  bulk/scRNA-seq P2RY6        (former Figure 2a-c)
#     Figure_S3  bulk/scRNA-seq calcium      (former Figure 3a-c)
#     Figure_S4  P2RY6 vs calcium            (fs1; former Figure S2, file Figure_S1)
#     Figure_S5  bulk/scRNA-seq PDGFRB       (former Figure 4a-c)
#     Figure_S6  P2RY6 vs PDGFRB             (fs2; former Figure S3, file Figure_S2)
#     Figure_S7  bulk/scRNA-seq ECM          (former Figure 5a-c)
#     Figure_S8  P2RY6 vs ECM                (fs3; former Figure S4, file Figure_S3)
#     Figure_S10 PT pyrimidine vs ECM        (f6a-f6e; former Figure 6)
#   Figures S1 and S9 come from S86 and S84. The panel objects keep their names
#   (f2a-f2e, ..., f6a-f6e, fs1a-fs3e), which S81, S85 and S89 use.
#   (**) Axes swapped relative to S50-S62 (P2RY6 on x, PDGFRB on y) so that
#       P2RY6 is on the x-axis in all supplementary correlation figures.
#
# Inputs (no heavy re-computation; Seurat objects are not needed):
#   S20_normalised.rds, S50_pseudobulk.rds, S60_pseudobulk.rds
#
# Differences to the panels saved by the source scripts:
#   - significance stars show Benjamini-Hochberg-adjusted p-values (two-sided
#     Wilcoxon rank-sum test, adjusted over the comparisons of each panel);
#     the source scripts showed unadjusted p-values
#   - correlation panels (Figures 6, S1-S3) show Spearman p-values that are
#     BH-adjusted over the five panels of each figure (Revision 1, R1.07)
#   - jitter uses a fixed seed (identical points in PNG and PDF, reproducible)
#   - long titles/subtitles are wrapped instead of being cut at the panel edge
#   - extra top margin per panel to make room for the panel label
#
# Outputs: 04-Outputs/S80_Manuscript_Figures/Figure_<n>.png / .pdf (manuscript
#          numbers: Figures 1-6 and Supplementary Figures S2-S8, S10)
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
s50_rds <- "../04-Outputs/S50_KPMP_scRNA/S50_pseudobulk.rds"
s60_rds <- "../04-Outputs/S60_KPMP_snRNA/S60_pseudobulk.rds"
sc_clin <- "../../KPMP_scRNA/OpenAccessClinicalData/20251205_OpenAccessClinicalData.csv"
sn_clin <- "../../KPMP_snRNA/6d3047c2-588f-4ad4-8eba-f4c46a1e5c14_20251205_OpenAccessClinicalData/20251205_OpenAccessClinicalData.csv"
out_dir <- "../04-Outputs/S80_Manuscript_Figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

pt_markers <- c("LRP2", "SLC5A12", "SLC22A6", "PRODH2", "SLC5A2", "SLC22A8", "SLC7A8",
                "SLC34A1", "SLC5A10", "SLC5A11", "SLC7A13")

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

# Revision 1 (AIReview1 D6): p-values of the correlation labels in the same
# format as Table S6 (S85 fmt_p): < 0.0001; 2 significant digits below 0.01;
# otherwise 2 decimals.
format_p_s6 <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-4) return("< 0.0001")
  if (p < 0.01) return(sprintf("%.2g", p))
  sprintf("%.2f", p)
}

# Revision 1 (AIReview1 D6): the same short group labels (Healthy, AKI, CKD)
# in the legends of the correlation figures as on the axes of Figures 1-5.
legend_short <- c("Healthy Reference Tissue" = "Healthy", "Healthy Reference" = "Healthy",
                  "Acute Kidney Injury" = "AKI", "Chronic Kidney Disease" = "CKD")
short_labels <- function(x) {
  x <- as.character(x)
  ifelse(x %in% names(legend_short), unname(legend_short[x]), x)
}

# Revision 1 (AIReview1 D6): headroom above the data for the correlation label
cor_y_expand <- ggplot2::expansion(mult = c(0.05, 0.22))

# NA-safe version (S50/S51/S52/S62); gives identical results to the
# version without NA handling on complete data.
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

# Correlation labels with Benjamini-Hochberg-adjusted p-values ----
# Revision 1 (R1.07): the Spearman p-values of the panels of one correlation
# figure are BH-adjusted together; the label layer of each panel is rewritten.
# The same family is used in S85_Correlation_Stats.R (Table S4).
bh_cor_panels <- function(panels) {
  is_lab <- function(l) inherits(l$geom, "GeomText") && !is.null(l$aes_params$label)
  cts <- lapply(panels, function(p) {
    x <- rlang::eval_tidy(p$mapping$x, p$data); y <- rlang::eval_tidy(p$mapping$y, p$data)
    spearman_ci(x, y)
  })
  p_adj <- p.adjust(vapply(cts, `[[`, numeric(1), "p"), method = "BH")
  p_txt <- vapply(p_adj, function(p) {
    f <- format_p_s6(p); if (startsWith(f, "<")) f else paste("=", f)
  }, character(1))
  Map(function(p, ct, pt) {
    i <- which(vapply(p$layers, is_lab, logical(1)))
    stopifnot(length(i) == 1)
    p$layers[[i]]$aes_params$label <- if (is.na(ct$ci_lo)) {
      sprintf("rho = %.2f\nBH-adjusted p %s", ct$rho, pt)
    } else {
      sprintf("rho = %.2f [%.2f, %.2f]\nBH-adjusted p %s", ct$rho, ct$ci_lo, ct$ci_hi, pt)
    }
    p
  }, panels, cts, p_txt)
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

# Assemble a 3 x 2 figure, label panels (a)-(e) and save PNG + PDF;
# NULL entries are empty positions (no panel label). An empty first position
# puts panel (a) alone in the centre of the top row.
save_figure <- function(panels, name, width, heights, ncol = 2) {
  empty  <- vapply(panels, is.null, logical(1))
  labels <- ifelse(empty, "", sprintf("(%s)", letters[pmax(cumsum(!empty), 1)]))
  panels <- lapply(panels, function(p) if (is.null(p)) NULL else tidy_panel(p))
  if (empty[1]) {
    # Empty top-left position: centre panel (a) in the top row, with the
    # same width as the panels below (A. Figurek, aesthetics)
    top <- plot_grid(NULL, panels[[2]], NULL, ncol = 3, rel_widths = c(1, 2, 1),
                     labels = c("", labels[2], ""),
                     label_size = 18, label_fontface = "bold")
    bottom <- plot_grid(plotlist = panels[-(1:2)], ncol = ncol,
                        labels = labels[-(1:2)],
                        label_size = 18, label_fontface = "bold",
                        rel_heights = heights[-1])
    fig <- plot_grid(top, bottom, ncol = 1,
                     rel_heights = c(heights[1], sum(heights[-1])))
  } else {
    fig <- plot_grid(plotlist = panels, ncol = ncol,
                     labels = labels,
                     label_size = 18, label_fontface = "bold",
                     rel_heights = heights)
  }
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

# . Bulk: Fan DKD cohort (S20; log2 quantile-normalized, no batch correction) ----
s20      <- readRDS(s20_rds)
meta_fan <- s20$meta %>%
  transmute(sample_id = sample,
            condition = factor(as.character(condition),
                               levels = c("Normal", "Early_DKD", "Advanced_DKD")),
            pyrim_score, calcium_score, ecm_score, P2RY6, PDGFRB) %>%
  # PT marker score as in S86 (KPMP atlas PT markers; proxy of PT content)
  mutate(pt_marker_score = colMeans(s20$expr_mat[pt_markers, sample_id]))

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
# PT filtered by each participant's own adjudicated category (S52 §1)
clin_sc <- read.csv(sc_clin)
pb_pt_s52 <- s50$pb_pt %>%
  mutate(adjudicated = clin_sc$Primary.Adjudicated.Category[
           match(SampleID, clin_sc$Participant.ID)]) %>%
  filter(!(as.character(EnrollmentCategory) == "Chronic Kidney Disease" &
           adjudicated != "Diabetic Kidney Disease")) %>%
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
# PT filtered by each participant's own adjudicated category (S62 §1)
clin_sn <- read.csv(sn_clin)
pb_pt_s62 <- s60$pb_pt %>%
  mutate(adjudicated = clin_sn$Primary.Adjudicated.Category[
           match(patient, clin_sn$Participant.ID)]) %>%
  filter(!(enrollment == "CKD" & adjudicated != "Diabetic Kidney Disease")) %>%
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

# . S20_FanDKD_BulkRNA.R (Fan DKD cohort) ----
fan_set <- list(
  condition_levels = c("Normal", "Early_DKD", "Advanced_DKD"),
  condition_labels = c("Normal", "Early DKD", "Adv. DKD"),
  condition_colors = c("Normal" = "#2166ac", "Early_DKD" = "#f4a582",
                       "Advanced_DKD" = "#d6604d"),
  comparisons_main = list(c("Normal", "Early_DKD"), c("Normal", "Advanced_DKD"))
)
fan_set$subtitle_n <- paste(sprintf("%s (n=%d)", fan_set$condition_labels,
                                    as.integer(table(meta_fan$condition))),
                            collapse = ", ")

fan_box <- function(y_col, y_label, title_str, df = meta_fan, S = fan_set) {
  ymax <- max(df[[y_col]], na.rm = TRUE)
  ymin <- min(df[[y_col]], na.rm = TRUE)
  rng  <- ymax - ymin
  ggplot(df, aes(x = condition, y = .data[[y_col]], fill = condition)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.45, width = 0.5) +
    geom_jitter(position = jitter_pos(0.12), size = 2.2, alpha = 0.85) +
    bh_signif(df, "condition", y_col, S$comparisons_main) +
    scale_fill_manual(values = S$condition_colors, guide = "none") +
    scale_x_discrete(labels = setNames(S$condition_labels, S$condition_levels)) +
    scale_y_continuous(limits = c(ymin - 0.05 * rng, ymax + 0.26 * rng)) +
    labs(x = NULL, y = y_label, title = title_str, subtitle = S$subtitle_n) +
    theme_box()
}

fan_scatter <- function(x_col, y_col, x_label, y_label, title_str,
                        df = meta_fan, S = fan_set) {
  ct <- spearman_ci(df[[x_col]], df[[y_col]])
  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]], color = condition)) +
    geom_point(size = 2.8, alpha = 0.9) +
    geom_smooth(method = "lm", se = FALSE, color = "black",
                linewidth = 0.8, inherit.aes = FALSE,
                aes(x = .data[[x_col]], y = .data[[y_col]])) +
    scale_color_manual(values = S$condition_colors, name = NULL,
                       labels = setNames(S$condition_labels, S$condition_levels)) +
    annotate("text", x = Inf, y = Inf,
             label = ct$label, hjust = 1.05, vjust = 1.2, size = 3.8) +
    scale_y_continuous(expand = cor_y_expand) +
    labs(x = x_label, y = y_label, title = title_str, subtitle = S$subtitle_n) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold"), legend.position = "right")
}

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
    scale_color_manual(values = colors, name = NULL, drop = drop, labels = short_labels) +
    scale_y_continuous(expand = cor_y_expand) +
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
                       labels = short_labels,
                       guide = guide_legend(ncol = 1, override.aes = list(size = 3))) +
    scale_y_continuous(expand = cor_y_expand) +
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

# (a) S86 bulk PT marker score (Fan DKD cohort; former Figure S5a)
f1a <- fan_box("pt_marker_score", "PT marker score (log2)",
               "Fan DKD bulk RNA-seq — PT marker score")
# (b) S20 (Fan DKD cohort)
f1b <- fan_box("pyrim_score", "Pyrimidine module score (log2)",
              "Fan DKD bulk RNA-seq — Pyrimidine score")
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
# FIGURE 2 AND SUPPLEMENTARY FIGURE S2 — P2RY6 ----
# ============================================================================#

# (a) S20 (Fan DKD cohort; centred in the top row)
f2a <- fan_box("P2RY6", "P2RY6 expression (log2)", "Fan DKD bulk RNA-seq — P2RY6")
# (b) S50 Fig_FB_P2RY6.png
#     (Figures.pptx used Fig_FB_P2RY6_pct.png, % P2RY6+ cells; replaced by
#      mean expression for consistency with panels (d)-(f), v5 revision)
f2b <- s50_fib_box("mean_P2RY6", "P2RY6 mean expression (pseudobulk)",
                   "KPMP scRNA — P2RY6 in stroma cells")
# (c) S52 Fig_FB_P2RY6.png
f2c <- s52_fib_box("mean_P2RY6", "P2RY6 mean expression (pseudobulk)",
                   "KPMP scRNA (DKD-focused) — P2RY6 in stroma")
# (d) S60 Fig_FB_P2RY6.png
f2d <- s60_fib_box("mean_P2RY6", "P2RY6 mean expression (pseudobulk)",
                   "KPMP snRNA — P2RY6 in fibroblasts")
# (e) S62 Fig_FB_P2RY6.png
f2e <- s62_fib_box("mean_P2RY6", "P2RY6 mean expression (pseudobulk)",
                   "KPMP snRNA (DKD-focused) — P2RY6 in fibroblasts")

# Revision 1 (A. Figurek, ms_rev1_v42_AF): main figure = snRNA-seq fibroblast
# panels (d, e -> a, b); bulk and scRNA-seq panels (a-c) -> Figure S2
save_figure(list(f2d, f2e), "Figure_2",
            width = 15, heights = 5.6)
save_figure(list(NULL, f2a, f2b, f2c), "Figure_S2",
            width = 15, heights = c(6.2, 5.6))

# ============================================================================#
# FIGURE 3 AND SUPPLEMENTARY FIGURE S3 — Calcium signaling score ----
# ============================================================================#

# (a) S20 (Fan DKD cohort; centred in the top row)
f3a <- fan_box("calcium_score", "Calcium signaling score (log2)",
              "Fan DKD bulk RNA-seq — Calcium score")
# (b) S50 Fig_FB_calcium.png
f3b <- s50_fib_box("mean_calcium_score", "Calcium signaling score (mean expression)",
                   "KPMP scRNA — Calcium score in stroma")
# (c) S52 Fig_FB_calcium.png
f3c <- s52_fib_box("mean_calcium_score", "Calcium signaling score (mean expression)",
                   "KPMP scRNA (DKD-focused) — Calcium score in stroma")
# (d) S60 Fig_FB_calcium.png
f3d <- s60_fib_box("mean_calcium_score", "Calcium signaling score (mean expression)",
                   "KPMP snRNA — Calcium score in fibroblasts")
# (e) S62 Fig_FB_calcium.png
f3e <- s62_fib_box("mean_calcium_score", "Calcium signaling score (mean expression)",
                   "KPMP snRNA (DKD-focused) — Calcium score in FBs")

# Revision 1 (A. Figurek, ms_rev1_v42_AF): main figure = snRNA-seq fibroblast
# panels (d, e -> a, b); bulk and scRNA-seq panels (a-c) -> Figure S3
save_figure(list(f3d, f3e), "Figure_3",
            width = 15, heights = 5.6)
save_figure(list(NULL, f3a, f3b, f3c), "Figure_S3",
            width = 15, heights = c(6.2, 5.6))

# ============================================================================#
# FIGURE 4 AND SUPPLEMENTARY FIGURE S5 — PDGFRB ----
# ============================================================================#

# (a) S20 (Fan DKD cohort; centred in the top row)
f4a <- fan_box("PDGFRB", "PDGFRB expression (log2)", "Fan DKD bulk RNA-seq — PDGFRB")
# (b) S50 Fig_FB_PDGFRB.png
f4b <- s50_fib_box("mean_PDGFRB", "PDGFRB mean expression (pseudobulk)",
                   "KPMP scRNA — PDGFRB in stroma cells")
# (c) S52 Fig_FB_PDGFRB.png
f4c <- s52_fib_box("mean_PDGFRB", "PDGFRB mean expression (pseudobulk)",
                   "KPMP scRNA (DKD-focused) — PDGFRB in stroma")
# (d) S60 Fig_FB_PDGFRB.png
f4d <- s60_fib_box("mean_PDGFRB", "PDGFRB mean expression (pseudobulk)",
                   "KPMP snRNA — PDGFRB in fibroblasts")
# (e) S62 Fig_FB_PDGFRB.png
f4e <- s62_fib_box("mean_PDGFRB", "PDGFRB mean expression (pseudobulk)",
                   "KPMP snRNA (DKD-focused) — PDGFRB in fibroblasts")

# Revision 1 (A. Figurek, ms_rev1_v42_AF): main figure = snRNA-seq fibroblast
# panels (d, e -> a, b); bulk and scRNA-seq panels (a-c) -> Figure S5
save_figure(list(f4d, f4e), "Figure_4",
            width = 15, heights = 5.6)
save_figure(list(NULL, f4a, f4b, f4c), "Figure_S5",
            width = 15, heights = c(6.2, 5.6))

# ============================================================================#
# FIGURE 5 AND SUPPLEMENTARY FIGURE S7 — ECM score ----
# ============================================================================#

# (a) S20 (Fan DKD cohort; centred in the top row)
f5a <- fan_box("ecm_score", "ECM / fibrotic score (log2)", "Fan DKD bulk RNA-seq — ECM score")
# (b) S50 Fig_FB_ECM.png
f5b <- s50_fib_box("mean_ecm_score", "ECM score (mean expression)",
                   "KPMP scRNA — ECM score in stroma")
# (c) S52 Fig_FB_ECM.png
f5c <- s52_fib_box("mean_ecm_score", "ECM score (mean expression)",
                   "KPMP scRNA (DKD-focused) — ECM score in stroma")
# (d) S60 Fig_FB_ECM.png
f5d <- s60_fib_box("mean_ecm_score", "ECM score (mean expression)",
                   "KPMP snRNA — ECM score in fibroblasts")
# (e) S62 Fig_FB_ECM.png
f5e <- s62_fib_box("mean_ecm_score", "ECM score (mean expression)",
                   "KPMP snRNA (DKD-focused) — ECM score in FBs")

# Revision 1 (A. Figurek, ms_rev1_v42_AF): main figure = snRNA-seq fibroblast
# panels (d, e -> a, b); bulk and scRNA-seq panels (a-c) -> Figure S7
save_figure(list(f5d, f5e), "Figure_5",
            width = 15, heights = 5.6)
save_figure(list(NULL, f5a, f5b, f5c), "Figure_S7",
            width = 15, heights = c(6.2, 5.6))

# ============================================================================#
# SUPPLEMENTARY FIGURE S10 (former Figure 6) — PT pyrimidine score vs ECM score ----
# ============================================================================#

# (a) S20 (Fan DKD cohort; centred in the top row)
f6a <- fan_scatter("pyrim_score", "ecm_score", "Pyrimidine score", "ECM score",
                  "Fan DKD bulk RNA-seq — Pyrimidine vs ECM")
# (b) S50 Fig_Cor_PT_Pyrimidine_vs_FB_ECM.png
f6b <- sc_cor(pb_cross_s50, "pyrim_score", "mean_ecm_score",
              "PT pyrimidine score (pseudobulk)", "Stroma ECM score (pseudobulk)",
              "KPMP scRNA — PT pyrimidine vs Stroma ECM",
              colors = s50_set$condition_colors)
# (c) S52 Fig_Cor_PT_Pyrimidine_vs_FB_ECM.png
f6c <- sc_cor(pb_cross_s52, "pyrim_score", "mean_ecm_score",
              "PT pyrimidine score (pseudobulk)", "Stroma ECM score (pseudobulk)",
              "KPMP scRNA (DKD-focused) — PT pyrimidine vs Stroma ECM",
              colors = s52_set$condition_colors, drop = FALSE)
# (d) S60 Fig_Cor_PT_Pyrimidine_vs_FB_ECM.png
f6d <- sc_cor(pb_cross_s60, "pyrim_score", "mean_ecm_score",
              "PT pyrimidine score (pseudobulk)", "FB ECM score (pseudobulk)",
              "KPMP snRNA — PT pyrimidine vs FB ECM",
              colors = s60_set$condition_colors, color_col = "enrollment",
              subtitle = sprintf("n = %d matched patients", nrow(pb_cross_s60)))
# (e) S62 Fig_Cor_PT_Pyrimidine_vs_FB_ECM.png
f6e <- sc_cor(pb_cross_s62, "pyrim_score", "mean_ecm_score",
              "PT pyrimidine score (pseudobulk)", "FB ECM score (pseudobulk)",
              "KPMP snRNA (DKD-focused) — PT pyrimidine vs FB ECM",
              colors = s62_set$condition_colors, color_col = "enrollment",
              subtitle = sprintf("n = %d matched patients", nrow(pb_cross_s62)),
              drop = FALSE)

# BH-adjusted p-values over the five panels (R1.07)
# Revision 1 (A. Figurek, ms_rev1_v42_AF): former Figure 6 -> Figure S10
f6_bh <- bh_cor_panels(list(f6a, f6b, f6c, f6d, f6e))
save_figure(c(list(NULL), f6_bh), "Figure_S10",
            width = 16, heights = c(5.6, 5.2, 5.2))

# ============================================================================#
# FIGURE 6 — PT pyrimidine vs FB ECM within AKI and within CKD (snRNA-seq) ----
# ============================================================================#
# Revision 1 (A. Figurek, ms_rev1_v42_AF): the within-condition correlations
# of Figure S10d (snRNA-seq, all conditions). As in Table S6 (S85), the
# p-values are BH-adjusted over all within-condition correlations (n >= 5) of
# the five panels of Figure S10; both panels use the same axis ranges.

within_bh <- function(panels, min_n = 5) {
  bind_rows(lapply(seq_along(panels), function(i) {
    p <- panels[[i]]; d <- p$data
    x <- rlang::eval_tidy(p$mapping$x, d); y <- rlang::eval_tidy(p$mapping$y, d)
    cond_col <- intersect(c("condition", "EnrollmentCategory", "enrollment"), names(d))[1]
    keep <- !is.na(x) & !is.na(y)
    dd <- data.frame(x = x[keep], y = y[keep], g = as.character(d[[cond_col]][keep]))
    bind_rows(lapply(split(dd, dd$g), function(s) {
      if (nrow(s) < min_n) return(NULL)
      ct <- spearman_ci(s$x, s$y)
      data.frame(panel = i, condition = s$g[1], n = ct$n, rho = ct$rho,
                 ci_lo = ct$ci_lo, ci_hi = ct$ci_hi, p = ct$p)
    }))
  })) %>% mutate(p_bh = p.adjust(p, method = "BH"))
}
w_s10 <- within_bh(list(f6a, f6b, f6c, f6d, f6e))

f6_within <- pb_cross_s60 %>% filter(enrollment %in% c("AKI", "CKD"))
f6_xlim <- range(f6_within$pyrim_score)
f6_ylim <- range(f6_within$mean_ecm_score)
f6_ylim <- f6_ylim + c(-0.05, 0.22) * diff(f6_ylim)

within_cor_panel <- function(cond) {
  df <- filter(f6_within, enrollment == cond)
  r  <- filter(w_s10, panel == 4, condition == cond)
  stopifnot(nrow(r) == 1, r$n == nrow(df))
  pt <- format_p_s6(r$p_bh); if (!startsWith(pt, "<")) pt <- paste("=", pt)
  ggplot(df, aes(x = pyrim_score, y = mean_ecm_score)) +
    geom_point(size = 2.8, alpha = 0.85, color = s60_set$condition_colors[[cond]]) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8) +
    annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.2, size = 3.8,
             label = sprintf("rho = %.2f [%.2f, %.2f]\nBH-adjusted p %s",
                             r$rho, r$ci_lo, r$ci_hi, pt)) +
    coord_cartesian(xlim = f6_xlim, ylim = f6_ylim) +
    labs(x = "PT pyrimidine score (pseudobulk)", y = "FB ECM score (pseudobulk)",
         title = sprintf("KPMP snRNA — %s: PT pyrimidine vs FB ECM", cond),
         subtitle = sprintf("n = %d matched patients", nrow(df))) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))
}
fig6_aki <- within_cor_panel("AKI")
fig6_ckd <- within_cor_panel("CKD")
save_figure(list(fig6_aki, fig6_ckd), "Figure_6", width = 13, heights = 5.6)

# ============================================================================#
# SUPPLEMENTARY FIGURE S4 (fs1) — P2RY6 vs Calcium score ----
# ============================================================================#

# (a) S20 (Fan DKD cohort; centred in the top row)
fs1a <- fan_scatter("P2RY6", "calcium_score", "P2RY6 (log2)", "Calcium score",
                   "Fan DKD bulk RNA-seq — P2RY6 vs Calcium")
# (b) S51 Fig_Cor_FB_P2RY6_vs_Calcium.png
fs1b <- s51_cor_adj("mean_P2RY6", "mean_calcium_score",
                    "Stroma P2RY6 (pseudobulk)", "Stroma calcium score (pseudobulk)",
                    "KPMP scRNA — Stroma P2RY6 vs Calcium (adjudicated)")
# (c) S52 Fig_Cor_FB_P2RY6_vs_Calcium.png
fs1c <- sc_cor(pb_fib_s52, "mean_P2RY6", "mean_calcium_score",
               "Stroma P2RY6 (pseudobulk)", "Stroma calcium score (pseudobulk)",
               "KPMP scRNA (DKD-focused) — Stroma P2RY6 vs Calcium",
               colors = s52_set$condition_colors, drop = FALSE,
               subtitle = sprintf("n = %d stroma samples", nrow(pb_fib_s52)))
# (d) S60 Fig_Cor_FB_P2RY6_vs_Calcium.png
fs1d <- sc_cor(pb_fib_s60, "mean_P2RY6", "mean_calcium_score",
               "FB P2RY6 (pseudobulk)", "FB calcium score (pseudobulk)",
               "KPMP snRNA — FB P2RY6 vs Calcium",
               colors = s60_set$condition_colors, color_col = "enrollment",
               subtitle = sprintf("n = %d fibroblast patients", nrow(pb_fib_s60)))
# (e) S62 Fig_Cor_FB_P2RY6_vs_Calcium.png
fs1e <- sc_cor(pb_fib_s62, "mean_P2RY6", "mean_calcium_score",
               "FB P2RY6 (pseudobulk)", "FB calcium score (pseudobulk)",
               "KPMP snRNA (DKD-focused) — FB P2RY6 vs Calcium",
               colors = s62_set$condition_colors, color_col = "enrollment",
               subtitle = sprintf("n = %d fibroblast patients", nrow(pb_fib_s62)),
               drop = FALSE)

# BH-adjusted p-values over the five panels (R1.07)
fs1_bh <- bh_cor_panels(list(fs1a, fs1b, fs1c, fs1d, fs1e))
save_figure(c(list(NULL), fs1_bh), "Figure_S4",
            width = 16, heights = c(5.6, 5.2, 5.2))

# ============================================================================#
# SUPPLEMENTARY FIGURE S6 (fs2) — P2RY6 vs PDGFRB ----
# ============================================================================#
# Axes swapped relative to the source panels (P2RY6 on x, PDGFRB on y), as in
# Figures S1 and S3. Spearman rho is symmetric; only the linear fit line changes.

# (a) S20 (Fan DKD cohort; centred in the top row)
fs2a <- fan_scatter("P2RY6", "PDGFRB", "P2RY6 (log2)", "PDGFRB (log2)",
                   "Fan DKD bulk RNA-seq — P2RY6 vs PDGFRB")
# (b) S51 Fig_Cor_FB_PDGFRB_vs_P2RY6.png
fs2b <- s51_cor_adj("mean_P2RY6", "mean_PDGFRB",
                    "Stroma P2RY6 (pseudobulk)", "Stroma PDGFRB (pseudobulk)",
                    "KPMP scRNA — Stroma P2RY6 vs PDGFRB (adjudicated)")
# (c) S52 Fig_Cor_FB_PDGFRB_vs_P2RY6.png
fs2c <- sc_cor(pb_fib_s52, "mean_P2RY6", "mean_PDGFRB",
               "Stroma P2RY6 (pseudobulk)", "Stroma PDGFRB (pseudobulk)",
               "KPMP scRNA (DKD-focused) — Stroma P2RY6 vs PDGFRB",
               colors = s52_set$condition_colors, drop = FALSE,
               subtitle = sprintf("n = %d stroma samples", nrow(pb_fib_s52)))
# (d) S60 Fig_Cor_FB_PDGFRB_vs_P2RY6.png
fs2d <- sc_cor(pb_fib_s60, "mean_P2RY6", "mean_PDGFRB",
               "FB P2RY6 (pseudobulk)", "FB PDGFRB (pseudobulk)",
               "KPMP snRNA — FB P2RY6 vs PDGFRB",
               colors = s60_set$condition_colors, color_col = "enrollment",
               subtitle = sprintf("n = %d fibroblast patients", nrow(pb_fib_s60)))
# (e) S62 Fig_Cor_FB_PDGFRB_vs_P2RY6.png
fs2e <- sc_cor(pb_fib_s62, "mean_P2RY6", "mean_PDGFRB",
               "FB P2RY6 (pseudobulk)", "FB PDGFRB (pseudobulk)",
               "KPMP snRNA (DKD-focused) — FB P2RY6 vs PDGFRB",
               colors = s62_set$condition_colors, color_col = "enrollment",
               subtitle = sprintf("n = %d fibroblast patients", nrow(pb_fib_s62)),
               drop = FALSE)

# BH-adjusted p-values over the five panels (R1.07)
fs2_bh <- bh_cor_panels(list(fs2a, fs2b, fs2c, fs2d, fs2e))
save_figure(c(list(NULL), fs2_bh), "Figure_S6",
            width = 16, heights = c(5.6, 5.2, 5.2))

# ============================================================================#
# SUPPLEMENTARY FIGURE S8 (fs3) — P2RY6 vs ECM score ----
# ============================================================================#
# Not present in Figures.pptx (slide duplicated Suppl. Figure 2); assembled
# from the P2RY6 vs ECM outputs following the layout of Figures S1-S2.

# (a) S20 (Fan DKD cohort; centred in the top row)
fs3a <- fan_scatter("P2RY6", "ecm_score", "P2RY6 (log2)", "ECM score",
                   "Fan DKD bulk RNA-seq — P2RY6 vs ECM")
# (b) S51 Fig_Cor_FB_P2RY6_vs_ECM.png
fs3b <- s51_cor_adj("mean_P2RY6", "mean_ecm_score",
                    "Stroma P2RY6 (pseudobulk)", "Stroma ECM score (pseudobulk)",
                    "KPMP scRNA — Stroma P2RY6 vs ECM (adjudicated)")
# (c) S52 Fig_Cor_FB_P2RY6_vs_ECM.png
fs3c <- sc_cor(pb_fib_s52, "mean_P2RY6", "mean_ecm_score",
               "Stroma P2RY6 (pseudobulk)", "Stroma ECM score (pseudobulk)",
               "KPMP scRNA (DKD-focused) — Stroma P2RY6 vs ECM",
               colors = s52_set$condition_colors, drop = FALSE,
               subtitle = sprintf("n = %d stroma samples", nrow(pb_fib_s52)))
# (d) S60 Fig_Cor_FB_P2RY6_vs_ECM.png
fs3d <- sc_cor(pb_fib_s60, "mean_P2RY6", "mean_ecm_score",
               "FB P2RY6 (pseudobulk)", "FB ECM score (pseudobulk)",
               "KPMP snRNA — FB P2RY6 vs ECM",
               colors = s60_set$condition_colors, color_col = "enrollment",
               subtitle = sprintf("n = %d fibroblast patients", nrow(pb_fib_s60)))
# (e) S62 Fig_Cor_FB_P2RY6_vs_ECM.png
fs3e <- sc_cor(pb_fib_s62, "mean_P2RY6", "mean_ecm_score",
               "FB P2RY6 (pseudobulk)", "FB ECM score (pseudobulk)",
               "KPMP snRNA (DKD-focused) — FB P2RY6 vs ECM",
               colors = s62_set$condition_colors, color_col = "enrollment",
               subtitle = sprintf("n = %d fibroblast patients", nrow(pb_fib_s62)),
               drop = FALSE)

# BH-adjusted p-values over the five panels (R1.07)
fs3_bh <- bh_cor_panels(list(fs3a, fs3b, fs3c, fs3d, fs3e))
save_figure(c(list(NULL), fs3_bh), "Figure_S8",
            width = 16, heights = c(5.6, 5.2, 5.2))

cat("\nS80 complete. Figures saved to:", out_dir, "\n")
