# ============================================================================#
# S60_KPMP_snRNA.R
# Tier 3 — KPMP single-nucleus RNA-seq (full fibroblast subtype resolution)
#
# Conditions : Healthy Reference (40) + AKI (33) + CKD (72) — DM-R excluded
# PT analysis : pyrimidine module score, pseudobulk by patient
# FB analysis : subclass.l1 == "FIB"  (FIB, M-FIB, pvFIB, infFIB, MYOF,
#               dFIB, dM-FIB) — combined + per-subtype pseudobulk
# Outputs    : 04-Outputs/S60_KPMP_snRNA/
#              S60_pseudobulk.rds  (used by S62 DKD-focused variant)
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

sn_file  <- "../../KPMP_snRNA/SingleNucleus_KPMP_Explorer_05182025.h5Seurat"
clin_file <- "../../KPMP_snRNA/6d3047c2-588f-4ad4-8eba-f4c46a1e5c14_20251205_OpenAccessClinicalData/20251205_OpenAccessClinicalData.csv"
out_dir  <- "../04-Outputs/S60_KPMP_snRNA"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

pyrim_genes   <- c("CAD", "DHODH", "UMPS", "UCK1", "CDA", "DCTPP1")
calcium_genes <- c("GNAQ", "ITPR2", "ITPR3", "ORAI1")
ecm_genes     <- c("COL1A1", "COL3A1", "FN1", "LUM", "DCN")

min_cells_per_patient <- 20

# Fibroblast subtypes (all 7 via subclass.l1 == "FIB")
fib_subtypes <- c("FIB", "M-FIB", "pvFIB", "infFIB", "MYOF", "dFIB", "dM-FIB")

conditions_keep  <- c("Healthy Reference", "AKI", "CKD")
condition_levels <- c("Healthy Reference", "AKI", "CKD")
condition_short  <- c("Healthy", "AKI", "CKD")

condition_colors <- c(
  "Healthy Reference" = "#2166ac",
  "AKI"              = "#e07b39",
  "CKD"              = "#756bb1"
)

comparisons_main <- list(
  c("Healthy Reference", "AKI"),
  c("Healthy Reference", "CKD")
)

dpi <- 300

# ============================================================================#
# HELPER FUNCTIONS ----
# ============================================================================#

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

make_subtitle <- function(pb, cond_col = "enrollment") {
  sprintf("n patients: %s",
    paste(mapply(function(lv, sh) {
      n <- sum(pb[[cond_col]] == lv, na.rm = TRUE)
      sprintf("%s=%d", sh, n)
    }, condition_levels, condition_short), collapse = ", "))
}

# ============================================================================#
# 1. LOAD OBJECT ONCE + ADD CLINICAL METADATA ----
# ============================================================================#

cat("Loading snRNA object once from h5Seurat...\n")
obj <- LoadH5Seurat(sn_file)

cat("Total cells:", ncol(obj), "| Patients:", length(unique(obj$patient)), "\n")

# Build clinical lookup: patient → enrollment, adjudicated
clin <- read.csv(clin_file)
names(clin)[names(clin) == "Participant.ID"]               <- "patient"
names(clin)[names(clin) == "Enrollment.Category"]          <- "enrollment"
names(clin)[names(clin) == "Primary.Adjudicated.Category"] <- "adjudicated"

obj$enrollment  <- clin$enrollment [match(obj$patient, clin$patient)]
obj$adjudicated <- clin$adjudicated[match(obj$patient, clin$patient)]
obj$enrollment  <- factor(obj$enrollment, levels = condition_levels)

cat("\nEnrollment category distribution (patients):\n")
pat_meta_all <- unique(obj@meta.data[, c("patient", "enrollment", "adjudicated")])
print(table(pat_meta_all$enrollment, useNA = "ifany"))

# Cell barcodes for each compartment (filtered to target conditions)
in_target <- obj$enrollment %in% conditions_keep
cells_pt  <- colnames(obj)[in_target & obj$subclass.l1 == "PT"]
cells_fib <- colnames(obj)[in_target & obj$subclass.l1 == "FIB"]

cat("\nTarget-condition cells — PT:", length(cells_pt), "| FIB:", length(cells_fib), "\n")

# Condition counts
pat_counts <- table(unique(
  obj@meta.data[in_target, c("patient", "enrollment")])$enrollment)
cat("Patients per condition:", paste(names(pat_counts), pat_counts,
                                     sep = "=", collapse = ", "), "\n")

# ============================================================================#
# 2. PT ANALYSIS — Pyrimidine module score ----
# ============================================================================#

cat("\n--- PT analysis ---\n")
cat("Subsetting PT cells from loaded object...\n")

pt <- subset(obj, cells = cells_pt)
pt <- NormalizeData(pt, normalization.method = "LogNormalize",
                    scale.factor = 10000, verbose = FALSE)

cat("PT cells:", ncol(pt), "| Patients:", length(unique(pt$patient)), "\n")

pyrim_found <- intersect(pyrim_genes, rownames(pt))
cat("Pyrimidine genes found:", paste(pyrim_found, collapse = ", "), "\n")

df_pt <- FetchData(pt, vars = c(pyrim_found, "patient", "enrollment"))
df_pt$pyrim_score <- module_mean(df_pt, pyrim_found)

# Pseudobulk by patient
pb_pt <- df_pt %>%
  group_by(patient, enrollment) %>%
  summarise(
    pyrim_score = mean(pyrim_score, na.rm = TRUE),
    across(all_of(pyrim_found), mean),
    n_cells = n(),
    .groups = "drop"
  ) %>%
  filter(n_cells >= min_cells_per_patient) %>%
  mutate(enrollment = factor(enrollment, levels = condition_levels))

cat("PT donors in pseudobulk:", nrow(pb_pt), "\n")
pw_pt <- pairwise.wilcox.test(pb_pt$pyrim_score, pb_pt$enrollment,
                               p.adjust.method = "BH")
cat("PT pyrimidine pairwise Wilcoxon (BH-FDR):\n"); print(pw_pt$p.value)

subtitle_pt <- make_subtitle(pb_pt)

# . Pyrimidine score boxplot ----
p_pt_score <- ggplot(pb_pt, aes(x = enrollment, y = pyrim_score,
                                 fill = enrollment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.45, width = 0.5) +
  geom_jitter(width = 0.12, size = 2.2, alpha = 0.85) +
  stat_compare_means(comparisons = comparisons_main, method = "wilcox.test",
                     label = "p.signif", hide.ns = FALSE) +
  scale_fill_manual(values = condition_colors, guide = "none") +
  scale_x_discrete(labels = setNames(condition_short, condition_levels)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.30))) +
  labs(x = NULL, y = "Pyrimidine module score (mean expression)",
       title = "KPMP snRNA — PT pyrimidine score (pseudobulk by patient)",
       subtitle = subtitle_pt) +
  theme_box()

ggsave(file.path(out_dir, "Fig_PT_pyrimidine_score.png"), p_pt_score,
       width = 6, height = 5, dpi = dpi)

# . Individual pyrimidine genes ----
df_pt_long <- pb_pt %>%
  select(patient, enrollment, n_cells, all_of(pyrim_found)) %>%
  pivot_longer(cols = all_of(pyrim_found), names_to = "gene", values_to = "expr")

p_pt_genes <- ggplot(df_pt_long, aes(x = enrollment, y = expr,
                                      fill = enrollment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.45, width = 0.5) +
  geom_jitter(width = 0.12, size = 1.4, alpha = 0.85) +
  stat_compare_means(comparisons = comparisons_main, method = "wilcox.test",
                     label = "p.signif", hide.ns = FALSE) +
  scale_fill_manual(values = condition_colors, guide = "none") +
  scale_x_discrete(labels = setNames(condition_short, condition_levels)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.55))) +
  facet_wrap(~gene, scales = "free_y", ncol = 3) +
  labs(x = NULL, y = "Mean expression (pseudobulk by patient)",
       title = "KPMP snRNA — Individual pyrimidine genes in PT",
       subtitle = subtitle_pt) +
  theme_box()

ggsave(file.path(out_dir, "Fig_PT_pyrimidine_genes.png"), p_pt_genes,
       width = 11, height = 7, dpi = dpi)

#rm(pt, df_pt, df_pt_long, p_pt_score, p_pt_genes); gc()

# ============================================================================#
# 3. FIBROBLAST ANALYSIS — subclass.l1 == "FIB" ----
# ============================================================================#

cat("\n--- Fibroblast analysis ---\n")
cat("Subsetting FIB cells from loaded object...\n")

fib <- subset(obj, cells = cells_fib)
fib <- NormalizeData(fib, normalization.method = "LogNormalize",
                     scale.factor = 10000, verbose = FALSE)

cat("FB cells:", ncol(fib), "| Patients:", length(unique(fib$patient)), "\n")
cat("FB subtype counts:\n")
print(sort(table(fib$subclass.l2), decreasing = TRUE))

all_fib_genes <- c("P2RY6", "PDGFRB", calcium_genes, ecm_genes)
fib_found     <- intersect(all_fib_genes, rownames(fib))
cat("FB genes found:", paste(fib_found, collapse = ", "), "\n")

df_fib <- FetchData(fib, vars = c(fib_found, "subclass.l2", "patient", "enrollment",
                                   "adjudicated"))

# Per-cell module scores
df_fib$calcium_score <- module_mean(df_fib, calcium_genes)
df_fib$ecm_score     <- module_mean(df_fib, ecm_genes)
if (!"PDGFRB" %in% fib_found) df_fib$PDGFRB <- NA_real_
if (!"P2RY6"  %in% fib_found) df_fib$P2RY6  <- NA_real_

# . 3a. Combined FB pseudobulk (all subtypes pooled) ----
pb_fib <- df_fib %>%
  group_by(patient, enrollment, adjudicated) %>%
  summarise(
    mean_P2RY6         = mean(.data$P2RY6, na.rm = TRUE),
    pct_P2RY6_pos = mean(.data$P2RY6 > 0, na.rm = TRUE) * 100,
    mean_PDGFRB        = mean(.data$PDGFRB, na.rm = TRUE),
    mean_calcium_score = mean(.data$calcium_score, na.rm = TRUE),
    mean_ecm_score     = mean(.data$ecm_score, na.rm = TRUE),
    n_cells       = n(),
    .groups       = "drop"
  ) %>%
  filter(n_cells >= min_cells_per_patient) %>%
  mutate(enrollment = factor(enrollment, levels = condition_levels))

cat("FB donors in pseudobulk:", nrow(pb_fib), "\n")

pw_fib_p2ry6 <- pairwise.wilcox.test(pb_fib$mean_P2RY6, pb_fib$enrollment,
                                       p.adjust.method = "BH")
cat("FB P2RY6 pairwise Wilcoxon (BH-FDR):\n"); print(pw_fib_p2ry6$p.value)

subtitle_fib <- make_subtitle(pb_fib)

plot_fib_box <- function(df, y_col, y_label, title_str) {
  ymax <- max(df[[y_col]], na.rm = TRUE)
  ymin <- min(df[[y_col]], na.rm = TRUE)
  rng  <- ymax - ymin
  ggplot(df, aes(x = enrollment, y = .data[[y_col]], fill = enrollment)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.45, width = 0.5) +
    geom_jitter(width = 0.12, size = 2.2, alpha = 0.85) +
    stat_compare_means(comparisons = comparisons_main, method = "wilcox.test",
                       label = "p.signif", hide.ns = FALSE) +
    scale_fill_manual(values = condition_colors, guide = "none") +
    scale_x_discrete(labels = setNames(condition_short, condition_levels)) +
    scale_y_continuous(limits = c(ymin - 0.05 * rng, ymax + 0.26 * rng)) +
    labs(x = NULL, y = y_label, title = title_str, subtitle = subtitle_fib) +
    theme_box()
}

p_p2ry6   <- plot_fib_box(pb_fib, "mean_P2RY6",         "P2RY6 mean expression (pseudobulk)",        "KPMP snRNA — P2RY6 in fibroblasts")
p_pct     <- plot_fib_box(pb_fib, "pct_P2RY6_pos", "% P2RY6+ cells per patient",               "KPMP snRNA — P2RY6 detection in fibroblasts")
p_calcium <- plot_fib_box(pb_fib, "mean_calcium_score", "Calcium signaling score (mean expression)", "KPMP snRNA — Calcium score in fibroblasts")
p_ecm     <- plot_fib_box(pb_fib, "mean_ecm_score",     "ECM score (mean expression)",               "KPMP snRNA — ECM score in fibroblasts")
p_pdgfrb  <- plot_fib_box(pb_fib, "mean_PDGFRB",        "PDGFRB mean expression (pseudobulk)",       "KPMP snRNA — PDGFRB in fibroblasts")

ggsave(file.path(out_dir, "Fig_FB_P2RY6.png"),       p_p2ry6,   width = 6, height = 5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_FB_P2RY6_pct.png"),   p_pct,     width = 6, height = 5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_FB_calcium.png"),      p_calcium, width = 6, height = 5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_FB_ECM.png"),          p_ecm,     width = 6, height = 5, dpi = dpi)
ggsave(file.path(out_dir, "Fig_FB_PDGFRB.png"),       p_pdgfrb,  width = 6, height = 5, dpi = dpi)

# . 3b. Per-subtype pseudobulk — P2RY6 ----
pb_fib_sub <- df_fib %>%
  group_by(patient, enrollment, subclass.l2) %>%
  summarise(
    mean_P2RY6   = mean(.data$P2RY6, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  ) %>%
  filter(n_cells >= min_cells_per_patient,
         subclass.l2 %in% fib_subtypes) %>%
  mutate(enrollment  = factor(enrollment, levels = condition_levels),
         subclass.l2 = factor(subclass.l2, levels = fib_subtypes))

p_sub_p2ry6 <- ggplot(pb_fib_sub, aes(x = enrollment, y = mean_P2RY6,
                                        fill = enrollment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.45, width = 0.5) +
  geom_jitter(width = 0.12, size = 1.4, alpha = 0.85) +
  stat_compare_means(comparisons = comparisons_main, method = "wilcox.test",
                     label = "p.signif", hide.ns = FALSE) +
  scale_fill_manual(values = condition_colors, guide = "none") +
  scale_x_discrete(labels = setNames(condition_short, condition_levels)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.55))) +
  facet_wrap(~subclass.l2, scales = "free_y", ncol = 4) +
  labs(x = NULL, y = "P2RY6 mean expression (pseudobulk)",
       title = "KPMP snRNA — P2RY6 per fibroblast subtype",
       subtitle = "Min 20 cells per patient per subtype") +
  theme_box()

ggsave(file.path(out_dir, "Fig_FB_P2RY6_by_subtype.png"), p_sub_p2ry6,
       width = 14, height = 8, dpi = dpi)

#rm(fib, df_fib, p_sub_p2ry6,
#   p_p2ry6, p_pct, p_calcium, p_ecm, p_pdgfrb); gc()

# Free the full object once both analyses are complete
#rm(obj); gc()

# ============================================================================#
# 4. CROSS-COMPARTMENT CORRELATION ----
# ============================================================================#

cat("\n--- Cross-compartment correlation ---\n")

pb_cross <- inner_join(
  pb_pt  %>% select(patient, enrollment, pyrim_score),
  pb_fib %>% select(patient, mean_P2RY6, pct_P2RY6_pos, mean_ecm_score),
  by = "patient"
)
cat("Matched patients for correlation:", nrow(pb_cross), "\n")

ct_p2ry6 <- spearman_ci(pb_cross$pyrim_score, pb_cross$mean_P2RY6)
ct_ecm   <- spearman_ci(pb_cross$pyrim_score, pb_cross$mean_ecm_score)
cat("PT pyrimidine vs FB P2RY6:", ct_p2ry6$label, "\n")
cat("PT pyrimidine vs FB ECM:  ", ct_ecm$label, "\n")

plot_cor <- function(df, x_col, y_col, x_lab, y_lab, title_str, ct) {
  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]],
                 color = enrollment)) +
    geom_point(size = 2.8, alpha = 0.85) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8,
                inherit.aes = FALSE, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    annotate("text", x = Inf, y = Inf,
             label = ct$label, hjust = 1.05, vjust = 1.2, size = 3.8) +
    scale_color_manual(values = condition_colors, name = NULL) +
    labs(x = x_lab, y = y_lab, title = title_str,
         subtitle = sprintf("n = %d matched patients", nrow(df))) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold"), legend.position = "right")
}

p_cor_p2ry6 <- plot_cor(pb_cross, "pyrim_score", "mean_P2RY6",
  "PT pyrimidine score (pseudobulk)", "FB P2RY6 (pseudobulk)",
  "KPMP snRNA — PT pyrimidine vs FB P2RY6", ct_p2ry6)
p_cor_ecm <- plot_cor(pb_cross, "pyrim_score", "mean_ecm_score",
  "PT pyrimidine score (pseudobulk)", "FB ECM score (pseudobulk)",
  "KPMP snRNA — PT pyrimidine vs FB ECM", ct_ecm)

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
  title_str = "KPMP snRNA — FB P2RY6 vs Calcium",
  ct        = ct_p2ry6_calcium
) + labs(subtitle = sprintf("n = %d fibroblast patients", ct_p2ry6_calcium$n))
ggsave(filename = file.path(out_dir, "Fig_Cor_FB_P2RY6_vs_Calcium.png"),
       plot = p_cor_p2ry6_calcium, width = 7, height = 4.5, dpi = dpi)

ct_p2ry6_ecm <- spearman_ci(x = pb_fib$mean_P2RY6, y = pb_fib$mean_ecm_score)
p_cor_p2ry6_ecm <- plot_cor(
  df        = pb_fib,
  x_col     = "mean_P2RY6",
  y_col     = "mean_ecm_score",
  x_lab     = "FB P2RY6 (pseudobulk)",
  y_lab     = "FB ECM score (pseudobulk)",
  title_str = "KPMP snRNA — FB P2RY6 vs ECM",
  ct        = ct_p2ry6_ecm
) + labs(subtitle = sprintf("n = %d fibroblast patients", ct_p2ry6_ecm$n))
ggsave(filename = file.path(out_dir, "Fig_Cor_FB_P2RY6_vs_ECM.png"),
       plot = p_cor_p2ry6_ecm, width = 7, height = 4.5, dpi = dpi)

ct_pdgfrb_p2ry6 <- spearman_ci(x = pb_fib$mean_PDGFRB, y = pb_fib$mean_P2RY6)
p_cor_pdgfrb_p2ry6 <- plot_cor(
  df        = pb_fib,
  x_col     = "mean_PDGFRB",
  y_col     = "mean_P2RY6",
  x_lab     = "FB PDGFRB (pseudobulk)",
  y_lab     = "FB P2RY6 (pseudobulk)",
  title_str = "KPMP snRNA — FB PDGFRB vs P2RY6",
  ct        = ct_pdgfrb_p2ry6
) + labs(subtitle = sprintf("n = %d fibroblast patients", ct_pdgfrb_p2ry6$n))
ggsave(filename = file.path(out_dir, "Fig_Cor_FB_PDGFRB_vs_P2RY6.png"),
       plot = p_cor_pdgfrb_p2ry6, width = 7, height = 4.5, dpi = dpi)

# ============================================================================#
# 5. SAVE PSEUDOBULK DATA FOR S62 ----
# ============================================================================#

saveRDS(
  list(
    pb_pt       = pb_pt,
    pb_fib      = pb_fib,
    pb_fib_sub  = pb_fib_sub,
    pb_cross    = pb_cross
  ),
  file.path(out_dir, "S60_pseudobulk.rds")
)

cat("\nS60 complete. Outputs saved to:", out_dir, "\n")
cat("PT patients:", nrow(pb_pt),
    "| FB patients:", nrow(pb_fib),
    "| Matched:", nrow(pb_cross), "\n")
