# ============================================================================#
# S86_PT_Abundance.R
# Proximal tubule (PT) abundance and PT cell states by condition, to interpret
# the pyrimidine results at the three levels (bulk, scRNA-seq, snRNA-seq).
# Reviewer comments R2.08, R2.09 (Life, life-4655118)
#
# Single-cell (per participant, KPMP annotation):
#   - PT cells/nuclei as % of all cells/nuclei of the participant
#   - reference PT (no altered state) as % of all cells/nuclei
#   - PT cell states as % of the participant's PT cells/nuclei:
#       reference (healthy PT), adaptive (aPT), failed repair (frPT),
#       cycling (cycPT), degenerative (dPT)
#     scRNA-seq: subclass.level2 (state.l1: reference / altered)
#     snRNA-seq: subclass.l3
#   - participants with < 20 PT cells/nuclei, i.e. not in the PT pyrimidine
#     pseudobulk (Figure 1b-e)
#   - all CKD and DKD-focused (CKD restricted to adjudicated DKD, as in
#     S52/S62)
# Bulk (Fan DKD cohort, S20): PT marker score = mean log2 expression of the
#   PT marker genes of the KPMP atlas (Lake et al. 2023, Extended Data
#   Fig. 2a): LRP2, SLC5A12, SLC22A6, PRODH2, SLC5A2, SLC22A8, SLC7A8,
#   SLC34A1, SLC5A10, SLC5A11, SLC7A13 (no PT label exists in bulk tissue, so
#   the marker score is used as a proxy of PT content).
#
# Only the cell metadata of the h5Seurat files are read (hdf5r); no expression
# matrices are loaded.
#
# Group tests: each disease group vs the control group (Healthy / Normal),
# two-sided Wilcoxon rank-sum test, BH-adjusted over the two comparisons of
# each measure and set (as in S80).
#
# Outputs: 04-Outputs/S86_PT_Abundance/
#   PT_abundance_per_participant.csv   per participant (single-cell)
#   Table_PT_abundance.csv/.md         median (IQR) by condition + BH p
#   Table_PT_pyrimidine_cor.csv        PT pyrimidine score vs PT abundance
#   Fig_PT_abundance.png               PT %, reference PT %, altered and
#                                      degenerative PT % by condition
#   Figure_S1.png/.pdf                 Supplementary Figure S1 (Revision 1;
#                                      Figure S5 before the renumbering of the
#                                      supplements): reference-state PT % of
#                                      all cells (scRNA-seq / snRNA-seq, all
#                                      CKD and DKD-focused); BH stars, disease
#                                      vs control. The bulk PT marker score
#                                      panel moved to Figure 1a (S80).
# ============================================================================#

rm(list = ls())

library(hdf5r)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsignif)
library(cowplot)

# ============================================================================#
# SETTINGS ----
# ============================================================================#

sc_file <- "../../KPMP_scRNA/KPMP_PREMIERE_SC_version2_ForExplorer_RemovedBatchEffect_Final2025.h5seurat"
sn_file <- "../../KPMP_snRNA/SingleNucleus_KPMP_Explorer_05182025.h5Seurat"
sc_clin <- "../../KPMP_scRNA/OpenAccessClinicalData/20251205_OpenAccessClinicalData.csv"
sn_clin <- "../../KPMP_snRNA/6d3047c2-588f-4ad4-8eba-f4c46a1e5c14_20251205_OpenAccessClinicalData/20251205_OpenAccessClinicalData.csv"
s20_rds <- "../04-Outputs/S20_FanDKD_BulkRNA/S20_normalised.rds"
s50_rds <- "../04-Outputs/S50_KPMP_scRNA/S50_pseudobulk.rds"
s60_rds <- "../04-Outputs/S60_KPMP_snRNA/S60_pseudobulk.rds"
out_dir <- "../04-Outputs/S86_PT_Abundance"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

min_cells <- 20   # PT pseudobulk threshold of S50/S60
groups    <- c("Healthy", "AKI", "CKD")
sc_lab    <- c("Healthy Reference Tissue" = "Healthy", "Acute Kidney Injury" = "AKI",
               "Chronic Kidney Disease" = "CKD")
sn_lab    <- c("Healthy Reference" = "Healthy", "AKI" = "AKI", "CKD" = "CKD")
pt_markers <- c("LRP2", "SLC5A12", "SLC22A6", "PRODH2", "SLC5A2", "SLC22A8", "SLC7A8",
                "SLC34A1", "SLC5A10", "SLC5A11", "SLC7A13")
group_colors <- c(Healthy = "#2166ac", AKI = "#e07b39", CKD = "#756bb1")

# ============================================================================#
# HELPER FUNCTIONS ----
# ============================================================================#

# one metadata column of an h5Seurat file (factors are stored as levels/values)
read_meta <- function(h, cols) {
  m <- h[["meta.data"]]
  as_tibble(lapply(setNames(cols, cols), function(col) {
    x <- m[[col]]
    if (inherits(x, "H5Group")) x[["levels"]][][x[["values"]][]] else x[]
  }))
}

# PT state from the KPMP subclass label (aPT*, frPT*, cycPT, dPT*, PT-S*)
pt_state <- function(lab) {
  case_when(startsWith(lab, "aPT")   ~ "adaptive",
            startsWith(lab, "frPT")  ~ "failed repair",
            startsWith(lab, "cycPT") ~ "cycling",
            startsWith(lab, "dPT")   ~ "degenerative",
            startsWith(lab, "PT-S")  ~ "reference",
            TRUE                     ~ NA_character_)
}

# per participant: total cells, PT cells, PT % and PT state %
pt_summary <- function(meta) {
  meta %>%
    group_by(id, group) %>%
    summarise(n_all        = n(),
              n_pt         = sum(is_pt),
              pct_pt       = 100 * n_pt / n_all,
              pct_ref_all  = 100 * sum(state == "reference",     na.rm = TRUE) / n_all,
              pct_ref      = 100 * sum(state == "reference",     na.rm = TRUE) / n_pt,
              pct_adaptive = 100 * sum(state == "adaptive",      na.rm = TRUE) / n_pt,
              pct_failed   = 100 * sum(state == "failed repair", na.rm = TRUE) / n_pt,
              pct_cycling  = 100 * sum(state == "cycling",       na.rm = TRUE) / n_pt,
              pct_degen    = 100 * sum(state == "degenerative",  na.rm = TRUE) / n_pt,
              .groups = "drop") %>%
    mutate(pct_altered = 100 - pct_ref,
           in_pb       = n_pt >= min_cells)
}

med_iqr <- function(v, digits = 1) {
  v <- v[!is.na(v)]
  q <- quantile(v, c(0.5, 0.25, 0.75), names = FALSE)
  sprintf(paste0("%.", digits, "f (%.", digits, "f–%.", digits, "f)"), q[1], q[2], q[3])
}

# Wilcoxon p-values of each disease group vs the control group (first level),
# BH-adjusted over the comparisons of one measure
pairwise_bh <- function(v, g) {
  lv <- levels(droplevels(g))
  pr <- lapply(lv[-1], function(x) c(lv[1], x))
  p  <- vapply(pr, function(x) wilcox.test(v[g == x[2]], v[g == x[1]])$p.value, numeric(1))
  setNames(p.adjust(p, method = "BH"), vapply(pr, function(x) paste(x[2], "vs", x[1]), ""))
}

format_p <- function(p) ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p))

# ============================================================================#
# 1. scRNA-seq METADATA ----
# ============================================================================#

h <- H5File$new(sc_file, mode = "r")
meta_sc <- read_meta(h, c("SampleID", "EnrollmentCategory", "cell_type", "subclass.level2"))
h$close_all()

clin_sc <- read.csv(sc_clin)
meta_sc <- meta_sc %>%
  transmute(id    = as.character(SampleID),
            group = unname(sc_lab[EnrollmentCategory]),
            is_pt = cell_type == "epithelial cell of proximal tubule",
            state = ifelse(is_pt, pt_state(subclass.level2), NA_character_)) %>%
  filter(!is.na(group))
stopifnot(!anyNA(meta_sc$state[meta_sc$is_pt]))
pt_sc <- pt_summary(meta_sc) %>%
  mutate(dataset = "scRNA-seq",
         adjudicated = clin_sc$Primary.Adjudicated.Category[match(id, clin_sc$Participant.ID)])

# ============================================================================#
# 2. snRNA-seq METADATA ----
# ============================================================================#

h <- H5File$new(sn_file, mode = "r")
meta_sn <- read_meta(h, c("patient", "subclass.l1", "subclass.l3", "percent.cortex"))
h$close_all()

clin_sn <- read.csv(sn_clin)
cortex_sn <- meta_sn %>% group_by(id = patient) %>%
  summarise(pct_cortex = mean(percent.cortex, na.rm = TRUE), .groups = "drop")
meta_sn <- meta_sn %>%
  transmute(id    = patient,
            group = unname(sn_lab[clin_sn$Enrollment.Category[match(patient, clin_sn$Participant.ID)]]),
            is_pt = subclass.l1 == "PT",
            state = ifelse(is_pt, pt_state(subclass.l3), NA_character_)) %>%
  filter(!is.na(group))
stopifnot(!anyNA(meta_sn$state[meta_sn$is_pt]))
pt_sn <- pt_summary(meta_sn) %>%
  left_join(cortex_sn, by = "id") %>%
  mutate(dataset = "snRNA-seq",
         adjudicated = clin_sn$Primary.Adjudicated.Category[match(id, clin_sn$Participant.ID)])

# ============================================================================#
# 3. CHECK AGAINST THE PT PSEUDOBULK OF FIGURE 1 ----
# ============================================================================#

# participants with >= 20 PT cells must be exactly the PT pseudobulk of S50/S60
s50 <- readRDS(s50_rds)
s60 <- readRDS(s60_rds)
stopifnot(setequal(pt_sc$id[pt_sc$in_pb], as.character(s50$pb_pt$SampleID)),
          setequal(pt_sn$id[pt_sn$in_pb], s60$pb_pt$patient))
stopifnot(all(pt_sc$n_pt[match(s50$pb_pt$SampleID, pt_sc$id)] == s50$pb_pt$n_cells),
          all(pt_sn$n_pt[match(s60$pb_pt$patient, pt_sn$id)] == s60$pb_pt$n_cells))

pt_all <- bind_rows(pt_sc, pt_sn) %>%
  mutate(group = factor(group, levels = groups))
pt_sets <- bind_rows(
  pt_all %>% mutate(set = "all CKD"),
  pt_all %>% filter(!(group == "CKD" & adjudicated != "Diabetic Kidney Disease")) %>%
    mutate(set = "DKD-focused")) %>%
  mutate(set = factor(set, levels = c("all CKD", "DKD-focused")))
write.csv(pt_all, file.path(out_dir, "PT_abundance_per_participant.csv"), row.names = FALSE)

# ============================================================================#
# 4. BULK: PT MARKER SCORE IN THE FAN DKD COHORT ----
# ============================================================================#

s20 <- readRDS(s20_rds)
stopifnot(all(pt_markers %in% rownames(s20$expr_mat)),
          identical(colnames(s20$expr_mat), s20$meta$sample))
fan <- s20$meta %>%
  transmute(id = sample,
            group = factor(as.character(condition), levels = c("Normal", "Early_DKD", "Advanced_DKD")),
            pt_marker_score = colMeans(s20$expr_mat[pt_markers, ]),
            pyrim_score, P2RY6)

# ============================================================================#
# 5. SUMMARY TABLE ----
# ============================================================================#

measures <- c(pct_pt = "PT, % of all cells/nuclei",
              pct_ref_all = "Reference-state PT, % of all cells/nuclei",
              pct_altered = "Altered PT (non-reference), % of PT",
              pct_adaptive = "Adaptive PT (aPT), % of PT",
              pct_failed = "Failed-repair PT (frPT), % of PT",
              pct_cycling = "Cycling PT (cycPT), % of PT",
              pct_degen = "Degenerative PT (dPT), % of PT",
              pct_cortex = "Cortex, % (snRNA-seq region estimate)")

sc_tab <- pt_sets %>%
  pivot_longer(all_of(names(measures)), names_to = "measure", values_to = "value") %>%
  filter(!is.na(value)) %>%
  group_by(dataset, set, measure) %>%
  group_modify(function(d, k) {
    p <- pairwise_bh(d$value, d$group)
    tibble(Healthy = med_iqr(d$value[d$group == "Healthy"]),
           AKI     = med_iqr(d$value[d$group == "AKI"]),
           CKD     = med_iqr(d$value[d$group == "CKD"]),
           p_AKI_vs_Healthy = p["AKI vs Healthy"],
           p_CKD_vs_Healthy = p["CKD vs Healthy"])
  }) %>%
  ungroup() %>%
  mutate(measure = factor(unname(measures[measure]), levels = measures)) %>%
  arrange(dataset, set, measure)

n_tab <- pt_sets %>%
  group_by(dataset, set, group) %>%
  summarise(n = n(), n_below = sum(!in_pb), .groups = "drop") %>%
  mutate(txt = sprintf("%d (%d < %d PT)", n, n_below, min_cells)) %>%
  select(dataset, set, group, txt) %>%
  pivot_wider(names_from = group, values_from = txt) %>%
  mutate(measure = "Participants (with < 20 PT cells/nuclei)")

p_fan <- pairwise_bh(fan$pt_marker_score, fan$group)
fan_tab <- tibble(
  dataset = "Bulk (Fan)", set = "DKD", measure = "PT marker score, log2",
  Healthy = med_iqr(fan$pt_marker_score[fan$group == "Normal"], 2),
  AKI     = med_iqr(fan$pt_marker_score[fan$group == "Early_DKD"], 2),
  CKD     = med_iqr(fan$pt_marker_score[fan$group == "Advanced_DKD"], 2),
  p_AKI_vs_Healthy = p_fan["Early_DKD vs Normal"],
  p_CKD_vs_Healthy = p_fan["Advanced_DKD vs Normal"])

tab <- bind_rows(fan_tab,
                 mutate(n_tab, set = as.character(set)),
                 mutate(sc_tab, set = as.character(set), measure = as.character(measure))) %>%
  mutate(across(starts_with("p_"), ~ ifelse(is.na(.x), "", format_p(.x))))
print(as.data.frame(tab))
write.csv(tab, file.path(out_dir, "Table_PT_abundance.csv"), row.names = FALSE)

md <- c("# PT abundance and PT cell states by condition",
        "",
        "Median (IQR) per participant; bulk columns are Normal / early DKD / advanced DKD. BH-adjusted two-sided Wilcoxon p.",
        "",
        paste0("| ", paste(names(tab), collapse = " | "), " |"),
        paste0(paste(rep("|---", ncol(tab)), collapse = ""), "|"),
        apply(tab, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
writeLines(md, file.path(out_dir, "Table_PT_abundance.md"))

# ============================================================================#
# 6. PT PYRIMIDINE SCORE VS PT ABUNDANCE ----
# ============================================================================#

# Spearman, pooled and within condition, on the PT pseudobulk of Figure 1b,d
# (scRNA-seq / snRNA-seq, all CKD) and on the Fan bulk samples
cor_one <- function(x, y) {
  ok <- complete.cases(x, y)
  if (sum(ok) < 5) return(tibble(n = sum(ok), rho = NA_real_, p = NA_real_))
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
  tibble(n = sum(ok), rho = unname(ct$estimate), p = ct$p.value)
}
pb_join <- bind_rows(
  s50$pb_pt %>% transmute(id = as.character(SampleID), pyrim_score, dataset = "scRNA-seq"),
  s60$pb_pt %>% transmute(id = patient, pyrim_score, dataset = "snRNA-seq")) %>%
  inner_join(pt_all, by = c("dataset", "id"))
cor_tab <- bind_rows(
  pb_join %>% mutate(group = "pooled"),
  pb_join %>% mutate(group = as.character(group))) %>%
  pivot_longer(c(pct_pt, pct_ref_all, pct_altered, pct_degen), names_to = "measure", values_to = "value") %>%
  group_by(dataset, measure, group) %>%
  group_modify(~ cor_one(.x$value, .x$pyrim_score)) %>%
  ungroup()
cor_fan <- bind_rows(
  fan %>% mutate(grp = "pooled"), fan %>% mutate(grp = as.character(group))) %>%
  group_by(group = grp) %>%
  group_modify(~ bind_rows(
    mutate(cor_one(.x$pt_marker_score, .x$pyrim_score), measure = "pt_marker_vs_pyrim"),
    mutate(cor_one(.x$pt_marker_score, .x$P2RY6),       measure = "pt_marker_vs_P2RY6"))) %>%
  ungroup() %>% mutate(dataset = "Bulk (Fan)")
cor_tab <- bind_rows(cor_tab, cor_fan) %>%
  group_by(dataset, measure) %>% mutate(p_BH = p.adjust(p, method = "BH")) %>% ungroup()
print(as.data.frame(cor_tab), digits = 3)
write.csv(cor_tab, file.path(out_dir, "Table_PT_pyrimidine_cor.csv"), row.names = FALSE)

# ============================================================================#
# 7. FIGURE ----
# ============================================================================#

fig_df <- pt_sets %>%
  filter(set == "all CKD") %>%
  pivot_longer(c(pct_pt, pct_ref_all, pct_altered, pct_degen), names_to = "measure", values_to = "value") %>%
  mutate(measure = factor(measures[measure],
                          levels = measures[c("pct_pt", "pct_ref_all", "pct_altered", "pct_degen")]))
p <- ggplot(fig_df, aes(group, value, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_point(aes(shape = in_pb), position = position_jitter(width = 0.15, seed = 2026),
             size = 1.2, alpha = 0.7) +
  scale_fill_manual(values = group_colors, guide = "none") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 4),
                     labels = c(`TRUE` = ">= 20 PT (in Figure 1)", `FALSE` = "< 20 PT (not in Figure 1)"),
                     name = NULL) +
  facet_grid(measure ~ dataset, scales = "free_y", labeller = label_wrap_gen(22)) +
  labs(x = NULL, y = "% per participant") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(out_dir, "Fig_PT_abundance.png"), p, width = 6.5, height = 9, dpi = 300)

# ============================================================================#
# 8. SUPPLEMENTARY FIGURE S1 ----
# ============================================================================#

# (a, b) scRNA-seq, (c, d) snRNA-seq, each all CKD and DKD-focused;
# reference-state PT cells only (total PT cells in Table S1). The bulk PT
# marker score (former panel a) is Figure 1a (S80).
# Each disease group vs the control group, BH-adjusted over the two
# comparisons of each panel (same p-values as Table S1).
p_stars <- function(p) {
  as.character(cut(p, c(-Inf, 1e-4, 1e-3, 1e-2, 0.05, Inf), c("****", "***", "**", "*", "ns")))
}
s5_box <- function(df, y_col, y_lab, title_str, colors, labels) {
  df  <- df[!is.na(df[[y_col]]), ]
  g   <- droplevels(df$group)
  lv  <- levels(g)
  cmp <- lapply(lv[-1], function(x) c(lv[1], x))
  p   <- pairwise_bh(df[[y_col]], g)
  ann <- p_stars(p[vapply(cmp, function(x) paste(x[2], "vs", x[1]), "")])
  rng <- range(df[[y_col]], na.rm = TRUE)
  sub <- paste(sprintf("%s (n=%d)", labels, as.integer(table(g))), collapse = ", ")
  ggplot(df, aes(x = group, y = .data[[y_col]], fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.45, width = 0.5) +
    geom_jitter(position = position_jitter(width = 0.12, seed = 2026), size = 1.6, alpha = 0.8) +
    geom_signif(comparisons = cmp, annotations = ann, step_increase = 0.12, tip_length = 0.03,
                size = 0.3, textsize = 3.88, vjust = 0) +
    scale_fill_manual(values = colors, guide = "none") +
    scale_x_discrete(labels = setNames(labels, lv)) +
    scale_y_continuous(limits = c(rng[1] - 0.05 * diff(rng), rng[2] + 0.42 * diff(rng))) +
    labs(x = NULL, y = y_lab, title = title_str, subtitle = sub) +
    theme_bw(base_size = 13) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1),
          plot.title = element_text(face = "bold"), plot.title.position = "plot",
          plot.margin = margin(t = 26, r = 8, b = 6, l = 8))
}

sc_panels <- list()
for (ds in c("scRNA-seq", "snRNA-seq")) {
  for (st in c("all CKD", "DKD-focused")) {
    df <- pt_sets %>% filter(dataset == ds, set == st)
    sc_panels[[length(sc_panels) + 1]] <- s5_box(
      df, "pct_ref_all",
      sprintf("Reference-state PT, %% of all %s", if (ds == "snRNA-seq") "nuclei" else "cells"),
      sprintf("KPMP %s%s\nReference-state PT cells", sub("-seq", "", ds),
              if (st == "DKD-focused") " (DKD-focused)" else ""),
      group_colors, c("Healthy", "AKI", if (st == "DKD-focused") "CKD (DKD)" else "CKD"))
  }
}
fig_s1 <- plot_grid(plotlist = sc_panels, ncol = 2, labels = sprintf("(%s)", letters[1:4]),
                    label_size = 18, label_fontface = "bold")
fig_s1 <- ggdraw(fig_s1) + theme(plot.background = element_rect(fill = "white", colour = NA))
ggsave(file.path(out_dir, "Figure_S1.png"), fig_s1, width = 15, height = 11.2, dpi = 600, bg = "white")
ggsave(file.path(out_dir, "Figure_S1.pdf"), fig_s1, width = 15, height = 11.2, device = cairo_pdf)
