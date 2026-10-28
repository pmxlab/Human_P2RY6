# =========================================================================== #
# S43_Combined_Sensitivity.R
# Sensitivity analysis of the combined bulk RNA-seq analysis (S40/S42)
# Reviewer comments R1.10, R2.02, R2.03 (Life, life-4655118)
#
# Problem: Normal samples exist only in Fan (log2 quantile-normalized) and
#   AKI/CKD samples only in KPMP (TMM log2 CPM). Cohort is confounded with
#   diagnosis, and S40 runs ComBat with a binary Normal/Diseased covariate,
#   which estimates the cohort shift from Fan DKD vs KPMP AKI+CKD
#   (Nygaard et al., Biostatistics 2016;17:29-39).
#
# Methods compared (same features, comparisons and tests as S40/S42):
#   1. ComBat_binary    : S40 as published (mod = ~ Normal vs Diseased).
#   2. Rank_in_sample   : no batch correction; module score = mean within-
#                         sample rank of the module genes / number of genes
#                         (singscore-style simple score); single genes =
#                         within-sample percentile rank. Invariant to the
#                         normalization of each cohort (quantile and TMM-CPM
#                         both preserve within-sample gene ranks).
#   3. ComBat_bridge_Adv: ComBat with DKD as a shared reference group:
#                         Fan Advanced DKD and KPMP CKD-enrolled participants
#                         with adjudicated DKD form one biological group; the
#                         cohort effect is estimated from this group only.
#   4. ComBat_bridge_All: as 3, but all Fan DKD (early + advanced) form the
#                         shared group with KPMP DKD.
#   Bridged methods assume that DKD tissue is comparable between cohorts;
#   CKD(DKD) vs Normal is then fixed by that assumption (not independent).
#
# Datasets: All (KPMP CKD = all CKD; as S40, Figures 1a-5a) and
#           DKD_focused (KPMP CKD = adjudicated DKD only; as S42, Figures 1b-5b).
#
# Effect size: Cliff's delta (group2 vs group1; positive = higher in group2),
#   comparable across methods with different scales.
# =========================================================================== #

rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(sva)

# =========================================================================== #
# SETTINGS ----
# =========================================================================== #

s20_rds <- "../04-Outputs/S20_FanDKD_BulkRNA/S20_normalised.rds"
s30_rds <- "../04-Outputs/S30_KPMP_BulkRNA/S30_normalised.rds"
out_dir <- "../04-Outputs/S43_Combined_Sensitivity"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

pyrim_genes   <- c("CAD", "DHODH", "UMPS", "UCK1", "CDA", "DCTPP1")
calcium_genes <- c("GNAQ", "ITPR2", "ITPR3", "ORAI1")
ecm_genes     <- c("COL1A1", "COL3A1", "FN1", "LUM", "DCN")

features <- c(pyrim_score = "Pyrimidine score", P2RY6 = "P2RY6",
              calcium_score = "Calcium score", PDGFRB = "PDGFRB",
              ecm_score = "ECM score")

methods <- c(ComBat_binary     = "ComBat, binary covariate (S40)",
             Rank_in_sample    = "Within-sample rank, no correction",
             ComBat_bridge_Adv = "ComBat, DKD bridge (Adv. DKD = KPMP DKD)",
             ComBat_bridge_All = "ComBat, DKD bridge (all DKD = KPMP DKD)")

condition_levels <- c("Normal", "Early_DKD", "Advanced_DKD", "AKI", "CKD")
condition_labels <- c("Normal", "Early DKD", "Adv. DKD", "AKI", "CKD")
condition_colors <- c(
  "Normal"       = "#2166ac",
  "Early_DKD"    = "#f4a582",
  "Advanced_DKD" = "#d6604d",
  "AKI"          = "#e07b39",
  "CKD"          = "#756bb1"
)

# Normal as reference (as S40) plus the within-KPMP contrast as a control
comparisons_main <- list(
  c("Normal", "Early_DKD"),
  c("Normal", "Advanced_DKD"),
  c("Normal", "AKI"),
  c("Normal", "CKD"),
  c("AKI", "CKD")
)

cor_pairs <- list(
  c("pyrim_score", "ecm_score"),
  c("pyrim_score", "P2RY6"),
  c("pyrim_score", "PDGFRB"),
  c("P2RY6", "calcium_score"),
  c("P2RY6", "PDGFRB"),
  c("P2RY6", "ecm_score")
)

dpi <- 300

# =========================================================================== #
# HELPER FUNCTIONS ----
# =========================================================================== #

signif_stars <- function(p) {
  ifelse(is.na(p), "n/a",
   ifelse(p < 0.001, "***",
    ifelse(p < 0.01, "**",
     ifelse(p < 0.05, "*", "ns"))))
}

module_score <- function(mat, genes) {
  found <- intersect(genes, rownames(mat))
  colMeans(mat[found, , drop = FALSE])
}

# Scores on an expression matrix (mean of log values, as S40)
score_expr <- function(mat) {
  data.frame(sample_id     = colnames(mat),
             pyrim_score   = module_score(mat, pyrim_genes),
             P2RY6         = as.numeric(mat["P2RY6", ]),
             calcium_score = module_score(mat, calcium_genes),
             PDGFRB        = as.numeric(mat["PDGFRB", ]),
             ecm_score     = module_score(mat, ecm_genes))
}

# Within-sample rank scores (0-1; higher = more highly expressed in the sample)
score_rank <- function(mat) {
  rk <- apply(mat, 2, rank, ties.method = "average") / nrow(mat)
  score_expr(rk)
}

# Wilcoxon + BH (per feature, across comparisons) + Cliff's delta
group_tests <- function(df, value_col) {
  res <- lapply(comparisons_main, function(cmp) {
    a <- df[[value_col]][df$condition == cmp[1]]
    b <- df[[value_col]][df$condition == cmp[2]]
    wt <- wilcox.test(b, a, exact = FALSE)
    data.frame(group1 = cmp[1], group2 = cmp[2],
               n_group1 = length(a), n_group2 = length(b),
               cliffs_delta = 2 * unname(wt$statistic) / (length(a) * length(b)) - 1,
               p_value = wt$p.value)
  })
  res <- bind_rows(res)
  res$p_adj  <- p.adjust(res$p_value, method = "BH")
  res$signif <- signif_stars(res$p_adj)
  res
}

pooled_cors <- function(df) {
  res <- bind_rows(lapply(cor_pairs, function(pr) {
    ct <- cor.test(df[[pr[1]]], df[[pr[2]]], method = "spearman", exact = FALSE)
    data.frame(x_var = pr[1], y_var = pr[2], rho = unname(ct$estimate),
               p = ct$p.value, n = nrow(df))
  }))
  res$p_adj  <- p.adjust(res$p, method = "BH")
  res$signif <- signif_stars(res$p_adj)
  res
}

run_combat <- function(mat, batch, mod) {
  ComBat(dat = mat, batch = batch, mod = mod,
         par.prior = TRUE, prior.plots = FALSE)
}

# =========================================================================== #
# 1. LOAD AND MERGE (as S40) ----
# =========================================================================== #

s20 <- readRDS(s20_rds)
s30 <- readRDS(s30_rds)

fan_meta <- s20$meta %>%
  transmute(sample_id = sample, condition = as.character(condition),
            cohort = "Fan", adjudicated = NA_character_)
kpmp_meta <- s30$meta %>%
  transmute(sample_id, condition = as.character(enrollment_category),
            cohort = "KPMP", adjudicated = primary_adjudicated_category)

common_genes <- intersect(rownames(s20$expr_mat), rownames(s30$expr_mat))
combined_mat <- cbind(s20$expr_mat[common_genes, ], s30$expr_mat[common_genes, ])

meta <- bind_rows(fan_meta, kpmp_meta)
meta$condition <- factor(meta$condition, levels = condition_levels)
stopifnot(identical(colnames(combined_mat), meta$sample_id))

meta$kpmp_dkd <- meta$cohort == "KPMP" & meta$condition == "CKD" &
  meta$adjudicated %in% "Diabetic Kidney Disease"

cat("Common genes:", length(common_genes), "\n")
print(table(meta$condition, meta$cohort))
cat("KPMP CKD with adjudicated DKD:", sum(meta$kpmp_dkd), "\n")

# =========================================================================== #
# 2. SCORES PER METHOD ----
# =========================================================================== #

batch <- meta$cohort

# . ComBat, binary covariate (S40) ----
mod_binary <- model.matrix(~ as.integer(meta$condition == "Normal"))
mat_binary <- run_combat(combined_mat, batch, mod_binary)

# . ComBat, DKD bridge ----
# Shared group: Fan DKD (advanced, or all) + KPMP CKD with adjudicated DKD.
# Remaining groups exist in one cohort only; the design is full rank because
# the cohort effect is identified by the shared group.
bridge_group <- function(fan_dkd_levels) {
  g <- as.character(meta$condition)
  g[meta$cohort == "Fan" & g %in% fan_dkd_levels] <- "DKD_shared"
  g[meta$kpmp_dkd] <- "DKD_shared"
  g[meta$cohort == "KPMP" & g == "CKD" & !meta$kpmp_dkd] <- "CKD_nonDKD"
  factor(g)
}
grp_adv <- bridge_group("Advanced_DKD")
grp_all <- bridge_group(c("Early_DKD", "Advanced_DKD"))
print(table(grp_adv, batch))
print(table(grp_all, batch))
mat_bridge_adv <- run_combat(combined_mat, batch, model.matrix(~ grp_adv))
mat_bridge_all <- run_combat(combined_mat, batch, model.matrix(~ grp_all))

scores <- bind_rows(
  score_expr(mat_binary)     %>% mutate(method = "ComBat_binary"),
  score_rank(combined_mat)   %>% mutate(method = "Rank_in_sample"),
  score_expr(mat_bridge_adv) %>% mutate(method = "ComBat_bridge_Adv"),
  score_expr(mat_bridge_all) %>% mutate(method = "ComBat_bridge_All")
) %>%
  left_join(meta, by = "sample_id")

# Datasets: all CKD (S40) and DKD-focused CKD (S42)
scores <- bind_rows(
  scores %>% mutate(dataset = "All"),
  scores %>% filter(condition != "CKD" | kpmp_dkd) %>%
    mutate(dataset = "DKD_focused")
)

# =========================================================================== #
# 3. GROUP COMPARISONS ----
# =========================================================================== #

tbl_groups <- scores %>%
  pivot_longer(names(features), names_to = "feature", values_to = "value") %>%
  group_by(dataset, method, feature) %>%
  group_modify(~ group_tests(.x, "value")) %>%
  ungroup() %>%
  mutate(feature = features[feature],
         cross_cohort = group1 == "Normal" & group2 %in% c("AKI", "CKD"))

# Check: ComBat_binary must reproduce the published S40 table
s40_tbl <- read.csv("../04-Outputs/S40_Combined_BulkRNA/Table_group_comparisons.csv")
chk <- tbl_groups %>%
  filter(dataset == "All", method == "ComBat_binary") %>%
  inner_join(s40_tbl, by = c(feature = "module", "group1", "group2"))
cat("S40 reproduced (max |p diff|):", max(abs(chk$p_value.x - chk$p_value.y)), "\n")
stopifnot(all(chk$signif.x == chk$signif.y))

write.csv(tbl_groups, file.path(out_dir, "Table_group_comparisons_sensitivity.csv"),
          row.names = FALSE)

# Compact wide table: "delta (stars)" per method
wide_groups <- tbl_groups %>%
  mutate(cell = sprintf("%+.2f %s", cliffs_delta, signif),
         comparison = paste(group2, "vs", group1)) %>%
  select(dataset, feature, comparison, method, cell) %>%
  pivot_wider(names_from = method, values_from = cell) %>%
  select(dataset, feature, comparison, all_of(names(methods)))
write.csv(wide_groups, file.path(out_dir, "Table_group_comparisons_summary.csv"),
          row.names = FALSE)

md <- c("# Combined bulk: group comparisons by method",
        "",
        "Cells: Cliff's delta (group vs Normal; AKI-CKD row: CKD vs AKI) and BH-adjusted significance (per feature, across the 5 comparisons).",
        "")
for (ds in unique(wide_groups$dataset)) {
  w <- wide_groups %>% filter(dataset == ds) %>% select(-dataset)
  md <- c(md, paste("##", ds), "",
          paste("|", paste(c("Feature", "Comparison", methods), collapse = " | "), "|"),
          paste(rep("|---", ncol(w)), collapse = ""), "|" ,
          apply(w, 1, function(r) paste("|", paste(r, collapse = " | "), "|")), "")
}
md <- md[md != "|"]
writeLines(md, file.path(out_dir, "Table_group_comparisons_summary.md"))

# =========================================================================== #
# 4. POOLED CORRELATIONS ----
# =========================================================================== #

tbl_cors <- scores %>%
  group_by(dataset, method) %>%
  group_modify(~ pooled_cors(.x)) %>%
  ungroup()
write.csv(tbl_cors, file.path(out_dir, "Table_pooled_correlations_sensitivity.csv"),
          row.names = FALSE)

# =========================================================================== #
# 5. FIGURES ----
# =========================================================================== #

# . Boxplots: rows = method, columns = feature ----
for (ds in c("All", "DKD_focused")) {
  df <- scores %>%
    filter(dataset == ds) %>%
    pivot_longer(names(features), names_to = "feature", values_to = "value") %>%
    mutate(feature = factor(features[feature], levels = features),
           method  = factor(methods[method], levels = methods))
  p <- ggplot(df, aes(condition, value, fill = condition)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(aes(shape = cohort), width = 0.15, size = 0.9, alpha = 0.7) +
    facet_grid(method ~ feature, scales = "free_y",
               labeller = labeller(method = label_wrap_gen(22))) +
    scale_fill_manual(values = condition_colors) +
    scale_shape_manual(values = c(Fan = 16, KPMP = 17)) +
    scale_x_discrete(labels = setNames(condition_labels, condition_levels)) +
    labs(x = NULL, y = "Score (method-specific scale)",
         title = paste("Combined bulk RNA-seq: sensitivity to batch handling -", ds)) +
    theme_bw(base_size = 10) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 35, hjust = 1),
          strip.text = element_text(face = "bold", size = 8))
  ggsave(file.path(out_dir, paste0("Fig_sensitivity_boxplots_", ds, ".png")),
         p, width = 13, height = 10, dpi = dpi)
}

# . Effect-size overview: Cliff's delta by method ----
df_eff <- tbl_groups %>%
  mutate(comparison = factor(paste(group2, "vs", group1),
                             levels = unique(paste(sapply(comparisons_main, `[`, 2),
                                                   "vs", sapply(comparisons_main, `[`, 1)))),
         feature = factor(feature, levels = features),
         method  = factor(methods[method], levels = methods),
         label   = ifelse(signif == "ns", "", signif))
p_eff <- ggplot(df_eff, aes(comparison, cliffs_delta, fill = method)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_col(position = position_dodge(0.85), width = 0.8) +
  geom_text(aes(label = label, y = cliffs_delta + sign(cliffs_delta) * 0.08),
            position = position_dodge(0.85), size = 2.6) +
  facet_grid(dataset ~ feature) +
  scale_fill_manual(values = c("#4d4d4d", "#1b9e77", "#7570b3", "#d95f02")) +
  coord_cartesian(ylim = c(-1.1, 1.1)) +
  labs(x = NULL, y = "Cliff's delta", fill = NULL,
       title = "Group differences by batch-handling method (BH-adjusted stars)") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 35, hjust = 1),
        strip.text = element_text(face = "bold"))
ggsave(file.path(out_dir, "Fig_sensitivity_effect_sizes.png"),
       p_eff, width = 14, height = 7.5, dpi = dpi)

saveRDS(list(scores = scores, groups = tbl_groups, cors = tbl_cors),
        file.path(out_dir, "S43_sensitivity.rds"))

cat("S43 done.\n")
