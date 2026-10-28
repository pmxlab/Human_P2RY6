# =========================================================================== #
# S84_Module_Genes.R
# Expression of the individual module genes in bulk, scRNA-seq and snRNA-seq
# Reviewer comments R1.09, R2.12 (Life, life-4655118)
#
# For each module score (pyrimidine, calcium signaling, ECM) one figure with
# three vertically stacked facets, each on its own y-axis:
#   1. Fan DKD bulk RNA-seq   (Normal, Early DKD, Advanced DKD)
#   2. KPMP scRNA-seq         (Healthy, AKI, CKD; pseudobulk per participant)
#   3. KPMP snRNA-seq         (Healthy, AKI, CKD; pseudobulk per participant)
# Within each condition, one colour-coded box per module gene; no individual
# points. Shows whether a single gene drives a module score.
#
# Compartments as in the module-score panels: pyrimidine in proximal tubule
# (PT); calcium and ECM in stroma (scRNA-seq) and fibroblasts (snRNA-seq).
# Single-cell data: all conditions (as panels c and e of Figures 1, 3, 5).
#
# Data:
#   S20_normalised.rds  Fan: log2 quantile-normalized values (Fan et al.)
#   S50_pseudobulk.rds  scRNA-seq: mean LogNormalize value per participant
#                       (pb_pt: pyrimidine genes; pb_fib: calcium/ECM genes)
#   S60_pseudobulk.rds  snRNA-seq: as S50 (pb_pt, pb_fib)
# Test: two-sided Wilcoxon rank-sum test of each disease group vs the
#   reference group of the dataset (Normal / Healthy; R defaults, as in S80);
#   BH-adjusted per gene and dataset over the two comparisons, as the
#   module-score panels. Effect size: Cliff's delta (positive = higher than
#   reference) and difference of group medians (expression units of the facet).
#
# Leave-one-gene-out (R1.09): each score recomputed without one gene at a
#   time (unweighted mean of the remaining genes, as the full score) and
#   re-tested with the comparisons of the score panels (Figures 1, 3, 5):
#   Fan (b), scRNA-seq (c, d) and snRNA-seq (e, f), i.e. all conditions and
#   DKD-focused (CKD restricted to adjudicated DKD, filters as in S80).
#   The full scores are checked against the stars drawn by S80.
#
# Outputs: 04-Outputs/S84_Module_Genes/
#   Fig_genes_<module>.png/.pdf, Figure_S5a-c.png (Supplementary Figure S5),
#   Table_gene_tests.csv/.md,
#   Table_LOGO_all.csv, Table_LOGO_summary.csv/.md
# =========================================================================== #

rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)

# =========================================================================== #
# SETTINGS ----
# =========================================================================== #

s20_rds <- "../04-Outputs/S20_FanDKD_BulkRNA/S20_normalised.rds"
s50_rds <- "../04-Outputs/S50_KPMP_scRNA/S50_pseudobulk.rds"
s60_rds <- "../04-Outputs/S60_KPMP_snRNA/S60_pseudobulk.rds"
sc_clin <- "../../KPMP_scRNA/OpenAccessClinicalData/20251205_OpenAccessClinicalData.csv"
sn_clin <- "../../KPMP_snRNA/6d3047c2-588f-4ad4-8eba-f4c46a1e5c14_20251205_OpenAccessClinicalData/20251205_OpenAccessClinicalData.csv"
s80_stars <- "../04-Outputs/S80_Manuscript_Figures/Check_significance_stars.csv"
out_dir <- "../04-Outputs/S84_Module_Genes"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

modules <- list(
  pyrimidine = list(genes = c("CAD", "DHODH", "UMPS", "UCK1", "CDA", "DCTPP1"),
                    title = "Pyrimidine module genes", sc_comp = "PT",
                    sn_comp = "PT"),
  calcium    = list(genes = c("GNAQ", "ITPR2", "ITPR3", "ORAI1"),
                    title = "Calcium signaling module genes", sc_comp = "stroma",
                    sn_comp = "fibroblasts"),
  ecm        = list(genes = c("COL1A1", "COL3A1", "FN1", "LUM", "DCN"),
                    title = "ECM module genes", sc_comp = "stroma",
                    sn_comp = "fibroblasts")
)

datasets <- c(fan = "Fan DKD bulk RNA-seq",
              sc  = "KPMP scRNA-seq",
              sn  = "KPMP snRNA-seq")

# Conditions of all datasets (each facet shows its own three)
condition_levels <- c("Normal", "Early DKD", "Advanced DKD", "Healthy", "AKI", "CKD")
reference        <- c(fan = "Normal", sc = "Healthy", sn = "Healthy")
disease          <- list(fan = c("Early DKD", "Advanced DKD"), sc = c("AKI", "CKD"),
                         sn = c("AKI", "CKD"))

sc_map <- c("Healthy Reference Tissue" = "Healthy", "Acute Kidney Injury" = "AKI",
            "Chronic Kidney Disease" = "CKD")
sn_map <- c("Healthy Reference" = "Healthy", AKI = "AKI", CKD = "CKD")
fan_map <- c(Normal = "Normal", Early_DKD = "Early DKD", Advanced_DKD = "Advanced DKD")

# Categorical gene colours in fixed slot order (dataviz reference palette)
gene_palette <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300")

dpi         <- 600
dodge_width <- 0.8

# =========================================================================== #
# HELPER FUNCTIONS ----
# =========================================================================== #

p_stars <- function(p) {
  as.character(cut(p, c(-Inf, 1e-4, 1e-3, 1e-2, 0.05, Inf),
                   c("****", "***", "**", "*", "ns")))
}

# Tests of one gene in one dataset: each disease group vs the reference
gene_tests <- function(df, ds) {
  bind_rows(lapply(disease[[ds]], function(g2) {
    a <- df$expr[df$condition == reference[[ds]]]
    b <- df$expr[df$condition == g2]
    wt <- suppressWarnings(wilcox.test(b, a))
    w  <- suppressWarnings(wilcox.test(b, a, exact = FALSE))$statistic
    data.frame(group1 = reference[[ds]], group2 = g2,
               n_group1 = length(a), n_group2 = length(b),
               median_group1 = median(a), median_group2 = median(b),
               median_diff = median(b) - median(a),
               cliffs_delta = 2 * unname(w) / (length(a) * length(b)) - 1,
               p_value = wt$p.value)
  })) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH"),
           signif = p_stars(p_adj))
}

# Long format: one row per participant x gene
to_long <- function(df, id_col, cond_col, cond_map, genes, ds, module) {
  stopifnot(all(genes %in% names(df)))
  df %>%
    transmute(sample_id = .data[[id_col]],
              condition = unname(cond_map[as.character(.data[[cond_col]])]),
              across(all_of(genes))) %>%
    pivot_longer(all_of(genes), names_to = "gene", values_to = "expr") %>%
    mutate(dataset = ds, module = module)
}

# =========================================================================== #
# 1. DATA ----
# =========================================================================== #

s20 <- readRDS(s20_rds)
s50 <- readRDS(s50_rds)
s60 <- readRDS(s60_rds)

fan_wide <- data.frame(s20$meta[, c("sample", "condition")],
                       t(s20$expr_mat[unlist(lapply(modules, `[[`, "genes")),
                                      s20$meta$sample]),
                       check.names = FALSE)

long <- bind_rows(lapply(names(modules), function(m) {
  g  <- modules[[m]]$genes
  sc <- if (modules[[m]]$sc_comp == "PT") s50$pb_pt else s50$pb_fib
  sn <- if (modules[[m]]$sn_comp == "PT") s60$pb_pt else s60$pb_fib
  bind_rows(
    to_long(fan_wide, "sample", "condition", fan_map, g, "fan", m),
    to_long(sc, "SampleID", "EnrollmentCategory", sc_map, g, "sc", m),
    to_long(sn, "patient", "enrollment", sn_map, g, "sn", m)
  )
})) %>%
  filter(!is.na(condition)) %>%
  mutate(condition = factor(condition, levels = condition_levels))

print(long %>% distinct(module, dataset, sample_id, condition) %>%
        count(module, dataset, condition) %>%
        pivot_wider(names_from = condition, values_from = n))

# =========================================================================== #
# 2. TESTS PER GENE ----
# =========================================================================== #

tests <- long %>%
  group_by(module, dataset, gene) %>%
  group_modify(~ gene_tests(.x, .y$dataset)) %>%
  ungroup() %>%
  mutate(module  = factor(module, levels = names(modules)),
         dataset = factor(dataset, levels = names(datasets)),
         gene    = factor(gene, levels = unlist(lapply(modules, `[[`, "genes")))) %>%
  arrange(module, dataset, gene, factor(group2, levels = condition_levels))

write.csv(tests %>% mutate(dataset = datasets[as.character(dataset)]),
          file.path(out_dir, "Table_gene_tests.csv"), row.names = FALSE)

md <- c("# Module genes vs reference group, per dataset",
        "",
        "Two-sided Wilcoxon rank-sum test; BH-adjusted per gene and dataset over the two comparisons.",
        "Cliff's delta > 0: higher than reference (Normal / Healthy).",
        "Median difference in the expression units of the dataset (Fan: log2; KPMP: ln(1 + normalized count), pseudobulk mean).",
        "",
        "| Module | Dataset | Gene | Comparison | Median diff | Cliff's delta | p | p (BH) | |",
        "|---|---|---|---|---|---|---|---|---|",
        sprintf("| %s | %s | %s | %s vs %s | %+.3g | %+.2f | %.2g | %.2g | %s |",
                tests$module, datasets[as.character(tests$dataset)], tests$gene,
                tests$group2, tests$group1, tests$median_diff, tests$cliffs_delta,
                tests$p_value, tests$p_adj, tests$signif))
writeLines(md, file.path(out_dir, "Table_gene_tests.md"))

# =========================================================================== #
# 3. FIGURES ----
# =========================================================================== #
# Stars above each disease-group box: that gene vs the reference group of the
# facet (BH per gene); "ns" is not drawn to keep the grouped layout readable.

facet_label <- function(mod) {
  m <- modules[[mod]]
  c(fan = "Fan DKD bulk RNA-seq (log2, quantile-normalized)",
    sc  = sprintf("KPMP scRNA-seq, %s (pseudobulk, mean log-normalized expression)",
                  m$sc_comp),
    sn  = sprintf("KPMP snRNA-seq, %s (pseudobulk, mean log-normalized expression)",
                  m$sn_comp))
}

gene_box <- function(mod) {
  genes <- modules[[mod]]$genes
  flab  <- facet_label(mod)
  df    <- long %>%
    filter(module == mod) %>%
    mutate(gene  = factor(gene, levels = genes),
           facet = factor(flab[dataset], levels = flab))

  # x labels with group sizes per facet
  ns <- df %>% distinct(dataset, sample_id, condition) %>% count(dataset, condition)
  x_lab <- setNames(sprintf("%s\n(n=%d)", ns$condition, ns$n),
                    paste(ns$dataset, ns$condition))
  df <- df %>% mutate(x = factor(paste(dataset, condition),
                                 levels = paste(rep(names(datasets), each = 6),
                                                condition_levels)))

  # star positions: upper whisker of each box + 3% of the facet range;
  # one label row per box (empty for reference and "ns"), so that the dodge
  # positions of the labels match those of the boxes
  tops <- df %>%
    group_by(dataset, facet, x, condition, gene) %>%
    summarise(top = boxplot.stats(expr)$stats[5], .groups = "drop") %>%
    group_by(dataset) %>%
    mutate(top = top + 0.03 * diff(range(df$expr[df$dataset == dataset[1]]))) %>%
    ungroup()
  lab <- tests %>%
    filter(module == mod) %>%
    transmute(dataset = as.character(dataset),
              condition = factor(group2, levels = condition_levels),
              gene = factor(as.character(gene), levels = genes),
              signif = ifelse(signif == "ns", "", signif)) %>%
    right_join(tops, by = c("dataset", "condition", "gene")) %>%
    mutate(signif = coalesce(signif, ""))

  ggplot(df, aes(x = x, y = expr, fill = gene)) +
    geom_boxplot(position = position_dodge(dodge_width), width = 0.7,
                 outliers = FALSE, linewidth = 0.35, colour = "grey25") +
    geom_text(data = lab, aes(y = top, label = signif, group = gene),
              position = position_dodge(dodge_width), vjust = 0, size = 3) +
    facet_wrap(~ facet, ncol = 1, scales = "free") +
    scale_fill_manual(values = setNames(gene_palette[seq_along(genes)], genes),
                      name = NULL) +
    scale_x_discrete(labels = x_lab) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.10))) +
    labs(x = NULL, y = "Gene expression", title = modules[[mod]]$title) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          plot.title.position = "plot",
          strip.text = element_text(face = "bold", hjust = 0),
          strip.background = element_rect(fill = "grey92", colour = NA),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "right")
}

for (m in names(modules)) {
  p <- gene_box(m)
  ggsave(file.path(out_dir, paste0("Fig_genes_", m, ".png")), p,
         width = 8, height = 11, dpi = dpi, bg = "white")
  ggsave(file.path(out_dir, paste0("Fig_genes_", m, ".pdf")), p,
         width = 8, height = 11, device = cairo_pdf)
  # Supplementary Figure S5 (a) pyrimidine, (b) calcium, (c) ECM
  ps5 <- p + labs(tag = sprintf("(%s)", letters[match(m, names(modules))])) +
    theme(plot.tag = element_text(face = "bold", size = 16))
  ggsave(file.path(out_dir, sprintf("Figure_S5%s.png", letters[match(m, names(modules))])),
         ps5, width = 8, height = 11, dpi = dpi, bg = "white")
}

# =========================================================================== #
# 4. LEAVE-ONE-GENE-OUT SENSITIVITY ANALYSIS ----
# =========================================================================== #

# . 4a. Per-gene data of the DKD-focused sets (filters as in S80) ----
# PT: each participant's own adjudicated category from the clinical file;
# stroma / fibroblasts: adjudicated column of the pseudobulk.
dkd_keep <- function(cond, adj) !(cond == "CKD" & (is.na(adj) | adj != "Diabetic Kidney Disease"))
clin_sc <- read.csv(sc_clin)
clin_sn <- read.csv(sn_clin)
adj_sc  <- setNames(clin_sc$Primary.Adjudicated.Category, clin_sc$Participant.ID)
adj_sn  <- setNames(clin_sn$Primary.Adjudicated.Category, clin_sn$Participant.ID)

long_dkd <- bind_rows(lapply(names(modules), function(m) {
  g  <- modules[[m]]$genes
  sc <- if (modules[[m]]$sc_comp == "PT") {
    s50$pb_pt %>% mutate(adjudicated = unname(adj_sc[SampleID]))
  } else s50$pb_fib
  sn <- if (modules[[m]]$sn_comp == "PT") {
    s60$pb_pt %>% mutate(adjudicated = unname(adj_sn[patient]))
  } else s60$pb_fib
  bind_rows(
    to_long(sc %>% filter(dkd_keep(unname(sc_map[as.character(EnrollmentCategory)]),
                                   adjudicated)),
            "SampleID", "EnrollmentCategory", sc_map, g, "sc", m),
    to_long(sn %>% filter(dkd_keep(unname(sn_map[as.character(enrollment)]),
                                   adjudicated)),
            "patient", "enrollment", sn_map, g, "sn", m)
  )
})) %>%
  mutate(condition = factor(condition, levels = condition_levels))

logo_data <- bind_rows(
  long %>% mutate(analysis = "All conditions"),
  long_dkd %>% mutate(analysis = "DKD-focused")
)
# Fan has no CKD without DKD: its comparisons are listed once
logo_data <- logo_data %>% filter(!(dataset == "fan" & analysis == "DKD-focused"))

# . 4b. Scores without one gene and tests ----
score_tests <- function(df, ds, genes_used) {
  sc <- df %>%
    filter(gene %in% genes_used) %>%
    group_by(sample_id, condition) %>%
    summarise(expr = mean(expr), .groups = "drop")
  gene_tests(sc, ds)
}

logo <- logo_data %>%
  group_by(module, analysis, dataset) %>%
  group_modify(function(df, key) {
    genes <- modules[[key$module]]$genes
    bind_rows(
      score_tests(df, key$dataset, genes) %>% mutate(left_out = "none (full score)"),
      bind_rows(lapply(genes, function(g)
        score_tests(df, key$dataset, setdiff(genes, g)) %>% mutate(left_out = g)))
    )
  }) %>%
  ungroup() %>%
  mutate(module  = factor(module, levels = names(modules)),
         dataset = factor(dataset, levels = names(datasets))) %>%
  arrange(module, dataset, analysis, factor(group2, levels = condition_levels))

# . 4c. Check: full scores reproduce the stars of the S80 score panels ----
# S80 panel letters: (c) scRNA all, (d) scRNA DKD, (e) snRNA all, (f) snRNA DKD
fig_no <- c(pyrimidine = 1, calcium = 3, ecm = 5)
s80 <- read.csv(s80_stars) %>%
  mutate(group2 = recode(group2, "Acute Kidney Injury" = "AKI",
                         "Chronic Kidney Disease" = "CKD"))
chk <- logo %>%
  filter(left_out == "none (full score)", dataset != "fan") %>%
  mutate(Panel = paste0(fig_no[as.character(module)],
                        ifelse(dataset == "sc",
                               ifelse(analysis == "All conditions", "c", "d"),
                               ifelse(analysis == "All conditions", "e", "f")))) %>%
  inner_join(s80 %>% select(Panel, group2, drawn), by = c("Panel", "group2"))
stopifnot(nrow(chk) == 24, all(chk$signif == chk$drawn))
cat("Full scores reproduce the S80 stars of", nrow(chk), "single-cell comparisons.\n")

write.csv(logo %>% mutate(dataset = datasets[as.character(dataset)]),
          file.path(out_dir, "Table_LOGO_all.csv"), row.names = FALSE)

# . 4d. Summary: one row per score x dataset x comparison ----
logo_summary <- logo %>%
  group_by(module, dataset, analysis, group1, group2) %>%
  summarise(
    full_delta   = cliffs_delta[left_out == "none (full score)"],
    full_p_adj   = p_adj[left_out == "none (full score)"],
    full_signif  = signif[left_out == "none (full score)"],
    logo_delta_min = min(cliffs_delta[left_out != "none (full score)"]),
    logo_delta_max = max(cliffs_delta[left_out != "none (full score)"]),
    logo_p_adj_max = max(p_adj[left_out != "none (full score)"]),
    logo_p_adj_min = min(p_adj[left_out != "none (full score)"]),
    # genes whose omission moves the BH-adjusted p across 0.05:
    # lost = significant full score becomes ns; gained = ns becomes significant
    lost_without = {
      lo <- left_out != "none (full score)"
      g  <- if (p_adj[!lo] < 0.05) left_out[lo][p_adj[lo] >= 0.05] else character()
      if (length(g)) paste(g, collapse = ", ") else "-"
    },
    gained_without = {
      lo <- left_out != "none (full score)"
      g  <- if (p_adj[!lo] >= 0.05) left_out[lo][p_adj[lo] < 0.05] else character()
      if (length(g)) paste(g, collapse = ", ") else "-"
    },
    .groups = "drop")

write.csv(logo_summary %>% mutate(dataset = datasets[as.character(dataset)]),
          file.path(out_dir, "Table_LOGO_summary.csv"), row.names = FALSE)

md <- c("# Leave-one-gene-out sensitivity of the module scores",
        "",
        "Each score recomputed without one gene at a time; same test as the score panels",
        "(two-sided Wilcoxon, BH over the two comparisons of a panel). Cliff's delta > 0:",
        "higher than reference. 'Lost without': significant full score becomes ns when",
        "this gene is omitted; 'Gained without': ns full score becomes significant.",
        "",
        "| Module | Dataset | Analysis | Comparison | Full score: delta (stars) | Leave-one-out: delta range | Leave-one-out: BH p range | Lost without | Gained without |",
        "|---|---|---|---|---|---|---|---|---|",
        sprintf("| %s | %s | %s | %s vs %s | %+.2f (%s) | %+.2f to %+.2f | %.2g to %.2g | %s | %s |",
                logo_summary$module, datasets[as.character(logo_summary$dataset)],
                logo_summary$analysis, logo_summary$group2, logo_summary$group1,
                logo_summary$full_delta, logo_summary$full_signif,
                logo_summary$logo_delta_min, logo_summary$logo_delta_max,
                logo_summary$logo_p_adj_min, logo_summary$logo_p_adj_max,
                logo_summary$lost_without, logo_summary$gained_without))
writeLines(md, file.path(out_dir, "Table_LOGO_summary.md"))

cat("S84 done.\n")
