# ============================================================================#
# S89_Group_Comparisons.R
# Exact BH-adjusted p-values and effect sizes for all group comparisons of
# Figures 1-5 and Supplementary Figures S3, S4, S6 and S8 (S2, S3, S5 and S7
# before ms_rev1_v48_AF; Supplementary Table S9; Table S8 before the supplements were renumbered in order of
# citation, Revision 1). Revision 1 (ms_rev1_v42_AF): panels (a)-(c) of the
# former Figures 2-5 are Supplementary Figures S2, S3, S5 and S7, and panels
# (d, e) are Figures 2-5 (a, b); column 'label' gives the manuscript panel.
# Reviewer comments R2.18, R3.05 (Life, life-4655118)
#
# The panels are rebuilt from S80_Manuscript_Figures.R (as in S81/S85, no
# files saved), so every statistic is computed on exactly the plotted data.
# Per panel and disease group vs the control group (Normal / Healthy):
#   n and median (IQR) of each group
#   Cliff's delta (P(disease > control) - P(disease < control))
#   two-sided Wilcoxon rank-sum p, BH-adjusted over the comparisons of the
#   panel - the same test and adjustment as the significance stars (S80
#   bh_test); the adjusted p-values are checked against the stars.
#
# Outputs: 04-Outputs/S89_Group_Comparisons/
#   Table_S8_group_comparisons.csv / .md
# ============================================================================#

rm(list = ls())

library(dplyr)
library(rlang)

# ============================================================================#
# SETTINGS ----
# ============================================================================#

s80_file <- "S80_Manuscript_Figures.R"
out_dir  <- "../04-Outputs/S89_Group_Comparisons"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Figure 1 has six panels: (a) bulk PT marker score, (b) bulk pyrimidine score,
# (c)-(f) KPMP; Figures 2-5 have five panels (a)-(e)
panel_ids <- c(paste0("f1", letters[1:6]), paste0("f", rep(2:5, each = 5), letters[1:5]))
measure_of <- c("1" = "Pyrimidine score (PT)", "1a" = "PT marker score", "2" = "P2RY6", "3" = "Calcium signaling score",
                "4" = "PDGFRB", "5" = "ECM score")
dataset_of <- c(a = "Fan DKD bulk RNA-seq", b = "KPMP scRNA-seq", c = "KPMP scRNA-seq, DKD-focused",
                d = "KPMP snRNA-seq", e = "KPMP snRNA-seq, DKD-focused")
cond_map <- c(Normal = "Normal", Early_DKD = "Early DKD", Advanced_DKD = "Advanced DKD",
              "Healthy Reference Tissue" = "Healthy", "Healthy Reference" = "Healthy",
              "Acute Kidney Injury" = "AKI", AKI = "AKI",
              "Chronic Kidney Disease" = "CKD", CKD = "CKD")

# ============================================================================#
# HELPER FUNCTIONS ----
# ============================================================================#

med_iqr <- function(v) {
  q <- quantile(v, c(0.5, 0.25, 0.75), names = FALSE)
  d <- if (abs(q[1]) >= 10) 2 else 3
  sprintf(paste0("%.", d, "f (%.", d, "f–%.", d, "f)"), q[1], q[2], q[3])
}

cliffs_delta <- function(x, y) {          # x = disease, y = control
  mean(sign(outer(x, y, "-")))
}

# ============================================================================#
# 1. BUILD THE S80 PANELS WITHOUT SAVING ----
# ============================================================================#

code <- readLines(s80_file)
code <- code[!grepl("^rm\\(list = ls\\(\\)\\)", code)]
first_fig <- grep("^# FIGURE 1 ", code)[1]
stopifnot(!is.na(first_fig))
code <- append(code, "save_figure <- function(...) invisible(NULL)", after = first_fig - 2)
env <- new.env()
eval(parse(text = code), envir = env)
stopifnot(all(panel_ids %in% ls(env)))

# ============================================================================#
# 2. COMPARISONS PER PANEL ----
# ============================================================================#

# manuscript panel labels (Revision 1, ms_rev1_v42_AF; renumbered after
# ms_rev1_v48_AF, where the module gene figure became Figure S2)
supp_of <- c("2" = "S3", "3" = "S4", "4" = "S6", "5" = "S8")
panel_label <- function(id) {
  fig <- substr(id, 2, 2); let <- substr(id, 3, 3)
  if (fig == "1") return(paste0("1", if (let == "a") "a" else chartr("bcdef", "abcde", let)))
  if (let %in% c("a", "b", "c")) paste0(supp_of[[fig]], let) else paste0(fig, chartr("de", "ab", let))
}

rows <- lapply(panel_ids, function(id) {
  p <- get(id, envir = env)
  d <- p$data
  x <- eval_tidy(p$mapping$x, d); y <- eval_tidy(p$mapping$y, d)
  keep <- !is.na(y)
  x <- droplevels(factor(x[keep])); y <- y[keep]
  lv <- levels(x)
  df <- data.frame(g = x, v = y)
  cmp <- lapply(lv[-1], function(l) c(lv[1], l))
  st <- env$bh_test(df, "g", "v", cmp)            # same test and adjustment as the stars
  fig <- substr(id, 2, 2); let <- substr(id, 3, 3)
  if (fig == "1") let <- if (let == "a") "a" else chartr("bcdef", "abcde", let)
  meas <- if (id == "f1a") measure_of[["1a"]] else measure_of[[fig]]
  bind_rows(lapply(seq_along(cmp), function(k) {
    ctl <- y[x == lv[1]]; dis <- y[x == cmp[[k]][2]]
    data.frame(panel = sub("^f", "", id), label = panel_label(id), measure = meas,
               dataset = dataset_of[[let]],
               control = cond_map[lv[1]], disease = cond_map[cmp[[k]][2]],
               n_control = length(ctl), n_disease = length(dis),
               median_control = med_iqr(ctl), median_disease = med_iqr(dis),
               cliffs_delta = cliffs_delta(dis, ctl),
               p = st$p[k], p_BH = st$p_adj[k], signif_BH = st$label[k],
               row.names = NULL)
  }))
})
tab <- bind_rows(rows)

# check against the stars drawn in the figures (S81)
chk <- read.csv("../04-Outputs/S80_Manuscript_Figures/Check_significance_stars.csv")
chk$disease <- unname(cond_map[chk$group2])
m <- merge(tab, chk, by.x = c("panel", "disease"), by.y = c("Panel", "disease"))
stopifnot(nrow(m) == nrow(tab), all(abs(m$p_BH.x - m$p_BH.y) <= 1e-3 * m$p_BH.y), all(m$signif_BH == m$drawn))   # CSV: 4 significant digits

print(tab, digits = 3)
write.csv(tab, file.path(out_dir, "Table_S8_group_comparisons.csv"), row.names = FALSE)

fmt_p <- function(p) ifelse(p < 1e-4, "< 0.0001", ifelse(p < 0.01, sprintf("%.4f", p), sprintf("%.3f", p)))
md_tab <- tab %>%
  transmute(Panel = label, Measure = measure, Dataset = dataset,
            Comparison = paste(disease, "vs", control),
            `n` = paste(n_disease, "/", n_control),
            `Median (IQR), disease` = median_disease, `Median (IQR), control` = median_control,
            `Cliff's delta` = sprintf("%.2f", cliffs_delta), `BH-adjusted p` = fmt_p(p_BH))
md <- c("# Table S9. Group comparisons of Figures 1-5 and Supplementary Figures S3, S4, S6 and S8", "",
        paste0("| ", paste(names(md_tab), collapse = " | "), " |"),
        paste0(paste(rep("|---", ncol(md_tab)), collapse = ""), "|"),
        apply(md_tab, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
writeLines(md, file.path(out_dir, "Table_S8_group_comparisons.md"))
