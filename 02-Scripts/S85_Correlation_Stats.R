# ============================================================================#
# S85_Correlation_Stats.R
# Correlation statistics for all correlation panels of the manuscript
# (Figure 6 and Supplementary Figures S1-S3): pooled, adjusted for condition
# and within condition, with multiple-testing correction.
# Reviewer comments R1.07, R2.04, R2.07 (Life, life-4655118)
#
# The panels are rebuilt from S80_Manuscript_Figures.R (as in S81, no files
# saved), so every statistic is computed on exactly the plotted data.
#
# Per panel:
#   pooled    : Spearman rho over all participants (as shown in the panel)
#   adjusted  : partial Spearman correlation adjusted for condition, i.e. the
#               Pearson correlation of the rank residuals of x and y after
#               linear regression on condition (df = n - 2 - (conditions - 1))
#   within    : Spearman rho within each condition (n >= 5)
# Multiple testing: Benjamini-Hochberg within each family, i.e. per figure
#   and analysis type (pooled: the five panels of one figure, as shown in the
#   figure labels of S80; adjusted: the five panels of one figure; within:
#   all within-condition correlations of one figure).
# 95% CIs: Fisher z-transformation (as in S80).
#
# Outputs: 04-Outputs/S85_Correlation_Stats/
#   Table_correlations_all.csv     long table (all statistics)
#   Table_S4_correlations.csv/.md  compact supplementary table (Table S4)
# ============================================================================#

rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(rlang)

# ============================================================================#
# SETTINGS ----
# ============================================================================#

s80_file <- "S80_Manuscript_Figures.R"
out_dir  <- "../04-Outputs/S85_Correlation_Stats"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

min_n_within <- 5

panel_ids <- c(paste0("f6", letters[1:5]), paste0("fs", rep(1:3, each = 5), letters[1:5]))

pairs <- c("6"  = "PT pyrimidine score vs ECM score",
           "S1" = "P2RY6 vs calcium score",
           "S2" = "P2RY6 vs PDGFRB",
           "S3" = "P2RY6 vs ECM score")

dataset_of <- c(a = "Fan DKD bulk RNA-seq", b = "KPMP scRNA-seq", c = "KPMP scRNA-seq",
                d = "KPMP snRNA-seq", e = "KPMP snRNA-seq")
analysis_of <- c(a = "All conditions", b = "All conditions", c = "DKD-focused",
                 d = "All conditions", e = "DKD-focused")

cond_map <- c(Normal = "Normal", Early_DKD = "Early DKD", Advanced_DKD = "Advanced DKD",
              "Healthy Reference Tissue" = "Healthy", "Healthy Reference" = "Healthy",
              "Acute Kidney Injury" = "AKI", AKI = "AKI",
              "Chronic Kidney Disease" = "CKD", CKD = "CKD")
cond_order <- c("Normal", "Early DKD", "Advanced DKD", "Healthy", "AKI", "CKD")

# ============================================================================#
# HELPER FUNCTIONS ----
# ============================================================================#

fisher_ci <- function(r, n, k = 0) {
  if (is.na(r) || abs(r) >= 1 || n - 3 - k <= 0) return(c(NA_real_, NA_real_))
  z <- atanh(r); se <- 1 / sqrt(n - 3 - k)
  tanh(c(z - 1.96 * se, z + 1.96 * se))
}

spearman <- function(x, y) {
  ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  r  <- unname(ct$estimate)
  ci <- fisher_ci(r, length(x))
  data.frame(n = length(x), rho = r, ci_lo = ci[1], ci_hi = ci[2], p = ct$p.value)
}

# Partial Spearman correlation adjusted for a categorical covariate
partial_spearman <- function(x, y, g) {
  g  <- droplevels(factor(g))
  k  <- nlevels(g) - 1
  rx <- resid(lm(rank(x) ~ g)); ry <- resid(lm(rank(y) ~ g))
  r  <- cor(rx, ry)
  n  <- length(x); df <- n - 2 - k
  tt <- r * sqrt(df / (1 - r^2))
  ci <- fisher_ci(r, n, k)
  data.frame(n = n, rho = r, ci_lo = ci[1], ci_hi = ci[2], p = 2 * pt(-abs(tt), df))
}

p_stars <- function(p) {
  as.character(cut(p, c(-Inf, 1e-4, 1e-3, 1e-2, 0.05, Inf),
                   c("****", "***", "**", "*", "ns")))
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

panel_data <- function(id) {
  p <- get(id, envir = env)
  d <- p$data
  x <- eval_tidy(p$mapping$x, d); y <- eval_tidy(p$mapping$y, d)
  cond_col <- intersect(c("condition", "EnrollmentCategory", "enrollment"), names(d))[1]
  keep <- !is.na(x) & !is.na(y)
  data.frame(x = x[keep], y = y[keep],
             condition = unname(cond_map[as.character(d[[cond_col]][keep])]))
}

# ============================================================================#
# 2. STATISTICS PER PANEL ----
# ============================================================================#

res <- bind_rows(lapply(panel_ids, function(id) {
  fig    <- sub("^f(s?)(\\d)([a-e])$", "\\1\\2", id); fig <- ifelse(grepl("^s", fig), toupper(fig), fig)
  letter <- sub(".*([a-e])$", "\\1", id)
  d <- panel_data(id)
  meta <- data.frame(figure = fig, panel = paste0(fig, letter), pair = pairs[[fig]],
                     dataset = dataset_of[[letter]], analysis = analysis_of[[letter]])
  pooled <- cbind(meta, type = "pooled", condition = "all", spearman(d$x, d$y))
  adj    <- cbind(meta, type = "adjusted for condition", condition = "all",
                  partial_spearman(d$x, d$y, d$condition))
  within <- bind_rows(lapply(split(d, d$condition), function(dd) {
    if (nrow(dd) < min_n_within) return(NULL)
    cbind(meta, type = "within condition", condition = dd$condition[1], spearman(dd$x, dd$y))
  }))
  bind_rows(pooled, adj, within)
}))

res <- res %>%
  group_by(figure, type) %>%
  mutate(p_BH = p.adjust(p, method = "BH")) %>%
  ungroup() %>%
  mutate(signif_BH = p_stars(p_BH),
         condition = factor(condition, levels = c("all", cond_order)),
         type = factor(type, levels = c("pooled", "adjusted for condition", "within condition"))) %>%
  arrange(factor(figure, levels = names(pairs)), panel, type, condition)

# check: pooled rho equals the value printed in the panel (S80 spearman_ci)
chk <- sapply(panel_ids, function(id) {
  d <- panel_data(id); env$spearman_ci(d$x, d$y)$rho
})
stopifnot(isTRUE(all.equal(unname(chk), res$rho[res$type == "pooled"])))

write.csv(res, file.path(out_dir, "Table_correlations_all.csv"), row.names = FALSE)

# ============================================================================#
# 3. COMPACT SUPPLEMENTARY TABLE ----
# ============================================================================#

fmt_p <- function(p) ifelse(p < 1e-4, "<0.0001", ifelse(p < 0.01, sprintf("%.2g", p), sprintf("%.2f", p)))
fmt_r <- function(r, lo, hi) ifelse(is.na(lo), sprintf("%.2f", r), sprintf("%.2f [%.2f, %.2f]", r, lo, hi))

tab <- res %>%
  mutate(stratum = ifelse(type == "within condition", as.character(condition), as.character(type)),
         rho_ci  = fmt_r(rho, ci_lo, ci_hi),
         p_txt   = fmt_p(p),
         pBH_txt = paste0(fmt_p(p_BH), ifelse(signif_BH == "ns", "", paste0(" (", signif_BH, ")")))) %>%
  select(Figure = panel, Correlation = pair, Dataset = dataset, Analysis = analysis,
         Stratum = stratum, n, `rho [95% CI]` = rho_ci, p = p_txt, `p (BH)` = pBH_txt)
write.csv(tab, file.path(out_dir, "Table_S4_correlations.csv"), row.names = FALSE)
writeLines(c("# Table S4. Correlations: pooled, adjusted for condition and within condition", "",
             paste("|", paste(names(tab), collapse = " | "), "|"),
             paste0(paste(rep("|---", ncol(tab)), collapse = ""), "|"),
             apply(tab, 1, function(r) paste("|", paste(r, collapse = " | "), "|"))),
           file.path(out_dir, "Table_S4_correlations.md"))

print(as.data.frame(res %>% filter(type != "within condition") %>%
                      select(panel, type, n, rho, p, p_BH, signif_BH)), digits = 3)
cat("\nS85 complete.\n")
