# ============================================================================#
# S81_Check_Figure_Ns.R
# Verifies that the participant numbers in Table 1 (S75) match the data
# plotted in every panel of the manuscript figures (S80), and that every
# significance star drawn in the figures equals the Benjamini-Hochberg-
# adjusted p-value recomputed here from the plotted data.
#
# The panels of S80 are rebuilt without saving any file; for every panel the
# number of plotted participants per condition is taken from the panel data
# and from the rendered point layer, and compared with the row of
# Table_analyzed_participants.csv that lists the panel.
#
# Inputs : S80_Manuscript_Figures.R (sourced, figures not saved)
#          04-Outputs/S75_Table1_Datasets/Table_analyzed_participants.csv
# Output : 04-Outputs/S75_Table1_Datasets/Check_figure_Ns.csv
#          04-Outputs/S80_Manuscript_Figures/Check_significance_stars.csv
#          Stops with an error if any panel does not match.
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
tab_file <- "../04-Outputs/S75_Table1_Datasets/Table_analyzed_participants.csv"
out_file <- "../04-Outputs/S75_Table1_Datasets/Check_figure_Ns.csv"

cond_cols <- c("Normal/Healthy", "Early DKD", "Advanced DKD", "AKI", "CKD")
cond_map <- c(Normal = "Normal/Healthy", Early_DKD = "Early DKD", Advanced_DKD = "Advanced DKD",
              "Healthy Reference Tissue" = "Normal/Healthy", "Healthy Reference" = "Normal/Healthy",
              "Acute Kidney Injury" = "AKI", AKI = "AKI",
              "Chronic Kidney Disease" = "CKD", CKD = "CKD")

# ============================================================================#
# 1. BUILD THE S80 PANELS WITHOUT SAVING ----
# ============================================================================#

code <- readLines(s80_file)
code <- code[!grepl("^rm\\(list = ls\\(\\)\\)", code)]
# replace save_figure() by a no-op just before the first figure section
first_fig <- grep("^# FIGURE 1 ", code)[1]
stopifnot(!is.na(first_fig))
code <- append(code, "save_figure <- function(...) invisible(NULL)", after = first_fig - 2)
env <- new.env()
eval(parse(text = code), envir = env)

panels <- mget(ls(env, pattern = "^fs?[0-9]+[a-f]$"), envir = env)
panel_id <- function(nm) {
  m <- regmatches(nm, regexec("^f(s?)([0-9]+)([a-f])$", nm))[[1]]
  paste0(if (m[2] == "s") "S" else "", m[3], m[4])
}

# ============================================================================#
# 2. COUNT PLOTTED PARTICIPANTS PER PANEL ----
# ============================================================================#

count_panel <- function(p) {
  d <- p$data
  x <- eval_tidy(p$mapping$x, d); y <- eval_tidy(p$mapping$y, d)
  cond_col <- intersect(c("condition", "EnrollmentCategory", "enrollment"), names(d))[1]
  keep <- !is.na(x) & !is.na(y)
  tab <- table(factor(unname(cond_map[as.character(d[[cond_col]][keep])]), levels = cond_cols))
  # rendered points (after scale limits): the point layer of the built plot
  b <- ggplot_build(p)
  pts <- which(vapply(p$layers, function(l) inherits(l$geom, "GeomPoint"), logical(1)))[1]
  ld <- b$data[[pts]]
  c(as.list(tab), Total = sum(tab), Rendered = sum(!is.na(ld$x) & !is.na(ld$y)))
}

fig <- bind_rows(lapply(names(panels), function(nm) {
  data.frame(Panel = panel_id(nm), count_panel(panels[[nm]]), check.names = FALSE)
}))

# ============================================================================#
# 3. COMPARE WITH TABLE 1 ----
# ============================================================================#

tab <- read.csv(tab_file, check.names = FALSE) %>%
  filter(`In figures` == "yes") %>%
  separate_rows(`Figure panels`, sep = ",\\s*") %>%
  rename(Panel = `Figure panels`)

cmp <- fig %>%
  inner_join(tab, by = "Panel", suffix = c(".figure", ".table")) %>%
  mutate(match = .data[["Normal/Healthy.figure"]] == .data[["Normal/Healthy.table"]] &
                 .data[["Early DKD.figure"]]     == .data[["Early DKD.table"]] &
                 .data[["Advanced DKD.figure"]]  == .data[["Advanced DKD.table"]] &
                 AKI.figure == AKI.table & CKD.figure == CKD.table &
                 Total.figure == Total.table & Rendered == Total.table) %>%
  arrange(Panel)

missing_in_table  <- setdiff(fig$Panel, tab$Panel)
missing_in_figure <- setdiff(tab$Panel, fig$Panel)

out <- cmp %>% select(Panel, Dataset, Analysis, Compartment,
                      ends_with(".figure"), Rendered, Total.table, match)
write.csv(out, out_file, row.names = FALSE)
print(as.data.frame(out %>% select(Panel, Dataset, Analysis, Compartment,
                                   Total.figure, Rendered, Total.table, match)))

cat("\nPanels checked:", nrow(cmp), "| mismatches:", sum(!cmp$match), "\n")
if (length(missing_in_table))  cat("Panels not listed in Table 1:", missing_in_table, "\n")
if (length(missing_in_figure)) cat("Table 1 panels not found in S80:", missing_in_figure, "\n")
stopifnot(all(cmp$match), !length(missing_in_table), !length(missing_in_figure))
cat("All figure panels match Table 1.\n")

# ============================================================================#
# 4. SIGNIFICANCE STARS: DRAWN vs. RECOMPUTED (BH within panel) ----
# ============================================================================#

stars <- function(p) as.character(cut(p, c(-Inf, 1e-4, 1e-3, 1e-2, 0.05, Inf),
                                      c("****", "***", "**", "*", "ns")))
sig <- bind_rows(lapply(names(panels), function(nm) {
  p  <- panels[[nm]]
  li <- which(vapply(p$layers, function(l) inherits(l$geom, "GeomSignif"), logical(1)))
  if (!length(li)) return(NULL)
  sd   <- suppressWarnings(ggplot_build(p))$data[[li]]
  ann  <- sd %>% distinct(group, annotation) %>%
    mutate(idx = as.integer(sub(".*-", "", group)))
  cmps <- p$layers[[li]]$stat_params$comparisons
  d <- p$data
  x <- as.character(eval_tidy(p$mapping$x, d)); y <- eval_tidy(p$mapping$y, d)
  pv <- sapply(cmps, function(cc) suppressWarnings(wilcox.test(y[x == cc[1]], y[x == cc[2]])$p.value))
  data.frame(Panel = panel_id(nm), group1 = sapply(cmps, `[`, 1), group2 = sapply(cmps, `[`, 2),
             p = signif(pv, 4), p_BH = signif(p.adjust(pv, "BH"), 4),
             recomputed = stars(p.adjust(pv, "BH")), unadjusted = stars(pv),
             drawn = ann$annotation[match(seq_along(cmps), ann$idx)])
})) %>% mutate(match = drawn == recomputed) %>% arrange(Panel)
write.csv(sig, "../04-Outputs/S80_Manuscript_Figures/Check_significance_stars.csv", row.names = FALSE)
cat("\nSignificance labels checked:", nrow(sig), "| mismatches:", sum(!sig$match), "\n")
changed <- sig %>% filter(unadjusted != recomputed)
cat("Labels that differ from the unadjusted p-values:\n")
print(as.data.frame(changed %>% select(Panel, group2, p, p_BH, unadjusted, drawn)), row.names = FALSE)
stopifnot(all(sig$match))
cat("S81 complete: all figure panels match Table 1 and all stars are BH-adjusted.\n")
