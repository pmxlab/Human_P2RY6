# ============================================================================#
# S82_P2RY6_Detection.R
# Supplementary Table S1: detection of P2RY6 in the stromal (KPMP scRNA-seq)
# and fibroblast (KPMP snRNA-seq) pseudobulk profiles shown in Figure 2b-e
# (reviewer comment R1.04).
#
# Per dataset and group: participants, cells/nuclei per participant, share
# of P2RY6-positive cells/nuclei per participant (UMI count > 0), participants
# without any P2RY6-positive cell/nucleus, and pooled P2RY6-positive cells/
# nuclei over all cells/nuclei of the group. Participant sets are those of
# Figure 2b-e (S80): all conditions from S50/S60; the DKD-focused panels share
# the Healthy and AKI groups and add CKD restricted to adjudicated Diabetic
# Kidney Disease as a separate row.
#
# Inputs : 04-Outputs/S50_KPMP_scRNA/S50_pseudobulk.rds
#          04-Outputs/S60_KPMP_snRNA/S60_pseudobulk.rds
# Output : 04-Outputs/S82_P2RY6_Detection/Table_S1_P2RY6_detection.csv
#          04-Outputs/S82_P2RY6_Detection/Table_S1_P2RY6_detection.md
# ============================================================================#

rm(list = ls())

library(dplyr)

# ============================================================================#
# SETTINGS ----
# ============================================================================#

s50_rds <- "../04-Outputs/S50_KPMP_scRNA/S50_pseudobulk.rds"
s60_rds <- "../04-Outputs/S60_KPMP_snRNA/S60_pseudobulk.rds"
out_dir <- "../04-Outputs/S82_P2RY6_Detection"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cond_map <- c("Healthy Reference Tissue" = "Healthy", "Healthy Reference" = "Healthy",
              "Acute Kidney Injury" = "AKI", "AKI" = "AKI",
              "Chronic Kidney Disease" = "CKD", "CKD" = "CKD")

# ============================================================================#
# 1. PARTICIPANT SETS OF FIGURE 2b-e ----
# ============================================================================#

s50 <- readRDS(s50_rds)
s60 <- readRDS(s60_rds)

sc <- s50$pb_fib %>%
  transmute(id = SampleID, condition = cond_map[as.character(EnrollmentCategory)],
            adjudicated, n_cells, pct_P2RY6_pos)
sn <- s60$pb_fib %>%
  transmute(id = patient, condition = cond_map[as.character(enrollment)],
            adjudicated, n_cells, pct_P2RY6_pos)
stopifnot(!anyNA(sc$condition), !anyNA(sn$condition))

dkd_only <- function(df) filter(df, !(condition == "CKD" & adjudicated != "Diabetic Kidney Disease"))

# Healthy and AKI are identical in the all-conditions and DKD-focused panels;
# only the CKD group differs, so CKD (DKD only) is added as a separate row.
sets <- list(
  list(data = sc,                                  dataset = "scRNA-seq", compartment = "Stroma",
       figure = c(Healthy = "2b, 2c", AKI = "2b, 2c", CKD = "2b")),
  list(data = dkd_only(sc) %>% filter(condition == "CKD") %>% mutate(condition = "CKD (DKD only)"),
       dataset = "scRNA-seq", compartment = "Stroma", figure = c("CKD (DKD only)" = "2c")),
  list(data = sn,                                  dataset = "snRNA-seq", compartment = "Fibroblasts",
       figure = c(Healthy = "2d, 2e", AKI = "2d, 2e", CKD = "2d")),
  list(data = dkd_only(sn) %>% filter(condition == "CKD") %>% mutate(condition = "CKD (DKD only)"),
       dataset = "snRNA-seq", compartment = "Fibroblasts", figure = c("CKD (DKD only)" = "2e"))
)

# ============================================================================#
# 2. DETECTION SUMMARY ----
# ============================================================================#

med_iqr <- function(v, digits) {
  q <- quantile(v, c(0.5, 0.25, 0.75), names = FALSE)
  sprintf(paste0("%.", digits, "f (%.", digits, "f–%.", digits, "f)"), q[1], q[2], q[3])
}

tab <- bind_rows(lapply(sets, function(s) {
  s$data %>%
    mutate(n_pos = round(pct_P2RY6_pos / 100 * n_cells),
           condition = factor(condition, levels = c("Healthy", "AKI", "CKD", "CKD (DKD only)"))) %>%
    group_by(condition) %>%
    summarise(participants     = n(),
              cells_med_iqr    = med_iqr(n_cells, 0),
              pct_pos_med_iqr  = med_iqr(pct_P2RY6_pos, 2),
              zero_pos         = sprintf("%d (%.0f%%)", sum(n_pos == 0), 100 * mean(n_pos == 0)),
              pooled_pos       = sprintf("%d / %d (%.2f%%)", sum(n_pos), sum(n_cells),
                                         100 * sum(n_pos) / sum(n_cells)),
              .groups = "drop") %>%
    mutate(condition = as.character(condition)) %>%
    mutate(dataset = s$dataset, compartment = s$compartment, figure = s$figure[condition], .before = 1)
}))

# consistency with the participant numbers in Table 1 / Figure 2
stopifnot(identical(tab$participants, c(36L, 15L, 41L, 15L, 36L, 31L, 63L, 27L)))

print(as.data.frame(tab))

# dependence of mean P2RY6 on cells/nuclei per participant (all conditions)
for (nm in c("2b", "2d")) {
  d <- if (nm == "2b") s50$pb_fib else s60$pb_fib
  ct <- cor.test(d$mean_P2RY6, d$n_cells, method = "spearman", exact = FALSE)
  cat(sprintf("Panel %s: Spearman rho(mean P2RY6, cells per participant) = %.2f, p = %.3g\n",
              nm, ct$estimate, ct$p.value))
}

# ============================================================================#
# 3. OUTPUT ----
# ============================================================================#

write.csv(tab, file.path(out_dir, "Table_S1_P2RY6_detection.csv"), row.names = FALSE)

hdr <- c("Dataset", "Compartment", "Figure", "Group", "Participants",
         "Cells/nuclei per participant, median (IQR)",
         "P2RY6+ cells/nuclei per participant, %, median (IQR)",
         "Participants without P2RY6+ cells/nuclei, n (%)",
         "P2RY6+ cells/nuclei, pooled n / total (%)")
md <- c(paste("|", paste(hdr, collapse = " | "), "|"),
        paste("|", paste(rep("---", length(hdr)), collapse = " | "), "|"),
        apply(tab, 1, function(r) paste("|", paste(r, collapse = " | "), "|")))
writeLines(md, file.path(out_dir, "Table_S1_P2RY6_detection.md"))
cat("Written:", file.path(out_dir, "Table_S1_P2RY6_detection.csv"), "\n")
