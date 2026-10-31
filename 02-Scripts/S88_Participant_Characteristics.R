# ============================================================================#
# S88_Participant_Characteristics.R
# Characteristics of the KPMP participants of the single-cell analyses by
# condition, as far as available in the KPMP open-access clinical data.
# Reviewer comment R2.17 (Life, life-4655118)
#
# Participants: all participants with a PT or stromal/fibroblast pseudobulk
# profile (S50 scRNA-seq, S60 snRNA-seq); condition as in the analyses;
# CKD (DKD): CKD participants with an adjudicated diagnosis of DKD.
# Characteristics (KPMP open-access clinical data, 2025-12-05): sex, age (only
# available in 10-year bins), baseline eGFR, diabetes and hypertension history.
# Revision 1 (AIReview1 B4, comments R2.03 / R2.17): tissue source (KPMP
# 'Sample Type', grouped into five categories) added.
# Each characteristic is summarized over the participants with a value and
# the number with a value is given; no tests (descriptive table).
#
# Outputs: 04-Outputs/S88_Participant_Characteristics/
#   Table_participant_characteristics.csv / .md
# ============================================================================#

rm(list = ls())

library(dplyr)
library(tidyr)

# ============================================================================#
# SETTINGS ----
# ============================================================================#

s50_rds <- "../04-Outputs/S50_KPMP_scRNA/S50_pseudobulk.rds"
s60_rds <- "../04-Outputs/S60_KPMP_snRNA/S60_pseudobulk.rds"
sc_clin <- "../../KPMP_scRNA/OpenAccessClinicalData/20251205_OpenAccessClinicalData.csv"
sn_clin <- "../../KPMP_snRNA/6d3047c2-588f-4ad4-8eba-f4c46a1e5c14_20251205_OpenAccessClinicalData/20251205_OpenAccessClinicalData.csv"
out_dir <- "../04-Outputs/S88_Participant_Characteristics"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

grp_lab <- c("Healthy Reference Tissue" = "Healthy", "Acute Kidney Injury" = "AKI",
             "Chronic Kidney Disease" = "CKD", "Healthy Reference" = "Healthy", "AKI" = "AKI",
             "CKD" = "CKD")
cols <- c("Healthy", "AKI", "CKD", "CKD (DKD)")
# KPMP 'Sample Type' -> tissue source categories
tissue_map <- c("Percutaneous Needle Biopsy"      = "needle biopsy (percutaneous or intra-operative)",
                "Intra-operative Needle Biopsy"   = "needle biopsy (percutaneous or intra-operative)",
                "Transplant Pre-perfusion Biopsy" = "transplant donor biopsy (pre-perfusion or surveillance)",
                "Transplant Surveillance Biopsy"  = "transplant donor biopsy (pre-perfusion or surveillance)",
                "Intra-operative Biopsy"          = "intra-operative biopsy",
                "Tumor Nephrectomy"               = "tumor nephrectomy",
                "Total Tumor Nephrectomy"         = "tumor nephrectomy",
                "Partial Tumor Nephrectomy"       = "tumor nephrectomy",
                "Deceased Donor Nephrectomy"      = "deceased-donor nephrectomy")
tissue_levels <- unique(unname(tissue_map))
missing_vals <- c("", "Unknown", "unknown", "Don't know", NA)

# ============================================================================#
# HELPER FUNCTIONS ----
# ============================================================================#

n_pct <- function(x, level) {
  x <- x[!is.na(x)]
  if (!length(x)) return("–")
  sprintf("%d (%.0f%%)", sum(x == level), 100 * mean(x == level))
}

med_iqr <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return("–")
  q <- quantile(x, c(0.5, 0.25, 0.75), names = FALSE)
  sprintf("%.0f (%.0f–%.0f)", q[1], q[2], q[3])
}

# age bins "20-29 Years" -> lower bound
age_low <- function(x) suppressWarnings(as.numeric(sub("-.*", "", x)))

characteristics <- function(d) {
  d <- d %>% mutate(across(c(sex, age, diab, htn, tissue), ~ ifelse(.x %in% missing_vals, NA, .x)))
  a <- age_low(d$age)
  tibble(
    measure = c("Participants, n", "With KPMP clinical data, n", "Female, n (%)",
                "Age < 40 years, n (%)", "Age 40–59 years, n (%)", "Age ≥ 60 years, n (%)",
                "Baseline eGFR, mL/min/1.73 m², median (IQR)", "Diabetes, n (%)",
                "Hypertension, n (%)", paste0("Tissue: ", tissue_levels, ", n (%)")),
    value = c(nrow(d), sum(d$in_clin), n_pct(d$sex, "Female"),
              n_pct(cut(a, c(-Inf, 40, 60, Inf), right = FALSE, labels = 1:3), 1),
              n_pct(cut(a, c(-Inf, 40, 60, Inf), right = FALSE, labels = 1:3), 2),
              n_pct(cut(a, c(-Inf, 40, 60, Inf), right = FALSE, labels = 1:3), 3),
              med_iqr(d$egfr), n_pct(d$diab, "Yes"), n_pct(d$htn, "Yes"),
              vapply(tissue_levels, function(l) n_pct(d$tissue, l), character(1))),
    n_with = c(NA, NA, sum(!is.na(d$sex)), rep(sum(!is.na(a)), 3), sum(!is.na(d$egfr)),
               sum(!is.na(d$diab)), sum(!is.na(d$htn)),
               rep(sum(!is.na(d$tissue)), length(tissue_levels))))
}

# ============================================================================#
# 1. PARTICIPANTS AND CLINICAL DATA ----
# ============================================================================#

s50 <- readRDS(s50_rds)
s60 <- readRDS(s60_rds)
parts <- bind_rows(
  bind_rows(s50$pb_pt, s50$pb_fib) %>%
    transmute(id = as.character(SampleID), group = as.character(EnrollmentCategory),
              dataset = "scRNA-seq"),
  bind_rows(s60$pb_pt, s60$pb_fib) %>%
    transmute(id = patient, group = as.character(enrollment), dataset = "snRNA-seq")) %>%
  distinct() %>%
  mutate(group = unname(grp_lab[group]))
stopifnot(!anyNA(parts$group), !any(duplicated(parts[, c("id", "dataset")])))

clin <- bind_rows(
  read.csv(sc_clin, check.names = FALSE) %>% mutate(dataset = "scRNA-seq"),
  read.csv(sn_clin, check.names = FALSE) %>% mutate(dataset = "snRNA-seq")) %>%
  transmute(dataset, id = `Participant ID`, adjudicated = `Primary Adjudicated Category`,
            sex = Sex, age = `Age (Years) (Binned)`, egfr = `Baseline eGFR (ml/min/1.73m2)`,
            diab = `Diabetes History`, htn = `Hypertension History`,
            tissue = unname(tissue_map[`Sample Type`]), found = TRUE)

dat <- parts %>%
  left_join(clin, by = c("dataset", "id")) %>%
  mutate(in_clin = !is.na(found))

# ============================================================================#
# 2. TABLE ----
# ============================================================================#

sets <- bind_rows(
  dat %>% mutate(col = group),
  dat %>% filter(group == "CKD", adjudicated %in% "Diabetic Kidney Disease") %>% mutate(col = "CKD (DKD)"))
tab <- sets %>%
  group_by(dataset, col) %>%
  group_modify(~ characteristics(.x)) %>%
  ungroup() %>%
  mutate(value = ifelse(is.na(n_with) | n_with == 0, value, sprintf("%s [%d]", value, n_with))) %>%
  select(-n_with) %>%
  pivot_wider(names_from = col, values_from = value) %>%
  select(dataset, measure, all_of(cols)) %>%
  mutate(measure = factor(measure, levels = unique(measure))) %>%
  arrange(dataset, measure)
print(as.data.frame(tab))
write.csv(tab, file.path(out_dir, "Table_participant_characteristics.csv"), row.names = FALSE)

md <- c("# Participant characteristics (KPMP open-access clinical data)", "",
        "[n] = participants with a value. Fan DKD cohort: no participant-level characteristics available.", "",
        paste0("| ", paste(names(tab), collapse = " | "), " |"),
        paste0(paste(rep("|---", ncol(tab)), collapse = ""), "|"),
        apply(tab, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
writeLines(md, file.path(out_dir, "Table_participant_characteristics.md"))
