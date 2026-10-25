# =============================================================================
# S70_DonorMembership.R
# Donor-level membership table across the four cohorts used in the revision:
#   Fan DKD bulk (S20), KPMP bulk (S30), KPMP scRNA (S50), KPMP snRNA (S60).
#
# Output: a CSV with one row per unique donor ID, columns indicating presence
# in each cohort and the per-cohort condition / adjudicated diagnosis.
#
# Used for: cohort overlap supplementary table in the manuscript revision.
# =============================================================================

rm(list = ls())

library(dplyr)
library(tidyr)
library(readr)

# =============================================================================
# SETTINGS ----
# =============================================================================

s20_rds <- "../04-Outputs/S20_FanDKD_BulkRNA/S20_normalised.rds"
s30_rds <- "../04-Outputs/S30_KPMP_BulkRNA/S30_normalised.rds"
s50_rds <- "../04-Outputs/S50_KPMP_scRNA/S50_pseudobulk.rds"
s60_rds <- "../04-Outputs/S60_KPMP_snRNA/S60_pseudobulk.rds"

out_dir <- "../04-Outputs/S70_DonorMembership"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# LOAD ----
# =============================================================================

s20 <- readRDS(s20_rds)   # list(expr_mat, meta) — Fan DKD
s30 <- readRDS(s30_rds)   # list(expr_mat, meta) — KPMP bulk
sc  <- readRDS(s50_rds)   # list(pb_pt, pb_fib, pb_fib_sub, pb_cross) — KPMP scRNA
sn  <- readRDS(s60_rds)   # list(pb_pt, pb_fib, pb_fib_sub, pb_cross) — KPMP snRNA

# =============================================================================
# . Per-cohort donor frames ----
# =============================================================================

fan <- s20$meta %>%
  transmute(
    donor_id = as.character(short_id),
    fan_condition = as.character(condition)
  ) %>%
  distinct()

kbulk <- s30$meta %>%
  transmute(
    donor_id = as.character(sample_id),
    kpmp_bulk_enrollment = as.character(enrollment_category),
    kpmp_bulk_adjudicated = as.character(primary_adjudicated_category)
  ) %>%
  distinct()

# scRNA: pull EnrollmentCategory from pb_pt; pull adjudicated from pb_fib
# (since adjudicated was joined into pb_fib, not pb_pt; see S50/S51 history).
sc_pt_ids  <- unique(as.character(sc$pb_pt$SampleID))
sc_fib_ids <- unique(as.character(sc$pb_fib$SampleID))

sc_enrollment <- sc$pb_pt %>%
  transmute(
    donor_id = as.character(SampleID),
    kpmp_scrna_enrollment = as.character(EnrollmentCategory)
  ) %>%
  distinct()

sc_adj_lookup <- if ("adjudicated" %in% colnames(sc$pb_fib)) {
  sc$pb_fib %>%
    transmute(
      donor_id = as.character(SampleID),
      kpmp_scrna_adjudicated = as.character(adjudicated)
    ) %>% distinct()
} else {
  tibble(donor_id = character(), kpmp_scrna_adjudicated = character())
}

scrna <- full_join(sc_enrollment, sc_adj_lookup, by = "donor_id") %>%
  mutate(
    in_kpmp_scrna_pt  = donor_id %in% sc_pt_ids,
    in_kpmp_scrna_fib = donor_id %in% sc_fib_ids
  )

# snRNA: same pattern with snRNA column names.
sn_pt_ids  <- unique(as.character(sn$pb_pt$patient))
sn_fib_ids <- unique(as.character(sn$pb_fib$patient))

sn_enrollment <- sn$pb_pt %>%
  transmute(
    donor_id = as.character(patient),
    kpmp_snrna_enrollment = as.character(enrollment)
  ) %>%
  distinct()

sn_adj_lookup <- if ("adjudicated" %in% colnames(sn$pb_fib)) {
  sn$pb_fib %>%
    transmute(
      donor_id = as.character(patient),
      kpmp_snrna_adjudicated = as.character(adjudicated)
    ) %>% distinct()
} else {
  tibble(donor_id = character(), kpmp_snrna_adjudicated = character())
}

snrna <- full_join(sn_enrollment, sn_adj_lookup, by = "donor_id") %>%
  mutate(
    in_kpmp_snrna_pt  = donor_id %in% sn_pt_ids,
    in_kpmp_snrna_fib = donor_id %in% sn_fib_ids
  )

# =============================================================================
# . Combine ----
# =============================================================================

all_ids <- sort(unique(c(
  fan$donor_id, kbulk$donor_id, scrna$donor_id, snrna$donor_id
)))

fan_donor_set     <- unique(fan$donor_id)
kbulk_donor_set   <- unique(kbulk$donor_id)
scrna_donor_set   <- unique(c(sc_pt_ids, sc_fib_ids))
snrna_donor_set   <- unique(c(sn_pt_ids, sn_fib_ids))

# Empty strings should be treated as missing in enrollment / adjudicated.
blank_to_na <- function(x) {
  x <- as.character(x)
  ifelse(is.na(x) | x == "", NA_character_, x)
}

# Each row tracks a donor's presence in four distinct DATASETS. Cohort and
# technology are jointly identified by the dataset name (Fan-bulk is bulk
# RNA-seq from the Fan DKD study; KPMP-bulk is bulk from KPMP; KPMP-scRNA
# and KPMP-snRNA are the cellular KPMP datasets).
membership <- tibble(donor_id = all_ids) %>%
  left_join(fan,   by = "donor_id") %>%
  left_join(kbulk, by = "donor_id") %>%
  left_join(scrna, by = "donor_id") %>%
  left_join(snrna, by = "donor_id") %>%
  mutate(
    in_fan_bulk        = donor_id %in% fan_donor_set,
    in_kpmp_bulk       = donor_id %in% kbulk_donor_set,
    in_kpmp_scrna      = donor_id %in% scrna_donor_set,
    in_kpmp_snrna      = donor_id %in% snrna_donor_set,
    in_kpmp_scrna_pt   = donor_id %in% sc_pt_ids,
    in_kpmp_scrna_fib  = donor_id %in% sc_fib_ids,
    in_kpmp_snrna_pt   = donor_id %in% sn_pt_ids,
    in_kpmp_snrna_fib  = donor_id %in% sn_fib_ids,
    n_datasets         = in_fan_bulk + in_kpmp_bulk + in_kpmp_scrna + in_kpmp_snrna,
    fan_condition           = blank_to_na(fan_condition),
    kpmp_bulk_enrollment    = blank_to_na(kpmp_bulk_enrollment),
    kpmp_bulk_adjudicated   = blank_to_na(kpmp_bulk_adjudicated),
    kpmp_scrna_enrollment   = blank_to_na(kpmp_scrna_enrollment),
    kpmp_scrna_adjudicated  = blank_to_na(kpmp_scrna_adjudicated),
    kpmp_snrna_enrollment   = blank_to_na(kpmp_snrna_enrollment),
    kpmp_snrna_adjudicated  = blank_to_na(kpmp_snrna_adjudicated),
    # Heuristic Lake-format flag: short numeric ID without hyphen (KPMP IDs
    # use the X-Y format). Used as a sanity check; refine against the Lake
    # GSE183276 donor list if a definitive flag is required.
    is_lake_format = grepl("^[0-9]+$", donor_id)
  )

# . Unified condition (best available source label, retained for traceability) ----
# Priority: Fan condition > KPMP scRNA adjudicated > snRNA adjudicated >
#           bulk adjudicated > any enrollment.
membership <- membership %>%
  mutate(
    condition_unified = case_when(
      !is.na(fan_condition)              ~ fan_condition,
      !is.na(kpmp_scrna_adjudicated)     ~ kpmp_scrna_adjudicated,
      !is.na(kpmp_snrna_adjudicated)     ~ kpmp_snrna_adjudicated,
      !is.na(kpmp_bulk_adjudicated)      ~ kpmp_bulk_adjudicated,
      !is.na(kpmp_scrna_enrollment)      ~ kpmp_scrna_enrollment,
      !is.na(kpmp_snrna_enrollment)      ~ kpmp_snrna_enrollment,
      !is.na(kpmp_bulk_enrollment)       ~ kpmp_bulk_enrollment,
      TRUE                                ~ NA_character_
    )
  )

# . Harmonised condition (common vocabulary across cohorts) ----
# "Normal" (Fan), "Healthy Reference Tissue" (KPMP scRNA), and "Healthy
# Reference" (KPMP snRNA) all denote the same healthy biological group; they
# are merged here. DKD is identified by Fan's Early/Advanced_DKD labels or by
# KPMP adjudicated == "Diabetic Kidney Disease". Non-DKD CKD is grouped as
# "CKD-other". KPMP enrollment "AKI" subsumes adjudicated subtypes such as
# Acute Tubular Injury and Acute Interstitial Nephritis.
healthy_lab <- c("Normal", "Healthy Reference", "Healthy Reference Tissue")
dkd_lab     <- c("Early_DKD", "Advanced_DKD", "Diabetic Kidney Disease")
aki_lab     <- c("AKI", "Acute Kidney Injury", "Acute Tubular Injury",
                 "Acute Interstitial Nephritis")
ckd_lab     <- c("CKD", "Chronic Kidney Disease", "Hypertensive Kidney Disease")

classify_one <- function(lab) {
  if (is.na(lab))                  return(NA_character_)
  if (lab %in% healthy_lab)        return("Healthy")
  if (lab %in% dkd_lab)            return("DKD")
  if (lab %in% aki_lab)            return("AKI")
  if (lab %in% ckd_lab)            return("CKD-other")
  return("Other")
}

# Combine cohort labels per donor — DKD takes precedence over generic CKD,
# so a donor flagged DKD anywhere is DKD; a donor flagged AKI anywhere
# (without a DKD flag) is AKI; etc.
harmonise <- function(...) {
  labs <- unique(stats::na.omit(vapply(c(...), classify_one,
                                       FUN.VALUE = character(1))))
  if (length(labs) == 0)            return(NA_character_)
  if ("DKD"       %in% labs)        return("DKD")
  if ("AKI"       %in% labs)        return("AKI")
  if ("CKD-other" %in% labs)        return("CKD-other")
  if ("Healthy"   %in% labs)        return("Healthy")
  return("Other")
}

membership$condition_harmonized <- mapply(
  harmonise,
  membership$fan_condition,
  membership$kpmp_bulk_adjudicated,
  membership$kpmp_bulk_enrollment,
  membership$kpmp_scrna_adjudicated,
  membership$kpmp_scrna_enrollment,
  membership$kpmp_snrna_adjudicated,
  membership$kpmp_snrna_enrollment,
  USE.NAMES = FALSE
)

# . Pairwise / triple / quadruple dataset intersections ----
# Each column is "is this donor present in BOTH (or all of) these datasets?"
# Names are fully qualified by dataset (Fan-bulk, KPMP-bulk, KPMP-scRNA,
# KPMP-snRNA) so cohort and technology are not conflated.
membership <- membership %>%
  mutate(
    # 2-way
    in_fan_bulk_x_kpmp_bulk    = in_fan_bulk   & in_kpmp_bulk,
    in_fan_bulk_x_kpmp_scrna   = in_fan_bulk   & in_kpmp_scrna,
    in_fan_bulk_x_kpmp_snrna   = in_fan_bulk   & in_kpmp_snrna,
    in_kpmp_bulk_x_kpmp_scrna  = in_kpmp_bulk  & in_kpmp_scrna,
    in_kpmp_bulk_x_kpmp_snrna  = in_kpmp_bulk  & in_kpmp_snrna,
    in_kpmp_scrna_x_kpmp_snrna = in_kpmp_scrna & in_kpmp_snrna,
    # 3-way
    in_fan_bulk_x_kpmp_bulk_x_kpmp_scrna  = in_fan_bulk  & in_kpmp_bulk  & in_kpmp_scrna,
    in_fan_bulk_x_kpmp_bulk_x_kpmp_snrna  = in_fan_bulk  & in_kpmp_bulk  & in_kpmp_snrna,
    in_fan_bulk_x_kpmp_scrna_x_kpmp_snrna = in_fan_bulk  & in_kpmp_scrna & in_kpmp_snrna,
    in_kpmp_bulk_x_kpmp_scrna_x_kpmp_snrna = in_kpmp_bulk & in_kpmp_scrna & in_kpmp_snrna,
    # 4-way
    in_all_four_datasets = in_fan_bulk & in_kpmp_bulk & in_kpmp_scrna & in_kpmp_snrna
  )

# Compact pattern label, e.g. "Fan-bulk" / "KPMP-bulk + KPMP-scRNA" /
# "KPMP-scRNA + KPMP-snRNA" / "KPMP-bulk + KPMP-scRNA + KPMP-snRNA".
build_pattern <- function(a, b, c, d) {
  parts <- c(
    if (isTRUE(a)) "Fan-bulk"   else NULL,
    if (isTRUE(b)) "KPMP-bulk"  else NULL,
    if (isTRUE(c)) "KPMP-scRNA" else NULL,
    if (isTRUE(d)) "KPMP-snRNA" else NULL
  )
  if (length(parts) == 0) return(NA_character_)
  paste(parts, collapse = " + ")
}
membership$membership_pattern <- mapply(
  build_pattern,
  membership$in_fan_bulk,
  membership$in_kpmp_bulk,
  membership$in_kpmp_scrna,
  membership$in_kpmp_snrna,
  USE.NAMES = FALSE
)

# Reorder columns for readability.
membership <- membership %>%
  select(
    donor_id,
    condition_harmonized,
    condition_unified,
    membership_pattern,
    n_datasets,
    in_fan_bulk, in_kpmp_bulk, in_kpmp_scrna, in_kpmp_snrna,
    in_fan_bulk_x_kpmp_bulk,
    in_fan_bulk_x_kpmp_scrna,
    in_fan_bulk_x_kpmp_snrna,
    in_kpmp_bulk_x_kpmp_scrna,
    in_kpmp_bulk_x_kpmp_snrna,
    in_kpmp_scrna_x_kpmp_snrna,
    in_fan_bulk_x_kpmp_bulk_x_kpmp_scrna,
    in_fan_bulk_x_kpmp_bulk_x_kpmp_snrna,
    in_fan_bulk_x_kpmp_scrna_x_kpmp_snrna,
    in_kpmp_bulk_x_kpmp_scrna_x_kpmp_snrna,
    in_all_four_datasets,
    in_kpmp_scrna_pt, in_kpmp_scrna_fib,
    in_kpmp_snrna_pt, in_kpmp_snrna_fib,
    is_lake_format,
    fan_condition,
    kpmp_bulk_enrollment, kpmp_bulk_adjudicated,
    kpmp_scrna_enrollment, kpmp_scrna_adjudicated,
    kpmp_snrna_enrollment, kpmp_snrna_adjudicated
  ) %>%
  arrange(desc(n_datasets), donor_id)

# =============================================================================
# WRITE ----
# =============================================================================

out_csv <- file.path(out_dir, "Table_donor_membership.csv")
write_csv(membership, out_csv)

cat("Wrote:", out_csv, "\n")
cat("Rows :", nrow(membership), "\n\n")

# =============================================================================
# SUMMARY ----
# =============================================================================

cat("=== Dataset sizes (n donors) ===\n")
cat(sprintf("  Fan-bulk    : %d\n", sum(membership$in_fan_bulk)))
cat(sprintf("  KPMP-bulk   : %d\n", sum(membership$in_kpmp_bulk)))
cat(sprintf("  KPMP-scRNA  : %d\n", sum(membership$in_kpmp_scrna)))
cat(sprintf("  KPMP-snRNA  : %d\n", sum(membership$in_kpmp_snrna)))

cat("\n=== Donors per dataset count ===\n")
print(table(membership$n_datasets, dnn = "n_datasets"))

cat("\n=== Pairwise dataset overlaps ===\n")
ov <- function(a, b) sum(membership[[a]] & membership[[b]])
cat(sprintf("  Fan-bulk   ∩ KPMP-bulk  : %d\n", ov("in_fan_bulk",   "in_kpmp_bulk")))
cat(sprintf("  Fan-bulk   ∩ KPMP-scRNA : %d\n", ov("in_fan_bulk",   "in_kpmp_scrna")))
cat(sprintf("  Fan-bulk   ∩ KPMP-snRNA : %d\n", ov("in_fan_bulk",   "in_kpmp_snrna")))
cat(sprintf("  KPMP-bulk  ∩ KPMP-scRNA : %d\n", ov("in_kpmp_bulk",  "in_kpmp_scrna")))
cat(sprintf("  KPMP-bulk  ∩ KPMP-snRNA : %d\n", ov("in_kpmp_bulk",  "in_kpmp_snrna")))
cat(sprintf("  KPMP-scRNA ∩ KPMP-snRNA : %d\n", ov("in_kpmp_scrna", "in_kpmp_snrna")))

cat("\n=== Triple overlap (KPMP-bulk ∩ KPMP-scRNA ∩ KPMP-snRNA) ===\n")
cat(sprintf("  n = %d\n",
            sum(membership$in_kpmp_bulk &
                membership$in_kpmp_scrna &
                membership$in_kpmp_snrna)))

cat("\n=== Lake-format donors (no-hyphen heuristic) ===\n")
cat(sprintf("  Lake-format AND in KPMP-scRNA : %d\n",
            sum(membership$is_lake_format & membership$in_kpmp_scrna)))
cat(sprintf("  Lake-format AND in KPMP-snRNA : %d\n",
            sum(membership$is_lake_format & membership$in_kpmp_snrna)))

cat("\n=== condition_unified distribution ===\n")
print(table(membership$condition_unified, useNA = "ifany"))

cat("\n=== condition_harmonized distribution ===\n")
print(table(membership$condition_harmonized, useNA = "ifany"))

cat("\n=== membership_pattern distribution ===\n")
print(sort(table(membership$membership_pattern), decreasing = TRUE))

cat("\n=== condition_harmonized × membership_pattern ===\n")
print(membership %>%
        count(condition_harmonized, membership_pattern) %>%
        pivot_wider(names_from = membership_pattern,
                    values_from = n, values_fill = 0) %>%
        arrange(desc(rowSums(across(where(is.numeric))))))

cat("\n=== Intersection counts (n donors) ===\n")
ic <- function(col) sum(membership[[col]])
cat(sprintf("  Fan-bulk   ∩ KPMP-bulk                : %d\n", ic("in_fan_bulk_x_kpmp_bulk")))
cat(sprintf("  Fan-bulk   ∩ KPMP-scRNA               : %d\n", ic("in_fan_bulk_x_kpmp_scrna")))
cat(sprintf("  Fan-bulk   ∩ KPMP-snRNA               : %d\n", ic("in_fan_bulk_x_kpmp_snrna")))
cat(sprintf("  KPMP-bulk  ∩ KPMP-scRNA               : %d\n", ic("in_kpmp_bulk_x_kpmp_scrna")))
cat(sprintf("  KPMP-bulk  ∩ KPMP-snRNA               : %d\n", ic("in_kpmp_bulk_x_kpmp_snrna")))
cat(sprintf("  KPMP-scRNA ∩ KPMP-snRNA               : %d\n", ic("in_kpmp_scrna_x_kpmp_snrna")))
cat(sprintf("  Fan-bulk  ∩ KPMP-bulk  ∩ KPMP-scRNA   : %d\n", ic("in_fan_bulk_x_kpmp_bulk_x_kpmp_scrna")))
cat(sprintf("  Fan-bulk  ∩ KPMP-bulk  ∩ KPMP-snRNA   : %d\n", ic("in_fan_bulk_x_kpmp_bulk_x_kpmp_snrna")))
cat(sprintf("  Fan-bulk  ∩ KPMP-scRNA ∩ KPMP-snRNA   : %d\n", ic("in_fan_bulk_x_kpmp_scrna_x_kpmp_snrna")))
cat(sprintf("  KPMP-bulk ∩ KPMP-scRNA ∩ KPMP-snRNA   : %d\n", ic("in_kpmp_bulk_x_kpmp_scrna_x_kpmp_snrna")))
cat(sprintf("  All four datasets                     : %d\n", ic("in_all_four_datasets")))
