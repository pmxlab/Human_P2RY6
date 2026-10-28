# ============================================================================#
# S75_Table1_Datasets.R
# Table 1 of the manuscript as a data.frame: for every dataset / analysis the
# original number of participants, each exclusion step and the analyzed
# numbers per condition.
#
# Terminology: "participant" = one human donor (Fan: patient / control,
# KPMP: Participant ID). All counts are participants unless stated otherwise.
# KPMP scRNA-seq: SampleID identifies the participant. For healthy reference
# donors the KPMPID column holds specimen/library IDs (e.g. 21_015 and 21-015_B
# are two specimens of donor 15; REF67-* are 8 sequencing runs of donor 67),
# so participants are counted by SampleID, the unit of the S50 pseudobulk.
#
# Inputs : S20/S30 rds, KPMP bulk metadata + count files (library sizes),
#          metadata columns of the KPMP scRNA/snRNA h5Seurat files (read with
#          hdf5r, no Seurat object is loaded), KPMP clinical data,
#          S50/S60 pseudobulk rds (cross-check).
# Outputs: 04-Outputs/S75_Table1_Datasets/
#          Table1_datasets.csv / .md        compact Table 1 (grid table, nested cells)
#          Table_analyzed_participants.*     one row per analysis set x compartment
#          Table_filtering_steps.csv / .md  long format: one row per step
# ============================================================================#

rm(list = ls())

library(dplyr)
library(tidyr)
library(hdf5r)
library(knitr)

# ============================================================================#
# SETTINGS ----
# ============================================================================#

s20_rds  <- "../04-Outputs/S20_FanDKD_BulkRNA/S20_normalised.rds"
s30_rds  <- "../04-Outputs/S30_KPMP_BulkRNA/S30_normalised.rds"
s50_rds  <- "../04-Outputs/S50_KPMP_scRNA/S50_pseudobulk.rds"
s60_rds  <- "../04-Outputs/S60_KPMP_snRNA/S60_pseudobulk.rds"
bulk_dir <- "../../KPMP_BulkRNA"
sc_file  <- "../../KPMP_scRNA/KPMP_PREMIERE_SC_version2_ForExplorer_RemovedBatchEffect_Final2025.h5seurat"
sc_clin  <- "../../KPMP_scRNA/OpenAccessClinicalData/20251205_OpenAccessClinicalData.csv"
sn_file  <- "../../KPMP_snRNA/SingleNucleus_KPMP_Explorer_05182025.h5Seurat"
sn_clin  <- "../../KPMP_snRNA/6d3047c2-588f-4ad4-8eba-f4c46a1e5c14_20251205_OpenAccessClinicalData/20251205_OpenAccessClinicalData.csv"
out_dir  <- "../04-Outputs/S75_Table1_Datasets"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

min_cells <- 20   # S50 / S60: minimum cells per pseudobulk unit and compartment

# ============================================================================#
# HELPERS ----
# ============================================================================#

# Read selected meta.data columns from an h5Seurat file (factors are stored
# as groups with 'levels' and 1-based 'values')
read_h5_meta <- function(f, cols) {
  h  <- H5File$new(f, mode = "r"); on.exit(h$close_all())
  md <- h[["meta.data"]]
  out <- lapply(cols, function(cl) {
    x <- md[[cl]]
    if (inherits(x, "H5Group")) x[["levels"]][][x[["values"]][]] else x[]
  })
  setNames(as.data.frame(out, stringsAsFactors = FALSE), cols)
}

# "Healthy 36, AKI 15, CKD 41" from a named count vector in a fixed order
fmt_counts <- function(x, labels) {
  x <- x[names(labels)]; x[is.na(x)] <- 0
  paste(sprintf("%s %d", labels, x), collapse = ", ")
}

steps <- list()   # long-format filtering table
add_step <- function(dataset, step, n, by_condition, note = "") {
  steps[[length(steps) + 1]] <<- data.frame(
    dataset = dataset, step = step, n_participants = n,
    by_condition = by_condition, note = note, stringsAsFactors = FALSE)
}

# ============================================================================#
# 1. FAN DKD BULK (GSE142025) ----
# ============================================================================#

fan <- readRDS(s20_rds)$meta
fan_lab <- c(Normal = "Normal", Early_DKD = "Early DKD", Advanced_DKD = "Advanced DKD")
fan_n   <- table(as.character(fan$condition))
add_step("Fan DKD bulk", "Deposited in GEO", nrow(fan), fmt_counts(fan_n, fan_lab),
         "one biopsy per participant")
add_step("Fan DKD bulk", "Analyzed", nrow(fan), fmt_counts(fan_n, fan_lab), "no exclusions")

# ============================================================================#
# 2. KPMP BULK (filtering steps of S30 §1-§2) ----
# ============================================================================#

atlas <- read.csv(file.path(bulk_dir, "0_ATLAS_PATIENT_FILE_METADATA.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
names(atlas) <- trimws(names(atlas))
cm <- atlas %>%
  filter(grepl("counts\\.csv$", `File Name`, ignore.case = TRUE),
         !grepl("raw-reads", `File Name`, ignore.case = TRUE)) %>%
  mutate(file_path = file.path(bulk_dir, `File Name`)) %>%
  filter(file.exists(file_path))
cm$total_counts <- vapply(cm$file_path, function(f) {
  d <- read.csv(f, header = TRUE, check.names = FALSE)
  if (ncol(d) < 2) d <- read.table(f, header = TRUE, sep = "\t", check.names = FALSE)
  sum(d[, 2], na.rm = TRUE)
}, numeric(1))

kb_lab <- c(AKI = "AKI", CKD = "CKD")
kb_cnt <- function(d) table(d$`Enrollment Category`[!duplicated(d$`Participant ID`)])
n_files <- nrow(cm)
add_step("KPMP bulk", "Participants with a count file", n_distinct(cm$`Participant ID`),
         fmt_counts(kb_cnt(cm), kb_lab),
         sprintf("%d count files; replicates reduced to the file with the largest library",
                 n_files))
kb <- cm %>% group_by(`Participant ID`) %>%
  slice_max(total_counts, n = 1, with_ties = FALSE) %>% ungroup()
lowlib <- kb %>% filter(total_counts < 1e6)
kb <- kb %>% filter(total_counts >= 1e6)
add_step("KPMP bulk", "Excluded: library < 1 million reads", nrow(lowlib),
         fmt_counts(kb_cnt(lowlib), kb_lab))
cbd <- kb %>% filter(`Primary Adjudicated Category` == "Cannot be determined")
kb  <- kb %>% filter(`Primary Adjudicated Category` != "Cannot be determined")
add_step("KPMP bulk", "Excluded: adjudicated diagnosis 'Cannot be determined'", nrow(cbd),
         fmt_counts(kb_cnt(cbd), kb_lab))
ckd_sub <- table(kb$`Primary Adjudicated Category`[kb$`Enrollment Category` == "CKD"])
ckd_sub_txt <- sprintf("CKD: DKD %d, HKD %d, other %d",
  sum(ckd_sub[names(ckd_sub) == "Diabetic Kidney Disease"]),
  sum(ckd_sub[names(ckd_sub) == "Hypertensive Kidney Disease"]),
  sum(ckd_sub[!names(ckd_sub) %in% c("Diabetic Kidney Disease", "Hypertensive Kidney Disease")]))
add_step("KPMP bulk", "Analyzed", nrow(kb), fmt_counts(kb_cnt(kb), kb_lab), ckd_sub_txt)
s30 <- readRDS(s30_rds)$meta
stopifnot(setequal(s30$sample_id, kb$`Participant ID`))   # matches S30 output

# ============================================================================#
# 3. COMBINED BULK (S40 / S42) ----
# ============================================================================#

n_ckd_nondkd <- sum(kb$`Enrollment Category` == "CKD" &
                    kb$`Primary Adjudicated Category` != "Diabetic Kidney Disease")
comb_n <- c(fan_n, kb_cnt(kb))
comb_lab <- c(fan_lab, AKI = "AKI", CKD = "CKD")
add_step("Combined bulk", "Analyzed, all conditions", sum(comb_n),
         fmt_counts(comb_n, comb_lab), "Fan DKD bulk + KPMP bulk; ComBat batch correction")
comb_dkd <- comb_n; comb_dkd["CKD"] <- comb_dkd["CKD"] - n_ckd_nondkd
add_step("Combined bulk", "Excluded for DKD-focused analysis: CKD without DKD",
         n_ckd_nondkd, fmt_counts(c(CKD = n_ckd_nondkd), c(CKD = "CKD")))
add_step("Combined bulk", "Analyzed, DKD-focused", sum(comb_dkd), fmt_counts(comb_dkd, comb_lab))

# ============================================================================#
# 4. KPMP scRNA-seq (S50 / S52) ----
# ============================================================================#

sc_lab <- c("Healthy Reference Tissue" = "Healthy", "Acute Kidney Injury" = "AKI",
            "Chronic Kidney Disease" = "CKD")
sc <- read_h5_meta(sc_file, c("SampleID", "KPMPID", "EnrollmentCategory", "cell_type", "class",
                               "subclass.level1"))
clin_sc <- read.csv(sc_clin)
sc$adjudicated <- clin_sc$Primary.Adjudicated.Category[match(sc$KPMPID, clin_sc$Participant.ID)]
sc$participant <- sc$SampleID   # see header: SampleID = participant

units_sc <- sc %>% distinct(SampleID, participant, EnrollmentCategory)
part_cnt <- function(u) table(u$EnrollmentCategory[!duplicated(u$participant)])
add_step("KPMP scRNA-seq", "Participants in the KPMP object", n_distinct(units_sc$participant),
         fmt_counts(part_cnt(units_sc), c(sc_lab, "Diabetes Mellitus - Resilient" = "DM-R")),
         sprintf("%d specimen/library IDs (KPMPID column) from %d participants",
                 n_distinct(sc$KPMPID), n_distinct(units_sc$participant)))
dmr <- units_sc %>% filter(!EnrollmentCategory %in% names(sc_lab))
add_step("KPMP scRNA-seq", "Excluded: enrollment category not studied",
         n_distinct(dmr$participant), fmt_counts(part_cnt(dmr), c("Diabetes Mellitus - Resilient" = "DM-R")))
units_sc <- units_sc %>% filter(EnrollmentCategory %in% names(sc_lab))

sc_compartment <- function(mask, label) {
  ncell <- sc[mask & sc$EnrollmentCategory %in% names(sc_lab), ] %>% count(SampleID)
  keep  <- ncell$SampleID[ncell$n >= min_cells]
  kept  <- units_sc %>% filter(SampleID %in% keep)
  drop  <- units_sc %>% filter(!SampleID %in% keep)
  smp   <- table(kept$EnrollmentCategory[!duplicated(kept$SampleID)])
  add_step("KPMP scRNA-seq", sprintf("Excluded: < %d %s cells", min_cells, label),
           n_distinct(drop$participant), fmt_counts(part_cnt(drop), sc_lab))
  add_step("KPMP scRNA-seq", sprintf("Analyzed, %s, all conditions", label),
           n_distinct(kept$participant), fmt_counts(part_cnt(kept), sc_lab))
  kept
}
sc_pt <- sc_compartment(sc$cell_type == "epithelial cell of proximal tubule", "PT")
sc_st <- sc_compartment(sc$class == "stroma cells", "stroma")

# DKD-focused: CKD samples restricted to adjudicated DKD (S52 filter on stroma table)
adj_sc <- sc %>% distinct(SampleID, adjudicated)
s50 <- readRDS(s50_rds)
pb_fib_dkd <- s50$pb_fib %>%
  filter(!(as.character(EnrollmentCategory) == "Chronic Kidney Disease" &
           adjudicated != "Diabetic Kidney Disease"))
keep_dkd_sc <- pb_fib_dkd$SampleID
sc_st_dkd <- sc_st %>% filter(SampleID %in% keep_dkd_sc)
# PT: filtered by each participant's own adjudicated category (SampleID is
# the Participant.ID for AKI/CKD), not by membership in the stroma set
sc_pt_dkd <- sc_pt %>%
  mutate(adjudicated = clin_sc$Primary.Adjudicated.Category[
           match(participant, clin_sc$Participant.ID)]) %>%
  filter(!(EnrollmentCategory == "Chronic Kidney Disease" &
           adjudicated != "Diabetic Kidney Disease"))
add_step("KPMP scRNA-seq", "Excluded for DKD-focused analysis: CKD without DKD (stroma)",
         n_distinct(sc_st$participant) - n_distinct(sc_st_dkd$participant),
         fmt_counts(part_cnt(anti_join(sc_st, sc_st_dkd, by = "participant")), sc_lab))
add_step("KPMP scRNA-seq", "Analyzed, stroma, DKD-focused", n_distinct(sc_st_dkd$participant),
         fmt_counts(part_cnt(sc_st_dkd), sc_lab))
add_step("KPMP scRNA-seq", "Excluded for DKD-focused analysis: CKD without DKD (PT)",
         n_distinct(sc_pt$participant) - n_distinct(sc_pt_dkd$participant),
         fmt_counts(part_cnt(anti_join(sc_pt, sc_pt_dkd, by = "participant")), sc_lab))
add_step("KPMP scRNA-seq", "Analyzed, PT, DKD-focused", n_distinct(sc_pt_dkd$participant),
         fmt_counts(part_cnt(sc_pt_dkd), sc_lab))
# cross-check against S50 pseudobulk tables
stopifnot(setequal(unique(sc_pt$SampleID), s50$pb_pt$SampleID),
          setequal(unique(sc_st$SampleID), s50$pb_fib$SampleID))

# ============================================================================#
# 5. KPMP snRNA-seq (S60 / S62) ----
# ============================================================================#

sn_lab <- c("Healthy Reference" = "Healthy", AKI = "AKI", CKD = "CKD")
sn <- read_h5_meta(sn_file, c("patient", "subclass.l1"))
clin_sn <- read.csv(sn_clin)
sn$enrollment  <- clin_sn$Enrollment.Category[match(sn$patient, clin_sn$Participant.ID)]
sn$enrollment[is.na(sn$enrollment) | sn$enrollment == ""] <- "No clinical record"
sn$adjudicated <- clin_sn$Primary.Adjudicated.Category[match(sn$patient, clin_sn$Participant.ID)]
pt_sn <- sn %>% distinct(patient, enrollment, adjudicated)
sn_cnt <- function(d) table(d$enrollment[!duplicated(d$patient)])
all_sn_lab <- c(sn_lab, setNames(setdiff(unique(pt_sn$enrollment), names(sn_lab)),
                                 setdiff(unique(pt_sn$enrollment), names(sn_lab))))
add_step("KPMP snRNA-seq", "Participants in the KPMP object", nrow(pt_sn),
         fmt_counts(sn_cnt(pt_sn), all_sn_lab), "pseudobulk per participant")
other_sn <- pt_sn %>% filter(!enrollment %in% names(sn_lab))
add_step("KPMP snRNA-seq", "Excluded: enrollment category not studied", nrow(other_sn),
         fmt_counts(sn_cnt(other_sn), all_sn_lab[!names(all_sn_lab) %in% names(sn_lab)]))
pt_sn <- pt_sn %>% filter(enrollment %in% names(sn_lab))

sn_compartment <- function(l1, label) {
  ncell <- sn %>% filter(subclass.l1 == l1, enrollment %in% names(sn_lab)) %>% count(patient)
  keep  <- ncell$patient[ncell$n >= min_cells]
  kept  <- pt_sn %>% filter(patient %in% keep)
  drop  <- pt_sn %>% filter(!patient %in% keep)
  add_step("KPMP snRNA-seq", sprintf("Excluded: < %d %s nuclei", min_cells, label),
           nrow(drop), fmt_counts(sn_cnt(drop), sn_lab))
  add_step("KPMP snRNA-seq", sprintf("Analyzed, %s, all conditions", label),
           nrow(kept), fmt_counts(sn_cnt(kept), sn_lab))
  kept
}
sn_pt <- sn_compartment("PT", "PT")
sn_fb <- sn_compartment("FIB", "fibroblast")
s60 <- readRDS(s60_rds)
stopifnot(setequal(sn_pt$patient, s60$pb_pt$patient),
          setequal(sn_fb$patient, s60$pb_fib$patient))

sn_fb_dkd <- sn_fb %>% filter(!(enrollment == "CKD" & adjudicated != "Diabetic Kidney Disease"))
sn_pt_dkd <- sn_pt %>% filter(!(enrollment == "CKD" & adjudicated != "Diabetic Kidney Disease"))
add_step("KPMP snRNA-seq", "Excluded for DKD-focused analysis: CKD without DKD (fibroblast)",
         nrow(sn_fb) - nrow(sn_fb_dkd), fmt_counts(sn_cnt(anti_join(sn_fb, sn_fb_dkd, by = "patient")), sn_lab))
add_step("KPMP snRNA-seq", "Analyzed, fibroblast, DKD-focused", nrow(sn_fb_dkd),
         fmt_counts(sn_cnt(sn_fb_dkd), sn_lab))
add_step("KPMP snRNA-seq", "Excluded for DKD-focused analysis: CKD without DKD (PT)",
         nrow(sn_pt) - nrow(sn_pt_dkd), fmt_counts(sn_cnt(anti_join(sn_pt, sn_pt_dkd, by = "patient")), sn_lab))
add_step("KPMP snRNA-seq", "Analyzed, PT, DKD-focused", nrow(sn_pt_dkd),
         fmt_counts(sn_cnt(sn_pt_dkd), sn_lab))

# ============================================================================#
# 6. CROSS-COMPARTMENT (PT and stroma/fibroblast available) ----
# ============================================================================#

cross <- function(pt, fb, idcol, lab, condcol) {
  k <- intersect(pt[[idcol]], fb[[idcol]])
  d <- pt[pt[[idcol]] %in% k, ]
  list(n = length(k), txt = fmt_counts(table(d[[condcol]][!duplicated(d[[idcol]])]), lab))
}
x1 <- cross(sc_pt %>% distinct(SampleID, .keep_all = TRUE), sc_st %>% distinct(SampleID, .keep_all = TRUE),
            "SampleID", sc_lab, "EnrollmentCategory")
x2 <- cross(sc_pt_dkd %>% distinct(SampleID, .keep_all = TRUE), sc_st_dkd %>% distinct(SampleID, .keep_all = TRUE),
            "SampleID", sc_lab, "EnrollmentCategory")
add_step("Cross-compartment scRNA-seq", "Analyzed, all conditions", x1$n, x1$txt)
add_step("Cross-compartment scRNA-seq", "Analyzed, DKD-focused", x2$n, x2$txt)
x3 <- cross(sn_pt, sn_fb, "patient", sn_lab, "enrollment")
x4 <- cross(sn_pt_dkd, sn_fb_dkd, "patient", sn_lab, "enrollment")
add_step("Cross-compartment snRNA-seq", "Analyzed, all conditions", x3$n, x3$txt)
add_step("Cross-compartment snRNA-seq", "Analyzed, DKD-focused", x4$n, x4$txt)
stopifnot(x1$n == nrow(s50$pb_cross), x3$n == nrow(s60$pb_cross))

steps_df <- bind_rows(steps)

# ============================================================================#
# 7. ANALYZED PARTICIPANTS PER ANALYSIS SET AND COMPARTMENT ----
# ============================================================================#
# One row per dataset x analysis set x compartment, with counts per condition
# and the figure panels that display these participants (verified against the
# figure data by S81_Check_Figure_Ns.R).

get <- function(ds, st) steps_df %>% filter(dataset == ds, step == st)
excl <- function(ds) {
  e <- steps_df %>% filter(dataset == ds, grepl("^Excluded", step), n_participants > 0)
  if (nrow(e) == 0) return("None")
  paste(sprintf("%s (n=%d)", sub("^Excluded( for)?:? ", "", e$step), e$n_participants), collapse = "; ")
}

cond_cols <- c("Normal/Healthy", "Early DKD", "Advanced DKD", "AKI", "CKD")
count_row <- function(dataset, analysis, compartment, cond, figures, in_figures = "yes",
                      note = "") {
  tab <- table(factor(cond, levels = cond_cols))
  data.frame(Dataset = dataset, Analysis = analysis, Compartment = compartment,
             as.list(tab), Total = sum(tab), `In figures` = in_figures,
             `Figure panels` = figures, Note = note,
             check.names = FALSE, stringsAsFactors = FALSE)
}
relabel <- function(x, map) unname(map[as.character(x)])
bulk_map <- c(Normal = "Normal/Healthy", Early_DKD = "Early DKD",
              Advanced_DKD = "Advanced DKD", AKI = "AKI", CKD = "CKD")
sc_map   <- c("Healthy Reference Tissue" = "Normal/Healthy",
              "Acute Kidney Injury" = "AKI", "Chronic Kidney Disease" = "CKD")
sn_map   <- c("Healthy Reference" = "Normal/Healthy", AKI = "AKI", CKD = "CKD")
uniq <- function(d, id, cond) d[!duplicated(d[[id]]), cond]

# . Bulk ----
fan_cond  <- as.character(fan$condition)
kb_cond   <- kb$`Enrollment Category`
kb_dkd    <- !(kb$`Enrollment Category` == "CKD" &
               kb$`Primary Adjudicated Category` != "Diabetic Kidney Disease")
comb_cond <- c(fan_cond, kb_cond)
comb_dkd_cond <- c(fan_cond, kb_cond[kb_dkd])
# Revision 1: panels (a) show the Fan DKD cohort only; the KPMP bulk data and
# the combined Fan + KPMP set are no longer analyzed (cohort and diagnosis
# confounded; see S43).
fan_fig <- "1a, 2a, 3a, 4a, 5a, 6a, S1a, S2a, S3a"

# . scRNA-seq: fibroblasts are not analyzed separately (too few cells) ----
sc_fib_n <- sc %>%
  filter(class == "stroma cells", subclass.level1 == "FIB",
         EnrollmentCategory %in% names(sc_lab)) %>%
  count(SampleID, EnrollmentCategory) %>% filter(n >= min_cells)
sc_stroma_comp <- sc %>% filter(class == "stroma cells") %>% count(subclass.level1) %>%
  mutate(pct = round(100 * n / sum(n)))
sc_stroma_txt <- paste(sprintf("%s %d%%", sc_stroma_comp$subclass.level1,
                               sc_stroma_comp$pct), collapse = ", ")
sc_fib_dkd <- sc_fib_n %>%
  left_join(sc %>% distinct(SampleID, adjudicated), by = "SampleID") %>%
  filter(!(EnrollmentCategory == "Chronic Kidney Disease" &
           (is.na(adjudicated) | adjudicated != "Diabetic Kidney Disease")))

# . snRNA-seq: a stroma compartment (FIB + VSM/P) is not analyzed ----
sn_stroma <- sn %>% filter(subclass.l1 %in% c("FIB", "VSM/P"),
                           enrollment %in% names(sn_lab)) %>%
  count(patient, enrollment, adjudicated) %>% filter(n >= min_cells)
sn_stroma_dkd <- sn_stroma %>%
  filter(!(enrollment == "CKD" & adjudicated != "Diabetic Kidney Disease"))

sc_cross     <- intersect(sc_pt$participant, sc_st$participant)
sc_cross_dkd <- intersect(sc_pt_dkd$participant, sc_st_dkd$participant)
sn_cross     <- intersect(sn_pt$patient, sn_fb$patient)
sn_cross_dkd <- intersect(sn_pt_dkd$patient, sn_fb_dkd$patient)
sc_u <- function(d) relabel(uniq(d, "participant", "EnrollmentCategory"), sc_map)
sn_u <- function(d) relabel(uniq(d, "patient", "enrollment"), sn_map)

analyzed <- bind_rows(
  count_row("Fan DKD bulk", "All conditions", "Whole tissue", relabel(fan_cond, bulk_map),
            fan_fig),
  count_row("KPMP bulk", "All conditions", "Whole tissue", relabel(kb_cond, bulk_map),
            "none", "no", note = paste(get("KPMP bulk", "Analyzed")$note,
              "| not analyzed: no healthy samples; AKI/CKD cannot be compared with Normal of another cohort")),
  count_row("Combined bulk", "All conditions", "Whole tissue", relabel(comb_cond, bulk_map),
            "none", "no", note = "not analyzed: cohort and diagnosis confounded (S43)"),
  count_row("Combined bulk", "DKD-focused", "Whole tissue", relabel(comb_dkd_cond, bulk_map),
            "none", "no", note = "not analyzed: cohort and diagnosis confounded (S43)"),
  count_row("KPMP scRNA-seq", "All conditions", "PT", sc_u(sc_pt), "1b"),
  count_row("KPMP scRNA-seq", "All conditions", "Stroma", sc_u(sc_st),
            "2b, 3b, 4b, 5b, S1b, S2b, S3b", note = paste("stroma class:", sc_stroma_txt)),
  count_row("KPMP scRNA-seq", "All conditions", "Fibroblasts",
            relabel(sc_fib_n$EnrollmentCategory, sc_map), "none", "no",
            note = sprintf("not analyzed: fewer than %d fibroblasts in most participants", min_cells)),
  count_row("KPMP scRNA-seq", "All conditions", "PT and stroma (matched)",
            sc_u(sc_pt %>% filter(participant %in% sc_cross)), "6b"),
  count_row("KPMP scRNA-seq", "DKD-focused", "PT", sc_u(sc_pt_dkd), "1c"),
  count_row("KPMP scRNA-seq", "DKD-focused", "Stroma", sc_u(sc_st_dkd),
            "2c, 3c, 4c, 5c, S1c, S2c, S3c"),
  count_row("KPMP scRNA-seq", "DKD-focused", "Fibroblasts",
            relabel(sc_fib_dkd$EnrollmentCategory, sc_map), "none", "no",
            note = sprintf("not analyzed: fewer than %d fibroblasts in most participants", min_cells)),
  count_row("KPMP scRNA-seq", "DKD-focused", "PT and stroma (matched)",
            sc_u(sc_pt_dkd %>% filter(participant %in% sc_cross_dkd)), "6c"),
  count_row("KPMP snRNA-seq", "All conditions", "PT", sn_u(sn_pt), "1d"),
  count_row("KPMP snRNA-seq", "All conditions", "Stroma", sn_u(sn_stroma), "none", "no",
            note = "not analyzed: fibroblast subclass analyzed instead; stroma = FIB + VSM/P"),
  count_row("KPMP snRNA-seq", "All conditions", "Fibroblasts", sn_u(sn_fb),
            "2d, 3d, 4d, 5d, S1d, S2d, S3d"),
  count_row("KPMP snRNA-seq", "All conditions", "PT and fibroblasts (matched)",
            sn_u(sn_pt %>% filter(patient %in% sn_cross)), "6d"),
  count_row("KPMP snRNA-seq", "DKD-focused", "PT", sn_u(sn_pt_dkd), "1e"),
  count_row("KPMP snRNA-seq", "DKD-focused", "Stroma", sn_u(sn_stroma_dkd), "none", "no",
            note = "not analyzed: fibroblast subclass analyzed instead; stroma = FIB + VSM/P"),
  count_row("KPMP snRNA-seq", "DKD-focused", "Fibroblasts", sn_u(sn_fb_dkd),
            "2e, 3e, 4e, 5e, S1e, S2e, S3e"),
  count_row("KPMP snRNA-seq", "DKD-focused", "PT and fibroblasts (matched)",
            sn_u(sn_pt_dkd %>% filter(patient %in% sn_cross_dkd)), "6e")
)
# consistency with the step table
stopifnot(
  analyzed$Total[analyzed$Dataset == "Combined bulk" & analyzed$Analysis == "All conditions"] ==
    get("Combined bulk", "Analyzed, all conditions")$n_participants,
  analyzed$Total[analyzed$Dataset == "KPMP scRNA-seq" & analyzed$Compartment == "PT and stroma (matched)" &
                   analyzed$Analysis == "All conditions"] == nrow(s50$pb_cross),
  analyzed$Total[analyzed$Dataset == "KPMP snRNA-seq" & analyzed$Compartment == "PT and fibroblasts (matched)" &
                   analyzed$Analysis == "All conditions"] == nrow(s60$pb_cross))

# ============================================================================#
# 8. COMPACT TABLE 1 ----
# ============================================================================#

fmt_row <- function(r) {
  v <- unlist(r[cond_cols]); v <- v[v > 0]
  lab <- sub("Normal/Healthy", if (grepl("bulk", r$Dataset)) "Normal" else "Healthy", names(v))
  sprintf("%s (n=%d)", paste(sprintf("%s %d", lab, v), collapse = ", "), r$Total)
}
analyzed_cell <- function(ds) {
  out <- character()
  for (an in c("All conditions", "DKD-focused")) {
    a <- analyzed %>% filter(Dataset == ds, Analysis == an)
    if (nrow(a) == 0) next
    out <- c(out, paste0(an, ":"))
    for (k in seq_len(nrow(a))) {
      r <- a[k, ]
      txt <- if (r$`In figures` == "yes" || grepl("bulk", ds)) fmt_row(r) else
        sprintf("not analyzed (%s); available: %s", sub("^not analyzed: ", "", r$Note), fmt_row(r))
      out <- c(out, sprintf("- %s: %s", r$Compartment, txt))
    }
  }
  paste(out, collapse = "\n")
}

proc <- c(
  "Fan DKD bulk" = paste(
    "By Fan et al.: read counts log2-transformed and quantile-normalized, transcripts with < 100",
    "reads across all samples removed; deposited in GEO in this form (raw counts not available).",
    "In this study: values used without further normalization. One kidney biopsy per participant."),
  "KPMP bulk" = paste(
    "In this study: one count file per participant (for replicates, the file with the largest",
    "library); participants with < 1 million mapped reads or an undetermined adjudicated diagnosis",
    "excluded. Raw counts merged on the genes present in all files; low-expression genes removed",
    "with filterByExpr, keeping all pathway and marker genes; TMM normalization and transformation",
    "to log2 counts per million (edgeR, prior count 1)."),
  "Combined bulk" = paste(
    "In this study: Fan DKD (log2 quantile-normalized) and KPMP (log2 CPM) matrices merged on",
    "the genes present in both; cohort batch effect, including the different normalization scales,",
    "removed with ComBat (sva), protecting a binary Normal vs. diseased covariate. Module scores",
    "computed on the corrected matrix. DKD-focused analysis: KPMP CKD restricted to participants",
    "with an adjudicated diagnosis of diabetic kidney disease."),
  "KPMP scRNA-seq" = paste(
    "By KPMP: integrated object including the Lake et al. (2023) atlas, with raw UMI counts,",
    "batch-corrected expression values and cell annotations.",
    "In this study: raw UMI counts log-normalized with Seurat LogNormalize (counts divided by the",
    "cell's total count, multiplied by 10,000, natural log of 1 + x); the batch-corrected values",
    "provided by KPMP were not used. Log-normalized values averaged per participant (pseudobulk),",
    sprintf("keeping participants with >= %d cells in the compartment.", min_cells),
    "PT: cell type 'epithelial cell of proximal tubule'. Stroma: cell class 'stroma cells'",
    sprintf("(%s).", sc_stroma_txt),
    "Matched: participants with both PT and stroma profiles.",
    "DKD-focused: CKD restricted to adjudicated diabetic kidney disease."),
  "KPMP snRNA-seq" = paste(
    "By KPMP: single-nucleus atlas with raw UMI counts and cell annotations.",
    "In this study: raw UMI counts log-normalized with Seurat LogNormalize (counts divided by the",
    "nucleus's total count, multiplied by 10,000, natural log of 1 + x; identical to the normalized",
    "values provided by KPMP). Log-normalized values averaged per participant (pseudobulk),",
    sprintf("keeping participants with >= %d nuclei in the compartment.", min_cells),
    "PT: subclass 'PT'. Fibroblasts: subclass 'FIB' (all fibroblast subtypes).",
    "Matched: participants with both PT and fibroblast profiles.",
    "DKD-focused: CKD restricted to adjudicated diabetic kidney disease.")
)

orig <- function(ds, st) { o <- get(ds, st); sprintf("%d (%s)", o$n_participants, o$by_condition) }
# Revision 1: KPMP bulk and the combined bulk set are not analyzed (S43) and
# are not listed.
table1 <- data.frame(
  Dataset  = c("Fan DKD (GSE142025)", "KPMP scRNA-seq (incl. Lake 2023)", "KPMP snRNA-seq"),
  Platform = c("Bulk RNA-seq", "scRNA-seq", "snRNA-seq"),
  `Original participants` = c(orig("Fan DKD bulk", "Deposited in GEO"),
                              orig("KPMP scRNA-seq", "Participants in the KPMP object"),
                              orig("KPMP snRNA-seq", "Participants in the KPMP object")),
  Exclusions = c(excl("Fan DKD bulk"), excl("KPMP scRNA-seq"), excl("KPMP snRNA-seq")),
  `Analyzed participants` = c(analyzed_cell("Fan DKD bulk"), analyzed_cell("KPMP scRNA-seq"),
                              analyzed_cell("KPMP snRNA-seq")),
  Processing = unname(proc[c("Fan DKD bulk", "KPMP scRNA-seq", "KPMP snRNA-seq")]),
  check.names = FALSE, stringsAsFactors = FALSE)
table1$Exclusions <- gsub("; ", ";\n", table1$Exclusions)

# ============================================================================#
# WRITE ----
# ============================================================================#
# Markdown: grid tables (pander) keep the line breaks inside cells and render
# in R Markdown to HTML, Word and PDF.

library(pander)
grid_md <- function(df, caption, widths) {
  panderOptions("table.split.table", Inf); panderOptions("keep.line.breaks", TRUE)
  c(caption, "", capture.output(pandoc.table(df, style = "grid", split.cells = widths,
                                             justify = "left")))
}

write.csv(table1,   file.path(out_dir, "Table1_datasets.csv"),       row.names = FALSE)
write.csv(analyzed, file.path(out_dir, "Table_analyzed_participants.csv"), row.names = FALSE)
write.csv(steps_df, file.path(out_dir, "Table_filtering_steps.csv"), row.names = FALSE)
# hard line breaks (pandoc: backslash at end of line) so that nested items
# stay on separate lines after rendering
table1_md <- table1
for (cl in c("Exclusions", "Analyzed participants"))
  table1_md[[cl]] <- gsub("\n", "\\\\\n", table1_md[[cl]])
table1_md$`Analyzed participants` <- gsub("(All conditions|DKD-focused):", "**\\1:**",
                                          table1_md$`Analyzed participants`)
writeLines(grid_md(table1_md, "Table 1. Summary of transcriptomic datasets used in this study.",
                   c(18, 12, 26, 32, 52, 48)),
           file.path(out_dir, "Table1_datasets.md"))
writeLines(c("Analyzed participants per analysis set and compartment.", "",
             kable(analyzed, format = "pipe")),
           file.path(out_dir, "Table_analyzed_participants.md"))
writeLines(c("Participant numbers at each filtering step.", "",
             kable(steps_df, format = "pipe")),
           file.path(out_dir, "Table_filtering_steps.md"))

print(kable(analyzed[, c("Dataset", "Analysis", "Compartment", cond_cols, "Total", "Figure panels")],
            format = "pipe"))
cat("\nS75 complete. Outputs saved to:", out_dir, "\n")
