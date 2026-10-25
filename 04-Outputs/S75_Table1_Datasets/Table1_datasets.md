Table 1. Summary of transcriptomic datasets used in this study.



+--------------------+--------------+----------------------------+----------------------------------+------------------------------------------------------+--------------------------------------------------+
| Dataset            | Platform     | Original participants      | Exclusions                       | Analyzed participants                                | Processing                                       |
+====================+==============+============================+==================================+======================================================+==================================================+
| Fan DKD            | Bulk RNA-seq | 36 (Normal 9, Early DKD 6, | None                             | **All conditions:**\                                 | By Fan et al.: read counts log2-transformed and  |
| (GSE142025)        |              | Advanced DKD 21)           |                                  | - Whole tissue: Normal 9, Early DKD 6, Advanced DKD  | quantile-normalized, transcripts with < 100      |
|                    |              |                            |                                  | 21 (n=36)                                            | reads across all samples removed; deposited in   |
|                    |              |                            |                                  |                                                      | GEO in this form (raw counts not available). In  |
|                    |              |                            |                                  |                                                      | this study: values used without further          |
|                    |              |                            |                                  |                                                      | normalization. One kidney biopsy per             |
|                    |              |                            |                                  |                                                      | participant.                                     |
+--------------------+--------------+----------------------------+----------------------------------+------------------------------------------------------+--------------------------------------------------+
| KPMP bulk RNA-seq  | Bulk RNA-seq | 27 (AKI 12, CKD 15)        | library < 1 million reads        | **All conditions:**\                                 | In this study: one count file per participant    |
|                    |              |                            | (n=1);\                          | - Whole tissue: AKI 11, CKD 13 (n=24)                | (for replicates, the file with the largest       |
|                    |              |                            | adjudicated diagnosis 'Cannot be |                                                      | library); participants with < 1 million mapped   |
|                    |              |                            | determined' (n=2)                |                                                      | reads or an undetermined adjudicated diagnosis   |
|                    |              |                            |                                  |                                                      | excluded. Raw counts merged on the genes present |
|                    |              |                            |                                  |                                                      | in all files; low-expression genes removed with  |
|                    |              |                            |                                  |                                                      | filterByExpr, keeping all pathway and marker     |
|                    |              |                            |                                  |                                                      | genes; TMM normalization and transformation to   |
|                    |              |                            |                                  |                                                      | log2 counts per million (edgeR, prior count 1).  |
+--------------------+--------------+----------------------------+----------------------------------+------------------------------------------------------+--------------------------------------------------+
| Combined bulk (Fan | Bulk RNA-seq | 60 (Normal 9, Early DKD 6, | All conditions: none;\           | **All conditions:**\                                 | In this study: Fan DKD (log2                     |
| + KPMP)            |              | Advanced DKD 21, AKI 11,   | DKD-focused: CKD without DKD     | - Whole tissue: Normal 9, Early DKD 6, Advanced DKD  | quantile-normalized) and KPMP (log2 CPM)         |
|                    |              | CKD 13)                    | (n=6)                            | 21, AKI 11, CKD 13 (n=60)\                           | matrices merged on the genes present in both;    |
|                    |              |                            |                                  | **DKD-focused:**\                                    | cohort batch effect, including the different     |
|                    |              |                            |                                  | - Whole tissue: Normal 9, Early DKD 6, Advanced DKD  | normalization scales, removed with ComBat (sva), |
|                    |              |                            |                                  | 21, AKI 11, CKD 7 (n=54)                             | protecting a binary Normal vs. diseased          |
|                    |              |                            |                                  |                                                      | covariate. Module scores computed on the         |
|                    |              |                            |                                  |                                                      | corrected matrix. DKD-focused analysis: KPMP CKD |
|                    |              |                            |                                  |                                                      | restricted to participants with an adjudicated   |
|                    |              |                            |                                  |                                                      | diagnosis of diabetic kidney disease.            |
+--------------------+--------------+----------------------------+----------------------------------+------------------------------------------------------+--------------------------------------------------+
| KPMP scRNA-seq     | scRNA-seq    | 111 (Healthy 37, AKI 20,   | enrollment category not studied  | **All conditions:**\                                 | By KPMP: integrated object including the Lake et |
| (incl. Lake 2023)  |              | CKD 51, DM-R 3)            | (n=3);\                          | - PT: Healthy 33, AKI 20, CKD 49 (n=102)\            | al. (2023) atlas, with raw UMI counts,           |
|                    |              |                            | < 20 PT cells (n=6);\            | - Stroma: Healthy 36, AKI 15, CKD 41 (n=92)\         | batch-corrected expression values and cell       |
|                    |              |                            | < 20 stroma cells (n=16);\       | - Fibroblasts: not analyzed (fewer than 20           | annotations. In this study: raw UMI counts       |
|                    |              |                            | DKD-focused analysis: CKD        | fibroblasts in most participants); available:        | log-normalized with Seurat LogNormalize (counts  |
|                    |              |                            | without DKD (stroma) (n=26);\    | Healthy 5, AKI 3, CKD 6 (n=14)\                      | divided by the cell's total count, multiplied by |
|                    |              |                            | DKD-focused analysis: CKD        | - PT and stroma (matched): Healthy 32, AKI 15, CKD   | 10,000, natural log of 1 + x); the               |
|                    |              |                            | without DKD (PT) (n=30)          | 39 (n=86)\                                           | batch-corrected values provided by KPMP were not |
|                    |              |                            |                                  | **DKD-focused:**\                                    | used. Log-normalized values averaged per         |
|                    |              |                            |                                  | - PT: Healthy 33, AKI 20, CKD 19 (n=72)\             | participant (pseudobulk), keeping participants   |
|                    |              |                            |                                  | - Stroma: Healthy 36, AKI 15, CKD 15 (n=66)\         | with >= 20 cells in the compartment. PT: cell    |
|                    |              |                            |                                  | - Fibroblasts: not analyzed (fewer than 20           | type 'epithelial cell of proximal tubule'.       |
|                    |              |                            |                                  | fibroblasts in most participants); available:        | Stroma: cell class 'stroma cells' (FIB 18%,      |
|                    |              |                            |                                  | Healthy 5, AKI 3, CKD 3 (n=11)\                      | VSMC/P 82%). Matched: participants with both PT  |
|                    |              |                            |                                  | - PT and stroma (matched): Healthy 32, AKI 15, CKD   | and stroma profiles. DKD-focused: CKD restricted |
|                    |              |                            |                                  | 15 (n=62)                                            | to adjudicated diabetic kidney disease.          |
+--------------------+--------------+----------------------------+----------------------------------+------------------------------------------------------+--------------------------------------------------+
| KPMP snRNA-seq     | snRNA-seq    | 157 (Healthy 40, AKI 33,   | enrollment category not studied  | **All conditions:**\                                 | By KPMP: single-nucleus atlas with raw UMI       |
|                    |              | CKD 72, DM-R 11, No        | (n=12);\                         | - PT: Healthy 38, AKI 33, CKD 71 (n=142)\            | counts and cell annotations. In this study: raw  |
|                    |              | clinical record 1)         | < 20 PT nuclei (n=3);\           | - Stroma: not analyzed (fibroblast subclass analyzed | UMI counts log-normalized with Seurat            |
|                    |              |                            | < 20 fibroblast nuclei (n=15);\  | instead; stroma = FIB + VSM/P); available: Healthy   | LogNormalize (counts divided by the nucleus's    |
|                    |              |                            | DKD-focused analysis: CKD        | 36, AKI 31, CKD 67 (n=134)\                          | total count, multiplied by 10,000, natural log   |
|                    |              |                            | without DKD (fibroblast)         | - Fibroblasts: Healthy 36, AKI 31, CKD 63 (n=130)\   | of 1 + x; identical to the normalized values     |
|                    |              |                            | (n=36);\                         | - PT and fibroblasts (matched): Healthy 35, AKI 31,  | provided by KPMP). Log-normalized values         |
|                    |              |                            | DKD-focused analysis: CKD        | CKD 62 (n=128)\                                      | averaged per participant (pseudobulk), keeping   |
|                    |              |                            | without DKD (PT) (n=43)          | **DKD-focused:**\                                    | participants with >= 20 nuclei in the            |
|                    |              |                            |                                  | - PT: Healthy 38, AKI 33, CKD 28 (n=99)\             | compartment. PT: subclass 'PT'. Fibroblasts:     |
|                    |              |                            |                                  | - Stroma: not analyzed (fibroblast subclass analyzed | subclass 'FIB' (all fibroblast subtypes).        |
|                    |              |                            |                                  | instead; stroma = FIB + VSM/P); available: Healthy   | Matched: participants with both PT and           |
|                    |              |                            |                                  | 36, AKI 31, CKD 28 (n=95)\                           | fibroblast profiles. DKD-focused: CKD restricted |
|                    |              |                            |                                  | - Fibroblasts: Healthy 36, AKI 31, CKD 27 (n=94)\    | to adjudicated diabetic kidney disease.          |
|                    |              |                            |                                  | - PT and fibroblasts (matched): Healthy 35, AKI 31,  |                                                  |
|                    |              |                            |                                  | CKD 26 (n=92)                                        |                                                  |
+--------------------+--------------+----------------------------+----------------------------------+------------------------------------------------------+--------------------------------------------------+

