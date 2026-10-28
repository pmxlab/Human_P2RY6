# Leave-one-gene-out sensitivity of the module scores

Each score recomputed without one gene at a time; same test as the score panels
(two-sided Wilcoxon, BH over the two comparisons of a panel). Cliff's delta > 0:
higher than reference. 'Lost without': significant full score becomes ns when
this gene is omitted; 'Gained without': ns full score becomes significant.

| Module | Dataset | Analysis | Comparison | Full score: delta (stars) | Leave-one-out: delta range | Leave-one-out: BH p range | Lost without | Gained without |
|---|---|---|---|---|---|---|---|---|
| pyrimidine | Fan DKD bulk RNA-seq | All conditions | Advanced DKD vs Normal | +0.38 (ns) | -0.24 to +0.99 | 5.6e-07 to 0.53 | - | CDA |
| pyrimidine | Fan DKD bulk RNA-seq | All conditions | Early DKD vs Normal | +1.00 (***) | +0.81 to +1.00 | 0.0004 to 0.015 | - | - |
| pyrimidine | KPMP scRNA-seq | All conditions | AKI vs Healthy | +0.81 (****) | +0.49 to +0.85 | 2e-08 to 0.0054 | - | - |
| pyrimidine | KPMP scRNA-seq | All conditions | CKD vs Healthy | +0.52 (****) | +0.35 to +0.57 | 5.7e-06 to 0.0072 | - | - |
| pyrimidine | KPMP scRNA-seq | DKD-focused | AKI vs Healthy | +0.81 (****) | +0.49 to +0.85 | 2e-08 to 0.0027 | - | - |
| pyrimidine | KPMP scRNA-seq | DKD-focused | CKD vs Healthy | +0.77 (****) | +0.61 to +0.83 | 4.9e-08 to 0.00035 | - | - |
| pyrimidine | KPMP snRNA-seq | All conditions | AKI vs Healthy | +0.48 (***) | +0.38 to +0.56 | 5.6e-05 to 0.012 | - | - |
| pyrimidine | KPMP snRNA-seq | All conditions | CKD vs Healthy | +0.13 (ns) | +0.10 to +0.19 | 0.11 to 0.4 | - | - |
| pyrimidine | KPMP snRNA-seq | DKD-focused | AKI vs Healthy | +0.48 (***) | +0.38 to +0.56 | 5.6e-05 to 0.012 | - | - |
| pyrimidine | KPMP snRNA-seq | DKD-focused | CKD vs Healthy | +0.25 (ns) | +0.18 to +0.33 | 0.022 to 0.22 | - | DHODH, UCK1 |
| calcium | Fan DKD bulk RNA-seq | All conditions | Advanced DKD vs Normal | +0.82 (***) | +0.23 to +0.83 | 0.00023 to 0.35 | ITPR3 | - |
| calcium | Fan DKD bulk RNA-seq | All conditions | Early DKD vs Normal | -0.30 (ns) | -0.81 to +0.07 | 0.015 to 0.86 | - | ITPR3 |
| calcium | KPMP scRNA-seq | All conditions | AKI vs Healthy | +0.80 (****) | +0.59 to +0.81 | 8.2e-07 to 0.00064 | - | - |
| calcium | KPMP scRNA-seq | All conditions | CKD vs Healthy | +0.69 (****) | +0.55 to +0.70 | 3.6e-08 to 3.7e-05 | - | - |
| calcium | KPMP scRNA-seq | DKD-focused | AKI vs Healthy | +0.80 (****) | +0.59 to +0.81 | 1.6e-06 to 0.00064 | - | - |
| calcium | KPMP scRNA-seq | DKD-focused | CKD vs Healthy | +0.76 (****) | +0.61 to +0.78 | 2.3e-06 to 0.00047 | - | - |
| calcium | KPMP snRNA-seq | All conditions | AKI vs Healthy | +0.33 (*) | +0.16 to +0.55 | 0.00015 to 0.26 | GNAQ, ITPR3 | - |
| calcium | KPMP snRNA-seq | All conditions | CKD vs Healthy | +0.23 (ns) | +0.14 to +0.37 | 0.0025 to 0.26 | - | ITPR2 |
| calcium | KPMP snRNA-seq | DKD-focused | AKI vs Healthy | +0.33 (*) | +0.16 to +0.55 | 0.00015 to 0.52 | GNAQ, ITPR3 | - |
| calcium | KPMP snRNA-seq | DKD-focused | CKD vs Healthy | +0.16 (ns) | +0.06 to +0.39 | 0.0082 to 0.71 | - | ITPR2 |
| ecm | Fan DKD bulk RNA-seq | All conditions | Advanced DKD vs Normal | +0.94 (****) | +0.92 to +0.96 | 3.4e-06 to 1.9e-05 | - | - |
| ecm | Fan DKD bulk RNA-seq | All conditions | Early DKD vs Normal | +0.22 (ns) | +0.15 to +0.30 | 0.39 to 0.69 | - | - |
| ecm | KPMP scRNA-seq | All conditions | AKI vs Healthy | +0.68 (****) | +0.66 to +0.72 | 1.8e-05 to 0.00011 | - | - |
| ecm | KPMP scRNA-seq | All conditions | CKD vs Healthy | +0.56 (****) | +0.55 to +0.58 | 1.4e-05 to 3.5e-05 | - | - |
| ecm | KPMP scRNA-seq | DKD-focused | AKI vs Healthy | +0.68 (***) | +0.66 to +0.72 | 3.6e-05 to 0.00021 | - | - |
| ecm | KPMP scRNA-seq | DKD-focused | CKD vs Healthy | +0.59 (***) | +0.57 to +0.61 | 0.00045 to 0.0011 | - | - |
| ecm | KPMP snRNA-seq | All conditions | AKI vs Healthy | +0.71 (****) | +0.66 to +0.76 | 1.8e-08 to 2.7e-06 | - | - |
| ecm | KPMP snRNA-seq | All conditions | CKD vs Healthy | +0.43 (***) | +0.34 to +0.47 | 0.00012 to 0.0046 | - | - |
| ecm | KPMP snRNA-seq | DKD-focused | AKI vs Healthy | +0.71 (****) | +0.66 to +0.76 | 1.8e-08 to 2.7e-06 | - | - |
| ecm | KPMP snRNA-seq | DKD-focused | CKD vs Healthy | +0.46 (**) | +0.39 to +0.50 | 0.0006 to 0.0075 | - | - |
