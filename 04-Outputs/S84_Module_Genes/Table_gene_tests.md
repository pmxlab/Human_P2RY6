# Module genes vs reference group, per dataset

Two-sided Wilcoxon rank-sum test; BH-adjusted per gene and dataset over the two comparisons.
Cliff's delta > 0: higher than reference (Normal / Healthy).
Median difference in the expression units of the dataset (Fan: log2; KPMP: ln(1 + normalized count), pseudobulk mean).

| Module | Dataset | Gene | Comparison | Median diff | Cliff's delta | p | p (BH) | |
|---|---|---|---|---|---|---|---|---|
| pyrimidine | Fan DKD bulk RNA-seq | CAD | Early DKD vs Normal | +0.409 | +1.00 | 0.0004 | 0.0004 | *** |
| pyrimidine | Fan DKD bulk RNA-seq | CAD | Advanced DKD vs Normal | +0.681 | +0.95 | 2.7e-06 | 5.3e-06 | **** |
| pyrimidine | Fan DKD bulk RNA-seq | DHODH | Early DKD vs Normal | +0.272 | +0.89 | 0.0028 | 0.0056 | ** |
| pyrimidine | Fan DKD bulk RNA-seq | DHODH | Advanced DKD vs Normal | -0.189 | -0.52 | 0.025 | 0.025 | * |
| pyrimidine | Fan DKD bulk RNA-seq | UMPS | Early DKD vs Normal | +0.366 | +0.93 | 0.0016 | 0.0032 | ** |
| pyrimidine | Fan DKD bulk RNA-seq | UMPS | Advanced DKD vs Normal | +0.117 | +0.49 | 0.039 | 0.039 | * |
| pyrimidine | Fan DKD bulk RNA-seq | UCK1 | Early DKD vs Normal | +0.203 | +0.85 | 0.0048 | 0.0096 | ** |
| pyrimidine | Fan DKD bulk RNA-seq | UCK1 | Advanced DKD vs Normal | -0.123 | -0.31 | 0.19 | 0.19 | ns |
| pyrimidine | Fan DKD bulk RNA-seq | CDA | Early DKD vs Normal | -0.565 | -0.33 | 0.33 | 0.33 | ns |
| pyrimidine | Fan DKD bulk RNA-seq | CDA | Advanced DKD vs Normal | -0.805 | -0.61 | 0.008 | 0.016 | * |
| pyrimidine | Fan DKD bulk RNA-seq | DCTPP1 | Early DKD vs Normal | +0.385 | +1.00 | 0.0004 | 0.0004 | *** |
| pyrimidine | Fan DKD bulk RNA-seq | DCTPP1 | Advanced DKD vs Normal | +0.602 | +0.89 | 1.9e-05 | 3.9e-05 | **** |
| pyrimidine | KPMP scRNA-seq | CAD | AKI vs Healthy | +0.0196 | +0.77 | 5.7e-07 | 5.7e-07 | **** |
| pyrimidine | KPMP scRNA-seq | CAD | CKD vs Healthy | +0.0202 | +0.73 | 2.9e-09 | 5.9e-09 | **** |
| pyrimidine | KPMP scRNA-seq | DHODH | AKI vs Healthy | -0.00077 | +0.03 | 0.85 | 0.85 | ns |
| pyrimidine | KPMP scRNA-seq | DHODH | CKD vs Healthy | -0.00341 | -0.10 | 0.43 | 0.85 | ns |
| pyrimidine | KPMP scRNA-seq | UMPS | AKI vs Healthy | +0.0112 | +0.40 | 0.014 | 0.028 | * |
| pyrimidine | KPMP scRNA-seq | UMPS | CKD vs Healthy | +0.0067 | +0.25 | 0.058 | 0.058 | ns |
| pyrimidine | KPMP scRNA-seq | UCK1 | AKI vs Healthy | +0.00428 | +0.06 | 0.74 | 0.74 | ns |
| pyrimidine | KPMP scRNA-seq | UCK1 | CKD vs Healthy | -0.00251 | -0.07 | 0.59 | 0.74 | ns |
| pyrimidine | KPMP scRNA-seq | CDA | AKI vs Healthy | -0.000405 | -0.04 | 0.8 | 0.8 | ns |
| pyrimidine | KPMP scRNA-seq | CDA | CKD vs Healthy | +0.00105 | +0.06 | 0.62 | 0.8 | ns |
| pyrimidine | KPMP scRNA-seq | DCTPP1 | AKI vs Healthy | +0.0985 | +0.89 | 1.2e-09 | 2.3e-09 | **** |
| pyrimidine | KPMP scRNA-seq | DCTPP1 | CKD vs Healthy | +0.0501 | +0.60 | 1.9e-06 | 1.9e-06 | **** |
| pyrimidine | KPMP snRNA-seq | CAD | AKI vs Healthy | +0.0197 | +0.38 | 0.0057 | 0.011 | * |
| pyrimidine | KPMP snRNA-seq | CAD | CKD vs Healthy | +0.00279 | +0.14 | 0.24 | 0.24 | ns |
| pyrimidine | KPMP snRNA-seq | DHODH | AKI vs Healthy | +0.00489 | +0.10 | 0.5 | 0.5 | ns |
| pyrimidine | KPMP snRNA-seq | DHODH | CKD vs Healthy | +0.00305 | +0.09 | 0.42 | 0.5 | ns |
| pyrimidine | KPMP snRNA-seq | UMPS | AKI vs Healthy | +0.0194 | +0.39 | 0.0049 | 0.0098 | ** |
| pyrimidine | KPMP snRNA-seq | UMPS | CKD vs Healthy | +0.0104 | +0.17 | 0.15 | 0.15 | ns |
| pyrimidine | KPMP snRNA-seq | UCK1 | AKI vs Healthy | -0.0028 | -0.14 | 0.33 | 0.66 | ns |
| pyrimidine | KPMP snRNA-seq | UCK1 | CKD vs Healthy | +0.000736 | +0.01 | 0.93 | 0.93 | ns |
| pyrimidine | KPMP snRNA-seq | CDA | AKI vs Healthy | -0.00233 | -0.28 | 0.041 | 0.062 | ns |
| pyrimidine | KPMP snRNA-seq | CDA | CKD vs Healthy | -0.00215 | -0.22 | 0.062 | 0.062 | ns |
| pyrimidine | KPMP snRNA-seq | DCTPP1 | AKI vs Healthy | +0.0176 | +0.64 | 1.1e-06 | 2.3e-06 | **** |
| pyrimidine | KPMP snRNA-seq | DCTPP1 | CKD vs Healthy | +0.00605 | +0.27 | 0.022 | 0.022 | * |
| calcium | Fan DKD bulk RNA-seq | GNAQ | Early DKD vs Normal | -0.256 | -1.00 | 0.0004 | 0.0008 | *** |
| calcium | Fan DKD bulk RNA-seq | GNAQ | Advanced DKD vs Normal | +0.0305 | +0.11 | 0.67 | 0.67 | ns |
| calcium | Fan DKD bulk RNA-seq | ITPR2 | Early DKD vs Normal | -0.16 | -0.63 | 0.05 | 0.05 | * |
| calcium | Fan DKD bulk RNA-seq | ITPR2 | Advanced DKD vs Normal | -0.175 | -0.88 | 0.00019 | 0.00038 | *** |
| calcium | Fan DKD bulk RNA-seq | ITPR3 | Early DKD vs Normal | +0.197 | +0.37 | 0.27 | 0.27 | ns |
| calcium | Fan DKD bulk RNA-seq | ITPR3 | Advanced DKD vs Normal | +1.29 | +0.85 | 6.7e-05 | 0.00013 | *** |
| calcium | Fan DKD bulk RNA-seq | ORAI1 | Early DKD vs Normal | -0.0133 | -0.07 | 0.86 | 0.86 | ns |
| calcium | Fan DKD bulk RNA-seq | ORAI1 | Advanced DKD vs Normal | +0.332 | +0.62 | 0.0069 | 0.014 | * |
| calcium | KPMP scRNA-seq | GNAQ | AKI vs Healthy | +0.152 | +0.57 | 0.00098 | 0.00098 | *** |
| calcium | KPMP scRNA-seq | GNAQ | CKD vs Healthy | +0.129 | +0.54 | 2.8e-05 | 5.5e-05 | **** |
| calcium | KPMP scRNA-seq | ITPR2 | AKI vs Healthy | +0.214 | +0.70 | 3.6e-05 | 3.6e-05 | **** |
| calcium | KPMP scRNA-seq | ITPR2 | CKD vs Healthy | +0.151 | +0.59 | 3.2e-06 | 6.5e-06 | **** |
| calcium | KPMP scRNA-seq | ITPR3 | AKI vs Healthy | +0.007 | +0.11 | 0.52 | 0.52 | ns |
| calcium | KPMP scRNA-seq | ITPR3 | CKD vs Healthy | +0.00836 | +0.09 | 0.51 | 0.52 | ns |
| calcium | KPMP scRNA-seq | ORAI1 | AKI vs Healthy | +0.0147 | +0.32 | 0.075 | 0.075 | ns |
| calcium | KPMP scRNA-seq | ORAI1 | CKD vs Healthy | +0.0244 | +0.33 | 0.014 | 0.028 | * |
| calcium | KPMP snRNA-seq | GNAQ | AKI vs Healthy | +0.244 | +0.49 | 0.00039 | 0.00079 | *** |
| calcium | KPMP snRNA-seq | GNAQ | CKD vs Healthy | +0.179 | +0.33 | 0.0064 | 0.0064 | ** |
| calcium | KPMP snRNA-seq | ITPR2 | AKI vs Healthy | +0.0736 | +0.11 | 0.43 | 0.43 | ns |
| calcium | KPMP snRNA-seq | ITPR2 | CKD vs Healthy | +0.0834 | +0.11 | 0.38 | 0.43 | ns |
| calcium | KPMP snRNA-seq | ITPR3 | AKI vs Healthy | +0.0407 | +0.50 | 0.00048 | 0.00097 | *** |
| calcium | KPMP snRNA-seq | ITPR3 | CKD vs Healthy | +0.00624 | +0.21 | 0.081 | 0.081 | ns |
| calcium | KPMP snRNA-seq | ORAI1 | AKI vs Healthy | +0.00261 | +0.12 | 0.41 | 0.41 | ns |
| calcium | KPMP snRNA-seq | ORAI1 | CKD vs Healthy | +0.00295 | +0.15 | 0.2 | 0.4 | ns |
| ecm | Fan DKD bulk RNA-seq | COL1A1 | Early DKD vs Normal | +0.207 | +0.30 | 0.39 | 0.39 | ns |
| ecm | Fan DKD bulk RNA-seq | COL1A1 | Advanced DKD vs Normal | +2.7 | +0.92 | 9.4e-06 | 1.9e-05 | **** |
| ecm | Fan DKD bulk RNA-seq | COL3A1 | Early DKD vs Normal | +0.483 | +0.44 | 0.18 | 0.18 | ns |
| ecm | Fan DKD bulk RNA-seq | COL3A1 | Advanced DKD vs Normal | +2.59 | +0.96 | 4.6e-05 | 9.3e-05 | **** |
| ecm | Fan DKD bulk RNA-seq | FN1 | Early DKD vs Normal | -0.0486 | +0.07 | 0.86 | 0.86 | ns |
| ecm | Fan DKD bulk RNA-seq | FN1 | Advanced DKD vs Normal | +1.61 | +0.89 | 0.00014 | 0.00029 | *** |
| ecm | Fan DKD bulk RNA-seq | LUM | Early DKD vs Normal | +0.0129 | +0.00 | 1 | 1 | ns |
| ecm | Fan DKD bulk RNA-seq | LUM | Advanced DKD vs Normal | +1.81 | +1.00 | 1.4e-07 | 2.8e-07 | **** |
| ecm | Fan DKD bulk RNA-seq | DCN | Early DKD vs Normal | +0.278 | +0.57 | 0.077 | 0.077 | ns |
| ecm | Fan DKD bulk RNA-seq | DCN | Advanced DKD vs Normal | +1.19 | +0.93 | 7.5e-05 | 0.00015 | *** |
| ecm | KPMP scRNA-seq | COL1A1 | AKI vs Healthy | +0.466 | +0.67 | 0.00019 | 0.00032 | *** |
| ecm | KPMP scRNA-seq | COL1A1 | CKD vs Healthy | +0.205 | +0.48 | 0.00032 | 0.00032 | *** |
| ecm | KPMP scRNA-seq | COL3A1 | AKI vs Healthy | +0.741 | +0.81 | 8.2e-07 | 1.6e-06 | **** |
| ecm | KPMP scRNA-seq | COL3A1 | CKD vs Healthy | +0.301 | +0.58 | 7.1e-06 | 7.1e-06 | **** |
| ecm | KPMP scRNA-seq | FN1 | AKI vs Healthy | +0.49 | +0.68 | 7e-05 | 0.00014 | *** |
| ecm | KPMP scRNA-seq | FN1 | CKD vs Healthy | +0.174 | +0.38 | 0.0042 | 0.0042 | ** |
| ecm | KPMP scRNA-seq | LUM | AKI vs Healthy | +0.227 | +0.56 | 0.0009 | 0.0009 | *** |
| ecm | KPMP scRNA-seq | LUM | CKD vs Healthy | +0.202 | +0.48 | 0.0002 | 0.00041 | *** |
| ecm | KPMP scRNA-seq | DCN | AKI vs Healthy | +0.484 | +0.59 | 0.00096 | 0.00096 | *** |
| ecm | KPMP scRNA-seq | DCN | CKD vs Healthy | +0.231 | +0.48 | 0.00032 | 0.00065 | *** |
| ecm | KPMP snRNA-seq | COL1A1 | AKI vs Healthy | +0.7 | +0.82 | 3.4e-10 | 6.8e-10 | **** |
| ecm | KPMP snRNA-seq | COL1A1 | CKD vs Healthy | +0.213 | +0.45 | 0.0002 | 0.0002 | *** |
| ecm | KPMP snRNA-seq | COL3A1 | AKI vs Healthy | +1.09 | +0.82 | 1.9e-10 | 3.8e-10 | **** |
| ecm | KPMP snRNA-seq | COL3A1 | CKD vs Healthy | +0.473 | +0.60 | 9.3e-07 | 9.3e-07 | **** |
| ecm | KPMP snRNA-seq | FN1 | AKI vs Healthy | +0.244 | +0.48 | 0.00057 | 0.0011 | ** |
| ecm | KPMP snRNA-seq | FN1 | CKD vs Healthy | -0.00934 | -0.12 | 0.32 | 0.32 | ns |
| ecm | KPMP snRNA-seq | LUM | AKI vs Healthy | +0.437 | +0.49 | 0.00046 | 0.00092 | *** |
| ecm | KPMP snRNA-seq | LUM | CKD vs Healthy | +0.198 | +0.23 | 0.059 | 0.059 | ns |
| ecm | KPMP snRNA-seq | DCN | AKI vs Healthy | +0.461 | +0.43 | 0.0024 | 0.0024 | ** |
| ecm | KPMP snRNA-seq | DCN | CKD vs Healthy | +0.623 | +0.39 | 0.0013 | 0.0024 | ** |
