# Combined bulk: group comparisons by method

Cells: Cliff's delta (group vs Normal; AKI-CKD row: CKD vs AKI) and BH-adjusted significance (per feature, across the 5 comparisons).

## All

| Feature | Comparison | ComBat, binary covariate (S40) | Within-sample rank, no correction | ComBat, DKD bridge (Adv. DKD = KPMP DKD) | ComBat, DKD bridge (all DKD = KPMP DKD) |
|---|---|---|---|---|---
| P2RY6 | Early_DKD vs Normal | +0.04 ns | +0.22 ns | +0.11 ns | +0.07 ns |
| P2RY6 | Advanced_DKD vs Normal | +0.74 ** | +0.83 ** | +0.71 ** | +0.76 ** |
| P2RY6 | AKI vs Normal | +0.80 ** | +0.66 * | +0.94 ** | +0.94 ** |
| P2RY6 | CKD vs Normal | +0.50 ns | +0.08 ns | +0.73 ** | +0.61 * |
| P2RY6 | CKD vs AKI | -0.52 ns | -0.52 ns | -0.64 * | -0.61 * |
| PDGFRB | Early_DKD vs Normal | -0.04 ns | -0.07 ns | -0.07 ns | -0.04 ns |
| PDGFRB | Advanced_DKD vs Normal | +0.77 ** | +0.76 ** | +0.77 ** | +0.77 ** |
| PDGFRB | AKI vs Normal | +0.56 ns | +0.72 * | +0.33 ns | +0.09 ns |
| PDGFRB | CKD vs Normal | +0.83 ** | +0.95 ** | +0.73 * | +0.57 ns |
| PDGFRB | CKD vs AKI | +0.48 ns | +0.55 * | +0.55 * | +0.51 ns |
| Calcium score | Early_DKD vs Normal | -0.52 ns | -0.07 ns | -0.26 ns | -0.44 ns |
| Calcium score | Advanced_DKD vs Normal | +0.79 ** | +0.83 ** | +0.78 ** | +0.81 ** |
| Calcium score | AKI vs Normal | +0.82 ** | +0.92 ** | +1.00 *** | +0.94 ** |
| Calcium score | CKD vs Normal | +0.78 ** | +0.81 ** | +0.95 *** | +0.86 ** |
| Calcium score | CKD vs AKI | -0.37 ns | -0.36 ns | -0.43 ns | -0.41 ns |
| ECM score | Early_DKD vs Normal | +0.11 ns | +0.15 ns | +0.22 ns | +0.22 ns |
| ECM score | Advanced_DKD vs Normal | +0.93 *** | +0.94 *** | +0.92 *** | +0.94 *** |
| ECM score | AKI vs Normal | +0.98 *** | +0.90 ** | +1.00 *** | +1.00 *** |
| ECM score | CKD vs Normal | +1.00 *** | +0.98 *** | +1.00 *** | +1.00 *** |
| ECM score | CKD vs AKI | -0.48 ns | -0.05 ns | -0.62 * | -0.55 * |
| Pyrimidine score | Early_DKD vs Normal | +0.89 * | +1.00 ** | +0.89 * | +0.93 * |
| Pyrimidine score | Advanced_DKD vs Normal | +0.16 ns | +0.64 ** | +0.26 ns | +0.22 ns |
| Pyrimidine score | AKI vs Normal | +0.64 * | -0.76 ** | +0.52 ns | +0.58 ns |
| Pyrimidine score | CKD vs Normal | +0.13 ns | -0.95 ** | -0.26 ns | -0.04 ns |
| Pyrimidine score | CKD vs AKI | -0.41 ns | -0.40 ns | -0.45 ns | -0.47 ns |

## DKD_focused

| Feature | Comparison | ComBat, binary covariate (S40) | Within-sample rank, no correction | ComBat, DKD bridge (Adv. DKD = KPMP DKD) | ComBat, DKD bridge (all DKD = KPMP DKD) |
|---|---|---|---|---|---
| P2RY6 | Early_DKD vs Normal | +0.04 ns | +0.22 ns | +0.11 ns | +0.07 ns |
| P2RY6 | Advanced_DKD vs Normal | +0.74 ** | +0.83 ** | +0.71 ** | +0.76 ** |
| P2RY6 | AKI vs Normal | +0.80 ** | +0.66 * | +0.94 ** | +0.94 ** |
| P2RY6 | CKD vs Normal | +0.62 ns | +0.14 ns | +0.84 ** | +0.75 * |
| P2RY6 | CKD vs AKI | -0.51 ns | -0.51 ns | -0.58 ns | -0.58 ns |
| PDGFRB | Early_DKD vs Normal | -0.04 ns | -0.07 ns | -0.07 ns | -0.04 ns |
| PDGFRB | Advanced_DKD vs Normal | +0.77 ** | +0.76 ** | +0.77 ** | +0.77 ** |
| PDGFRB | AKI vs Normal | +0.56 ns | +0.72 * | +0.33 ns | +0.09 ns |
| PDGFRB | CKD vs Normal | +0.94 ** | +1.00 ** | +0.90 ** | +0.78 * |
| PDGFRB | CKD vs AKI | +0.66 * | +0.74 * | +0.79 * | +0.71 * |
| Calcium score | Early_DKD vs Normal | -0.52 ns | -0.07 ns | -0.26 ns | -0.44 ns |
| Calcium score | Advanced_DKD vs Normal | +0.79 ** | +0.83 ** | +0.78 ** | +0.81 ** |
| Calcium score | AKI vs Normal | +0.82 ** | +0.92 ** | +1.00 *** | +0.94 ** |
| Calcium score | CKD vs Normal | +0.68 * | +0.84 ** | +0.90 ** | +0.84 ** |
| Calcium score | CKD vs AKI | -0.40 ns | -0.25 ns | -0.45 ns | -0.45 ns |
| ECM score | Early_DKD vs Normal | +0.11 ns | +0.15 ns | +0.22 ns | +0.22 ns |
| ECM score | Advanced_DKD vs Normal | +0.93 *** | +0.94 *** | +0.92 *** | +0.94 *** |
| ECM score | AKI vs Normal | +0.98 *** | +0.90 ** | +1.00 *** | +1.00 *** |
| ECM score | CKD vs Normal | +1.00 ** | +1.00 ** | +1.00 ** | +1.00 ** |
| ECM score | CKD vs AKI | -0.43 ns | +0.19 ns | -0.58 ns | -0.51 ns |
| Pyrimidine score | Early_DKD vs Normal | +0.89 * | +1.00 ** | +0.89 * | +0.93 * |
| Pyrimidine score | Advanced_DKD vs Normal | +0.16 ns | +0.64 ** | +0.26 ns | +0.22 ns |
| Pyrimidine score | AKI vs Normal | +0.64 * | -0.76 ** | +0.52 ns | +0.58 ns |
| Pyrimidine score | CKD vs Normal | +0.65 ns | -0.90 ** | +0.21 ns | +0.52 ns |
| Pyrimidine score | CKD vs AKI | -0.06 ns | -0.14 ns | -0.09 ns | -0.12 ns |

