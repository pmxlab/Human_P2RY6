# Human_P2RY6

Analysis code for the P2RY6 / pyrimidine-metabolism study in human kidney and
diabetic kidney disease (DKD), using public bulk RNA-seq (Fan et al.) and the
KPMP Kidney Tissue Atlas single-cell (scRNA-seq) and single-nucleus
(snRNA-seq) data.

## Repository layout

```
Human_P2RY6/
├── 02-Scripts/      R analysis scripts (S20 … S90) and the driver S00_RUNALL.R
└── 04-Outputs/      Figures, tables and intermediate .rds files, one folder per script
```

## R and package requirements

The analysis was run with R 4.4.1. The packages and versions used:

| Package    | Version     | Source                                   |
|------------|-------------|------------------------------------------|
| Seurat     | 5.2.1       | CRAN                                     |
| SeuratDisk | 0.0.0.9021  | GitHub: `mojaveazure/seurat-disk`        |
| hdf5r      | 1.3.12      | CRAN (needs the system HDF5 library)     |
| Matrix     | 1.7.0       | CRAN                                     |
| edgeR      | 4.4.2       | Bioconductor                             |
| sva        | 3.54.0      | Bioconductor                             |
| dplyr      | 1.1.4       | CRAN                                     |
| tidyr      | 1.3.1       | CRAN                                     |
| readr      | 2.1.5       | CRAN                                     |
| stringr    | 1.5.1       | CRAN                                     |
| rlang      | 1.1.6       | CRAN                                     |
| ggplot2    | 3.5.2       | CRAN                                     |
| ggpubr     | 0.6.0       | CRAN                                     |
| ggsignif   | 0.6.4       | CRAN                                     |
| cowplot    | 1.1.3       | CRAN                                     |
| knitr      | 1.50        | CRAN                                     |
| pander     | 0.6.6       | CRAN                                     |

Installation:

```r
install.packages(c("Seurat", "hdf5r", "Matrix", "dplyr", "tidyr", "readr",
                   "stringr", "rlang", "ggplot2", "ggpubr", "ggsignif",
                   "cowplot", "knitr", "pander", "remotes", "BiocManager"))
BiocManager::install(c("edgeR", "sva"))
remotes::install_github("mojaveazure/seurat-disk")
```

`hdf5r` needs the HDF5 system library (e.g. `sudo apt install libhdf5-dev` on
Debian/Ubuntu, `brew install hdf5` on macOS).

The two KPMP single-cell objects are 7–8 GB each on disk. The scripts that load
them (S50, S53, S60, S75, S83, S86, S87, S90) need a machine with ample RAM
(≥ 64 GB recommended).

## Raw data (not included)

The raw data files are **not** part of this repository. Download them from the
sources below and place them at the paths given. All paths are relative to
`02-Scripts/`, so the data folders sit **next to** the `Human_P2RY6/` project
folder:

```
<parent folder>/
├── Human_P2RY6/          this repository
├── Fan_DKD_2020/
├── KPMP_BulkRNA/
├── KPMP_scRNA/
└── KPMP_snRNA/
```

### 1. Fan et al. bulk RNA-seq (GEO GSE142025)

Fan Y, et al. *Comparison of kidney transcriptomic profiles of early and
advanced diabetic nephropathy reveals potential new mechanisms for disease
progression.* Diabetes 2019;68(12):2301–2314.

- Download: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142025>
  (supplementary file `GSE142025_RAW.tar`)
- Extract the tar archive so that the 36 per-sample files
  (`GSM42177xx_*.txt.gz`) are in:

| Expected path                           | Used by |
|-----------------------------------------|---------|
| `../../Fan_DKD_2020/GSE142025_RAW/`     | S20     |

### 2. KPMP Kidney Tissue Atlas

Kidney Precision Medicine Project (KPMP), Kidney Tissue Atlas:

- Data repository: <https://atlas.kpmp.org/repository>
- Explorer: <https://atlas.kpmp.org/explorer>

Download the following files from the repository:

| Data | Expected path | Used by |
|------|---------------|---------|
| Bulk RNA-seq count files (`*_counts.csv` / `*_Counts.csv`) and the file manifest `0_ATLAS_PATIENT_FILE_METADATA.csv` | `../../KPMP_BulkRNA/` (all files in one folder) | S30, S75 |
| scRNA-seq object | `../../KPMP_scRNA/KPMP_PREMIERE_SC_version2_ForExplorer_RemovedBatchEffect_Final2025.h5seurat` | S50, S53, S75, S83, S86, S87 |
| Open-access clinical data, scRNA copy | `../../KPMP_scRNA/OpenAccessClinicalData/20251205_OpenAccessClinicalData.csv` | S50, S52, S53, S75, S80, S84, S86, S87, S88 |
| snRNA-seq object | `../../KPMP_snRNA/SingleNucleus_KPMP_Explorer_05182025.h5Seurat` | S60, S75, S83, S86, S87, S90 |
| Open-access clinical data, snRNA copy (unzipped download, folder name kept as-is) | `../../KPMP_snRNA/6d3047c2-588f-4ad4-8eba-f4c46a1e5c14_20251205_OpenAccessClinicalData/20251205_OpenAccessClinicalData.csv` | S60, S62, S75, S80, S83, S84, S86, S87, S88, S90 |

The two clinical-data paths hold the same file (open-access clinical data
release of 2025-12-05). If a newer release is downloaded, update the file
names at the top of the scripts listed above.

## Running the analysis

Set the working directory to `02-Scripts/` (all paths are relative to it) and
run the scripts in order:

```r
setwd("02-Scripts")
source("S00_RUNALL.R")
```

Each script writes its results to `04-Outputs/<script name>/`. Later scripts
read the `.rds` files written by earlier ones (e.g. S40 reads S20 and S30), so
keep the order of `S00_RUNALL.R`. S80 assembles the manuscript figures from
these outputs.
