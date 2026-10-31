# Human_P2RY6

Analysis code for the P2RY6 / pyrimidine-metabolism study in human kidney and
diabetic kidney disease (DKD), using public bulk RNA-seq (Fan et al.) and the
KPMP Kidney Tissue Atlas single-cell (scRNA-seq) and single-nucleus
(snRNA-seq) data.

## Repository layout

```
Human_P2RY6/
├── 02-Scripts/      R analysis scripts (S20 … S90) and the driver S00_RUNALL.R
├── 04-Outputs/      Figures, tables and intermediate .rds files, one folder per script
├── renv.lock        Exact versions of R and all R packages used
├── renv/            renv infrastructure (activate.R, settings)
└── .Rprofile        Activates renv when R is started in the project folder
```

Scripts that were used in earlier versions of the analysis but are no longer
part of it (e.g. the KPMP bulk RNA-seq and combined bulk analyses, S30–S43)
have been removed; they remain available in the git history.

## Setting up the R environment

The R package environment is managed with [renv](https://rstudio.github.io/renv/).
The analysis was run with **R 4.4.1** and Bioconductor 3.20; `renv.lock`
records the exact version and source of every package (CRAN, Bioconductor and
GitHub).

**Before running any script, restore the environment once:**

1. Install R 4.4.x and the system libraries needed to build the packages,
   in particular HDF5 for `hdf5r` (Debian/Ubuntu:
   `sudo apt install libhdf5-dev libcurl4-openssl-dev libssl-dev libxml2-dev`;
   macOS: `brew install hdf5`).
2. Start R in the project root folder (`Human_P2RY6/`, or open
   `Human_P2RY6.Rproj` in RStudio). The `.Rprofile` bootstraps renv
   automatically.
3. Run:

   ```r
   renv::restore()
   ```

   This installs all packages listed in `renv.lock` into a project-private
   library. Installation of Seurat and its dependencies can take a while.

The main packages used by the scripts:

| Package    | Version     | Source                              |
|------------|-------------|-------------------------------------|
| Seurat     | 5.2.1       | CRAN                                |
| SeuratDisk | 0.0.0.9021  | GitHub: `mojaveazure/seurat-disk`   |
| hdf5r      | 1.3.12      | CRAN                                |
| Matrix     | 1.7.0       | CRAN                                |
| edgeR      | 4.4.2       | Bioconductor                        |
| dplyr      | 1.1.4       | CRAN                                |
| tidyr      | 1.3.1       | CRAN                                |
| readr      | 2.1.5       | CRAN                                |
| stringr    | 1.5.1       | CRAN                                |
| rlang      | 1.1.6       | CRAN                                |
| ggplot2    | 3.5.2       | CRAN                                |
| ggpubr     | 0.6.0       | CRAN                                |
| ggsignif   | 0.6.4       | CRAN                                |
| cowplot    | 1.1.3       | CRAN                                |
| knitr      | 1.50        | CRAN                                |
| pander     | 0.6.6       | CRAN                                |

If a package is added to a script, record it with `renv::snapshot()` and
commit the updated `renv.lock`.

## Raw data (not included)

The raw data files are **not** part of this repository. Download them from the
sources below and place them at the paths given. All paths are relative to
`02-Scripts/`, so the data folders sit **next to** the `Human_P2RY6/` project
folder:

```
<parent folder>/
├── Human_P2RY6/          this repository
├── Fan_DKD_2020/
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

The analysis uses the processed, open-access Seurat objects (format
h5Seurat, workflow "Aggregated Clustered Data"), not raw FASTQ/BAM files
(which are controlled access). Both objects belong to the KPMP dataset with
DOI [10.48698/16dd-vj20](https://doi.org/10.48698/16dd-vj20). To list them in
the repository, filter by **Data Format = h5Seurat**
(<https://atlas.kpmp.org/repository/?size=n_1000_n&filters%5B0%5D%5Bfield%5D=data_format&filters%5B0%5D%5Bvalues%5D%5B0%5D=h5Seurat&filters%5B0%5D%5Btype%5D=any>)
and pick the two files below. Older Explorer versions with similar names
are listed there too; check the file name and size.

If this link stops working, set the filter by hand (labels as of October 2026):

1. Open <https://atlas.kpmp.org/repository>.
2. In the filter panel on the left (open it with the arrow button if it is
   collapsed), expand **Dataset metadata**.
3. Expand **Data Format** and tick **h5Seurat**. Alternatively, expand
   **DOIs** and tick **10.48698/16dd-vj20**, the dataset that contains the
   two files used here.
4. Optionally narrow down further with **Experimental Strategy**
   (Single-cell RNA-Seq / Single-nucleus RNA-Seq).
5. Download a file with the download button in its row. The files are open
   access, so no data access request is needed.

The files used in this analysis:

| Experimental strategy | File name in the repository | Size (bytes) | DOI |
|-----------------------|-----------------------------|--------------|-----|
| Single-cell RNA-Seq   | `KPMP_PREMIERE_SC_version2_ForExplorer_RemovedBatchEffect_Final2025.h5seurat` | 7,302,897,490 | 10.48698/16dd-vj20 |
| Single-nucleus RNA-Seq | `SingleNucleus_KPMP_Explorer_05182025.h5Seurat` | 8,013,374,548 | 10.48698/16dd-vj20 |

Downloaded files can carry a file-ID prefix (e.g.
`b8732f23-4372-4ea1-b374-4779b318381e_SingleNucleus_KPMP_Explorer_05182025.h5Seurat`);
remove the prefix so that the names match the paths below. The open-access
clinical data (release 2025-12-05, `20251205_OpenAccessClinicalData`) are
downloaded from the same repository as a zip file.

| Data | Expected path | Used by |
|------|---------------|---------|
| scRNA-seq object | `../../KPMP_scRNA/KPMP_PREMIERE_SC_version2_ForExplorer_RemovedBatchEffect_Final2025.h5seurat` | S50, S53, S75, S83, S86, S87 |
| Open-access clinical data, scRNA copy | `../../KPMP_scRNA/OpenAccessClinicalData/20251205_OpenAccessClinicalData.csv` | S50, S52, S53, S75, S80, S84, S86, S87, S88 |
| snRNA-seq object | `../../KPMP_snRNA/SingleNucleus_KPMP_Explorer_05182025.h5Seurat` | S60, S75, S83, S86, S87, S90 |
| Open-access clinical data, snRNA copy (unzipped download, folder name kept as-is) | `../../KPMP_snRNA/6d3047c2-588f-4ad4-8eba-f4c46a1e5c14_20251205_OpenAccessClinicalData/20251205_OpenAccessClinicalData.csv` | S60, S62, S75, S80, S83, S84, S86, S87, S88, S90 |

The two clinical-data paths hold the same file (open-access clinical data
release of 2025-12-05). If a newer release is downloaded, update the file
names at the top of the scripts listed above.

## Running the analysis

After `renv::restore()`, set the working directory to `02-Scripts/` (all paths
are relative to it) and run the scripts in order:

```r
setwd("02-Scripts")
source("S00_RUNALL.R")
```

Each script writes its results to `04-Outputs/<script name>/`. Later scripts
read the `.rds` files written by earlier ones (e.g. S52 reads S50, S62 reads
S60), so keep the order of `S00_RUNALL.R`. S80 assembles the manuscript
figures from these outputs; S81 checks the participant numbers and
significance labels of every figure panel.

### Memory requirements

Several scripts are **memory heavy** and may crash with out-of-memory errors
(R aborting with `cannot allocate vector of size …`, or the R session being
killed by the operating system) on machines with too little RAM. The KPMP
scRNA-seq and snRNA-seq objects are about 7 GB and 7.5 GB on disk and take
considerably more memory once loaded into Seurat. This affects:

- **S50, S53** – load the scRNA-seq object into Seurat;
- **S60** – loads the snRNA-seq object into Seurat;
- **S83** – loads both objects (raw counts only), one after the other;
- **S90** – loads the snRNA-seq object (raw counts only).

We recommend ≥ 64 GB of RAM, running these scripts in a fresh R session one at
a time (`Rscript S50_KPMP_scRNA_loadOnce.R` from within `02-Scripts/`) rather
than via `S00_RUNALL.R`, and closing other memory-intensive programs. The
remaining scripts work on the per-participant summaries in `04-Outputs/` (S75,
S86 and S87 read only metadata columns of the `.h5Seurat` files via `hdf5r`)
and run on an ordinary workstation.
