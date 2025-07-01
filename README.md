<p align="center">
  <img src="ScanneR_hex.png" width="180"/>
</p>

# ScanneR: QC Pipeline for Single-Cell RNA-Seq with Seurat

**ScanneR** is an R-based command-line pipeline for quality control of single-cell RNA-seq data using the Seurat framework. It computes mitochondrial and ribosomal gene percentages, suggests filtering thresholds, and produces violin and scatter plots.

## 📦 Features

- Accepts named Seurat object **lists** via `.rds` or `.RData`  
- Calculates per-cell mitochondrial gene percentage (`percent.mt`) and ribosomal gene percentage (`percent.ribo`) using user-defined regex patterns  
- Suggests `nFeature_RNA` thresholds based on **mean ± 3 SD** and a lower bound  
- Generates:
  - Violin plots of QC metrics (`nFeature_RNA`, `nCount_RNA`, `percent.mt`, `percent.ribo`)
  - Violin plots with dashed filter cutoffs overlaid  
  - Scatter plots of QC metrics  
- Saves intermediate Seurat objects and thresholds table  
- Timestamped output folder to avoid overwrites  

## 🚀 Quick Start

### 🔧 Requirements

```r
install.packages(c(
  "tidyverse",
  "Seurat",
  "optparse",
  "magrittr",
  "patchwork",
  "ggplot2"
))
```

### 🖥️ Usage

```bash
Rscript ScanneR.R \
  --seu_list path/to/your_seurat_list.rds \
  --output_dir results/ \
  --mt_pattern "^mt-" \
  --ribo_pattern "^(rpl|rps)" \
  --plot_group_by orig.ident \
  --mt_cutoff 10 \
  --lower_cutoff 200
```
Run `Rscript ScanneR.R --help` to see all options.

### 📝 Parameters

- `--seu_list`  
  **(Required)** Path to a .rds or .RData file containing a named list of Seurat objects..

- `--output_dir`  
  Directory where output files will be saved. Default is `./`.

- `--mt_pattern`  
  Regex pattern to identify mitochondrial genes. Default is `^mt-`.

- `--ribo_pattern`  
  Regex pattern to identify ribosomal protein genes. Default is `^(rpl|rps)`.

- `--plot_group_by`  
  Metadata column used to group cells in the QC plots. Default is `orig.ident`.

- `--dot_topN_wilcox`  
  Number of top markers per cluster to show in the Wilcoxon-based dot plot. Default is `5`.

- `--mt_cutoff`  
  Percent mitochondrial cutoff for filtering.. Default is `10`.

- `--lower_cutoff`  
  Lower cutoff for nFeature_RNA threshold. Default is `200`.

## 📂 Output

All files are saved in a timestamped subdirectory under `--output_dir` (e.g. `ScanneR_20250701_153045`).

- RDS files
  - seu_list_with_mt_ribo.rds — your input list annotated with %mt and %ribo
	- seu_merged_for_plots.rds — merged Seurat object for plotting
  - qc_thresholds.rds    — data frame of lower and upper nFeature_RNA per group
- QC plots
  - plot_qc.png          — violin plots of QC metrics
	- plot_qc_filter.png   — violin plots with dashed threshold lines
	- plot_sca.png         — scatter plots of QC metrics

## 📌 Notes

- If multiple named lists are present in an .RData, the first one is selected (with a warning).
- Unnamed lists will be auto-named Sample1, Sample2, … before merging.

