suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(optparse)))
suppressMessages(suppressWarnings(library(magrittr)))
suppressMessages(suppressWarnings(library(patchwork)))
suppressMessages(suppressWarnings(library(ggplot2)))

## ===== define a function of the analysis ===== ##
ScanneR <- function(seu.list,
                    output_dir,
                    mt_pattern,
                    ribo_pattern,
                    plot_group_by,
                    mt_cutoff,
                    lower_cutoff) {
  ## ----- handle input seurat obj list ----- ##
  # User may input .rds or RData of the seu.list, we will read it in
  if (grepl("\\.rds$", seu.list)) {
    seu.list <- readRDS(seu.list)
  } else if (grepl("\\.RData$", seu.list)) {
    temp_env <- new.env()
    load(seu.list, envir = temp_env)
    
    # Look for the first list in the loaded environment
    seurat_vars <- ls(temp_env)
    seu.lists <- Filter(function(x)
      inherits(temp_env[[x]], "list"), seurat_vars)
    
    if (length(seu.lists) == 0) {
      stop("No Seurat object list found in the RData file.")
    } else if (length(seu.lists) > 1) {
      warning("Multiple Seurat object lists found. Using the first one: ",
              seu.lists[1])
    }
    
    seu.list <- temp_env[[seu.lists[1]]]
    message("Using Seurat object from .RData: ", seu.lists[1])
  } else {
    stop("The input file must be either .rds or .RData format.")
  }
  
  # Check if the input is a list of Seurat object
  # must be a list
  if (!is.list(seu.list)) {
    stop("seu.list must be a list.")
  }
  
  # every element of the list must be a Seurat object
  is_seurat <- vapply(seu.list, inherits, logical(1), what = "Seurat")#vapply(X, FUN, FUN.VALUE, ...)
  if (!all(is_seurat)) {
    bad <- which(!is_seurat)
    stop(
      sprintf(
        "All elements of seu.list must be Seurat objects.  Offending indices: %s",
        paste(bad, collapse = ", ")
      )
    )
  }
  
  ## ----- handle output directory ----- ##
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  output_dir <- file.path(output_dir, paste0("ScanneR_", timestamp))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  ## ----- calculate percentage of mitochondrial genes and ribosomal genes ----- ##
  ##Ribosomal genes also tend to be very highly represented, and can vary between cell types, so it can be instructive to see how prevalent they are in the data.
  ##These are ribosomal protein genes rather than the actual rRNA, so they are more a measure of the translational activity of the cell rather than the cleanliness of the polyA selection.
  for (i in seq_along(seu.list)) {
    seu.list[[i]][["percent.mt"]] <-
      PercentageFeatureSet(seu.list[[i]], pattern = mt_pattern)
    seu.list[[i]][["percent.ribo"]] <-
      PercentageFeatureSet(seu.list[[i]], pattern = ribo_pattern)
  }
  saveRDS(seu.list, file = file.path(output_dir, "seu_list_with_mt_ribo.rds"))
  message(
    "Seurat objects with mitochondrial and ribosomal gene percentages saved as 'seu_list_with_mt_ribo.rds' in the output directory: ",
    output_dir
  )
  
  ## ----- merge for visualization ----- ##
  if (is.null(names(seu.list))) {
    names(seu.list) <- paste0("Sample", seq_along(seu.list))
  }
  seu_m <- merge(seu.list[[1]],
                 y = seu.list[2:length(seu.list)],
                 add.cell.ids = names(seu.list))
  saveRDS(seu_m, file = file.path(output_dir, "seu_merged_for_plots.rds"))#save this just in case the y.max needs to be adjusted later
  message("Merged Seurat object for visualization purpose saved as 'seu_merged_for_plots.rds'")
  
  pv1 <- VlnPlot(
    seu_m,
    features = c("nFeature_RNA"),
    pt.size = 0,
    ncol = 1,
    group.by = plot_group_by
  ) +
    NoLegend() +
    xlab(label = "")
  pv2 <- VlnPlot(
    seu_m,
    features = c("nCount_RNA"),
    pt.size = 0,
    ncol = 1,
    group.by = plot_group_by
  ) +
    NoLegend() +
    xlab(label = "")
  pv3 <- VlnPlot(
    seu_m,
    features = c("percent.mt"),
    pt.size = 0,
    ncol = 1,
    group.by = plot_group_by
  ) +
    NoLegend() +
    xlab(label = "")
  pv4 <- VlnPlot(
    seu_m,
    features = c("percent.ribo"),
    pt.size = 0,
    ncol = 1,
    group.by = plot_group_by
  ) +
    NoLegend() +
    xlab(label = "")
  
  n_samples <- length(seu.list)
  p <- wrap_plots(pv1, pv3, pv2, pv4, ncol = 2)
  png(
    filename = file.path(output_dir, "plot_qc.png"),
    width = 200 + 60 * n_samples,
    height = 800
  )
  print(p)
  dev.off()
  
  message("QC plots saved as 'plot_qc.png' in the output directory: ",
          output_dir)
  
  ## ----- detection based filtering suggestion ----- ##
  df <- FetchData(seu_m, vars = c("nFeature_RNA", plot_group_by))
  
  thresholds <- df %>%
    group_by(!!as.name(plot_group_by)) %>%
    summarise(
      mean = mean(nFeature_RNA, na.rm = TRUE),
      sd = sd(nFeature_RNA, na.rm = TRUE),
      upper = round((mean + 3 * sd), digits = 0),
      lower = lower_cutoff
    ) %>%
    mutate(x = as.integer(factor(!!as.name(plot_group_by))))
  # create x = 1, 2, 3, … in the order of the factor levels
  
  p1 <- VlnPlot(seu_m,
                features = "nFeature_RNA",
                group.by = plot_group_by,
                pt.size = 0) +
    geom_segment(
      data = thresholds,
      aes(
        x = x - 0.4,
        xend = x + 0.4,
        y = lower,
        yend = lower
      ),
      inherit.aes = FALSE,
      color = "tomato",
      linetype = "dashed"
    ) +
    geom_segment(
      data = thresholds,
      aes(
        x = x - 0.4,
        xend = x + 0.4,
        y = upper,
        yend = upper
      ),
      inherit.aes = FALSE,
      color = "tomato",
      linetype = "dashed"
    ) +
    NoLegend() +
    labs(x = "")
  saveRDS(thresholds, file = file.path(output_dir, "qc_thresholds.rds"))
  message(
    "Suggested QC thresholds saved as 'qc_thresholds.rds' in the output directory: ",
    output_dir
  )
  
  ##mt based filtering
  p2 <- VlnPlot(
    seu_m,
    features = "percent.mt",
    group.by = plot_group_by,
    pt.size = 0,
    ncol = 1
  ) +
    xlab(label = "") +
    geom_hline(yintercept = mt_cutoff,
               color = "tomato",
               linetype = "dashed") +
    NoLegend()
  
  p3 <- wrap_plots(p1, p2, ncol = 1)
  png(
    filename = file.path(output_dir, "plot_qc_filter.png"),
    width = 200 + 45 * n_samples,
    height = 800
  )
  print(p3)
  dev.off()
  
  ## ----- Visualization of feature scatter plot ----- ##
  plot_list <- list()
  for (i in seq_along(seu.list)) {
    plot_list[[i]] <-
      FeatureScatter(seu.list[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt")  + NoLegend() +
      FeatureScatter(seu.list[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")  + NoLegend() +
      FeatureScatter(seu.list[[i]], feature1 = "percent.ribo", feature2 = "percent.mt") + NoLegend()
    plot_list[[i]] <- plot_list[[i]] + labs(caption = names(seu.list)[i])
  }
  
  p4 <- wrap_plots(plot_list, ncol = 1)
  png(
    filename = file.path(output_dir, "plot_sca.png"),
    width = 300 * 3,
    height = 300 * n_samples
  )
  print(p4)
  dev.off()
  
  message("Feature scatter plots saved as 'plot_sca.png' in the output directory: ",
          output_dir)
}

## ===== define options for the script ===== ##
description_text <- '
This script takes one named list of Seurat objects (or file paths to .rds/.RData), calculates mitochondrial and ribosomal gene percentages per sample, and then merges them for QC visualization. It produces violin plots of gene counts, feature counts, percent.mt and percent.ribo; suggests filtering thresholds based on mean±3×SD of nFeature_RNA; and generates scatter plots of QC metrics. All intermediate Seurat objects, thresholds tables, and composite figures are saved in a timestamped subdirectory.

ScanneR expects a list of Seurat objects without QC filtering. It writes out:
  • A Seurat list with percent.mt/ribo calculation in the metadata   
  • A merged Seurat object for plotting, in case y.max to be adjusted later (`seu_merged_for_plots.rds`)  
  • QC violin plots (`plot_qc.png`)  
  • QC filter violin plots (`plot_qc_filter.png`)  
  • QC scatter plots (`plot_sca.png`)  
  • A table of suggested nFeature_RNA thresholds (`qc_thresholds.rds`)

Usage:
  Rscript ScanneR.R \\
    --seu_list <Path to a seurat object list .rds|.RData> \\ 
    --output_dir <output_directory> (default: current directory) \\ 
    --mt_pattern <mito_gene_pattern> (default: "^mt-") \\ 
    --ribo_pattern <ribo_gene_pattern> (default: "^(rpl|rps)") \\ 
    --plot_group_by <metadata_column> (default: "orig.ident") \\ 
    --mt_cutoff <numeric_value> (default: 10)
'
  
option_list <- list(
  make_option(
    c("--seu_list"),
    type = "character",
    help = "Path to a .rds or .RData file containing a named list of Seurat objects."
  ),
  make_option(
    c("--output_dir"),
    type = "character",
    default = ".",
    help = "Directory to save output files. Default is current directory."
  ),
  make_option(
    c("--mt_pattern"),
    type = "character",
    default = "^mt-",
    help = "Pattern to identify mitochondrial genes. Default is '^mt-'."
  ),
  make_option(
    c("--ribo_pattern"),
    type = "character",
    default = "^(rpl|rps)",
    help = "Pattern to identify ribosomal genes. Default is '^(rpl|rps)'."
  ),
  make_option(
    c("--plot_group_by"),
    type = "character",
    default = "orig.ident",
    help = "Metadata column to group by for plotting. Default is 'orig.ident'."
  ),
  make_option(
    c("--mt_cutoff"),
    type = "numeric",
    default = 10,
    help = "Cutoff for mitochondrial gene percentage. Default is 10."
  ),
  make_option(
    c("--lower_cutoff"),
    type = "numeric",
    default = 200,
    help = "Lower cutoff for nFeature_RNA. Default is 200."
  )
)

opt_parser <- OptionParser(option_list = option_list, description = description_text)
opt <- parse_args(opt_parser)

## ===== check the input parameters ===== ##
ScanneR(
  seu.list = opt$seu_list,
  output_dir = opt$output_dir,
  mt_pattern = opt$mt_pattern,
  ribo_pattern = opt$ribo_pattern,
  plot_group_by = opt$plot_group_by,
  mt_cutoff = opt$mt_cutoff,
  lower_cutoff = opt$lower_cutoff
)