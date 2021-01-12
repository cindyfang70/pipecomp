suppressPackageStartupMessages({
  library(scater) # BioConductor
  library(SingleCellExperiment) # BioConductor
  library(DropletUtils) # BioConductor
  library(tidyverse) # CRAN
  library(here) # CRAN
  library(DT) # CRAN
  library(pheatmap) # CRAN
  library(here)
  library(Matrix)
  library(DESeq2)
})
knitr::opts_chunk$set(echo = TRUE,
                      cache = TRUE)

files <- system("ls test-files", intern = TRUE)
print(files)


for (file in files){
  # reading in the data
  file_address <- paste("test-files", file, sep="/")
  print(file_address)
  list_exp <- paste("ls", file_address, sep = " ")
  print(list_exp)
  exp_files <- system(list_exp, intern = TRUE)
  data_dir <- paste(file_address, exp_files[1], sep = "/")
  print(data_dir)
  data_col_dir <- paste(file_address, exp_files[2], sep = "/")
  data_row_dir <- paste(file_address, exp_files[3], sep = "/")
  cluster_data_dir <-paste(file_address, exp_files[4], sep="/")
  data_design <- paste(file_address, exp_files[5], sep = "/")
  
  data <- readMM(data_dir)
  data_col <- read_tsv(data_col_dir, col_names = FALSE)
  data_row <- read_tsv(data_row_dir, col_names = FALSE)
  cluster_data <- read_tsv(cluster_data_dir, col_names = TRUE)
  expdesign <- read_tsv(data_design)
  expdesign = semi_join(expdesign[,1], data_col[,1], by = c("Assay" = "X1"))

  # creating the SCE
  data <- as.matrix(data)
  sce <- SingleCellExperiment(list(counts = data), colData = as.vector(expdesign))
  colnames(sce) <- as.vector(data_col$X1)
  rownames(sce) <- as.vector(data_row$X1)

  # formatting the clustering data
  sel_k <- which(cluster_data$sel.K)
  clusters <- cluster_data[sel_k, -c(1:2)]
  clusters <- as.vector(as.matrix(clusters))
  cell_names <- names(cluster_data)[-c(1,2)]
  names(clusters) <- cell_names
  sce$phenoid <- as.vector(clusters[colnames(sce)])

  # computing QC metrics and metadata
  sce <- addPerCellQC(sce)
  sce <- addPerFeatureQC(sce)
  sce <- add_meta(sce)

  # renaming some colData to match what's used in pipecomp
  pct_Mt <- sce$pct_Mt
  sce$pct_counts_Mt <- pct_Mt
  pct_top_20 <- sce$percent_top_20
  sce$pct_counts_in_top_20_features <- pct_top_20

  # saving the SCE as an RDS
  sce_filename <- paste(file_address, ".RDS", sep="")
  saveRDS(sce, sce_filename)

 }


