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
  library(pipeComp)
  library(data.table)
})
knitr::opts_chunk$set(echo = TRUE,
                      cache = TRUE)

source(system.file("extdata", "scrna_alternatives.R", package="pipeComp"))

pip_def <- scrna_pipeline(pipeClass = "sce")
print(pip_def)
alternatives <- list(
  doubletmethod=c("none"),
  filt=c("filt.lenient", "filt.stringent"),
  norm=c("norm.seurat", "norm.sctransform", "norm.scran"),
  sel=c("sel.vst"),
  selnb=2000,
  dr=c("seurat.pca"),
  clustmethod=c("clust.seurat"),
  dims=c(10, 15, 20, 30),
  resolution=c(0.01, 0.1, 0.2, 0.3, 0.5, 0.8, 1, 1.2, 2)
)

print(checkPipelinePackages(alternatives, scrna_pipeline()))
datasets <- list.files(path="test-files", pattern = "(.)*.RDS")
datasets <- paste("test-files", datasets, sep="/")
datasets <- c("E-MTAB-9221" = datasets)
print(datasets)
print(typeof(datasets))
sce <- readRDS(datasets)
print(class(sce))
print(sce[[f[1]]])
print(head(sce$phenoid))
print(colnames(sce))
# dat_list <- lapply(datasets, function (x) readRDS(x))
res <- runPipeline(datasets, alternatives, pip_def, nthreads=1, debug = TRUE)
print(res)



