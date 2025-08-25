library(Seurat)
library(presto)
library(ggplot2)
setwd("/home/main/current_work/tigit/")

## Initial object provided by Nicole

seu <- readRDS("merged_celltypes-final.rds")

levels(seu@meta.data[["Tissue"]]) <- c("Lung", "Spleen")
seu@meta.data[["Condition"]] <- factor(seu@meta.data$condition)
levels(seu@meta.data[["Condition"]]) <- c("LCMV", "Naive")
seu@meta.data[["Condition"]] <- factor(seu@meta.data$Condition, levels=c("Naive", "LCMV"))
seu@meta.data[["Genotype"]] <- factor(seu@meta.data$Genotype, levels=c("WT", "KO"))

process.obj <- function(seu) {
  library(dplyr)
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu, selection.method = "mvp")
  hvg <- VariableFeatures(object=seu)
  seu <- seu %>% ScaleData() %>% RunPCA(features = hvg)
  seu <- FindNeighbors(seu, dims = 1:16, verbose = FALSE)
  seu <- FindClusters(seu, resolution = 0.5, verbose = FALSE)
  seu <- RunUMAP(seu, dims = 1:16, verbose = FALSE)
  seu
}
seu <- process.obj(seu)

saveRDS(seu, "dataset-with-celltypes.rds")

