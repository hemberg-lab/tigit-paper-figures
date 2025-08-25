setwd("~/current_work/tigit")
library(Seurat)
seu.orig = readRDS("merged-with-pub.rds")

seu = subset(seu.orig, dataset == "public" | (CellType == "Treg" & Genotype == "WT"))
seu@meta.data$dataset[seu@meta.data$dataset == "public"] = "Delacher"

library(harmony)
process.obj <- function(seu) {
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu, selection.method = "mvp")
  hvg <- VariableFeatures(object=seu)
  seu <- seu %>% ScaleData(features=hvg) %>% RunPCA(features = hvg)
  seu
}

seu.treg <- process.obj(seu)

seu.treg <- RunHarmony(seu.treg, c("dataset", "sample_id"), plot = T, dims.use = 1:20)
seu.treg <- FindNeighbors(seu.treg, reduction="harmony", verbose = FALSE)
seu.treg <- FindClusters(seu.treg, resolution = 0.5, verbose = FALSE)
seu.treg <- RunUMAP(seu.treg, reduction = "harmony",dims=1:20)

saveRDS(seu.treg, "tigit-delachler-treg.rds")
