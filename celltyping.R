setwd("~/current_work/tigit/")
library(Seurat)
library(ggplot2)
library(gtools)
seu = readRDS("dataset-with-celltypes.rds")

seu.sub <- subset(seu, seurat_clusters %in% c("14"))
seu.sub@meta.data$coarse <- seu.sub@meta.data$seurat_clusters
seu.sub <- FindNeighbors(seu.sub, reduction="harmony", dims = 1:20)
seu.sub <- FindClusters(seu.sub, resolution = 0.5)

seu@meta.data$sub14 = "other"
seu@meta.data[colnames(seu.sub),]$sub14 = seu.sub@meta.data$seurat_clusters
seu@meta.data[colnames(seu.sub),]$sub14 = seu.sub@meta.data$seurat_clusters

seu@meta.data$clustersFinal = as.character(seu@meta.data$clustersFinal)
seu@meta.data$clustersFinal[seu@meta.data$sub14 %in% c("2")] = "14_CD4"
seu@meta.data$clustersFinal[seu@meta.data$sub14 %in% c("3")] = "14_CD8"
seu@meta.data$clustersFinal[seu@meta.data$sub14 %in% c("4")] = "14_Treg"
seu@meta.data$clustersFinal[seu@meta.data$sub14 %in% c("1")] = "14_CD4"
seu@meta.data$clustersFinal = factor(seu@meta.data$clustersFinal, levels= mixedsort(unique(seu@meta.data$clustersFinal)))

cluster.mappings = list("Naive CD4 T cells" =c("0", "11", "12"),
                        "Naive CD8 T cells" = "1",
                        "Effector Tregs" = "2",
                        "Effector CD8 T cells" = c("3", "14", "14_CD8"),
                        "Effector CD4 T cells" = c("14_CD4"),
                        "Naive Tregs" = "4",
                        "Stem-like exhausted CD4 T cells" = "5",
                        "Stem-like exhausted CD8 T cells" = "10",
                        "Exhausted CD4 T cells" = "6_CD4",
                        "Exhausted CD8 T cells" = "6_CD8",
                        "Proliferating CD4 T cells" = "7_CD4",
                        "Proliferating CD8 T cells" = "7_CD8",
                        "Central Memory CD8 T cells" = "8",
                        "Proliferating Tregs" = "9",
                        "Tox+ Tregs" = c("13", "14_Treg")
                        )

seu@meta.data$celltype = "None"

for(name in names(cluster.mappings)) {
  seu@meta.data$celltype[seu@meta.data$clustersFinal %in% cluster.mappings[[name]]] = name
}

new.levels = c(
  "Proliferating CD8 T cells",
  "Central Memory CD8 T cells",
  "Effector CD8 T cells",
  "Exhausted CD8 T cells",
  "Stem-like exhausted CD8 T cells",
  "Naive CD8 T cells",
  "Naive Tregs",  
  "Effector Tregs",
  "Proliferating Tregs",
  "Tox+ Tregs",  
  "Proliferating CD4 T cells",
  "Exhausted CD4 T cells",
  "Effector CD4 T cells",
  "Stem-like exhausted CD4 T cells",
  "Naive CD4 T cells"
  )

seu@meta.data$celltype = factor(seu@meta.data$celltype, levels  = new.levels)

seu@meta.data$ct_short = seu@meta.data$celltype

levels(seu@meta.data[["ct_short"]]) <- c(
  "Ki67+ CD8",
  "CM CD8",
  "Eff CD8",
  "Exh CD8",
  "Stem-like exh CD8",
  "Naive CD8",
  "Naive Tregs",  
  "Eff Tregs",
  "Ki67+ Tregs",
  "Tox+ Tregs",  
  "Ki67+ CD4",
  "Exh CD4",
  "Eff CD4",
  "Stem-like exh CD4",
  "Naive CD4"  
  )



write.csv(with(seu@meta.data, table(ct_short, celltype)), "ct_short-celltype-mappings.csv")

levels(seu@meta.data[["Tissue"]]) <- c("Lung", "Spleen")
seu@meta.data[["Condition"]] <- factor(seu@meta.data$condition)
levels(seu@meta.data[["Condition"]]) <- c("LCMV", "Naive")
seu@meta.data[["Condition"]] <- factor(seu@meta.data$Condition, levels=c("Naive", "LCMV"))
seu@meta.data[["Genotype"]] <- factor(seu@meta.data$Genotype, levels=c("WT", "KO"))

seu@meta.data$CellType = "none"
seu@meta.data$CellType[grepl("CD4", seu@meta.data$celltype)] = "CD4"
seu@meta.data$CellType[grepl("CD8", seu@meta.data$celltype)] = "CD8"
seu@meta.data$CellType[grepl("Treg", seu@meta.data$celltype)] = "Treg"
seu@meta.data$CellType = factor(seu@meta.data$CellType, levels= c("CD4", "CD8", "Treg"))
Idents(seu) <- seu@meta.data$celltype

saveRDS(seu, "dataset.rds")
