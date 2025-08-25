setwd("~/current_work/tigit")
library(Seurat)
seu.treg <- readRDS("tigit-delachler-treg.rds")
get.tissue <-function(x){strsplit(x,"_")[[1]][3]}
mask = seu.treg@meta.data$dataset == "Delacher"

seu.treg@meta.data$tissue = paste0(seu.treg@meta.data$Tissue,"_",seu.treg@meta.data$Condition,"_", seu.treg@meta.data$dataset)
seu.treg@meta.data$tissue[mask] = paste0(unlist(lapply(colnames(seu.treg)[mask], get.tissue)), "_Delacher")



## table(seu.treg@meta.data$tissue)
## lung_INF_tigit   lung_STS_tigit
## spleen_INF_tigit spleen_STS_tigit

delacher.tissue = c("Fat_Delacher", "Liver_Delacher", "Lung_Delacher", "Skin_Delacher")

delacher.slo = c("Blood_Delacher", "BM_Delacher", "LN_Delacher", "Spleen_Delacher")

seu.treg@meta.data$tissue_show =  seu.treg@meta.data$tissue
seu.treg@meta.data$tissue_show[seu.treg@meta.data$tissue %in% delacher.slo] = "SLO_Delacher"
seu.treg@meta.data$tissue_show[seu.treg@meta.data$tissue %in% delacher.tissue] = "Tissue_Delacher"
seu.treg@meta.data$tissue_show = factor(seu.treg@meta.data$tissue_show, levels = c("SLO_Delacher", "spleen_STS_tigit", "spleen_INF_tigit", "Tissue_Delacher", "lung_STS_tigit", "lung_INF_tigit"))
levels(seu.treg@meta.data$tissue_show) = c("SLO (Delacher)", "Spleen Naive (Tigit)", "Spleen LCMV (Tigit)", "Tissue (Delacher)", "Lung Naive (Tigit)", "Lung LCMV (Tigit)")


## DimPlot(seu.treg,group.by="celltype",split.by="tissue", ncol=6)
## ggsave("Treg-tissue-celltype.png", width = 20, height=8)

## DimPlot(seu.treg,group.by="seurat_clusters",split.by="tissue", ncol=6)
## ggsave("Treg-tissue-seurat_clusters.png", width = 20, height=8)

setwd("~/current_work/tigit")
library(Seurat)
seu.treg <- readRDS("tigit-delachler-treg.rds")
get.tissue <-function(x){strsplit(x,"_")[[1]][3]}
mask = seu.treg@meta.data$dataset == "Delacher"

seu.treg@meta.data$tissue = paste0(seu.treg@meta.data$Tissue,"_",seu.treg@meta.data$Condition,"_", seu.treg@meta.data$dataset)
seu.treg@meta.data$tissue[mask] = paste0(unlist(lapply(colnames(seu.treg)[mask], get.tissue)), "_Delacher")



## table(seu.treg@meta.data$tissue)
## lung_INF_tigit   lung_STS_tigit
## spleen_INF_tigit spleen_STS_tigit

delacher.tissue = c("Fat_Delacher", "Liver_Delacher", "Lung_Delacher", "Skin_Delacher")

delacher.slo = c("Blood_Delacher", "BM_Delacher", "LN_Delacher", "Spleen_Delacher")

seu.treg@meta.data$tissue_show =  seu.treg@meta.data$tissue
seu.treg@meta.data$tissue_show[seu.treg@meta.data$tissue %in% delacher.slo] = "SLO_Delacher"
seu.treg@meta.data$tissue_show[seu.treg@meta.data$tissue %in% delacher.tissue] = "Tissue_Delacher"
seu.treg@meta.data$tissue_show = factor(seu.treg@meta.data$tissue_show, levels = c("SLO_Delacher", "spleen_STS_tigit", "spleen_INF_tigit", "Tissue_Delacher", "lung_STS_tigit", "lung_INF_tigit"))
levels(seu.treg@meta.data$tissue_show) = c("SLO (Delacher)", "Spleen Naive (Tigit)", "Spleen LCMV (Tigit)", "Tissue (Delacher)", "Lung Naive (Tigit)", "Lung LCMV (Tigit)")


## DimPlot(seu.treg,group.by="celltype",split.by="tissue", ncol=6)
## ggsave("Treg-tissue-celltype.png", width = 20, height=8)

## DimPlot(seu.treg,group.by="seurat_clusters",split.by="tissue", ncol=6)
## ggsave("Treg-tissue-seurat_clusters.png", width = 20, height=8)
table(seu.treg@meta.data$tissue_show)

library(ggplot2)
library(dplyr)
cov = "tissue_show"
df_counts = seu.treg@meta.data %>% group_by(!!rlang::sym(cov))%>% slice_sample(n=290) %>% group_by(!!rlang::sym(cov), seurat_clusters) %>% summarize(total_cells = n())
ggplot(df_counts, aes(x = !!rlang::sym(cov), y = total_cells, fill = seurat_clusters))+ xlab("") +
  geom_bar(stat="identity", position="fill")+ ylab("Cell Proportion") + guides(fill=guide_legend(ncol=2)) + coord_flip()
ggsave(paste0("bar-",cov, ".pdf"), width=10, height = 2.5)


DimPlot(seu.treg, group.by="seurat_clusters", split.by="tissue_show", ncol=3, pt.size=0.0) + theming +
  ggtitle("") + xlab("UMAP1") + ylab("UMAP2")
ggsave("Treg-tissue_show-seurat_clusters.pdf", width = 10, height=6)

DimPlot(seu.treg, group.by="celltype", split.by="tissue_show", ncol=3, pt.size=0.0) + theming +
  ggtitle("") + xlab("UMAP1") + ylab("UMAP2")
ggsave("Treg-tissue_show-ct_short.pdf", width = 10, height=6)


genes = c("Tigit", "Prdm1", "Il10", "Foxp3", "Nfil3", "Batf", "Gata3", "Il1rl1", "Klrg1", "Pdcd1", "Tnfrsf4", "Pparg")

library(viridis)


X.old = as.matrix(seu.treg@assays$RNA@layers$data)
factors = rowMaxs(X.old)
X = sweep(X.old, 1, factors, "/")
seu.treg@assays$RNA@layers$data = X

genes = c("Tigit", "Prdm1", "Il10", "Foxp3", "Nfil3", "Batf", "Gata3",
          "Il1rl1", "Klrg1", "Pdcd1", "Tnfrsf4", "Pparg")

for (ti in c("tissue", "tissue_show")) {  
  DoHeatmap(seu.treg, genes, slot = "data", group.by=ti,size = 3.5, disp.max=0.8) + scale_fill_viridis(name = "Normalized\nExpression") # scale_fill_gradient('expression',  low = "white", high = "firebrick3")
  ggsave(paste0("Treg-",ti,"-markers.pdf"), width = 10, height=5)
}


## DoHeatmap(seu.treg, genes, slot = "data", group.by="tissue")
## ggsave("Treg-tissue-markers.pdf", width = 10, height=6)

## DoHeatmap(seu.treg, genes, slot = "data", col.low = "white", color.mid = "pink", color.hi="red", group.by="tissue_show")

## ggsave("Treg-tissue_show-markers.pdf", width = 10, height=6)

all(Idents(seu.treg)== seu.treg@meta.data$seurat_clusters)

Idents(seu.treg) = seu.treg@meta.data$tissue_show
mm = FindAllMarkers(seu.treg)

library(ggplot2)
library(dplyr)
cov = "celltype"
df_counts = seu.treg@meta.data %>% filter(dataset== "tigit")%>% group_by(!!rlang::sym(cov))%>% slice_sample(n=1000) %>% group_by(!!rlang::sym(cov), seurat_clusters) %>% summarize(total_cells = n())
ggplot(df_counts, aes(x = seurat_clusters, y = total_cells, fill = !!rlang::sym(cov)))+
  geom_bar(stat="identity", position="fill") + ylab("Cell Proportion")
ggsave(paste0("bar",cov, "tigit-only.png"), width=5, height = 3)


