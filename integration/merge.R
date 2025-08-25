library(Seurat)
setwd("~/current_work/tigit/")
base.dir = "./integration/"
file.prefix = paste0(base.dir, lapply(list.files(base.dir, pattern = "*h5"), function (x){paste0(strsplit(x,"_")[[1]][1:4],sep="_", collapse = "")}))
names(file.prefix) = lapply(list.files(base.dir, pattern = "*h5"), function (x){paste0(strsplit(x,"_")[[1]][1:4],sep="_", collapse = "")})
## strplit(files, "_")

all.mat = lapply(file.prefix, function(f) {message(f);ReadMtx(paste0(f,"matrix.mtx.gz"), paste0(f,"barcodes.tsv.gz"), paste0(f,"genes.tsv.gz"))})

cells = lapply(all.mat, function(x){ x[,(colSums(x>0) > 500)] })

cells.X = lapply(names(cells), function(n){
  x <- cells[[n]]
  colnames(x) = paste0(n,colnames(x))
  x
})

genes = rownames(cells.X[[1]])

stopifnot(all(as.logical(lapply(cells.X, function(x){all(genes == rownames(x))}))))
X = do.call(cbind,cells.X)


seu = readRDS("/home/main/current_work/tigit/dataset-with-celltypes.rds")

batch.prefix = unlist(lapply(colnames(X), function (x){strsplit(x,"_")[[1]][1]}))
seu.pub = CreateSeuratObject(X)
seu.pub@meta.data$sample_id = batch.prefix
seu.pub@meta.data$dataset = "public"

## Fuck Seurat
manual.merge.features = intersect(rownames(seu.pub), rownames(seu))
seu.sub = seu[manual.merge.features,]
seu = NULL
gc()
seu.sub.pub = seu.pub[manual.merge.features,]
saveRDS(seu.sub.pub,"pub.rds")

library(harmony)
library(dplyr)
library(Matrix)

seu.sub@meta.data$dataset = "tigit"
seu.sub.pub@meta.data$dataset = "public"
## Did I tell you seurat 5 sucks?
stopifnot(all(rownames(seu.sub@assays$RNA@counts) == rownames(seu.sub.pub@assays$RNA@features)))

X.new.sub = seu.sub.pub@assays$RNA@layers$counts
colnames(X.new.sub) = colnames(seu.sub.pub)
rownames(X.new.sub) = rownames(seu.sub.pub)

X = cbind(seu.sub@assays$RNA@counts, X.new.sub)
## genes = rownames(X)
## rownames(X) = paste0("G",rownames(X))

df = bind_rows(seu.sub@meta.data, seu.sub.pub@meta.data)
rownames(df) = colnames(X)
seu.merge = CreateSeuratObject(X, meta.data = df)

process.obj <- function(seu) {
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu, selection.method = "mvp")
  hvg <- VariableFeatures(object=seu)
  seu <- seu %>% ScaleData() %>% RunPCA(features = hvg)
  ## seu <- FindNeighbors(seu, dims = 1:16, verbose = FALSE)
  ## seu <- FindClusters(seu, resolution = 0.5, verbose = FALSE)
  ## seu <- RunUMAP(seu, dims = 1:16, verbose = FALSE)
  seu
}
seu.merge <- process.obj(seu.merge)

library(harmony)
seu.all <- RunHarmony(seu.merge,c("sample_id"),plot=T)
seu.all <- FindNeighbors(seu.all, reduction="harmony", verbose = FALSE)
seu.all <- FindClusters(seu.all, resolution = 0.5, verbose = FALSE)
seu.all <- RunUMAP(seu.all, reduction = "harmony", dims=1:50)

DimPlot(seu.all, group="orig.ident", raster = TRUE)
ggsave("integration-pub.png")

library(ggplot2)
seu.all@meta.data$CellType = df$celltype
DimPlot(seu.all, group="CellType")
ggsave("integration-pub-celltype.png")

library(ggplot2)
DimPlot(seu.all, group="dataset")
ggsave("integration-pub-dataset.png")

head(seu.all@meta.data)
