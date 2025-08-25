setwd("~/current_work/tigit/")
library(Seurat)
library(ggplot2)
library(gtools)
seu = readRDS("dataset.rds")


repair.signature = list(
  patsalos = c("Igf1", "Gdf15", "Tgfbr2", "Tgfbr1", "Slc40a1", "Tpp1", "Adgre1",
               "Ms4a7", "Siglec1", "Itgax", "Aif1", "Mertk", "Folr2", "Gas6",
               "C1qc", "C1qa", "C1qb", "Cd74", "H2-Aa", "H2-Ab1", "H2-Eb1", "Snx5",
               "Pld3", "Dnmt3a", "Hpgds", "Pla2g15", "Serpinb6a", "Acsbg1", "Fcrls",
               "Timp2", "Ctsb", "Apoe", "Gpnmb", "Grn", "Selenop", "Pltp", "Trem2",
               "Lipa", "Cxcl16", "Acp5", "Chil3", "Chil4", "Sell", "Cd177",
               "Itgb7", "Gsr", "Ace", "Ifitm6", "Hp", "Ccr2", "Plac8", "Ly6c2",
               "Ly6c1", "Serpinb10"),

  werner = c("Tgfb1", "Tgfb2", "Tgfb3", "Ccl2", "Cxcl2", "Cxcl15", "Il6", "Il10",
             "Ccl3", "Vegfa", "Pgf", "Fgf2", "Angpt1", "Angpt2", "Hgf", "Ccn1",
             "Csf2", "Cxcl10", "Pdgfb", "Ccn2", "Igf1", "Igf2", "Ngf", "Inhba",
             "Inhbb", "Inhbc", "Inhbe", "Fgf7", "Fgf10", "Egf", "Tgfa", "Hbegf",
             "Glyr1", "Lep"),

  ## These two are more recent, they are focused on T-cells
  belkaid = c("Furin", "Mmp10", "Areg", "Fgf16", "Fgf18"), # Bigger Weight, maybe combine this with burzyn
  
  ## Bigger Weight here
  burzyn = c("Itgae", "Klrg1", "Ccr2", "Il10", "Havcr2", "Il1rl1", "Ikzf2", "Nrp1",
             "Il23r", "Rorc", "Bhlhe40", "Chd7", "Vcl", "Rgs2",
             "Dusp1", "Klf10", "Ankrd6", "Pcsk1", "Sytl2", "Plp2", "Cck",
             "Fam185a", "Agpat2", "Il12rb1", "Ccr3", "Cd80", "Disp3",
             "Rhoc", "Mapkapk3", "Bcl2l1", "Niban2", "Tnf", "Plxdc1", "Areg"),
  
  burzynplus = c("Itgae", "Klrg1", "Ccr2", "Il10", "Havcr2", "Il1rl1",
                 "Ikzf2", "Nrp1", "Il23r", "Rorc", "Bhlhe40", "Chd7",
                 "Vcl", "Rgs2", "Dusp1", "Klf10", "Ankrd6", "Pcsk1",
                 "Sytl2", "Plp2", "Agpat2", "Il12rb1", "Ccr3", ## Remove Cck, not expressed
                 "Cd80", "Disp3", "Rhoc", "Mapkapk3", "Bcl2l1", "Tnf",
                 "Plxdc1", "Areg", "Hif1a", "Prdm1", "Tgfb1", "Nfil3",
                 "Batf", "Tigit")
)

## TFs from burzyn

TFs_burzynplus = c("Ikzf2", "Hif1a", "Batf", "Prdm1", "Rhoc",
                   "Bhlhe40", "Klf10", "Rorc")

genes = as.character(unlist(repair.signature))
genes = unique(genes[(genes %in% rownames(seu))])


ct_colors = c("#FFCC99", "#F0A0FF", "#FFA405", "#FFA8BB", "#993F00", "#C20088", ## Cd8 cols
              "#426600", "#2BCE48", "#94FFB5", "#9DCC00", ## Treg cols
              "#4C005C", "#0075DC", "#740AFF", "#003380", "#CCCCFF") ## Cd4 cols

feature.palette = c("lightgrey", "#E30000")


theming = list(FontSize(main=10), NoAxes())
featureplot.theming = theme(legend.key.width = unit(0.11, 'cm'), legend.key.height = unit(0.3, 'cm'),  legend.text = element_text(size=7))

sig = repair.signature$burzynplus
## sig = sig[1:length(sig)-1] ## Remove tigit, it should be last in order

seu = ScaleData(seu, features = sig, do.center=FALSE)
repair.sig.mat = as.data.frame(seu@assays$RNA@scale.data)
seu@meta.data$repair_score = colSums(repair.sig.mat[rownames(repair.sig.mat) != "Tigit",])


genes = c("Sell", "Cd44", "Il7r", "Ccr7", "Mki67", "Tcf7", "Slamf6",
          "Tbx21", "Gata3", "Rorc", "Bcl6", "Tox", "Cxcr3", "Cxcr5",
          "Ifng", "Il21", "Il10", "Pdcd1", "Havcr2", "Lag3",
          "Tigit", "Cd160", "Gzma", "Il2", "Ccr6" , "Irf4","Ccr4",
          "Fas", "Cd8a", "Cd4", "Il2ra", "Foxp3")

celltype.marker.genes = genes

for (cond in c("Genotype", "Condition", "Tissue") ){  
  fn.path = paste0("figures/umap-celltypes-labels_", cond)
  message(fn.path)
  DimPlot(seu, cols = ct_colors, split.by = cond, raster = TRUE, pt.size = 1.2) + theming
  ggsave(paste0(fn.path, ".pdf"), width = 11.4, height = 5.46)
  ggsave(paste0(fn.path, ".png"), width = 11.4, height = 5.46)  
}



featureplot.theming = theme(
  legend.position="right",
  legend.justification="top",
  legend.key.width = unit(0.12, 'cm'),
  legend.key.height = unit(0.25, 'cm'),
  legend.text = element_text(size=7),
  plot.title = element_text(size=10),
  legend.margin=margin(0,0,0,0),
  legend.box.spacing = unit(0, "pt")
  )

genes = celltype.marker.genes

p = FeaturePlot(subset(seu, Condition == "Naive"), genes,
                order = TRUE,
                cols = feature.palette,
                pt.size = 0.05,
                ncol = 4,
                raster = FALSE) & update_geom_defaults("point", aes(stroke = 0)) & theming & featureplot.theming

ggsave("figures/featurePlots-MarkerGenes-Naive.png", width = 5.0, height = 9)



p = FeaturePlot(subset(seu, Condition != "Naive"), genes,
                order = TRUE,
                cols = feature.palette,
                pt.size = 0.05,
                ncol = 4,
                raster = FALSE) & update_geom_defaults("point", aes(stroke = 0)) & theming & featureplot.theming

ggsave("figures/featurePlots-MarkerGenes-Infected.png", width = 5.0, height = 9)

library(dplyr)

plot.cell.numbers <- function(df, gb = "CellType") {
  df.sizes = df %>% group_by(!!sym(gb), Condition, Genotype, Tissue) %>% summarise(cells = n())
  ggplot(df.sizes, aes(y=!!sym(gb), x=cells, fill = Condition)) +
  geom_bar(stat = "identity", position= position_dodge(width=0.7, preserve = "single"), width= 0.5) +  scale_x_continuous(trans = 'log10') +  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  annotation_logticks(sides = "b",
                      base = 10,
                      short = unit(0.01, "cm"),
                      mid = unit(0.06, "cm"),                      
                      long = unit(0.12, "cm")) +
  facet_grid(Genotype~Tissue) + xlab("Cells") + ylab("")
}

plot.cell.numbers(seu[[]])

ggsave("figures/cells-per-sample-all-sum.pdf", width = 6, height = 4)

plot.cell.numbers(seu[[]], gb = "ct_short")

ggsave("figures/cells-per-sample-all-ct_short-sum.pdf", width = 6, height = 9)




p = FeaturePlot(subset(seu, Genotype == "WT"), c("Tigit","repair_score"),
                order = TRUE,                
                pt.size = 0.5,
                ## blend.threshold = 0,
                min.cutoff="q01",
                max.cutoff="q99",
                raster = FALSE,
                blend = TRUE) & update_geom_defaults("point", aes(stroke = 0))
p
ggsave("figures/tigit-repair-signature-coexpression2.png", width = 15, height = 4)


p = FeaturePlot(subset(seu, Genotype == "WT"), c("Tigit","Areg"),
                order = TRUE,                
                pt.size = 0.5,
                ## blend.threshold = 0,
                min.cutoff="q10",
                max.cutoff="q99",
                raster = FALSE,
                blend = TRUE) & update_geom_defaults("point", aes(stroke = 0))
p
ggsave("figures/tigit-areg-signature-coexpression2.png", width = 15, height = 4)
library(pheatmap)
library(reshape2)
library(RColorBrewer)
seu.wt = subset(seu, Genotype == "WT")
for(ct in unique(seu.wt@meta.data$CellType)) {
  message(ct)
  seu.sub = subset(seu.wt, CellType == ct)
  
  seu.sub.sub = seu.sub
  fn = paste0("figures/heatmap-dendro_wt-", ct, "-merged.pdf")
  message(fn)
  gene.wide = as.data.frame(t(seu.sub.sub@assays$RNA@scale.data))    
  stopifnot(all(colnames(gene.wide) %in% sig) & length(sig) ==  ncol(gene.wide))
  X = cor(gene.wide)
  X[!is.finite(X)] = 0
  pdf(fn)
  print(pheatmap(X,
                 color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                 main = ct,
                   scale = "none",
                   cluster_rows = T,
                   cluster_cols = T,
                   breaks = pracma::linspace(-0.15, 0.15, 101)
                   ))
  dev.off()
}

library(ggpubr)
library(dplyr)
fig.height=3
fig.width = 6
expr.vector = function(x){
  y = as.numeric(x)
  y
}

ct.covar = "CellType"

scatter.genes = list(tigit.areg = c("Tigit", "Areg"),
                     areg.prdm1 = c("Areg", "Prdm1"),
                     hif1a.areg = c("Hif1a", "Areg"),
                     tigit.prdm1 = c("Tigit","Prdm1"),
                     tigit.hif1a = c("Tigit","Hif1a")
                     )
seu.active = seu

font.size.annotation = 4
df = seu.active[[]]
## df$ExpCond = with(df, paste(Condition, Tissue))

df$technical.replicate = sample(1:10,nrow(df), replace = TRUE)

for (g in c("Prdm1", "Tigit", "Areg", "Hif1a")) {  
  df[[paste0(g,"_norm")]] = expr.vector(seu.active@assays$RNA@data[g,])
}

for (g in scatter.genes) {  
  x = paste0(g[1], "_norm")
  y = paste0(g[2], "_norm")
  df.pseudo = df %>% group_by(!!rlang::sym(ct.covar), Condition, Genotype, technical.replicate) %>%
    summarize(x = mean(!!rlang::sym(x)), y = mean(!!rlang::sym(y)))
  ggplot(as.data.frame(df.pseudo),
         aes(x = x, y = y, color = !!rlang::sym(ct.covar), shape = Genotype))  + geom_point(size=2,stroke=0.4) + scale_shape(solid = FALSE) + ylim(c(0,max(df.pseudo$y)*1.45)) + xlab(g[1]) + ylab(g[2]) +
    facet_wrap(vars(Condition)) + 
    stat_cor(aes(x = x, y = y, col=!!rlang::sym(ct.covar)), method = "pearson", inherit.aes = FALSE, show.legend=FALSE, size = font.size.annotation) + guides(color = guide_legend(override.aes = list(size=3.5))) +
    theme_classic()
  ggsave(paste0("figures/scatter_pseudobulk-", g[1], "-", g[2], ".pdf"), width = fig.width, height = fig.height)  
  ## ggplot(as.data.frame(df.pseudo),
  ##        aes(x = x, y = y, color = !!rlang::sym(ct.covar), shape = Genotype)) + geom_point() + xlab(g[1]) + ylab(g[2]) +
  ##   facet_wrap(vars(Condition), ncol=1) +
  ##   stat_cor(aes(x = x, y = y, col=!!rlang::sym(ct.covar)), method = "pearson", inherit.aes = FALSE, show.legend = FALSE, size = font.size.annotation) +
  ##   theme_classic()
  ## ggsave(paste0("figures/scatter_pseudobulk-", g[1], "-", g[2], "_t.pdf"), width = fig.height + 1, height = 6)
}



for (g in c("Prdm1", "Tigit", "Areg")) {  
  x = paste0(g, "_norm")
  y = "repair_score"
  library(dplyr)
  df.pseudo = df %>% group_by(!!rlang::sym(ct.covar), Condition, Genotype, technical.replicate) %>% summarize(x = mean(!!rlang::sym(x)), y = mean(!!rlang::sym(y)))
  ggplot(as.data.frame(df.pseudo), aes(x = x, y = y, col = !!rlang::sym(ct.covar), shape = Genotype)) + geom_point(size=2,stroke=0.4) + scale_shape(solid = FALSE) + ylim(c(min(df.pseudo$y), max(df.pseudo$y)*1.45)) + xlab(g) + ylab("Repair Score") + facet_wrap(~ Condition) + stat_cor(aes(x = x, y = y, col=!!rlang::sym(ct.covar)), method = "pearson", inherit.aes = FALSE, show.legend=FALSE, size = font.size.annotation) + theme_classic() 
  ggsave(paste0("figures/scatter_pseudobulk-",g, "-", "repair_score.pdf"), width = fig.width, height = fig.height)
  ## ggplot(as.data.frame(df.pseudo), aes(x = x, y = y, col = !!rlang::sym(ct.covar), shape = Genotype)) + geom_point() + xlab(g) + ylab("Repair Score") + facet_wrap(~ Condition, ncol=1) + stat_cor(aes(x = x, y = y, col=!!rlang::sym(ct.covar)), method = "pearson", inherit.aes = FALSE, show.legend=FALSE, size = font.size.annotation) + theme_classic()
  ## ggsave(paste0("figures/scatter_pseudobulk-",g, "-", "repair_score_t.pdf"), height = 7, width = fig.height+1)
}


genes = c("Sell", "Cd44", "Il7r", "Ccr7", "Mki67", "Tcf7", "Slamf6",
          "Tbx21", "Gata3", "Rorc", "Bcl6", "Tox", "Cxcr3", "Cxcr5",
          "Ifng", "Il21", "Il10", "Pdcd1", "Havcr2", "Lag3",
          "Tigit", "Cd160", "Gzma", "Il2", "Ccr6" , "Irf4","Ccr4",
          "Fas", "Cd8a", "Cd4", "Il2ra", "Foxp3")

seu.wt = subset(seu, Genotype == "WT")

p = DotPlot(
  seu.wt,
  genes,  
  group.by = "ct_short"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + FontSize(10) + xlab("") + ylab("")

ggsave("figures/dotplot-marker-genes.pdf", width = 12, height = 4)


p = DotPlot(
  seu.wt,
  genes,  
  group.by = "clustersFinal"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + FontSize(10) + xlab("") + ylab("")

ggsave("figures/dotplot-marker-genes-clusters.pdf", width = 12, height = 4)


seu.wt = subset(seu, Genotype == "WT")
## Idents(seu) = seu.wt@meta.data$ct_short
## seu.wt@meta.data$ct_short =  as.character(seu.wt@meta.data$ct_short)
DotPlot(seu.wt, sig, group.by = "ct_short") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + FontSize(10) + xlab("") + ylab("")
ggsave("figures/dotplot-repair_signature.pdf", width = 12)


seu.wt = subset(seu, Genotype == "WT")
## Idents(seu) = seu.wt@meta.data$ct_short
## seu.wt@meta.data$ct_short =  as.character(seu.wt@meta.data$ct_short)
DotPlot(seu.wt, c("Tigit", TFs_burzynplus), group.by = "celltype")
ggsave("figures/dotplot-tfs.pdf", width = 12)

p = FeaturePlot(subset(seu, Genotype == "WT"), c("repair_score"),
                order = TRUE,                
                pt.size = 0.5,
                raster = FALSE
                ) + theming
p

p = FeaturePlot(subset(seu, Genotype == "WT"), c("Tigit","repair_score"),
                order = TRUE,                
                pt.size = 0.5,
                ## blend.threshold = 0,
                min.cutoff="q01",
                max.cutoff="q99",
                raster = FALSE,
                blend = TRUE) & update_geom_defaults("point", aes(stroke = 0))
p
ggsave("figures/tigit-repair-signature-coexpression2.png", width = 15, height = 4)


p = FeaturePlot(subset(seu, Genotype == "WT"), c("Tigit","Areg"),
                order = TRUE,                
                pt.size = 0.5,
                ## blend.threshold = 0,
                min.cutoff="q10",
                max.cutoff="q99",
                raster = FALSE,
                blend = TRUE) & update_geom_defaults("point", aes(stroke = 0))
p
ggsave("figures/tigit-areg-signature-coexpression2.png", width = 15, height = 4)
