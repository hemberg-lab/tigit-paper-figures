## The co-inhibitory receptor TIGIT promotes tissue protective functions in T cells

Panetti et al. 2025

Scripts for the analysis of the manuscript of the single-cell RNAseq data.

## Experimental design

TIGIT KO and wild-type mouse T cells were sorted and collected from the lung and spleen in steady state and LCMV conditions. Briefly, cells were sorted into three sub-populations and subsequently pooled for scRNAseq. The individual subpopulations were CD3+/CD8+, CD3+/CD4+/Foxp3-, and CD3+/CD4+/Foxp3+.

## Dataset location

The datasets can be found in (https://zenodo.org/records/14041419)[this] Zenodo archive

tigit-full-dataset.rds: contains the full processed dataset ~198645 cells and 25906 genes formatted as a Seurat (v5) object


## Dataset samples

There are several 10x runs per condition (technical replicates)

KO_Lung_LCMV   KO_Lung_Naive  KO_Spleen_LCMV KO_Spleen_Naive    WT_Lung_LCMV
           4             2             6             5             4
WT_Lung_Naive  WT_Spleen_LCMV WT_Spleen_Naive
           3             6             4

#### Copyright

Nikos Patikas 2025

This code is released under  GPL v2 license
