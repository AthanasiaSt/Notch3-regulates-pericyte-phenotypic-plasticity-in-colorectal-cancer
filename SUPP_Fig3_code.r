#--------- Code for SUPP Fig. 4
#libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(readxl)

set.seed(17)

#working directory
setwd('/home/astavropoulou/Notch3_PAPER/')

#load the processed seurat file:
seurat_obj <- readRDS('analysis/8_Vega_et_al_analysis/Vega_allCells.rds')

#SUPP Fig4.a, Fig4.b
pdf('SUPP_Fig4_a_b_umaps.pdf', width=8, height = 8)
DimPlot(seurat_obj, reduction = "umap", group.by = c("celltypes_general"))
DimPlot(seurat_obj, reduction = "umap", group.by = c("group"),  cols = c("#B90E0A", "lightgrey"))
dev.off()

#SUPP Fig4.c
pdf('SUPP_Fig4_c_marker_genes_dotplot.pdf', width = 8, height = 3)
genes<- c("Epcam", "Klf5", "Dsp","Cd52",'Cd53', "Igkc","Dpt", "Col1a1", "Pdgfra", "Plvap", "Pecam1", "Kdr","Rgs5", "Notch3", "Ndufa4l2","Hba-a1", "Hbb-bt", "Hba-a2")
DotPlot(seurat_obj, features =genes , group.by = 'celltypes_general', dot.scale = 8,cluster.idents=F) + FontSize(10) + theme(axis.text.x= element_text(angle=90, hjust=1))
dev.off()

# SUPP Fig4.d
pdf('SUPP_Fig4_Notch3_featureplot.pdf', width = 8, height = 8)

FeaturePlot(
  seurat_obj,
  features = c('Notch3'),
  blend = F, order = TRUE, 
  blend.threshold = 0.3, reduction = 'umap',
  pt.size = 0.5, combine = T, min.cutoff = 'q9', ncol = 1)

dev.off()


#-------SUPP fig 4.e
pdf('SUPP_Fig3_e_vlnplot.pdf', width=5, height =4)
VlnPlot(seurat_obj, split.by = 'group', group.by = 'celltypes_general',cols = c("#B90E0A", "lightgrey"), features = c('Notch3'),ncol = 1, pt.size = 0)
dev.off()


pdf('SUPP_Fig4_Notch3_featureplot.pdf', width = 8, height = 8)

FeaturePlot(
  seurat_obj,
  features = c('Notch3'),
  blend = F, order = TRUE, 
  blend.threshold = 0.3, reduction = 'umap',
  pt.size = 0.5, combine = T, min.cutoff = 'q9', ncol = 1)

dev.off()
