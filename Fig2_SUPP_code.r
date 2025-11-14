#--------- Code for Sup Fig. 2
#libraries
library(Seurat)
library(monocle3)

set.seed(17)

#working directory
setwd('/home/astavropoulou/Notch3_PAPER/')

#read the processed seurat file with the Human data:
seurat_obj <- readRDS('analysis/5_Human_datasets/stroma/seurat_harmony_Qi_Lee_Pelka_Khaliq_Stromal_filtered.rds')

#-------------------change the names of the clusters 
new.cluster.ids <- c("Fibroblasts",'Fibroblasts','Fibroblasts','CAFs','SMCs','BECs','LECs','Pericytes')
names(new.cluster.ids) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj, new.cluster.ids)
seurat_obj$cell_types_2<- Idents(seurat_obj)

seurat_obj$cell_types_2 <- factor(seurat_obj$cell_types_2, levels = c('Fibroblasts','CAFs','SMCs','BECs','LECs','Pericytes'))


pdf('Supp_Fig2_a_b_umaps_human.pdf', width=10, height = 4)
DimPlot(seurat_obj, reduction = "umap", group.by = c("disease_state",'orig.ident'))
dev.off()

pdf('Supp_Fig2_c_dotplot_human.pdf', width=8, height = 3)
genes<- c("COL1A1, DPT, PDGFRA, PROCR, FAP, INHBA, THBS2, ACTA2, MYH11, CNN1, KDR, PLVAP, PTPRB, ERG, LYVE1, CCL21, PROX1, RGS5, NOTCH3, PDGFRB")
genes<-strsplit(genes, ", ")[[1]]
DotPlot(seurat_obj, features = unique(genes), group.by = 'cell_types_2', dot.scale = 8,cluster.idents=F) + FontSize(10) + theme(axis.text.x= element_text(angle=90, hjust=1))
dev.off()
