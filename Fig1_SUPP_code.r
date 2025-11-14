#--------- Code for Supplementary Fig. 1
#libraries
library(Seurat)
library(ggplot2)

set.seed(17)

#working directory
setwd('/home/astavropoulou/Notch3_PAPER/')

#load the processed seurat file with the Human data:
seurat_obj <- readRDS('analysis/3_subclustering/BECs/BECs_subclustering.rds')
new.cluster.ids <- c('CapECs',"VenECs",'ArtECs')
names(new.cluster.ids) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj, new.cluster.ids)
seurat_obj$cell_types<- Idents(seurat_obj)


pdf('Supp_Fig1_a_umaps_BECs.pdf', width=5, height = 5)
DimPlot(seurat_obj, group.by = 'cell_types', reduction = "umap")
dev.off()

pdf('Supp_Fig1_b_dotplot_BECs.pdf', width=5, height = 6)
genes<- c("Thrsp, Ccdc85a, Rgcc, Piezo2, Scgb3a1, Rbp7, Ramp3, Hecw2, Cd300lg, Ubd, Ackr1, Cd14, Fjx1, Adgrg6, Cxcl9, Ctsh, Selp, Cxcl1, Lbp, Gja5, Mal, Nebl, Bmx, S100a4, Pcsk5, Sema3g, Enpp6, Sdcbp2, Slc45a4")
genes<-strsplit(genes, ", ")[[1]]
DotPlot(seurat_obj, features = rev(unique(genes)), group.by = 'cell_types', dot.scale = 8,cluster.idents=F) + FontSize(10) + theme(axis.text.x= element_text(angle=90, hjust=1)) +   coord_flip() 
dev.off()

pdf('Supp_Fig1_c_vlnplots_BECs.pdf', width=7, height = 10)
VlnPlot(seurat_obj, alpha = 0.1, group.by = 'cell_types', features = c('Rgcc','Ackr1','Sema3g'),ncol = 1, pt.size = 1, combine = T)
dev.off()

