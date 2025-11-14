#--------- Code for Main Fig. 6
#libraries
library(Seurat)
library(monocle3)

set.seed(17)

#working directory
setwd('/home/astavropoulou/Notch3_PAPER/')

# load Pericyte subclustering
seurat_obj <- readRDS('analysis/3_subclustering/Pericytes/Pericyte_subclustering.rds')
seurat_obj$cell_types <- factor(seurat_obj$cell_types, levels = c('ecmPCs','myPCs','iPCs','endoPCs'))

# Fig6.a
pdf('Fig6_a_umap.pdf', width=5, height = 5)
DimPlot(seurat_obj, reduction = "umap", group.by = c("cell_types"))
dev.off()

# Fig6.b
pdf('Fig6_b_dotplot.pdf', width=5, height = 6)
genes<- c("Rgs5, Notch3, Pdgfrb, Ndufa4l2, Mcam, Heyl, Col3a1, Vtn, Mmp11, Tagln, Acta2, Myh11, Cnn1, Pln, Cxcl1, Il6, Mmp10, Saa3, Pecam1, Plvap, Cdh5, Kdr")
genes<-strsplit(genes, ", ")[[1]]
DotPlot(seurat_obj, features = rev(unique(genes)), group.by = 'cell_types', dot.scale = 8,cluster.idents=F) +  FontSize(10) + theme(axis.text.x= element_text(angle=45, hjust=1)) +   coord_flip()
dev.off()

# explanation on myPCs, why they are not SMCs 
seurat_obj <- readRDS('analysis/3_subclustering/Pericytes/Pericyte_subclustering.rds')
seurat_obj$cell_types <- factor(seurat_obj$cell_types, levels = c('ecmPCs','myPCs','iPCs','endoPCs'))

seurat_obj2<- readRDS('analysis/2_integration/Cnt_AOM_DSS_colon_allCells_Harmony_origIdent.rds')
seurat_obj2<-seurat_obj2[,seurat_obj2$cell_types_2 %in% c('SMCs','BECs')]

# merge Seurat objects
combined <- merge(seurat_obj, y = seurat_obj2, add.cell.ids = c("Pericytes", "SMCs_BECs"))
combined$cell_types <- factor(combined$cell_types, levels = c('ecmPCs','myPCs','iPCs','endoPCs','SMCs','BECs'))

pdf('Explain_fig6b_myPCs_SMCs_BECs.pdf', width=8, height = 6)
DotPlot(
    combined,
    features = rev(unique(genes)),
    group.by = "cell_types",   # or any metadata variable
    dot.scale = 8) + coord_flip()
dev.off()


# Fig6.e
pdf('Fig6_e_vlnplot.pdf', width=5, height =8)
VlnPlot(seurat_obj, alpha = 0.5, group.by = 'cell_types', features = c('Notch3','Heyl'),ncol = 1, pt.size = 1, combine = T)
dev.off()

#Fig6.f

monocle_object <- load_monocle_objects(directory_path='analysis/5_monocle/Pericytes_subclustering_monocle_umap/AOM_Pericytes_subclustering_monocle/')

pdf('Fig6_f_monocle_trajectory.pdf', width = 5, height = 5)

plot_cells(monocle_object,
           color_cells_by = "cell_types", cell_size = 1,
           graph_label_size=3,
           label_cell_groups=F,
           show_trajectory_graph = TRUE)# + scale_color_manual(values = colours_celltypes)

plot_cells(monocle_object,
           color_cells_by = "pseudotime", cell_size = 1,
           graph_label_size=3,
           label_cell_groups=F,
           show_trajectory_graph = TRUE, )


dev.off()

#Fig6.g
pdf('Fig6_g_monocle_umaps_genes.pdf', width = 5, height=5)
plot_cells(monocle_object, genes=c("Notch3","Heyl",'Rgs5','Acta2','Col3a1'),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE,cell_size = 1)

dev.off()

# Fig6.h
seurat_obj_becs<- readRDS('analysis/3_subclustering/BECs/BECs_subclustering.rds')

pdf('Fig6_h_pericytes_BECs_featureplot.pdf', width=4, height =8)
FeaturePlot(
  seurat_obj,
  features = c("Jag1", "Dll4"),
  blend = F, order = TRUE, 
  blend.threshold = 0.3, reduction = 'umap',
  pt.size = 0.5, combine = T, min.cutoff = 'q9', ncol = 1)

FeaturePlot(
  seurat_obj_becs,
  features = c("Jag1", "Dll4"),
  blend = F, order = TRUE, 
  blend.threshold = 0.3, reduction = 'umap',
  pt.size = 0.5, combine = T, min.cutoff = 'q9', ncol = 1)
dev.off()

#Fig6.I
pdf('Fig6_i_pericytes_Becs_vlnplot.pdf', width=4, height =8)
VlnPlot(seurat_obj, alpha = 0.1, group.by = 'cell_types', features = c('Jag1','Dll4'),ncol = 1, pt.size = 1, combine = T)
VlnPlot(seurat_obj_becs, alpha = 0.1, group.by = 'cell_types', features = c('Jag1','Dll4'),ncol = 1, pt.size = 1, combine = T)
dev.off()

