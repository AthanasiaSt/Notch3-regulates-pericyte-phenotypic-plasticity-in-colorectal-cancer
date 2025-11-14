#--------- Code for Main Fig. 7
#libraries
library(Seurat)
library(monocle3)

set.seed(17)

#working directory
setwd('/home/astavropoulou/Notch3_PAPER/')

# load Pericyte subclustering
seurat_obj <- readRDS('analysis/5_Human_datasets/Pericyte_analysis/resolution_1/Final/seurat_harmony_Qi_Lee_Pelka_Khaliq_Pericytes_subclustering_res1_subcategories_final.rds')

#-------------------change the names of the clusters 
new.cluster.ids <- c("myPC",'ecmPC_T','ecmPC_N','intPC','prolifPC')
names(new.cluster.ids) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj, new.cluster.ids)
seurat_obj$subcategories_detailed<- Idents(seurat_obj)


# Fig7.a_b
pdf('Fig7_a_b_umap.pdf', width=5, height = 5)
DimPlot(seurat_obj, reduction = "umap", group.by = c("subcategories_detailed"))
DimPlot(seurat_obj, reduction = "umap", group.by = c("disease_state"))
dev.off()

# Fig7.c
pdf('Fig7_c_vlnplot.pdf', width=4, height =20)
genes<- c("RGS5, Col3a1, Mmp11, Acta2, Myh11, Cnn1, Mki67, Top2a, Notch3, Heyl")
genes<-strsplit(genes, ", ")[[1]]
VlnPlot(seurat_obj, alpha = 0.1, group.by = 'subcategories_detailed', features = toupper(genes),ncol = 1, pt.size = 0, combine = T)
dev.off()

# Fig7.d
library(UCell)

signatures<- list()
signatures$myPCs<- c('MYH11','TAGLN','PLN','CNN1')
signatures$ecmPCs<- c('RGS5','MMP11','COL3A1')

seurat_obj <- AddModuleScore_UCell(seurat_obj,  features = signatures)

pdf('Fig7_d_Ucell_signatures.pdf', width = 5, height = 5)
FeaturePlot(seurat_obj, reduction = 'umap', features = 'myPCs_UCell' ,  ncol = 1, pt.size = 0.5, combine = T)  & NoAxes()
FeaturePlot(seurat_obj, reduction = 'umap', features = 'ecmPCs_UCell' ,  ncol = 1, pt.size = 0.5, combine = T)  & NoAxes()

dev.off()

#Fig7.h

monocle_object <- load_monocle_objects(directory_path='analysis/5_Human_datasets/Pericyte_analysis/Monocle/Human_Pericytes_subclustering_monocle/')

pdf('Fig7_h_monocle_trajectory.pdf', width = 5, height = 5)

plot_cells(monocle_object,
           color_cells_by = "subcategories_detailed", cell_size = 1,
           graph_label_size=3,
           label_cell_groups=F,
           show_trajectory_graph = TRUE)# + scale_color_manual(values = colours_celltypes)

plot_cells(monocle_object,
           color_cells_by = "pseudotime", cell_size = 1,
           graph_label_size=3,
           label_cell_groups=F,
           show_trajectory_graph = TRUE, )


dev.off()

#Fig7.h
pdf('Fig7_i_monocle_umaps_genes.pdf', width = 5, height=5)
plot_cells(monocle_object, genes=c("RGS5","ACTA2",'NOTCH3','COL3A1','CNN1','HEYL'),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE,cell_size = 1)

dev.off()


# Fig7.e, f, g
seurat_obj_all_cells<- readRDS('analysis/5_Human_datasets/stroma/seurat_harmony_Qi_Lee_Pelka_Khaliq_Stromal_filtered.rds')

pdf('Fig7_e_pericytes_inflammatory.pdf', width=6, height =8)
FeaturePlot(
  seurat_obj,
  features = c("IL6",'CCL8','CXCL1', 'HEYL'),
  blend = F, order = TRUE, 
  blend.threshold = 0.3, reduction = 'umap',
  pt.size = 0.5, combine = T, min.cutoff = 'q9', ncol = 2)

FeaturePlot(
  seurat_obj_all_cells,
  features = c("PLVAP", "RGS5"),
  blend = F, order = TRUE, 
  blend.threshold = 0.3, reduction = 'umap',
  pt.size = 0.5, combine = T, min.cutoff = 'q9', ncol = 1)
dev.off()

