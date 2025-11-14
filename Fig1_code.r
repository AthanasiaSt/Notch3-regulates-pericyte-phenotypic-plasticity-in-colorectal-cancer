#--------- Code for Main Fig. 1
#libraries
library(Seurat)
library(ggplot2)

set.seed(17)

#working directory
setwd('/home/astavropoulou/Notch3_PAPER/')

#read the processed seurat file:
seurat_obj <- readRDS('analysis/2_integration/Cnt_AOM_DSS_colon_allCells_Harmony_origIdent.rds')

#change the levels of the seurat obbject
seurat_obj$cell_types_2 <- factor(seurat_obj$cell_types_2, levels = c('Fibroblasts','Pericytes','SMCs','BECs','LECs','Proliferating'))

# Fig1.a, Fig1.b
pdf('Fig1_a_b_umaps.pdf', width=5, height = 5)
DimPlot(seurat_obj, reduction = "umap", group.by = c("cell_types_2"))
DimPlot(seurat_obj, reduction = "umap", group.by = c("orig.ident"),  cols = c("#B90E0A", "lightgrey"))
dev.off()

# Fig1.c
pdf('Fig1_c_marker_genes_dotplot.pdf', width = 8, height = 3)
genes<- c("Col1a1, Dpt, Pdgfra, Rgs5, Notch3, Ndufa4l2, Acta2, Myh11, Cnn1, Plvap, Pecam1, Kdr, Lyve1, Reln, Ccl21a, Birc5, Mki67, Top2a")
genes<-strsplit(genes, ", ")[[1]]
DotPlot(seurat_obj, features = unique(genes), group.by = 'cell_types_2', dot.scale = 8,cluster.idents=F) + FontSize(10) + theme(axis.text.x= element_text(angle=90, hjust=1))
dev.off()


#-----------------------------Figure 1.d
df<-as.data.frame(table(seurat_obj$orig.ident,seurat_obj$cell_types_2)/rowSums(table(seurat_obj$orig.ident,seurat_obj$cell_types_2)))

pdf('Fig1_d_Percentage_celltypes.pdf', width=4, height = 2)
ggplot(df,                  # Stacked barplot using ggplot2
       aes(x = Freq,
           y = Var1,
           fill = Var2)) + 
  labs( 
    y="samples", x='percentage')+
  
  geom_bar(stat = "identity")  + theme(text = element_text(size=3),axis.text  = element_text(size=3)) 

dev.off()

#-----------------------------Figure 1.g

#preprocessing again
subset<-seurat_obj[,seurat_obj$cell_types_2  == 'Proliferating']
subset <- RunPCA(subset)
subset <- FindNeighbors(subset, dims = 1:20)
subset <- FindClusters(subset, resolution = 0.3)  # Adjust resolution for finer subclustering
subset <- RunUMAP(subset, dims = 1:20)

pdf('Fig1_g_Proliferating_featureplot.pdf', width = 3, height = 3)
FeaturePlot(
  subset,
  features = c("Pecam1"),
  blend = F, order = TRUE, 
  blend.threshold = 0.3, reduction = 'umap',
  pt.size = 0.5, combine = T, min.cutoff = 'q9', ncol = 1)
FeaturePlot(
  subset,
  features = c("Rgs5"),
  blend = F, order = TRUE, 
  blend.threshold = 0.3, reduction = 'umap',
  pt.size = 0.5, combine = T, min.cutoff = 'q9', ncol = 1)
FeaturePlot(
  subset,
  features = c("Pdgfra"),
  blend = F, order = TRUE, 
  blend.threshold = 0.3, reduction = 'umap',
  pt.size = 0.5, combine = T, min.cutoff = 'q9', ncol = 1)
dev.off()


#Fig1.i 
#read the processed seurat file with the Human data:
seurat_obj <- readRDS('analysis/5_Human_datasets/stroma/seurat_harmony_Qi_Lee_Pelka_Khaliq_Stromal_filtered.rds')

#-------------------change the names of the clusters 
new.cluster.ids <- c("Fibroblasts",'Fibroblasts','Fibroblasts','CAFs','SMCs','BECs','LECs','Pericytes')
names(new.cluster.ids) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj, new.cluster.ids)
seurat_obj$cell_types_2<- Idents(seurat_obj)

seurat_obj$cell_types_2 <- factor(seurat_obj$cell_types_2, levels = c('Fibroblasts','Pericytes','SMCs','BECs','LECs','CAFs'))

my_colors <- c(
  'Fibroblasts' = '#F8766D',  # red
  'Pericytes'   = '#C49A00',  # mustard
  'SMCs'        = '#3C8C00',  # darker green
  'BECs'        = '#33C6E8',  # lighter blue
  'LECs'        = '#0091C5',  # darker cyan
  'CAFs'        = '#7B5FE0'   # darker purple
)


pdf('Fig1_i_umaps_human.pdf', width=6, height = 6)
DimPlot(seurat_obj, reduction = "umap", group.by = c("cell_types_2"), cols = my_colors) 
dev.off()


#-----------------------------Figure 1.j
df<-as.data.frame(table(seurat_obj$disease_state,seurat_obj$cell_types_2)/rowSums(table(seurat_obj$disease_state,seurat_obj$cell_types_2)))
df<- df[df$Var1 != 'Border',]
pdf('Fig1_j_Percentage_celltypes.pdf', width=3, height = 1.5)
ggplot(df,                  # Stacked barplot using ggplot2
       aes(x = Freq,
           y = Var1,
           fill = Var2)) + 
  labs( 
    y="samples", x='percentage')+
  
  geom_bar(stat = "identity") +scale_fill_manual(values = my_colors)  + theme(text = element_text(size=3),axis.text  = element_text(size=3)) 

dev.off()

