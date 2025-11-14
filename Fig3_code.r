#--------- Code for Main Fig. 3
#libraries
library(Seurat)
library(ggplot2)
library(CellChat)
library(patchwork)
library(pheatmap)

set.seed(17)

#working directory
setwd('/home/astavropoulou/Notch3_PAPER/')

#load cellchat object 
cellchat <- readRDS('analysis/4_cellchat/AOM/cell_chat_pval_0_05_logfc_0_2_pct_0_1_AOM_all_cells.rds')
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

#------Fig3_a
pdf('Fig3_a_outgoing_ingoing_signalling_pathways.pdf', width = 15, height = 15)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", width = 5, height = 15)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", width = 5, height = 15)
ht1+ht2
dev.off()

#------Fig3_b_c
# Chord diagram
pathways.show <- c("NOTCH") 
pdf('Fig3_b_c_Signalling_NOTCH_circle.pdf', width = 5, height = 5)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
netAnalysis_contribution(cellchat, signaling = pathways.show)
dev.off()

#------Fig3_d
pdf('Fig3_e_vlnplots_NOTCH_expression.pdf', width = 3, height = 3)
plotGeneExpression(cellchat, signaling = "NOTCH", enriched.only = TRUE, type = "violin")
dev.off()


#------Fig3_e
pdf('Fig3_e_Bubbleplot_Pericytes_BECs_AOM.pdf', width = 4, height = 4)
netVisual_bubble(cellchat, sources.use = 'Pericytes', targets.use = c('Pericytes','BECs','LECs','Proliferating'), signaling = c("NOTCH"), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use =  'BECs', targets.use = c('Pericytes','BECs','LECs','Proliferating'), signaling = c("NOTCH"), remove.isolate = FALSE)
dev.off()

#------Fig3_f

pairLR.notch <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
pdf('Fig3_f_chordPlots.pdf', width = 5, height = 5)
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = pairLR.notch[3,], layout = "chord")
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = pairLR.notch[7,], layout = "chord")
dev.off()


