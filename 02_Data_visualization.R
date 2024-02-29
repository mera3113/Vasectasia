### Written by Minjun Kim, M.S., Ph.D student in Dept. of Human Genetics @ McGill Univ.
### minjun.kim@mail.mcgill.ca for any inquiries

## Figures For Publcation
##################### 01. Environment Setup ##############################################
####### Desktop ######
wd = "E:/Lab data/Projects_McGill/Janus_Rak_GBM/Janus_EGFR_scRNA/" 
setwd(wd)
Janus_dir3 = paste0(getwd(), "/10x_scRNA_RAK1434/")
#save.image(paste0(Janus_dir3,"/.RData"))
load(paste0(Janus_dir3,"/.RData"))
library(Seurat); library(scater); library(ensembldb); library(DropletUtils); library(biomaRt)
library(tidyverse); library(dbscan); library(scran); library(uwot); library(dittoSeq)

#######################################################################

seurat_EGFR = readRDS(file.path(Janus_dir3, "seurat_EGFR.rds")) # Whole data object w/o filtering
seurat_EGFR_human = readRDS(file.path(Janus_dir3, "seurat_EGFR_human.rds")) # human MSCs
seurat_EGFR_mouse_merged = readRDS(file.path(Janus_dir3, "seurat_EGFR_mouse_fil.rds")) # total murine stroma
seurat_EGFR_merge1 = readRDS(file.path(Janus_dir3, "seurat_EGFR_mouse_fil.rds")) # Endo + Peri
seurat_EGFR_Endo = readRDS(file.path(Janus_dir3, "seurat_EGFR_Endo.rds")) # Endo
seurat_EGFR_Fibro = readRDS(file.path(Janus_dir3, "seurat_EGFR_Fibro.rds")) # Peri
seurat_EGFR_Myeloid = readRDS(file.path(Janus_dir3, "seurat_EGFR_Myeloid.rds")) # Myeloid
seurat_EGFR_Brain = readRDS(file.path(Janus_dir3, "seurat_EGFR_Brain.rds")) # OPC + OLG
seurat_EGFR_OPC = readRDS(file.path(Janus_dir3, "seurat_EGFR_OPC.rds")) # OPC

saveRDS(seurat_EGFR_human, file.path(Janus_dir3, "seurat_EGFR_human.rds"))
saveRDS(seurat_EGFR_mouse_merged, file.path(Janus_dir3, "seurat_EGFR_mouse_merged.rds"))
saveRDS(seurat_EGFR_merge1, file.path(Janus_dir3, "seurat_EGFR_merge1.rds"))
saveRDS(seurat_EGFR_Endo, file.path(Janus_dir3, "seurat_EGFR_Endo.rds"))
saveRDS(seurat_EGFR_Fibro, file.path(Janus_dir3, "seurat_EGFR_Fibro.rds"))
saveRDS(seurat_EGFR_Myeloid, file.path(Janus_dir3, "seurat_EGFR_Myeloid.rds"))
saveRDS(seurat_EGFR_Brain, file.path(Janus_dir3, "seurat_EGFR_Brain.rds"))
saveRDS(seurat_EGFR_OPC, file.path(Janus_dir3, "seurat_EGFR_OPC.rds"))

rm(seurat_EGFR)
rm(seurat_EGFR_human)
rm(seurat_EGFR_mouse_merged)
rm(seurat_EGFR_merge1)
rm(seurat_EGFR_Endo)
rm(seurat_EGFR_Fibro)
rm(seurat_EGFR_Myeloid)
rm(seurat_EGFR_Brain)
rm(seurat_EGFR_OPC)

##

EGFR_color = list(EGFR_status = c(WT = '#FF0000', KO = '#1A12FF'))
Idents(seurat_EGFR_mouse_merged) = seurat_EGFR_mouse_merged$Cell_type2
seurat_EGFR_mouse_merged = RenameIdents(seurat_EGFR_mouse_merged, 'Fibroblast' = 'Pericyte')
seurat_EGFR_mouse_merged$Cell_type2 =Idents(seurat_EGFR_mouse_merged)

seurat_EGFR_mouse_merged$EGFR_status = factor(seurat_EGFR_mouse_merged$EGFR_status, levels = c('WT', 'KO'))
seurat_EGFR_mouse_merged$Cell_type3 = factor(seurat_EGFR_mouse_merged$Cell_type3, 
                                             levels = c('Microglia', 'Microglia_Proliferating', 'Microglia_SPP1+', 'Mono/Mac_1',
                                                        'Mono/Mac_2', 'Mono/Mac_3', 'OLG', 'OPC_1', 'OPC_2', 'Pericyte_1', 'Pericyte_2',
                                                        'Endo_Angiogenic', 'Endo_Migrating', 'Endo_Permeable', 'Endo_Proliferative'))

# 800 * 600
DimPlot(seurat_EGFR_mouse_merged, reduction = "tsne", pt.size=2, label = T, group.by = 'Cell_type2', repel = T)
# # 700 * 600 
DimPlot(seurat_EGFR_mouse_merged, reduction = "tsne", pt.size=2, label = T, group.by = 'EGFR_status',cols = c('#FF0000', '#1A12FF'), shuffle = T, seed = '1234')
# 1400 * 600
DimPlot(seurat_EGFR_mouse_merged, reduction = "tsne", pt.size=2, label = T, group.by = 'Cell_type3', split.by = 'EGFR_status')
DimPlot(seurat_EGFR_mouse_merged, reduction = "tsne", pt.size=2, label = T, group.by = 'Cell_type3', repel = T) + ggtitle('Cell types')
DimPlot(seurat_EGFR_mouse_merged[,seurat_EGFR_mouse_merged$EGFR_status =='WT'], reduction = "tsne", pt.size=2, label = T, group.by = 'Cell_type3', repel = T) + ggtitle('WT')
DimPlot(seurat_EGFR_mouse_merged[,seurat_EGFR_mouse_merged$EGFR_status =='KO'], reduction = "tsne", pt.size=2, label = T, group.by = 'Cell_type3', repel = T)+ ggtitle('KO')

Idents(seurat_EGFR_mouse_merged) = seurat_EGFR_mouse_merged$Cell_type3
seurat_EGFR_mouse_merged.markers <- FindAllMarkers(seurat_EGFR_mouse_merged, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25, features = rownames(seurat_EGFR_mouse_merged)[!grepl('^MT-|^mt-|^RPS|^RPL|^Rps|^Rpl', rownames(seurat_EGFR_mouse_merged))])
seurat_EGFR_mouse_merged.markers$selectivity = seurat_EGFR_mouse_merged.markers$pct.1 - seurat_EGFR_mouse_merged.markers$pct.2
seurat_EGFR_mouse_merged.markers[seurat_EGFR_mouse_merged.markers$cluster=='6',]
top5 = seurat_EGFR_mouse_merged.markers %>% group_by(cluster) %>% top_n(n=10, wt = selectivity)

# 1000 * 1000
dittoHeatmap(seurat_EGFR_mouse_merged, main = 'All cells', annot.colors = c(dittoColors()), 
             c(Janus_genesets2[[1]], Janus_genesets2[[2]], Janus_genesets2[[3]], Janus_genesets2[[7]], Janus_genesets2[[8]], Janus_genesets2[[5]], Janus_genesets2[[6]]),  
             scaled.to.max = T, complex = T, use_raster = T, cluster_cols =F, cluster_rows =T, slot = 'data',
             annot.by = c('Cell_type3', 'EGFR_status'), order.by = c('Cell_type3'), 
             annotation_colors = EGFR_color,
             heatmap.colors.max.scaled = colorRampPalette(c('blue', 'white', 'red'))(50),  
             heatmap.colors = colorRampPalette(c('blue', 'white', 'red'))(50))

# 200 * 400
dittoPlot(seurat_EGFR_mouse_merged, main = '', plots = c('vlnplot'), c("Myeloid_comm_sig"), group.by = 'Cell_type2', xlab = '', ylab = 'Myeloid_signature',legend.show = T)
dittoPlot(seurat_EGFR_mouse_merged, main = '', plots = c('vlnplot'), c("Fibroblast_sig"), group.by = 'Cell_type2', xlab = '', ylab = 'Pericyte_signature',legend.show = F)
dittoPlot(seurat_EGFR_mouse_merged, main = '', plots = c('vlnplot'), c("Endothelial_sig"), group.by = 'Cell_type2', xlab = '', ylab = 'Endothelial_signature',legend.show = F)
dittoPlot(seurat_EGFR_mouse_merged, main = '', plots = c('vlnplot'), c("Oligodendrocyte_sig"), group.by = 'Cell_type2', xlab = '', ylab = 'OPC_signature',legend.show = F)
# 300 * 400
dittoPlot(seurat_EGFR_mouse_merged, main = '', plots = c('vlnplot'), c("Myelin_sig"), group.by = 'Cell_type2', xlab ='', ylab = 'Oigodendrocyte_signature', legend.show = F)


#FeaturePlot(seurat_EGFR_mouse_merged, features = c("Myeloid_comm_sig", "Fibroblast_sig", "Endothelial_sig", 
#                                                'Oligodendrocyte_sig', 'Myelin_sig'), min.cutoff = "q9", max.cutoff = 3, reduction = 'tsne', pt.size = 2, ncol = 5 )
#dittoPlot(seurat_EGFR_mouse_merged,
#          plots = c('jitter', 'vlnplot'), "Myeloid_comm_sig",
#          #adjustment = 'z-score',
#          group.by = 'Cell_type2')
#dittoPlot(seurat_EGFR_mouse_merged[,seurat_EGFR_mouse_merged$Species_Doublets == 'Mouse'],
#          plots = c('jitter', 'ridgeplot'), "EGFR_TPM",
#          adjustment = 'z-score',
#          group.by = 'Cell_type2', split.by = 'EGFR_status')
######################################################################
## Vascular structure
# 650 * 600
DimPlot(seurat_EGFR_merge1, reduction = "tsne", pt.size=2, label = T, group.by = 'EGFR_status',cols = c('#FF0000', '#1A12FF'), shuffle = T, seed = '1234')
# 1400 * 600
DimPlot(seurat_EGFR_merge1, reduction = "tsne", pt.size=2, label = T, group.by = 'Cell_type3' , split.by = 'EGFR_status')

# 600 * 800
DotPlot(seurat_EGFR_Endo, features = Endothelial_genes, dot.scale = 8, cluster.idents = F)+ 
    coord_flip() +
    RotatedAxis() +  
    scale_colour_gradient2(low = "#54A2FF", mid = "#FFFEC3", high = "Tomato", midpoint = 0)
# 500 * 800
DotPlot(seurat_EGFR_Fibro, features = Fibroblast_genes, dot.scale = 8, cluster.idents = F)+ 
    coord_flip() +
    RotatedAxis() +  
    scale_colour_gradient2(low = "#54A2FF", mid = "#FFFEC3", high = "Tomato", midpoint = 0)

# 400 * 400
dittoBarPlot(seurat_EGFR_Endo, 'EGFR_status', group.by = 'Cell_type3', main = '', color.panel = c('#1A12FF', '#FF0000'), colors = c(1,2))
# 250 * 350
dittoBarPlot(seurat_EGFR_Fibro, 'EGFR_status', group.by = 'Cell_type3', main = '', color.panel = c('#1A12FF', '#FF0000'), colors = c(1,2), )

#######################
## Myeloid
# 650 * 600
DimPlot(seurat_EGFR_Myeloid, reduction = "umap", pt.size=2, label = T, group.by = 'EGFR_status',cols = c('#FF0000', '#1A12FF'), shuffle = T, seed = '1234')
# 1400 * 600
DimPlot(seurat_EGFR_Myeloid, reduction = "umap", pt.size=2, label = T, group.by = 'Cell_type3' , split.by = 'EGFR_status')

# 800 * 600
DimPlot(seurat_EGFR_Myeloid, reduction = "umap", pt.size=2, label = T, group.by = 'Cell_type3', repel = F)
# 650 * 600
DimPlot(seurat_EGFR_Myeloid, reduction = "umap", pt.size=2, label = T, group.by = 'EGFR_status',cols = c('#1A12FF', '#FF0000'), shuffle = T, seed = '1234', repel = T)
# 650 * 600
DimPlot(seurat_EGFR_Myeloid[,seurat_EGFR_Myeloid$EGFR_status =='WT'], reduction = "umap", pt.size=2, label = T, group.by = 'Cell_type3', repel = T) + ggtitle('WT')
DimPlot(seurat_EGFR_Myeloid[,seurat_EGFR_Myeloid$EGFR_status =='KO'], reduction = "umap", pt.size=2, label = T, group.by = 'Cell_type3', repel = T)+ ggtitle('KO')


# 600 * 800
DotPlot(seurat_EGFR_Myeloid, features = Myeloid_genes, dot.scale = 8, cluster.idents = F)+ 
    coord_flip() +
    RotatedAxis() +  
    scale_colour_gradient2(low = "#54A2FF", mid = "#FFFEC3", high = "Tomato", midpoint = 0.5)

# 450 * 400
seurat_EGFR_Myeloid$EGFR_status = factor(seurat_EGFR_Myeloid$EGFR_status, levels = c('WT', 'KO'))
dittoBarPlot(seurat_EGFR_Myeloid, 'EGFR_status', group.by = 'Cell_type3', main = '', color.panel = c('#1A12FF', '#FF0000'), colors = c(1,2))


seurat_EGFR_Brain$EGFR_status = factor(seurat_EGFR_Brain$EGFR_status, levels = c('WT', 'KO'))
## Brain
# 700 * 600
DimPlot(seurat_EGFR_Brain, reduction = "umap", pt.size=2, label = T, group.by = 'Cell_type3', repel = T)
# 650 * 600
DimPlot(seurat_EGFR_Brain, reduction = "umap", pt.size=2, label = T, group.by = 'EGFR_status',cols = c('#FF0000', '#1A12FF'), shuffle = T, seed = '1234')

# 650 * 600
DimPlot(seurat_EGFR_Brain[,seurat_EGFR_Brain$EGFR_status =='WT'], reduction = "umap", pt.size=2, label = T, group.by = 'Cell_type3', repel = T) + ggtitle('WT')
DimPlot(seurat_EGFR_Brain[,seurat_EGFR_Brain$EGFR_status =='KO'], reduction = "umap", pt.size=2, label = T, group.by = 'Cell_type3', repel = T)+ ggtitle('KO')


# 1400 * 600
DimPlot(seurat_EGFR_Brain, reduction = "umap", pt.size=2, label = T, group.by = 'Cell_type3' , split.by = 'EGFR_status')
seurat_EGFR_Brain$Cell_type3 = factor(seurat_EGFR_Brain$Cell_type3, levels = c('OPC_1', 'OPC_2', 'OLG'))
Idents(seurat_EGFR_Brain) = seurat_EGFR_Brain$Cell_type3
# 550 * 800
DotPlot(seurat_EGFR_Brain, 
        features = c(Brain_genes, seurat_EGFR_OPC.markers[seurat_EGFR_OPC.markers$cluster=='1',]$gene[1:8]), , dot.scale = 8, cluster.idents = F)+ 
    coord_flip() +
    RotatedAxis() +  
    scale_colour_gradient2(low = "#54A2FF", mid = "#FFFEC3", high = "Tomato", midpoint = 0)

# 350 * 400
dittoBarPlot(seurat_EGFR_Brain, 'EGFR_status', group.by = 'Cell_type3', main = '', color.panel = c('#1A12FF', '#FF0000'), colors = c(1,2))

dittoPlot(seurat_EGFR_Brain,
          plots = c('jitter', 'ridgeplot'), var =  c("C4b"),
                                                     #'Kcnj10', 'Fth1', 'Ptprz', 'Fabp7', 'Ppp1r14b', 'Rtn1', 'Marcks', 'Basp1', 'Tubb2b', 'Gap43'),
          adjustment = 'z-score',
          group.by = 'Cell_type3')
dittoPlot(seurat_EGFR_Brain,
          plots = c('jitter', 'ridgeplot'), var =  c("Fth1"),
          #'Kcnj10', 'Fth1', 'Ptprz', 'Fabp7', 'Ppp1r14b', 'Rtn1', 'Marcks', 'Basp1', 'Tubb2b', 'Gap43'),
          adjustment = 'z-score',
          group.by = 'Cell_type3')
dittoPlot(seurat_EGFR_Brain,
          plots = c('jitter', 'ridgeplot'), var =  c("Ppp1r14b"),
          #'Kcnj10', 'Fth1', 'Ptprz', 'Fabp7', 'Ppp1r14b', 'Rtn1', 'Marcks', 'Basp1', 'Tubb2b', 'Gap43'),
          adjustment = 'z-score',
          group.by = 'Cell_type3')
dittoPlot(seurat_EGFR_Brain,
          plots = c('jitter', 'ridgeplot'), var =  c("Marcks"),
          #'Kcnj10', 'Fth1', 'Ptprz', 'Fabp7', 'Ppp1r14b', 'Rtn1', 'Marcks', 'Basp1', 'Tubb2b', 'Gap43'),
          adjustment = 'z-score',
          group.by = 'Cell_type3')

DimPlot(seurat_EGFR_OPC, reduction = "tsne", pt.size=2, label = T, group.by = 'Cell_type3', repel = T)


DotPlot(seurat_EGFR_OPC, features = unique(top5$gene), dot.scale = 8, cluster.idents = F)+ 
    coord_flip() +
    RotatedAxis() +  
    scale_colour_gradient2(low = "#54A2FF", mid = "#FFFEC3", high = "Tomato", midpoint = 0)


#####

# 650 * 500
FeatureScatter(seurat_EGFR[,seurat_EGFR$sample.name=='WT'], feature1 = 'frac_hg_genes', feature2 = 'WT_sig',  pt.size = 1.5, group.by = 'Species_Doublets')
FeatureScatter(seurat_EGFR[,seurat_EGFR$sample.name=='KO'], feature1 = 'frac_hg_genes', feature2 = 'KO_sig',  pt.size = 1.5, group.by = 'Species_Doublets')

# 800 * 1200
dittoHeatmap(seurat_EGFR, main = 'Expression of DEGs for human MSCs in murine stroma', annot.colors = c(dittoColors()), 
             c(WT_markers[WT_markers$pct.1>0.9,'gene'][1:50],
               KO_markers[KO_markers$pct.1>0.9,'gene'][1:50]),  
             scaled.to.max = T, complex = T, use_raster = T, cluster_cols =T, cluster_rows =T, slot = 'data',
             annot.by = c('frac_hg_genes', 'sample.name'), #order.by = c('Cell_type3'), 
             annotation_colors = EGFR_color,
             heatmap.colors.max.scaled = colorRampPalette(c('blue', 'white', 'red'))(50),  
             heatmap.colors = colorRampPalette(c('blue', 'white', 'red'))(50))

#
seurat_EGFR_merge1$EGFR_TPM = seurat_EGFR[,colnames(seurat_EGFR_merge1)]$EGFR_TPM
seurat_EGFR_merge1$EGFR = seurat_EGFR[,colnames(seurat_EGFR_merge1)]$EGFR
seurat_EGFR_merge1$Species_Doublets = seurat_EGFR[,colnames(seurat_EGFR_merge1)]$Species_Doublets

hist(seurat_EGFR_merge1$EGFR_TPM, breaks =100)

# 2000 * 1000
FeaturePlot(object = seurat_EGFR_merge1, reduction = 'tsne', features = c("EGFR_TPM", 'Pecam1'), min.cutoff = 0, max.cutoff = 2,  pt.size = 2, split.by = 'EGFR_status', blend = T, blend.threshold = 0.5 )
FeaturePlot(object = seurat_EGFR_merge1, reduction = 'tsne', features = c("EGFR_TPM", 'Cd34'), min.cutoff = 0, max.cutoff = 2,  pt.size = 2, split.by = 'EGFR_status', blend = T, blend.threshold = 0.5 )


FeaturePlot(object = seurat_EGFR_merge1, reduction = 'tsne', features = c('Cd34'), min.cutoff = 0, max.cutoff = 2,  pt.size = 2, blend = F)

FeaturePlot(object = seurat_EGFR_Endo, reduction = 'umap', features = c('EGFR', "Pecam1"),  min.cutoff = 0.2, max.cutoff = 1, pt.size = 1.5, blend = T)



####### Human
# 700 * 600
DimPlot(seurat_EGFR_human, reduction = "umap", pt.size=2, label = T, group.by = 'seurat_clusters', repel = F)
# 650 * 600
DimPlot(seurat_EGFR_human, reduction = "umap", pt.size=2, label = T, group.by = 'EGFR_status',cols = c('#FF0000', '#1A12FF'), shuffle = T, seed = '1234')
# 650 * 600
DimPlot(seurat_EGFR_human[,seurat_EGFR_human$EGFR_status =='WT'], reduction = "umap", pt.size=2, label = T, group.by = 'Cell_type2', repel = T) + ggtitle('WT')
DimPlot(seurat_EGFR_human[,seurat_EGFR_human$EGFR_status =='KO'], reduction = "umap", pt.size=2, label = T, group.by = 'Cell_type2', repel = T)+ ggtitle('KO')

seurat_EGFR_human$Cell_type2 = factor(seurat_EGFR_human$Cell_type2, levels = c('WT_1', 'WT_2', 'WT_3', 'WT_4', 'KO_1', 'KO_2', 'KO_3', 'KO_4', 'KO_5'))
seurat_EGFR_human.markers[seurat_EGFR_human.markers$cluster=='',]
top5 = seurat_EGFR_human.markers %>% group_by(cluster) %>% top_n(n=5, wt =selectivity)
top5[grep('0|4|6|7',top5$cluster),]$gene 
top5[grep('2|3|5|8',top5$cluster),]$gene
top5[grep('1',top5$cluster),]$gene
# 800 * 800
dittoHeatmap(seurat_EGFR_human, main = 'Expression of DEGs for human MSCs', annot.colors = c(dittoColors()), 
             unique(c(WT_markers[WT_markers$pct.1>0.9,'gene'][1:5], top5[grep('0|4|6|7',top5$cluster),]$gene, 
               KO_markers[KO_markers$pct.1>0.9,'gene'][1:5], top5[grep('2|3|5|8',top5$cluster),]$gene,
               top5[grep('1',top5$cluster),]$gene)),  
             scaled.to.max = T, complex = T, use_raster = T, cluster_cols =F, cluster_rows =T, slot = 'data',
             annot.by = c('Cell_type2', 'EGFR_status'), order.by = c('EGFR_status', 'Cell_type2'), 
             annotation_colors = EGFR_color,
             heatmap.colors.max.scaled = colorRampPalette(c('blue', 'white', 'red'))(50),  
             heatmap.colors = colorRampPalette(c('blue', 'white', 'red'))(50))
