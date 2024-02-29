### Written by Minjun Kim, M.S., Ph.D student in Dept. of Human Genetics @ McGill Univ.
### minjun.kim@mail.mcgill.ca for any inquiries

##################### 01. Environment Setup ##############################################
####### Desktop ######
wd = "E:/Lab data/Projects_McGill/Janus_Rak_GBM/Janus_EGFR_scRNA/" 
setwd(wd)
Janus_dir3 = paste0(getwd(), "/10x_scRNA_RAK1434/")
save.image(paste0(Janus_dir3,"/.RData"))
load(paste0(Janus_dir3,"/.RData"))
#######################################################################
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))
library(Seurat); library(scater); library(ensembldb); library(DropletUtils); library(biomaRt)
library(tidyverse); library(dbscan); library(scran); library(uwot); library(dittoSeq)

######################################################## 
## Load data
Janus_raw3 = paste0(Janus_dir3, list.files(Janus_dir3), '/filtered_feature_bc_matrix/')
## For Seurat
EGFR_KO_m = Read10X(data.dir = Janus_raw3[2])
EGFR_WT_m = Read10X(data.dir = Janus_raw3[4])
## For SCE
EGFR_KO_m <- read10xCounts(Janus_raw3[2])
EGFR_WT_m <- read10xCounts(Janus_raw3[4])

## Feature process by species
head(EGFR_KO_m)
rowData(EGFR_KO_m) <- rowData(EGFR_KO_m) %>%
    as_tibble() %>%
    separate(ID, c("genome", "ID"), "_") %>%
    separate(Symbol, c("genome2", "Symbol"), "_") %>%
    #select( -NA., -genome2) %>%
    remove_rownames()   
table(as.tibble(rowData(EGFR_KO_m))$genome2 == as.tibble(rowData(EGFR_KO_m))$genome)
rowData(EGFR_KO_m)$genome2 = NULL
rownames(EGFR_KO_m) <- rowData(EGFR_KO_m)$ID

species_counts <- 
    assay(EGFR_KO_m) %>%
    as.matrix()  %>%
    as_tibble()  %>%
    group_by(genome=rowData(EGFR_KO_m)$genome) %>%
    summarize_all(sum) 

species_counts <-  
    species_counts %>%
    select(-genome) %>%
    t() %>%
    as_tibble() %>%
    setNames(species_counts$genome) %>%
    mutate(frac_hg_genes = hg19/(hg19+mm10))

ggplot(species_counts)+
    geom_point( aes(hg19,
                    mm10,
                    colour=frac_hg_genes))

colData(EGFR_KO_m) <- cbind(colData(EGFR_KO_m), species_counts) 

EGFR_KO_m_mtx <- 
    t(assay(EGFR_KO_m) %>%
          as.matrix)

## Make duplicated genes unique
rowData(EGFR_KO_m)[32739:60736,]$Symbol[rowData(EGFR_KO_m)[32739:60736,]$Symbol %in% rowData(EGFR_KO_m)[1:32738,]$Symbol] = paste0(rowData(EGFR_KO_m)[32739:60736,]$Symbol[rowData(EGFR_KO_m)[32739:60736,]$Symbol %in% rowData(EGFR_KO_m)[1:32738,]$Symbol], '_mouse')
grep('mt-', rowData(EGFR_KO_m)[32738:60736,]$Symbol, value = T)


## Feature process by species
rowData(EGFR_WT_m) <- rowData(EGFR_WT_m) %>%
    as_tibble() %>%
    separate(ID, c("genome", "ID"), "_") %>%
    separate(Symbol, c("genome2", "Symbol"), "_") %>%
    #select( -NA., -genome2) %>%
    remove_rownames()   
table(as.tibble(rowData(EGFR_WT_m))$genome2 == as.tibble(rowData(EGFR_WT_m))$genome)
rowData(EGFR_WT_m)$genome2 = NULL
rownames(EGFR_WT_m) <- rowData(EGFR_WT_m)$ID

species_counts <- 
    assay(EGFR_WT_m) %>%
    as.matrix()  %>%
    as_tibble()  %>%
    group_by(genome=rowData(EGFR_WT_m)$genome) %>%
    summarize_all(sum) 

species_counts <-  
    species_counts %>%
    select(-genome) %>%
    t() %>%
    as_tibble() %>%
    setNames(species_counts$genome) %>%
    mutate(frac_hg_genes = hg19/(hg19+mm10))

species_counts 

ggplot(species_counts)+
    geom_point( aes(hg19,
                    mm10,
                    colour=frac_hg_genes))

colData(EGFR_WT_m) <- cbind(colData(EGFR_WT_m), species_counts) 

table(rowData(EGFR_WT_m)[1:32738,]$Symbol %in% rowData(EGFR_WT_m)[32739:60736,]$Symbol) ## 19 genes are relapsed in Symbol
rowData(EGFR_WT_m)[32739:60736,]$Symbol[rowData(EGFR_WT_m)[32739:60736,]$Symbol %in% rowData(EGFR_WT_m)[1:32738,]$Symbol] = paste0(rowData(EGFR_WT_m)[32739:60736,]$Symbol[rowData(EGFR_WT_m)[32739:60736,]$Symbol %in% rowData(EGFR_WT_m)[1:32738,]$Symbol], '_mouse')


##### Merge data
length(rowData(EGFR_KO_m)$Symbol); length(unique(rowData(EGFR_KO_m)$Symbol)) # some duplications, but its fine to be rownames
length(rowData(EGFR_WT_m)$Symbol); length(unique(rowData(EGFR_WT_m)$Symbol))
summary(factor(rowData(EGFR_WT_m)$Symbol))

rownames(EGFR_KO_m) = rowData(EGFR_KO_m)$Symbol
rownames(EGFR_WT_m) = rowData(EGFR_WT_m)$Symbol
EGFR_KO_m = sceToSeurat(EGFR_KO_m)
EGFR_WT_m = sceToSeurat(EGFR_WT_m)

hist(breaks=100, EGFR_KO_m_hg19$frac_hg_genes)
hist(breaks=100, EGFR_WT_m_hg19$frac_hg_genes)

seurat_EGFR = merge(EGFR_KO_m, y=c(EGFR_WT_m), add.cell.ids=c('KO', 'WT'), project='EGFR', merge.data=TRUE)
seurat_EGFR$sample.name = ifelse(grepl('KO', colnames(seurat_EGFR)), 'KO', 'WT')
hist(seurat_EGFR$nCount_RNA[seurat_EGFR$nCount_RNA<100000], breaks=100)
hist(seurat_EGFR$nFeature_RNA[seurat_EGFR$nFeature_RNA<15000], breaks=100)
hist(seurat_EGFR$percent.mt, breaks=100)
hist(seurat_EGFR$percent.mt_mouse, breaks=100)
FeatureScatter(seurat_EGFR[,seurat_EGFR$hg19 < 5000 & seurat_EGFR$mm10 <1000], feature1 ='hg19', feature2 = 'mm10' )
plot2 = FeatureScatter(seurat_EGFR[,!grepl('Human', seurat_EGFR$Cell_type1)], feature1 = 'hg19', feature2 = 'mm10')

select.cells = CellSelector(plot=plot2)
Idents(seurat_EGFR, cells=select.cells) = 'Mouse'
Idents(seurat_EGFR, cells=select.cells) = 'Human'
seurat_EGFR$Cell_type1 = Idents(seurat_EGFR) # Total 29764 cells before filtering

plot2 = FeatureScatter(seurat_EGFR[,grepl('Mouse', seurat_EGFR$Cell_type1) & seurat_EGFR$mm10 < 2000], feature1 = 'mm10', feature2 = 'percent.mt_mouse')
plot2 = FeatureScatter(seurat_EGFR[,grepl('Mouse', seurat_EGFR$Cell_type1) & seurat_EGFR$mm10 < 4000 & seurat_EGFR$percent.mt_mouse <0.4], feature1 = 'mm10', feature2 = 'percent.mt_mouse')
plot2 = FeatureScatter(seurat_EGFR[,grepl('Human', seurat_EGFR$Cell_type1) & seurat_EGFR$hg19 < 12000 & seurat_EGFR$percent.mt <0.2], feature1 = 'hg19', feature2 = 'percent.mt')

select.cells = CellSelector(plot=plot2)
Idents(seurat_EGFR, cells=select.cells) = 'Mouse_bad'
Idents(seurat_EGFR, cells=select.cells) = 'Human_bad'
seurat_EGFR$Cell_type1 = Idents(seurat_EGFR) # Total 29764 cells before filtering

FeatureScatter(seurat_EGFR[,grepl('Mouse', seurat_EGFR$Cell_type1) & seurat_EGFR$mm10 < 8000 & seurat_EGFR$percent.mt_mouse <0.4], feature1 = 'mm10', feature2 = 'percent.mt_mouse')
FeatureScatter(seurat_EGFR[,grepl('Human', seurat_EGFR$Cell_type1) & seurat_EGFR$hg19 < 80000], feature1 = 'hg19', feature2 = 'percent.mt')
table(seurat_EGFR$Cell_type1) # Human_bad = 761, Mouse_bad = 16757, Human = 10234, Mouse = 2012 cells

seurat_EGFR_human <- subset(seurat_EGFR, subset = Cell_type1 == 'Human') # 10995 human cells before filtering
seurat_EGFR_mouse <- subset(seurat_EGFR, subset = Cell_type1 == 'Mouse') # 2247 mouse cells before filtering
saveRDS(seurat_EGFR, file.path(Janus_dir3, "seurat_EGFR.rds")) ## Object 1. All cells with both hu & mu genes, unfiltered
rm(seurat_EGFR)
summary(factor(colnames(seurat_EGFR_mouse)))
seurat_EGFR_human = ProcessSeurat_pdx2(seurat_EGFR_human, scale = 2, 0.5, 1)
seurat_EGFR_mouse = ProcessSeurat_pdx2(seurat_EGFR_mouse, 1, 1, 2)
saveRDS(seurat_EGFR_human, file.path(Janus_dir3, "seurat_EGFR_human.rds")) ## Object 2. All human cells with only human genes
saveRDS(seurat_EGFR_mouse, file.path(Janus_dir3, "seurat_EGFR_mouse.rds")) ## Obejct 3. All murine cells with only murine genes
seurat_EGFR_human = readRDS(file.path(Janus_dir3, "seurat_EGFR_human.rds"))
seurat_EGFR_mouse = readRDS(file.path(Janus_dir3, "seurat_EGFR_mouse.rds"))

summary(seurat_EGFR_human$hg19)
summary(seurat_EGFR_human$percent.mt)
summary(seurat_EGFR_mouse$mm10)
summary(seurat_EGFR_mouse$percent.mt_mouse)

summary(seurat_EGFR_mouse$mm10)
min(seurat_EGFR_human$hg19)
min(seurat_EGFR_human$hg19)

###################################

seurat_EGFR_human.markers <- FindAllMarkers(seurat_EGFR_human, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25, features = rownames(seurat_EGFR_human)[!grepl('^MT-|^mt-|^RPS|^RPL|^Rps|^Rpl', rownames(seurat_EGFR_human))])
seurat_EGFR_human.markers$selectivity = seurat_EGFR_human.markers$pct.1 - seurat_EGFR_human.markers$pct.2
seurat_EGFR_human.markers[seurat_EGFR_human.markers$cluster=='',]
top5 = seurat_EGFR_human.markers %>% group_by(cluster) %>% top_n(n=10, wt = p_val_adj)
write.csv(top5, 'GBM_top20_selectivity.csv')
write.csv(top5, 'GBM_top20_log2FC.csv')
write.csv(top5, 'GBM_top20_adj-pvalue.csv')


seurat_EGFR_human$EGFR_status = factor(seurat_EGFR_human$sample.name, levels = c('WT', 'KO'))
seurat_EGFR_human = RenameIdents(seurat_EGFR_human, '0' = 'WT_1', '4' = 'WT_2', '6' = 'WT_3', '7' = 'WT_4', '2' = 'KO_1', '3' = 'KO_2', '5' = 'KO_3', '8' = 'KO_4', '1' = 'KO_5')
seurat_EGFR_human$Cell_type2 =Idents(seurat_EGFR_human)

saveRDS(seurat_EGFR_human, file.path(Janus_dir3, "seurat_EGFR_human.rds")) ## Object 2.1 annotated human cell object
seurat_EGFR_human = readRDS(file.path(Janus_dir3, "seurat_EGFR_human.rds"))


###################

seurat_EGFR_mouse.markers <- FindAllMarkers(seurat_EGFR_mouse, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25, features = rownames(seurat_EGFR_mouse)[!grepl('^MT-|^mt-|^RPS|^RPL|^Rps|^Rpl', rownames(seurat_EGFR_mouse))])
seurat_EGFR_mouse.markers$selectivity = seurat_EGFR_mouse.markers$pct.1 - seurat_EGFR_mouse.markers$pct.2
seurat_EGFR_mouse.markers[seurat_EGFR_mouse.markers$cluster=='6',]

top5 = seurat_EGFR_mouse.markers %>% group_by(cluster) %>% top_n(n=20, wt = selectivity)

## top marker genes in each cell type in the murine stroma
Janus_genesets = list(
    Myeloid_markers = unique(c(top5[top5$cluster=='8',]$gene,
                               top5[top5$cluster=='0',]$gene,
                               top5[top5$cluster=='5',]$gene,
                               top5[top5$cluster=='7',]$gene,
                               top5[top5$cluster=='4',]$gene,
                               top5[top5$cluster=='11',]$gene,
                               top5[top5$cluster=='10',]$gene)),
    Myelin_markers = top5[top5$cluster=='6',]$gene,
    OligoDendro_markers = top5[top5$cluster=='3',]$gene,
    RBC_markers = top5[top5$cluster=='9',]$gene,
    Fibro_markers = top5[top5$cluster=='1',]$gene,
    Endo_markers = top5[top5$cluster=='2',]$gene)
## Curated top markers for each cell type
Janus_genesets2 = list(
    Myeloid_comm = c('Tyrobp', 'Ctss', 'Fcer1g'),
    Macrophage = c('Aif1', 'Hexb', 'C1qc', 'C1qb', 'C1qa', 'Trem2', 'Fcgr1'),
    Monocyte = c('Cd74', 'H2-Eb1', 'H2-Aa', 'H2-Ab1', 'Lsp1', 'Napsa', 'Cytip', 'Gpr171', 'Plbd1', 'Chil3', 'Ccr2', 'S100a4'),
    RBC = c("Hbb-bt", "Hba-a2", "Alas2", "Hba-a1","Hbb-bs"),
    Fibroblast = c("Col1a1", "Rgs5", "Vtn", "Col3a1", "Bgn", "Higd1b", "Cd248", "Pdgfrb", "Des", "Rarres2", "Fstl1", "Col1a2", "Acta2"),
    Endothelial = c("Cldn5", "Ramp2", "Egfl7", "Ptprb", "Adgrl4", "Ctla2a", "Tmem252", "Pecam1", "Scgb3a1", "Ly6c1", "Itm2a", "Kdr", "Cd34"),
    OligoDendro = c("Ptprz1","Cspg5", "Lhfpl3", "Fabp7",  "Scrg1", "Olig1", "Gap43", "Gng3", "Gpm6a", "Bex2", "Olig2"),
    Myelin = c("Klk6", "Ermn", "Mag", "Mog", "Nkx6-2", "Cldn11", "Stmn4", "Ugt8a", "Aplp1", "Tubb4a", "Plp1", "Serpina3n", "Fez1")
)
names(Janus_genesets2[1])
dittoHeatmap(seurat_EGFR_mouse, main = 'All cells', annot.colors = c(dittoColors()), 
             c(Janus_genesets2[[1]], Janus_genesets2[[7]], Janus_genesets2[[8]], Janus_genesets2[[5]], Janus_genesets2[[6]])[isGene(c(Janus_genesets2[[1]], Janus_genesets2[[7]], Janus_genesets2[[8]], Janus_genesets2[[5]], Janus_genesets2[[6]]),seurat_EGFR_mouse)],
             scaled.to.max = T, complex = T, use_raster = T, cluster_cols =T, cluster_rows =F, slot = 'data', 
             annot.by = c('Cell_type2'), #order.by = 'Cell_type4',
             heatmap.colors.max.scaled = colorRampPalette(c('blue', 'white', 'red'))(50),  
             heatmap.colors = colorRampPalette(c('blue', 'white', 'red'))(50))
dittoHeatmap(seurat_EGFR_mouse, main = 'All cells', annot.colors = c(dittoColors()), 
             c(Janus_genesets2[[1]], Janus_genesets2[[2]], Janus_genesets2[[3]], Janus_genesets2[[7]], Janus_genesets2[[8]], Janus_genesets2[[5]], Janus_genesets2[[6]]),
             scaled.to.max = T, complex = T, use_raster = T, cluster_cols =T, cluster_rows =F, slot = 'data', 
             annot.by = c('Cell_type2'), #order.by = 'Cell_type4',
             heatmap.colors.max.scaled = colorRampPalette(c('blue', 'white', 'red'))(50),  
             heatmap.colors = colorRampPalette(c('blue', 'white', 'red'))(50))

seurat_EGFR_mouse$Myeloid_comm_sig = log(((colMeans(seurat_EGFR_mouse[['RNA']]@counts[Janus_genesets2[[1]],])/colMeans(seurat_EGFR_mouse@assays$RNA@counts)))+1)
seurat_EGFR_mouse$Macrophage_sig = log(((colMeans(seurat_EGFR_mouse[['RNA']]@counts[Janus_genesets2[[2]],])/colMeans(seurat_EGFR_mouse@assays$RNA@counts)))+1)
seurat_EGFR_mouse$Monocyte_sig = log(((colMeans(seurat_EGFR_mouse[['RNA']]@counts[Janus_genesets2[[3]],])/colMeans(seurat_EGFR_mouse@assays$RNA@counts)))+1)
seurat_EGFR_mouse$RBC_sig = log(((colMeans(seurat_EGFR_mouse[['RNA']]@counts[Janus_genesets2[[4]],])/colMeans(seurat_EGFR_mouse@assays$RNA@counts)))+1)
seurat_EGFR_mouse$Fibroblast_sig = log(((colMeans(seurat_EGFR_mouse[['RNA']]@counts[Janus_genesets2[[5]],])/colMeans(seurat_EGFR_mouse@assays$RNA@counts)))+1)
seurat_EGFR_mouse$Endothelial_sig = log(((colMeans(seurat_EGFR_mouse[['RNA']]@counts[Janus_genesets2[[6]],])/colMeans(seurat_EGFR_mouse@assays$RNA@counts)))+1)
seurat_EGFR_mouse$Oligodendrocyte_sig = log(((colMeans(seurat_EGFR_mouse[['RNA']]@counts[Janus_genesets2[[7]],])/colMeans(seurat_EGFR_mouse@assays$RNA@counts)))+1)
seurat_EGFR_mouse$Myelin_sig = log(((colMeans(seurat_EGFR_mouse[['RNA']]@counts[Janus_genesets2[[8]],])/colMeans(seurat_EGFR_mouse@assays$RNA@counts)))+1)

############## Doublet filtering, filter out cells expressing two or more lineage-specific gene sets
Idents(seurat_EGFR_mouse) = seurat_EGFR_mouse$Cell_type2
plot2 = DimPlot(seurat_EGFR_mouse, reduction = "umap", pt.size=1.5, label = T, group.by = 'Cell_type2')
plot2 = FeatureScatter(seurat_EGFR_mouse[,seurat_EGFR_mouse$Cell_type2!='Doublets'], feature1 = 'Myeloid_comm', feature2 = 'RBC',  pt.size = 1.5)
plot2 = FeatureScatter(seurat_EGFR_mouse[,seurat_EGFR_mouse$Cell_type2!='Doublets'], feature1 = 'Myeloid_comm_sig', feature2 = 'Myelin_sig',  pt.size = 1.5)
plot2 = FeatureScatter(seurat_EGFR_mouse[,seurat_EGFR_mouse$Cell_type2!='Doublets'], feature1 = 'Myeloid_comm', feature2 = 'Oligodendrocyte',  pt.size = 1.5)
plot2 = FeatureScatter(seurat_EGFR_mouse[,seurat_EGFR_mouse$Cell_type2!='Doublets'], feature1 = 'Myeloid_comm_sig', feature2 = 'Endothelial_sig',  pt.size = 1.5)
plot2 = FeatureScatter(seurat_EGFR_mouse[,seurat_EGFR_mouse$Cell_type2!='Doublets'], feature1 = 'Myeloid_comm', feature2 = 'Fibroblast',  pt.size = 1.5)
plot2 = FeatureScatter(seurat_EGFR_mouse[,seurat_EGFR_mouse$Cell_type2!='Doublets'], feature1 = 'Fibroblast_sig', feature2 = 'Endothelial_sig',  pt.size = 1.5)
plot2 = FeatureScatter(seurat_EGFR_mouse[,seurat_EGFR_mouse$Cell_type2!='Doublets'], feature1 = 'Fibroblast', feature2 = 'Myelin',  pt.size = 1.5)
plot2 = FeatureScatter(seurat_EGFR_mouse[,seurat_EGFR_mouse$Cell_type2!='Doublets'], feature1 = 'Fibroblast', feature2 = 'Oligodendrocyte',  pt.size = 1.5)
plot2 = FeatureScatter(seurat_EGFR_mouse[,seurat_EGFR_mouse$Cell_type2!='Doublets'], feature1 = 'Fibroblast', feature2 = 'RBC',  pt.size = 1.5)
plot2 = FeatureScatter(seurat_EGFR_mouse[,seurat_EGFR_mouse$Cell_type2!='Doublets'], feature1 = 'Endothelial_sig', feature2 = 'Myelin',  pt.size = 1.5)
plot2 = FeatureScatter(seurat_EGFR_mouse[,seurat_EGFR_mouse$Cell_type2!='Doublets'], feature1 = 'Endothelial', feature2 = 'Oligodendrocyte',  pt.size = 1.5)
plot2 = FeatureScatter(seurat_EGFR_mouse[,seurat_EGFR_mouse$Cell_type2=='Doublets'], feature1 = 'Endothelial', feature2 = 'Oligodendrocyte',  pt.size = 1.5)
plot2 = FeatureScatter(seurat_EGFR_mouse[,seurat_EGFR_mouse$Cell_type2!='Doublets'], feature1 = 'Monocyte_sig', feature2 = 'Oligodendrocyte_sig',  pt.size = 1.5)
plot2 = FeatureScatter(seurat_EGFR_mouse[,seurat_EGFR_mouse$Cell_type2!='Doublets'], feature1 = 'Myelin_sig', feature2 = 'Oligodendrocyte_sig',  pt.size = 1.5)
plot2 = FeatureScatter(seurat_EGFR_mouse[,seurat_EGFR_mouse$Cell_type2!='Doublets'], feature1 = 'Lsp1', feature2 = 'Oligodendrocyte_sig',  pt.size = 1.5)
plot2 = FeatureScatter(seurat_EGFR_mouse[,seurat_EGFR_mouse$Cell_type2!='Doublets'], feature1 = 'Lsp1', feature2 = 'Fibroblast_sig',  pt.size = 1.5)
plot2 = FeatureScatter(seurat_EGFR_mouse[,seurat_EGFR_mouse$Cell_type2!='Doublets'], feature1 = 'Lsp1', feature2 = 'Myeloid_comm_sig',  pt.size = 1.5)
plot2 = FeatureScatter(seurat_EGFR_mouse[,seurat_EGFR_mouse$Cell_type2=='Unknown'], feature1 = 'Lsp1', feature2 = 'Myeloid_comm_sig',  pt.size = 1.5)
plot2 = FeatureScatter(seurat_EGFR_mouse[,seurat_EGFR_mouse$Cell_type2=='Unknown'], feature1 = 'Monocyte_sig', feature2 = 'Myelin_sig',  pt.size = 1.5)
plot2 = FeatureScatter(seurat_EGFR_mouse[,seurat_EGFR_mouse$Cell_type2=='Unknown'], feature1 = 'Monocyte_sig', feature2 = 'Myelin_sig',  pt.size = 1.5)
plot2 = FeatureScatter(seurat_EGFR_mouse[,seurat_EGFR_mouse$Cell_type2!='Doublets'], feature1 = 'Myelin_sig', feature2 = 'Oligodendrocyte_sig',  pt.size = 1.5)

select.cells = CellSelector(plot=plot2)
Idents(seurat_EGFR_mouse, cells=select.cells) = 'Doublets'
seurat_EGFR_mouse$Cell_type2 = Idents(seurat_EGFR_mouse)

table(seurat_EGFR_mouse$Cell_type1) 

seurat_EGFR_mouse_fil <- subset(seurat_EGFR_mouse, subset = Cell_type2 != 'Doublets') # 2247 mouse cells before filtering
saveRDS(seurat_EGFR_mouse, file.path(Janus_dir3, "seurat_EGFR_mouse.rds"))
####################################################################################
seurat_EGFR_mouse_fil = ProcessSeurat_pdx2(seurat_EGFR_mouse_fil, res = 1.5, scale = 1,species = 2)
seurat_EGFR_mouse_fil$EGFR_status = factor(seurat_EGFR_mouse_fil$sample.name, levels = c('WT', 'KO')) 
saveRDS(seurat_EGFR_mouse_fil, file.path(Janus_dir3, "seurat_EGFR_mouse_fil.rds")) ## Object 4. All murine cells after doublet filtering

################################# 
plot2 = DimPlot(seurat_EGFR_mouse_fil, reduction = "umap", pt.size=1.5, label = T)
plot2 = DimPlot(seurat_EGFR_mouse_fil, reduction = "umap", pt.size=1.5, label = T, group.by = 'Cell_type2')
plot2 = FeaturePlot(seurat_EGFR_mouse_fil, features = 'mm10', reduction = "umap", pt.size=1.5, label = T, max.cutoff = 1000)
plot2 = FeaturePlot(seurat_EGFR_mouse_fil, features = 'Myeloid_comm_sig', reduction = "umap", pt.size=1.5, label = T)
plot2 = FeaturePlot(seurat_EGFR_mouse_fil, features = 'Monocyte_sig', reduction = "umap", pt.size=1, label = T)
plot2 = FeaturePlot(seurat_EGFR_mouse_fil, features = 'Ptprz1', reduction = "umap", pt.size=1, label = T)
plot2 = FeaturePlot(seurat_EGFR_mouse_fil, features = 'Fibroblast_sig', reduction = "umap", pt.size=1, label = T)
plot2 = FeaturePlot(seurat_EGFR_mouse_fil, features = 'Endothelial_sig', reduction = "umap", pt.size=1, label = T)
plot2 = FeaturePlot(seurat_EGFR_mouse_fil, features = 'Klk6', reduction = "umap", pt.size=1, label = T)
plot2 = FeaturePlot(seurat_EGFR_mouse_fil, features = 'Ermn', reduction = "umap", pt.size=1, label = T)
plot2 = FeaturePlot(seurat_EGFR_mouse_fil, features = 'C1qa', reduction = "umap", pt.size=1, label = T)
plot2 = FeaturePlot(seurat_EGFR_mouse_fil, features = 'Plac8', reduction = "umap", pt.size=1, label = T)
plot2 = FeaturePlot(seurat_EGFR_mouse_fil, features = 'Cd74', reduction = "umap", pt.size=1, label = T)
plot2 = FeaturePlot(seurat_EGFR_mouse_fil, features = 'Ms4a7', reduction = "umap", pt.size=1, label = T)
plot2 = FeaturePlot(seurat_EGFR_mouse_fil, features = 'Lgals3', reduction = "umap", pt.size=1, label = T)
plot2 = FeaturePlot(seurat_EGFR_mouse_fil, features = 'Lyz2', reduction = "umap", pt.size=1, label = T)
plot2 = FeaturePlot(seurat_EGFR_mouse_fil, features = 'Col1a1', reduction = "umap", pt.size=1, label = T)
plot2 = FeaturePlot(seurat_EGFR_mouse_fil, features = 'Rgs5', reduction = "umap", pt.size=1, label = T)
plot2 = FeaturePlot(seurat_EGFR_mouse_fil, features = 'Cldn5', reduction = "umap", pt.size=1, label = T)
plot2 = FeaturePlot(seurat_EGFR_mouse_fil, features = 'Ramp2', reduction = "umap", pt.size=1, label = T)
plot2 = FeaturePlot(seurat_EGFR_mouse_fil, features = 'Ctla2a', reduction = "umap", pt.size=1, label = T)
plot2 = FeaturePlot(seurat_EGFR_mouse_fil, features = 'S100a4', reduction = "umap", pt.size=1, label = T)
plot2 = FeaturePlot(seurat_EGFR_mouse_fil, features = 'Aldh1l1', reduction = "umap", pt.size=1, label = T)
plot2 = FeaturePlot(seurat_EGFR_mouse_fil, features = 'Aldoc', reduction = "umap", pt.size=1, label = T)
plot2 = FeaturePlot(seurat_EGFR_mouse_fil, features = 'S100b', reduction = "umap", pt.size=1, label = T)
plot2 = FeaturePlot(seurat_EGFR_mouse_fil, features = 'Itgam', reduction = "umap", pt.size=1, label = T)
plot2 = FeaturePlot(seurat_EGFR_mouse_fil, features = 'Sparc', reduction = "umap", pt.size=1, label = T)
plot2 = FeaturePlot(seurat_EGFR_mouse_fil, features = 'Sparcl1', reduction = "umap", pt.size=1, label = T)

select.cells = CellSelector(plot=plot2)

Idents(seurat_EGFR_mouse_fil, cells=select.cells) = 'Endothelial'
Idents(seurat_EGFR_mouse_fil, cells=select.cells) = 'Pericyte'
Idents(seurat_EGFR_mouse_fil, cells=select.cells) = 'OPG'
Idents(seurat_EGFR_mouse_fil, cells=select.cells) = 'OLG'
Idents(seurat_EGFR_mouse_fil, cells=select.cells) = 'Myeloid'

seurat_EGFR_mouse_fil$Cell_type2 = Idents(seurat_EGFR_mouse_fil)
Idents(seurat_EGFR_mouse_fil) = seurat_EGFR_mouse_fil$RNA_snn_res.1.5
Idents(seurat_EGFR_mouse_fil) = seurat_EGFR_mouse_fil$Cell_type2

seurat_EGFR_mouse_fil.markers <- FindAllMarkers(seurat_EGFR_mouse_fil, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25, features = rownames(seurat_EGFR_mouse_fil)[!grepl('^MT-|^mt-|^RPS|^RPL|^Rps|^Rpl', rownames(seurat_EGFR_mouse_fil))])
seurat_EGFR_mouse_fil.markers$selectivity = seurat_EGFR_mouse_fil.markers$pct.1 - seurat_EGFR_mouse_fil.markers$pct.2
seurat_EGFR_mouse_fil.markers[seurat_EGFR_mouse_fil.markers$cluster=='6',]

top5 = seurat_EGFR_mouse_fil.markers %>% group_by(cluster) %>% top_n(n=20, wt = selectivity)
dittoHeatmap(seurat_EGFR_mouse_fil, main = 'All cells', annot.colors = c(dittoColors()), 
             c(Janus_genesets2[[1]], Janus_genesets2[[2]], Janus_genesets2[[3]], Janus_genesets2[[7]], Janus_genesets2[[8]], Janus_genesets2[[5]], Janus_genesets2[[6]]),
             scaled.to.max = T, complex = T, use_raster = T, cluster_cols =T, cluster_rows =F, slot = 'data', 
             annot.by = c('Cell_type2'), #order.by = 'Cell_type4',
             heatmap.colors.max.scaled = colorRampPalette(c('blue', 'white', 'red'))(50),  
             heatmap.colors = colorRampPalette(c('blue', 'white', 'red'))(50))


FeaturePlot(seurat_EGFR_mouse_fil, features = c("Myeloid_comm_sig", "Fibroblast_sig", "Endothelial_sig", 
                                                'Oligodendrocyte_sig', 'Myelin_sig'), min.cutoff = "q9", max.cutoff = 5)
dittoPlot(seurat_EGFR_mouse_fil,
          plots = c('jitter', 'vlnplot'), "Myeloid_comm_sig",
          #adjustment = 'z-score',
          group.by = 'Cell_type2')
dittoPlot(seurat_EGFR_mouse_fil,
          plots = c('jitter', 'ridgeplot'), "Olig1",
          adjustment = 'z-score',
          group.by = 'Cell_type2')

DimPlot(seurat_EGFR_mouse_fil, reduction = 'tsne', label = TRUE, pt.size = 1.5, group.by = 'seurat_clusters')
DimPlot(seurat_EGFR_mouse_fil, reduction = 'tsne', label = TRUE, pt.size = 1.5, group.by = 'Cell_type2')
DimPlot(seurat_EGFR_mouse_fil, label = TRUE, pt.size = 1.5, group.by = 'sample.name')

Idents(seurat_EGFR_mouse_fil) <- factor(Idents(seurat_EGFR_mouse_fil))
markers.to.plot <- c("Tyrobp", "Ctss", "Fcer1g", 
                     "Cd74",   "H2-Eb1", "H2-Aa","H2-Ab1", "Lsp1","Napsa","Cytip", "Gpr171", "Plbd1","Chil3","Ccr2","S100a4",
                     "Col1a1", "Rgs5", "Vtn", "Bgn", "Higd1b", "Cd248", "Pdgfrb", "Des", "Rarres2", "Fstl1", "Col1a2", "Acta2",
                     "Cldn5", "Ramp2", "Egfl7", "Ptprb", "Adgrl4", "Ctla2a", "Tmem252", "Pecam1", "Scgb3a1", "Ly6c1", "Itm2a", "Kdr", "Cd34",
                     "Ptprz1","Cspg5", "Lhfpl3", "Fabp7",  "Scrg1", "Olig1", "Gap43", "Gng3", "Gpm6a", "Bex2", "Olig2",
                     "Klk6", "Ermn", "Mag", "Mog", "Nkx6-2", "Cldn11", "Stmn4", "Ugt8a", "Aplp1", "Tubb4a", "Plp1", "Serpina3n", "Fez1")


DotPlot(seurat_EGFR_mouse_fil, features = c(markers.to.plot, 'P2ry12'), dot.scale = 8, cluster.idents = T)+ 
    coord_flip() +
    RotatedAxis() +  
    scale_colour_gradient2(low = "#54A2FF", mid = "#FFFEC3", high = "Tomato", midpoint = 1)

## Pericytes are labeled as fibroblast here
seurat_EGFR_Endo <- subset(seurat_EGFR_mouse_fil, subset = Cell_type2 == 'Endothelial') 
seurat_EGFR_Fibro <- subset(seurat_EGFR_mouse_fil, subset = Cell_type2 == 'Fibroblast') ## Pericytes
seurat_EGFR_Myeloid <- subset(seurat_EGFR_mouse_fil, subset = Cell_type2 == 'Myeloid') 
seurat_EGFR_OLG <- subset(seurat_EGFR_mouse_fil, subset = Cell_type2 == 'OLG') 
seurat_EGFR_OPC <- subset(seurat_EGFR_mouse_fil, subset = Cell_type2 == 'OPC') 

seurat_EGFR_Brain = merge(seurat_EGFR_OLG, y=c(seurat_EGFR_OPC), project='Janus', merge.data=TRUE)
seurat_EGFR_Brain = ProcessSeurat_pdx2(seurat_EGFR_Brain, 1, 1, species = 2)

seurat_EGFR_merge1 = merge(seurat_EGFR_Endo, y=c(seurat_EGFR_Fibro), project='Janus', merge.data=TRUE)
seurat_EGFR_merge1$EGFR_status = factor(seurat_EGFR_merge1$sample.name, levels = c('WT', 'KO'))
seurat_EGFR_Endo$EGFR_status = factor(seurat_EGFR_Endo$sample.name, levels = c('WT', 'KO'))
seurat_EGFR_Fibro$EGFR_status = factor(seurat_EGFR_Fibro$sample.name, levels = c('WT', 'KO'))
seurat_EGFR_Brain$EGFR_status = factor(seurat_EGFR_Brain$sample.name, levels = c('WT', 'KO'))
seurat_EGFR_Myeloid$EGFR_status = factor(seurat_EGFR_Myeloid$sample.name, levels = c('WT', 'KO'))
seurat_EGFR_OPC$EGFR_status = factor(seurat_EGFR_OPC$sample.name, levels = c('WT', 'KO'))
seurat_EGFR_mouse_fil$EGFR_status = factor(seurat_EGFR_mouse_fil$sample.name, levels = c('WT', 'KO'))
seurat_EGFR_human$EGFR_status = factor(seurat_EGFR_human$sample.name, levels = c('WT', 'KO'))

seurat_EGFR_mouse_merged = merge(seurat_EGFR_merge1, y=c(seurat_EGFR_Brain, seurat_EGFR_Myeloid), project='Janus', merge.data=TRUE)
seurat_EGFR_mouse_merged = ProcessSeurat_pdx2(seurat_EGFR_mouse_merged, 1, 1, species = 2)
saveRDS(seurat_EGFR_mouse_merged, file.path(Janus_dir3, "seurat_EGFR_mouse_merged.rds")) ## Object 4.1 Murine cells withv cell type annotations


seurat_EGFR_merge1 <- subset(seurat_EGFR_mouse_fil, subset = Cell_type2 == 'Endothelial' | Cell_type2 == 'Fibroblast') 
seurat_EGFR_merge1 = ProcessSeurat_pdx2(seurat_EGFR_merge1, 1, 1, species = 2)

seurat_EGFR_merge1.markers <- FindAllMarkers(seurat_EGFR_merge1, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25, features = rownames(seurat_EGFR_merge1)[!grepl('^MT-|^mt-|^RPS|^RPL|^Rps|^Rpl', rownames(seurat_EGFR_merge1))])
seurat_EGFR_merge1.markers$selectivity = seurat_EGFR_merge1.markers$pct.1 - seurat_EGFR_merge1.markers$pct.2
seurat_EGFR_merge1.markers[seurat_EGFR_merge1.markers$cluster=='5',] #

top5 = seurat_EGFR_merge1.markers %>% group_by(cluster) %>% top_n(n=20, wt = selectivity)
write.csv(top5, 'GBM_Vascular_top20_selectivity.csv')
top5 = seurat_EGFR_merge1.markers %>% group_by(cluster) %>% top_n(n=20, wt = avg_log2FC)
write.csv(top5, 'GBM_Vascular_top20_FC.csv')
top5 = seurat_EGFR_merge1.markers %>% group_by(cluster) %>% top_n(n=20, wt = p_val_adj)
write.csv(top5, 'GBM_Vascular_top20_adj-pval.csv')

top5$gene
top5$p_val_adj
merge1_genes = c("Cx3cr1", "Sparc", "Cd81", "Olfml3", "Ccl12",
                 "C1qb","Cxcl13","Hexb", "Ms4a7","Lgals1","Msrb1","Lgals3","Fxyd5",
                 "Hp","Plac8","S100a6","Ccr2","S100a4","Ahnak","Crip1",
                 "H2-Eb1","H2-Aa", "Plbd1","Cytip","Cd74","Lsp1",
                 "Spp1","Fam20c","Ctla2b","Fabp5","Dab2", 'Cd68')
dittoHeatmap(seurat_EGFR_merge1, main = 'Vascular cells', annot.colors = c(dittoColors()), 
             c(unique(top5$gene)),
             scaled.to.max = T, complex = T, use_raster = T, cluster_cols =F, cluster_rows =T, slot = 'data', 
             annot.by = c('Cell_type2', 'seurat_clusters', 'sample.name'), order.by = 'seurat_clusters',
             heatmap.colors.max.scaled = colorRampPalette(c('blue', 'white', 'red'))(50),  
             heatmap.colors = colorRampPalette(c('blue', 'white', 'red'))(50))

seurat_EGFR_merge1$RNA_snn_res.1 = factor(seurat_EGFR_merge1$RNA_snn_res.1, levels = c('1', '2','3','6', '0', '4', '5'))
Idents(seurat_EGFR_merge1) = seurat_EGFR_merge1$RNA_snn_res.1

####
seurat_EGFR_OPC = ProcessSeurat_pdx2(seurat_EGFR_OPC, 0.5, 1, species = 2)
seurat_EGFR_OPC.markers <- FindAllMarkers(seurat_EGFR_OPC, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25, features = rownames(seurat_EGFR_OPC)[!grepl('^MT-|^mt-|^RPS|^RPL|^Rps|^Rpl', rownames(seurat_EGFR_OPC))])
seurat_EGFR_OPC.markers$selectivity = seurat_EGFR_OPC.markers$pct.1 - seurat_EGFR_OPC.markers$pct.2
seurat_EGFR_OPC.markers[seurat_EGFR_OPC.markers$cluster=='1',]
seurat_EGFR_OPC.markers = seurat_EGFR_OPC.markers[order(seurat_EGFR_OPC.markers$avg_log2FC, decreasing =T),]
top5 = seurat_EGFR_OPC.markers %>% group_by(cluster) %>% top_n(n=10, wt = selectivity)

seurat_EGFR_OPC <- RenameIdents(seurat_EGFR_OPC, '0' = 'OPC_1', '1' = 'OPC_2')
seurat_EGFR_OPC$Cell_type3 = Idents(seurat_EGFR_OPC)

DotPlot(seurat_EGFR_OPC, features = unique(top5$gene), dot.scale = 8, cluster.idents = F)+ 
    coord_flip() +
    RotatedAxis() +  
    scale_colour_gradient2(low = "#54A2FF", mid = "#FFFEC3", high = "Tomato", midpoint = 0)

#### OLG
seurat_EGFR_OLG = ProcessSeurat_pdx2(seurat_EGFR_OLG, 0.5, 1, species = 2)
seurat_EGFR_OLG.markers <- FindAllMarkers(seurat_EGFR_OLG, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25, features = rownames(seurat_EGFR_OLG)[!grepl('^MT-|^mt-|^RPS|^RPL|^Rps|^Rpl', rownames(seurat_EGFR_OLG))])
seurat_EGFR_OLG.markers$selectivity = seurat_EGFR_OLG.markers$pct.1 - seurat_EGFR_OLG.markers$pct.2
seurat_EGFR_OLG.markers[seurat_EGFR_OLG.markers$cluster=='1',]
top5 = seurat_EGFR_OLG.markers %>% group_by(cluster) %>% top_n(n=10, wt = selectivity)

seurat_EGFR_OLG <- RenameIdents(seurat_EGFR_OLG, '0' = 'OLG')
seurat_EGFR_OLG$Cell_type3 = Idents(seurat_EGFR_OLG)

DotPlot(seurat_EGFR_OLG, features = unique(top5$gene), dot.scale = 8, cluster.idents = F)+ 
    coord_flip() +
    RotatedAxis() +  
    scale_colour_gradient2(low = "#54A2FF", mid = "#FFFEC3", high = "Tomato", midpoint = 0)



####################

seurat_EGFR_Endo = ProcessSeurat_pdx2(seurat_EGFR_Endo, 0.5, 1, species = 2)
saveRDS(seurat_EGFR_Endo, file.path(Janus_dir3, "seurat_EGFR_Endo.rds"))
rm(seurat_EGFR_Endo)
seurat_EGFR_Endo = readRDS(file.path(Janus_dir3, "seurat_EGFR_Endo.rds"))

seurat_EGFR_Endo.markers <- FindAllMarkers(seurat_EGFR_Endo, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25, features = rownames(seurat_EGFR_Endo)[!grepl('^MT-|^mt-|^RPS|^RPL|^Rps|^Rpl', rownames(seurat_EGFR_Endo))])
seurat_EGFR_Endo.markers$selectivity = seurat_EGFR_Endo.markers$pct.1 - seurat_EGFR_Endo.markers$pct.2
write.csv(seurat_EGFR_Endo.markers, 'seurat_EGFR_Endo.markers.csv')

top5 = seurat_EGFR_mouse.markers %>% group_by(cluster) %>% top_n(n=20, wt = selectivity)

#######################

seurat_EGFR_Fibro = ProcessSeurat_pdx2(seurat_EGFR_Fibro, 0.5, 1, 2)
saveRDS(seurat_EGFR_Fibro, file.path(Janus_dir3, "seurat_EGFR_Fibro.rds"))
rm(seurat_EGFR_Fibro)
seurat_EGFR_Fibro = readRDS(file.path(Janus_dir3, "seurat_EGFR_Fibro.rds"))

seurat_EGFR_Myeloid = ProcessSeurat_pdx2(seurat_EGFR_Myeloid, 1.5, 2, 2)
saveRDS(seurat_EGFR_Myeloid, file.path(Janus_dir3, "seurat_EGFR_Myeloid.rds"))
rm(seurat_EGFR_Myeloid)
seurat_EGFR_Myeloid = readRDS(file.path(Janus_dir3, "seurat_EGFR_Myeloid.rds"))

seurat_EGFR_Brain = ProcessSeurat_pdx2(seurat_EGFR_Brain, 0.5, 1, 2)
saveRDS(seurat_EGFR_Brain, file.path(Janus_dir3, "seurat_EGFR_Brain.rds"))
rm(seurat_EGFR_Brain)
seurat_EGFR_Brain = readRDS(file.path(Janus_dir3, "seurat_EGFR_Brain.rds"))

###################################

seurat_EGFR_Myeloid.markers <- FindAllMarkers(seurat_EGFR_Myeloid, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25, features = rownames(seurat_EGFR_Myeloid)[!grepl('^MT-|^mt-|^RPS|^RPL|^Rps|^Rpl', rownames(seurat_EGFR_Myeloid))])
seurat_EGFR_Myeloid.markers$selectivity = seurat_EGFR_Myeloid.markers$pct.1 - seurat_EGFR_Myeloid.markers$pct.2
seurat_EGFR_Myeloid.markers[seurat_EGFR_Myeloid.markers$cluster=='1',]

top5 = seurat_EGFR_Myeloid.markers %>% group_by(cluster) %>% top_n(n=20, wt = selectivity)
write.csv(top5, 'GBM_Myeloid_top20_selectivity.csv')
top5 = seurat_EGFR_Myeloid.markers %>% group_by(cluster) %>% top_n(n=10, wt = avg_log2FC)
write.csv(top5, 'GBM_Myeloid_top20_FC.csv')
top5 = seurat_EGFR_Myeloid.markers %>% group_by(cluster) %>% top_n(n=20, wt = p_val_adj)
write.csv(top5, 'GBM_Myeloid_top20_adj-pval.csv')

top5$gene

Myeloid_genes = c("Cx3cr1", "Sparc", "Cd81", "Olfml3", "Ccl12",
                  "C1qb","Cxcl13","Hexb", "Ms4a7","Lgals1","Msrb1","Lgals3","Fxyd5",
                  "Hp","Plac8","S100a6","Ccr2","S100a4","Ahnak","Crip1",
                  "H2-Eb1","H2-Aa", "Plbd1","Cytip","Cd74","Lsp1",
                  "Spp1","Fam20c","Ctla2b","Fabp5","Dab2", 'Cd68')
seurat_EGFR_Myeloid$RNA_snn_res.0.5 = factor(seurat_EGFR_Myeloid$RNA_snn_res.0.5, levels = c('1', '0','3','6', '2', '5', '4'))
seurat_EGFR_Myeloid$seurat_clusters = seurat_EGFR_Myeloid$RNA_snn_res.0.5

Idents(seurat_EGFR_Myeloid) = seurat_EGFR_Myeloid$RNA_snn_res.0.5
seurat_EGFR_Myeloid <- RenameIdents(seurat_EGFR_Myeloid, '0' = "Microglia", `1` = "Microglia", `3` = "Microglia_Proliferating",
                                    `6` = "Microglia_SPP1+", '2' = 'Mono/Mac_1',  `4` = 'Mono/Mac_2', '5' = 'Mono/Mac_3')
seurat_EGFR_Myeloid$Cell_type3 = Idents(seurat_EGFR_Myeloid)

DotPlot(seurat_EGFR_Myeloid, features = Myeloid_genes, dot.scale = 8, cluster.idents = F)+ 
    coord_flip() +
    RotatedAxis() +  
    scale_colour_gradient2(low = "#54A2FF", mid = "#FFFEC3", high = "Tomato", midpoint = 0.5)

######################################################################### # Pericytes

seurat_EGFR_Fibro.markers <- FindAllMarkers(seurat_EGFR_Fibro, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25, features = rownames(seurat_EGFR_Fibro)[!grepl('^MT-|^mt-|^RPS|^RPL|^Rps|^Rpl', rownames(seurat_EGFR_Fibro))])
seurat_EGFR_Fibro.markers$selectivity = seurat_EGFR_Fibro.markers$pct.1-seurat_EGFR_Fibro.markers$pct.2
seurat_EGFR_Fibro.markers[seurat_EGFR_Fibro.markers$cluster=='1',]

top5 = seurat_EGFR_Fibro.markers %>% group_by(cluster) %>% top_n(n=15, wt = selectivity)
top5$gene
top5[top5$cluster=='1',]
# Pericytes genes
Fibroblast_genes = c("AY036118", "Tmem176b", "Apoe", "Crip2", "Il34", "Col6a2", "Gpx3", "Ypel3",
                     "Gm9493", "Xist", "Erdr1", "Nme1", "Gm2000", "Gm9843", "Ddx21", 
                     "Wdr89", "Cycs", "Uqcr11","Bola2","Eif2s1","Uqcc2","Zfp91","Ndufv3","Hdgf")

seurat_EGFR_Fibro <- RenameIdents(seurat_EGFR_Fibro, '0' = 'Pericyte_1', '1' = 'Pericyte_2')
seurat_EGFR_Fibro$Cell_type3 = Idents(seurat_EGFR_Fibro)

DotPlot(seurat_EGFR_Fibro, features = Fibroblast_genes, dot.scale = 8, cluster.idents = F)+ 
    coord_flip() +
    RotatedAxis() +  
    scale_colour_gradient2(low = "#54A2FF", mid = "#FFFEC3", high = "Tomato", midpoint = 0)

#####################################################################################

seurat_EGFR_Endo.markers <- FindAllMarkers(seurat_EGFR_Endo, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25, features = rownames(seurat_EGFR_Endo)[!grepl('^MT-|^mt-|^RPS|^RPL|^Rps|^Rpl', rownames(seurat_EGFR_Endo))])
seurat_EGFR_Endo.markers$selectivity = seurat_EGFR_Endo.markers$pct.1 - seurat_EGFR_Endo.markers$pct.2
seurat_EGFR_Endo.markers[seurat_EGFR_Endo.markers$cluster=='3',]

seurat_EGFR_Endo.markers <- FindAllMarkers(seurat_EGFR_Endo, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25, features = rownames(seurat_EGFR_Endo)[!grepl('^MT-|^mt-|^RPS|^RPL|^Rps|^Rpl', rownames(seurat_EGFR_Endo))])
seurat_EGFR_Endo.markers$selectivity = seurat_EGFR_Endo.markers$pct.1-seurat_EGFR_Endo.markers$pct.2
seurat_EGFR_Endo.markers[seurat_EGFR_Endo.markers$cluster=='1',]

top5 = seurat_EGFR_Endo.markers %>% group_by(cluster) %>% top_n(n=10, wt = selectivity)
top5$gene
top5[top5$cluster=='1',]
Endothelial_genes = c("Pmepa1","Sec11c","Tmem167","G3bp2","Mif", "Ywhah",
                      "Cycs","Igfbp3",
                      "Lrg1","Col18a1","Ifitm2","Cxcl12","Sema3c","Ndrg1","Podxl",
                      "Cenpf", "Top2a" ,"Mki67","Dut","Smc2")

seurat_EGFR_Endo$RNA_snn_res.0.5 = factor(seurat_EGFR_Endo$RNA_snn_res.0.5, levels = c('3','0', '1', '2'))
Idents(seurat_EGFR_Endo) = seurat_EGFR_Endo$RNA_snn_res.0.5

seurat_EGFR_Endo <- RenameIdents(seurat_EGFR_Endo, '3' = 'Endo_Proliferative', '0' = 'Endo_Permeable', '1' = 'Endo_Angiogenic', '2'= 'Endo_Migrating')
seurat_EGFR_Endo$Cell_type3 = Idents(seurat_EGFR_Endo)

DotPlot(seurat_EGFR_Endo, features = Endothelial_genes, dot.scale = 8, cluster.idents = F)+ 
    coord_flip() +
    RotatedAxis() +  
    scale_colour_gradient2(low = "#54A2FF", mid = "#FFFEC3", high = "Tomato", midpoint = 0)

###################################################

seurat_EGFR_Brain.markers <- FindAllMarkers(seurat_EGFR_Brain, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25, features = rownames(seurat_EGFR_Brain)[!grepl('^MT-|^mt-|^RPS|^RPL|^Rps|^Rpl', rownames(seurat_EGFR_Brain))])
seurat_EGFR_Brain.markers$selectivity = seurat_EGFR_Brain.markers$pct.1 - seurat_EGFR_Brain.markers$pct.2
seurat_EGFR_Brain.markers[seurat_EGFR_Brain.markers$cluster=='3',]

top5 = seurat_EGFR_Brain.markers %>% group_by(cluster) %>% top_n(n=20, wt = selectivity)
write.csv(top5, 'GBM_MyelinSheath_top20_selectivity.csv')
top5 = seurat_EGFR_Brain.markers %>% group_by(cluster) %>% top_n(n=20, wt = avg_log2FC)
write.csv(top5, 'GBM_MyelinSheath_top20_FC.csv')
top5 = seurat_EGFR_Brain.markers %>% group_by(cluster) %>% top_n(n=20, wt = p_val_adj)
write.csv(top5, 'GBM_MyelinSheath_top20_adj-pval.csv')


top5$gene
Brain_genes = top5$gene 
Brain_genes = c("Ptprz1", "Fabp7","Cspg5","Ppp1r14b" ,"Rtn1", "Marcks","Basp1","Tubb2b",
                "Gap43","Pdgfra","Ermn","Mal","Car2","Trf",
                "Klk6","Apod","Mag","Plp1","Nkx6-2","Ndrg1"
)

seurat_EGFR_Brain$RNA_snn_res.0.5 = factor(seurat_EGFR_Brain$RNA_snn_res.0.5, levels = c('0', '1'))
Idents(seurat_EGFR_Brain) = seurat_EGFR_Brain$RNA_snn_res.0.5

DotPlot(seurat_EGFR_Brain, features = Brain_genes, dot.scale = 8, cluster.idents = F)+ 
    coord_flip() +
    RotatedAxis() +  
    scale_colour_gradient2(low = "#54A2FF", mid = "#FFFEC3", high = "Tomato", midpoint = 0)

Idents(seurat_EGFR_Endo) = seurat_EGFR_Endo$sample.name
seurat_EGFR_Endo.markers2 <- FindAllMarkers(seurat_EGFR_Endo, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
seurat_EGFR_Endo.markers[seurat_EGFR_Endo.markers$cluster=='3',]
seurat_EGFR_Endo.markers2[seurat_EGFR_Endo.markers2$cluster=='WT',]
seurat_EGFR_Endo.markers2[seurat_EGFR_Endo.markers2$cluster=='KO',]


####### Investigate human genes in murine cells 
seurat_EGFR = readRDS(file.path(Janus_dir3, "seurat_EGFR.rds"))
colnames(seurat_EGFR_mouse_merged)

## Select the murine cells in the object with both hg & mm genes
seurat_EGFR = seurat_EGFR[,colnames(seurat_EGFR_mouse_merged)]
seurat_EGFR$frac_hg_genes

## Get gene signatures for human cells by EGFR status due to the significanct difference
Idents(seurat_EGFR_human) = seurat_EGFR_human$EGFR_status
seurat_EGFR_human.markers2 <- FindAllMarkers(seurat_EGFR_human, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25, features = rownames(seurat_EGFR_human)[!grepl('^MT-|^mt-|^RPS|^RPL|^Rps|^Rpl', rownames(seurat_EGFR_human))])
seurat_EGFR_human.markers2$selectivity = seurat_EGFR_human.markers2$pct.1 - seurat_EGFR_human.markers2$pct.2
seurat_EGFR_human.markers2
WT_markers = seurat_EGFR_human.markers2[seurat_EGFR_human.markers2$cluster == 'WT',]
KO_markers = seurat_EGFR_human.markers2[seurat_EGFR_human.markers2$cluster == 'KO',]
WT_markers = WT_markers[order(WT_markers$avg_log2FC, decreasing = T),]
KO_markers = KO_markers[order(KO_markers$avg_log2FC, decreasing = T),]
WT_markers[WT_markers$pct.1>0.9,'gene']
KO_markers[KO_markers$pct.1>0.9,'gene']

seurat_EGFR$WT_sig = (colMeans(seurat_EGFR@assays$RNA@counts[WT_markers[WT_markers$pct.1>0.9,'gene'][1:50],]) / colMeans(seurat_EGFR@assays$RNA@counts))
seurat_EGFR$KO_sig = (colMeans(seurat_EGFR@assays$RNA@counts[KO_markers[KO_markers$pct.1>0.9,'gene'][1:50],]) / colMeans(seurat_EGFR@assays$RNA@counts))
hist(seurat_EGFR[,seurat_EGFR$sample.name =='WT']$WT_sig, breaks = 100, main = 'WT_sig in EGFR_WT stroma')
hist(seurat_EGFR[,seurat_EGFR$sample.name =='KO']$KO_sig, breaks = 100, main = 'KO_sig in EGFR_KO stroma')

dittoHeatmap(seurat_EGFR, main = 'Expression of DEGs for human MSCs in murine stroma', annot.colors = c(dittoColors()), 
             c(WT_markers[WT_markers$pct.1>0.9,'gene'][1:50],
               KO_markers[KO_markers$pct.1>0.9,'gene'][1:50]),  
             scaled.to.max = T, complex = T, use_raster = T, cluster_cols =T, cluster_rows =T, slot = 'data',
             annot.by = c('frac_hg_genes', 'sample.name'), #order.by = c('Cell_type3'), 
             annotation_colors = EGFR_color,
             heatmap.colors.max.scaled = colorRampPalette(c('blue', 'white', 'red'))(50),  
             heatmap.colors = colorRampPalette(c('blue', 'white', 'red'))(50))


## Filter out human-murine doublets to avoid false positive human gene expression in murine cells
plot2 = FeatureScatter(seurat_EGFR[,seurat_EGFR$sample.name=='WT'], feature1 = 'frac_hg_genes', feature2 = 'WT_sig',  pt.size = 1.5)
plot2 = FeatureScatter(seurat_EGFR[,seurat_EGFR$sample.name=='KO'], feature1 = 'frac_hg_genes', feature2 = 'KO_sig',  pt.size = 1.5)

# 650 * 500
FeatureScatter(seurat_EGFR[,seurat_EGFR$sample.name=='WT'], feature1 = 'frac_hg_genes', feature2 = 'WT_sig',  pt.size = 1.5, group.by = 'Species_Doublets')
FeatureScatter(seurat_EGFR[,seurat_EGFR$sample.name=='KO'], feature1 = 'frac_hg_genes', feature2 = 'KO_sig',  pt.size = 1.5, group.by = 'Species_Doublets')

seurat_EGFR$EGFR_TPM = (seurat_EGFR@assays$RNA@counts['EGFR',]/seurat_EGFR$nCount_RNA)*1000000
hist(seurat_EGFR$EGFR_TPM, breaks = 100)
seurat_EGFR$EGFR = seurat_EGFR@assays$RNA@data['EGFR',] 

select.cells = CellSelector(plot=plot2)
Idents(seurat_EGFR, cells=select.cells) = 'WT_Doublets'
Idents(seurat_EGFR, cells=select.cells) = 'KO_Doublets'
Idents(seurat_EGFR, cells=select.cells) = 'Oligodendrocytes'
seurat_EGFR$Species_Doublets = Idents(seurat_EGFR)
table(colnames(seurat_EGFR) == colnames(seurat_EGFR_mouse_merged))
table(seurat_EGFR$Species_Doublets)

seurat_EGFR_mouse_merged$Species_Doublets = seurat_EGFR$Species_Doublets
seurat_EGFR_mouse_merged$EGFR_TPM = seurat_EGFR$EGFR_TPM
seurat_EGFR_mouse_merged$EGFR = seurat_EGFR$EGFR
seurat_EGFR_mouse_merged$WT_sig = seurat_EGFR$WT_sig
seurat_EGFR_mouse_merged$KO_sig = seurat_EGFR$KO_sig

saveRDS(seurat_EGFR_mouse_merged, file.path(Janus_dir3, "seurat_EGFR_mouse_merged.rds")) ## Object 4.1 Murine cells withv cell type annotations
