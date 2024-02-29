### Written by Minjun Kim, M.S., Ph.D student in Dept. of Human Genetics @ McGill Univ.
### minjun.kim@mail.mcgill.ca for any inquiries

### 00. Setting up functions
library(Seurat); library(scater); library(ensembldb); library(DropletUtils); library(biomaRt)
library(tidyverse); library(dbscan); library(scran); library(uwot); library(dittoSeq)


sceToSeurat = function(sce){
    colnames(sce) = colData(sce)[,'Barcode']
    temp_counts = assay(sce, "counts")
    colnames(temp_counts) = colData(sce)$Barcode
    meta = data.frame(colData(sce))
    rownames(meta) = meta$Barcode
    seu = CreateSeuratObject(counts=temp_counts, meta.data=meta)
    seu[["percent.mt"]] = colSums(as.matrix(seu@assays$RNA@counts)[grep('^MT-', rownames(seu)),])/seu$hg19
    seu[["percent.mt_mouse"]] = colSums(as.matrix(seu@assays$RNA@counts)[grep('^mt-', rownames(seu)),])/seu$mm10
    rm(sce)
    return(seu)
}

BasicSeurat_pdx = function(obj_10x){
    seu = CreateSeuratObject(counts=obj_10x)
    seu[["percent.mt"]] = PercentageFeatureSet(object = seu, pattern = "^hg19_MT-")
    seu[["percent.rb"]] = PercentageFeatureSet(object = seu, pattern = "^hg19_RP[SL]")
    seu[["percent.mt_m"]] = PercentageFeatureSet(object = seu, pattern = "^mm10_MT-")
    seu[["percent.rb_m"]] = PercentageFeatureSet(object = seu, pattern = "^mm10_RP[SL]")
    seu$GenePerRNA = seu$nFeature_RNA/seu$nCount_RNA
    seu <- subset(seu, subset = nCount_RNA > 5000) 
    seu <- subset(seu, subset = nFeature_RNA > 1000) # 650? 
    seu <- subset(seu, subset = percent.mt > 1) 
    seu <- subset(seu, subset = percent.mt < 40) 
    rm(obj_10x)
    return(seu)
}

ProcessSeurat_v3 = function(seu, res, scale){
    if(scale == '1'){
        seu = NormalizeData(object = seu, normalization.method = "LogNormalize", scale.factor = 10000)
        seu = FindVariableFeatures(seu, selection.method = "vst", assay = 'RNA', nfeatures = 2100)
        seu@assays$RNA@var.features = seu@assays$RNA@var.features[!grepl('^MT-|^RP[SL]' ,seu@assays$RNA@var.features)]
        seu=ScaleData(seu, vars.to.regress=c("nCount_RNA", "percent.mt"), features = rownames(seu))
    } else if(scale == '2'){
        print('normalization and scale will be skipped')
    }
    seu = FindVariableFeatures(seu, selection.method = "vst", assay = 'RNA', nfeatures = 2100)
    seu@assays$RNA@var.features = seu@assays$RNA@var.features[!grepl('^MT-|^RP[SL]' ,seu@assays$RNA@var.features)]
    #seu@assays$RNA@var.features = seu@assays$RNA@var.features[!seu@assays$RNA@var.features %in% excluded_genes]
    seu = Seurat::RunPCA(seu, features = VariableFeatures(object = seu), npcs = 100, vars.to.regress=c("nCount_RNA", "percent.mt"))
    seu = JackStraw(seu, num.replicate = 100, dims = 100)
    seu = ScoreJackStraw(seu, dims = 1:100)
    dim = c(1:100)[ifelse(seu@reductions$pca@jackstraw$overall.p.values[,2] < 0.05, T, F)]
    seu = FindNeighbors(seu, dims = dim, random.seed = 1234)
    seu = FindClusters(seu, random.seed = 1234, resolution=res)
    seu = Seurat::RunTSNE(object = seu, seed.use=1234, dims=dim, nthreads = 8)
    seu = Seurat::RunUMAP(object = seu, seed.use=1234, dims=dim, nthreads = 8)
    return(seu)
}

ProcessSeurat_pdx= function(seu, res, scale){
    if(scale == '1'){
        seu = NormalizeData(object = seu, normalization.method = "LogNormalize", scale.factor = 10000)
        seu = FindVariableFeatures(seu, selection.method = "vst", assay = 'RNA', nfeatures = 2100)
        seu=ScaleData(seu, features = rownames(seu))
    } else if(scale == '2'){
        print('normalization and scale will be skipped')
    }
    seu = FindVariableFeatures(seu, selection.method = "vst", assay = 'RNA', nfeatures = 2100)
    #seu@assays$RNA@var.features = seu@assays$RNA@var.features[!seu@assays$RNA@var.features %in% excluded_genes]
    seu = Seurat::RunPCA(seu, features = VariableFeatures(object = seu), npcs = 100)
    seu = JackStraw(seu, num.replicate = 100, dims = 100)
    seu = ScoreJackStraw(seu, dims = 1:100)
    dim = c(1:100)[ifelse(seu@reductions$pca@jackstraw$overall.p.values[,2] < 0.05, T, F)]
    seu = FindNeighbors(seu, dims = dim, random.seed = 1234)
    seu = FindClusters(seu, random.seed = 1234, resolution=res)
    seu = Seurat::RunTSNE(object = seu, seed.use=1234, dims=dim, nthreads = 8)
    seu = Seurat::RunUMAP(object = seu, seed.use=1234, dims=dim, nthreads = 8)
    return(seu)
}

ProcessSeurat_pdx2= function(seu, res, scale, species){
    if(species == '1'){
        print('Human genes will be considered only')
        seu = seu[1:32738,]
        if(scale == '1'){
            seu = NormalizeData(object = seu, normalization.method = "LogNormalize", scale.factor = 10000)
            seu = FindVariableFeatures(seu, selection.method = "vst", assay = 'RNA', nfeatures = 2100)
            seu=ScaleData(seu, features = rownames(seu))
        } else if(scale == '2'){
            print('normalization and scale will be skipped')
        }
        seu = FindVariableFeatures(seu, selection.method = "vst", assay = 'RNA', nfeatures = 2100)
        #seu@assays$RNA@var.features = seu@assays$RNA@var.features[!seu@assays$RNA@var.features %in% excluded_genes]
        seu = Seurat::RunPCA(seu, features = VariableFeatures(object = seu), npcs = 100)
        seu = JackStraw(seu, num.replicate = 100, dims = 100)
        seu = ScoreJackStraw(seu, dims = 1:100)
        dim = c(1:100)[ifelse(seu@reductions$pca@jackstraw$overall.p.values[,2] < 0.05, T, F)]
        seu = FindNeighbors(seu, dims = dim, random.seed = 1234)
        seu = FindClusters(seu, random.seed = 1234, resolution=res)
        seu = Seurat::RunTSNE(object = seu, seed.use=1234, dims=dim, nthreads = 8, check_duplicates = F)
        seu = Seurat::RunUMAP(object = seu, seed.use=1234, dims=dim, nthreads = 8, check_duplicates = F)
        return(seu)
        
    } else if(species == '2'){
        print('Mouse genes will be considered only')
        seu = seu[32739:60736,]
        if(scale == '1'){
            seu = NormalizeData(object = seu, normalization.method = "LogNormalize", scale.factor = 10000)
            seu = FindVariableFeatures(seu, selection.method = "vst", assay = 'RNA', nfeatures = 2100)
            seu=ScaleData(seu, features = rownames(seu))
        } else if(scale == '2'){
            print('normalization and scale will be skipped')
        }
        seu = FindVariableFeatures(seu, selection.method = "vst", assay = 'RNA', nfeatures = 2100)
        #seu@assays$RNA@var.features = seu@assays$RNA@var.features[!seu@assays$RNA@var.features %in% excluded_genes]
        seu = Seurat::RunPCA(seu, features = VariableFeatures(object = seu), npcs = 100)
        seu = JackStraw(seu, num.replicate = 100, dims = 100)
        seu = ScoreJackStraw(seu, dims = 1:100)
        dim = c(1:100)[ifelse(seu@reductions$pca@jackstraw$overall.p.values[,2] < 0.05, T, F)]
        seu = FindNeighbors(seu, dims = dim, random.seed = 1234)
        seu = FindClusters(seu, random.seed = 1234, resolution=res)
        seu = Seurat::RunTSNE(object = seu, seed.use=1234, dims=dim, nthreads = 8, check_duplicates = F)
        seu = Seurat::RunUMAP(object = seu, seed.use=1234, dims=dim, nthreads = 8, check_duplicates = F)
        return(seu)
    }
}
