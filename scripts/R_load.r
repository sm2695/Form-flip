# Script to load files saved from Anndata object (.h5ad) in scanpy. 

library(Matrix)
save_dir<-'./data_for_R' #set  directory
counts<-readMM(paste0(save_dir,'/counts.mtx'))
dim(counts)
cellMeta<-read.csv(paste0(save_dir,'/counts_cellMeta.csv'))
head(cellMeta)
geneMeta<-read.csv(paste0(save_dir,'/counts_geneMeta.csv'))
dim(geneMeta)
head(geneMeta)
### Set the rownames and colnames
rownames(counts)<-cellMeta$Barcode
colnames(counts)<-geneMeta$GeneName

#Make a seurat object

seo <- CreateSeuratObject(counts = t(counts), project = "min", min.cells = 3, min.features = 200)
### Set the meta data
seo@meta.data<-cbind(cellMeta,seo@meta.data)
rownames(seo@meta.data)<-colnames(seo)
### Normalize the data
seo <- NormalizeData(seo)


## Add precimputed dimensionality reductions save in cell metatdata

seo@reductions$umap2<-seo@reductions$umap
head(seo@reductions$umap@cell.embeddings)
scanpy_umap<-seo@meta.data[,c('UMAP1','UMAP2')]
colnames(scanpy_umap)<-c('UMAP_1','UMAP_2')
seo@reductions$umap@cell.embeddings<-as.matrix(scanpy_umap)

saveRDS(seo, "Converted_seuratobj.rds")
