library(Matrix)
library(Seurat)

#where are the files
input = ''
#where to save
dir=''

#adapted from 'https://medium.com/@daimin0514/how-to-convert-singlecellexperiment-to-anndata-8ec678b3a99e'
#Add wnn coordinates as needed
convertAnnData2Seurat <- function(
    input_dir='',
    save_dir='',
    save_prefix='',
    project_name='SEO',
    verbose=FALSE
){

    counts<-readMM(paste0(input_dir,'/counts.mtx'))
    cellMeta<-read.csv(paste0(input_dir,'/counts_cellMeta.csv'))
    geneMeta<-read.csv(paste0(input_dir,'/counts_geneMeta.csv'))
    ### Set the rownames and colnames
    rownames(counts)<-cellMeta$Barcode
    colnames(counts)<-geneMeta$GeneName
    
    
    seo <- CreateSeuratObject(counts = t(counts), project = "P1"
                              # , min.cells = 3, min.features = 200
                             )
    seo@meta.data<-cbind(cellMeta,seo@meta.data)
    rownames(seo@meta.data)<-colnames(seo)
    seo <- NormalizeData(seo, verbose = verbose)

    seo[["percent.mt"]] <- PercentageFeatureSet(seo, pattern = "^Mt-")

    seo <- FindVariableFeatures(seo, selection.method = "vst", nfeatures = 3000, verbose = verbose)

    seo <- ScaleData(seo, verbose = verbose)

    seo <- RunPCA(seo, features = VariableFeatures(object = seo), verbose = verbose)

    seo <- RunUMAP(seo, dims = 1:30, verbose = verbose)

    seo@reductions$umap2<-seo@reductions$umap
    head(seo@reductions$umap@cell.embeddings)
    scanpy_umap<-seo@meta.data[,c('UMAP1','UMAP2')]
    colnames(scanpy_umap)<-c('UMAP_1','UMAP_2')
    seo@reductions$umap@cell.embeddings<-as.matrix(scanpy_umap)
    
    
    print("Saving the Seurat object as .rds file:")
    if(nchar(save_dir)>0){
    dir.create(save_dir, showWarnings = FALSE)
        saveRDS(seo, paste0(save_dir,'/',save_prefix,'.rds'))
    } else{
        saveRDS(seo, paste0('./',save_prefix,'.rds'))
    
    } 
    print("Done")
   
    return(seo)
   
}

seo<-convertAnnData2Seurat(
    input_dir=input,
    save_dir=dir,
    save_prefix='seo',
    project_name='SEO'
)
