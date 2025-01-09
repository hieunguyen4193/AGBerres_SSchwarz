gc()
rm(list = ls())

#####----------------------------------------------------------------------#####
##### preparation
#####----------------------------------------------------------------------#####
path.to.storage <- "/media/hieunguyen/HD01/storage"

path.to.main.project.src <- "/home/hieunguyen/CRC1382/src_2023/AGBerres_SSchwarz"
source(file.path(path.to.main.project.src, "00_import_libraries.R"))
source(file.path(path.to.main.project.src, "00_helper_functions.R"))

path.to.downloaded.Rdata <- file.path(path.to.storage, "AGBerres_Maier_et_al_Nature_pub", "human_dc.rd")
load(path.to.downloaded.Rdata)

path.to.outdir <- "/home/hieunguyen/CRC1382/outdir/AGBerres"
path.to.main.output <- file.path(path.to.outdir, "human_dc_Maier_et_al_data")
dir.create(path.to.main.output, showWarnings = FALSE, recursive = TRUE)

path.to.03.output <- file.path(path.to.main.output, "03_output")
dir.create(path.to.03.output, showWarnings = FALSE, recursive = TRUE)

#####----------------------------------------------------------------------#####
##### create the seurat object
#####----------------------------------------------------------------------#####
count.data <- human_dc$filtered_umitab
adt.obj <- human_dc$adt_matrix_by_sample

MINCELLS <- 0
MINGENES <- 0
PROJECT <- "Maier_data"
s.obj <- CreateSeuratObject(counts = count.data, 
                            min.cells = MINCELLS, 
                            min.features = MINGENES, 
                            project = PROJECT)

s.obj[["percent.mt"]] <- PercentageFeatureSet(s.obj, pattern = "^mt-|^MT-")
s.obj[["percent.ribo"]] <- PercentageFeatureSet(s.obj, pattern = "^Rpl|^Rps|^RPL|^RPS")

adtdf <- adt.obj[[names(adt.obj)[[1]]]] %>%  as.data.frame() %>% rownames_to_column("hashtag")
for (group in names(adt.obj)[2:length(names(adt.obj))]){
  tmpdf <- adt.obj[[group]] %>% as.data.frame()  %>% rownames_to_column("hashtag")
  adtdf <- merge(adtdf, tmpdf, by.x = "hashtag", by.y = "hashtag", all.x = TRUE, all.y = TRUE)
}
adtdf <- adtdf %>% column_to_rownames("hashtag") %>% t() %>% as.data.frame() %>% rownames_to_column("barcode")
meta.data <- s.obj@meta.data  %>% rownames_to_column("barcode")  

colnames(adtdf) <- unlist(lapply(colnames(adtdf), function(x){
  if (x == "barcode"){
    return(x)
  } else {
    return(sprintf("ht_%s", x))
  }
  
}))

meta.data <- merge(meta.data, adtdf, by.x = "barcode", by.y = "barcode", all.x = TRUE) %>%
  column_to_rownames("barcode")
adtdf <- adtdf %>% column_to_rownames("barcode")
meta.data <- meta.data[row.names(s.obj@meta.data), ]
s.obj <- AddMetaData(object = s.obj, metadata = meta.data[, colnames(adtdf)], col.name = colnames(adtdf))

#####----------------------------------------------------------------------#####
##### Add cell annotation to the main seurat object
#####----------------------------------------------------------------------#####
annotationdf <- human_dc[["cell_to_annot"]] %>% as.data.frame() %>% rownames_to_column("barcode")
colnames(annotationdf) <- c("barcode", "celltype")

meta.data <- s.obj@meta.data %>% rownames_to_column("barcode")

meta.data <- merge(meta.data, annotationdf, by.x = "barcode", by.y = "barcode", all.x = TRUE) %>%
  column_to_rownames("barcode")
meta.data <- meta.data[row.names(s.obj@meta.data),]

s.obj <- AddMetaData(object = s.obj, metadata = meta.data$celltype, col.name = "celltype")

#####----------------------------------------------------------------------#####
##### Running clustering and dimensional reduction
#####----------------------------------------------------------------------#####
if (file.exists(file.path(path.to.03.output, "Maier_data_seurat_object.rds")) == FALSE){
  num.PCA <- 25
  num.PC.used.in.UMAP <- 25
  num.PC.used.in.Clustering <- 25
  my_random_seed <- 42
  cluster.resolution <- 0.5
  s.obj <- NormalizeData(s.obj) # ---> use Log Normalized
  s.obj <- FindVariableFeatures(s.obj, selection.method = "vst")
  s.obj <- ScaleData(s.obj, features = rownames(s.obj))
  
  chosen.assay <- "RNA"
  DefaultAssay(s.obj) <- chosen.assay
  
  s.obj <- RunPCA(s.obj, npcs = num.PCA, verbose = FALSE, reduction.name=sprintf("%s_PCA", chosen.assay))
  
  s.obj <- RunUMAP(s.obj, reduction = sprintf("%s_PCA", chosen.assay), 
                   dims = 1:num.PC.used.in.UMAP, reduction.name=sprintf("%s_UMAP", chosen.assay),
                   seed.use = my_random_seed, umap.method = "uwot")
  s.obj <- RunTSNE(s.obj, sprintf("%s_PCA", chosen.assay), 
                   dims = 1:num.PC.used.in.UMAP, reduction.name=sprintf("%s_TSNE", chosen.assay),
                   seed.use = my_random_seed)    
  s.obj <- FindNeighbors(s.obj, reduction = sprintf("%s_PCA", chosen.assay), dims = 1:num.PC.used.in.Clustering)
  
  s.obj <- FindClusters(s.obj, 
                        resolution = cluster.resolution, random.seed = 0)
  
  saveRDS(object = s.obj, file.path(path.to.03.output, "Maier_data_seurat_object.rds"))
} else {
  s.obj <- readRDS(file.path(path.to.03.output, "Maier_data_seurat_object.rds"))
}

