
reference = readRDS("~/Documents/starting_kit_phase1/data/reference_pdac.rds")
mixes = readRDS("~/Documents/starting_kit_phase1/data/mixes1_insilicodirichletEMFA_pdac.rds")



library(Seurat)
sc = list()
for (i in 1:length(reference$ref_scRNA)) {
  sc[[i]] = CreateSeuratObject(as.matrix(reference$ref_scRNA[[i]]$counts), metadata = reference$ref_scRNA[[i]]$metadata)
  sc[[i]]@meta.data$cell.annot = reference$ref_scRNA[[i]]$metadata$cell_type
  Idents(sc[[i]]) = sc[[i]]@meta.data$cell.annot
}



library(SeuratData)
datasets <- c("Peng", "Baron", "Raghavan")
seuratobject <- merge(sc[[1]], y = c(sc[[2]], sc[[3]]), add.cell.ids = datasets, project = "seuratobject")

pattern <- paste(datasets, collapse = "|")
library(stringr)
seuratobject@meta.data$Dataset <- str_extract(colnames(seuratobject), pattern)


seuratobject <- FindVariableFeatures(seuratobject, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seuratobject)
seuratobject <- ScaleData(seuratobject, features = all.genes)
seuratobject <- RunPCA(seuratobject, features = VariableFeatures(object = seuratobject))
ElbowPlot(seuratobject)
seuratobject <- FindNeighbors(seuratobject, dims = 1:20)
seuratobject <- FindClusters(seuratobject, resolution = 0.5)
seuratobject <- RunUMAP(seuratobject, dims = 1:10)
Idents(seuratobject) = seuratobject@meta.data$Dataset
DimPlot(seuratobject, reduction = "umap")


obj.list <- SplitObject(seuratobject, split.by = 'Dataset')
for(i in 1:length(obj.list)){
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
}


features <- SelectIntegrationFeatures(object.list = obj.list)

anchors <- FindIntegrationAnchors(object.list = obj.list,
                                  anchor.features = features)

seurat.integrated <- IntegrateData(anchorset = anchors)

seurat.integrated <- ScaleData(object = seurat.integrated)
seurat.integrated <- RunPCA(object = seurat.integrated)
seurat.integrated <- RunUMAP(object = seurat.integrated, dims = 1:50)

Idents(seurat.integrated) = seurat.integrated@meta.data$Dataset
DimPlot(seurat.integrated, reduction = "umap")



saveRDS(seurat.integrated, "scRNAseq_integrated.rds")

