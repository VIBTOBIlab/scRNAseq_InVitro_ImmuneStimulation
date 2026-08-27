library(R.utils)
Args<-commandArgs(trailingOnly = TRUE, asValues = TRUE)
string<-"##################\nInput Arguments:\n##################\n"
for (i in names(Args)){
	string<-paste0(string, i, " : ", Args[i], "\n")}
string<-paste0(string,"#################\n")
cat(string)

# Correct path for libraries
.libPaths(Args$LibFolder)

# Load libraries
library(Seurat)
library(SeuratWrappers)
library(tidyverse)
library(patchwork)

# Expand memory if necessary
#options(future.globals.maxSize = 1e9)
options(future.globals.maxSize = 356000 * 1024^2)
options(Seurat.object.assay.version = Args$SeuratVersion)

packageVersion("Seurat")

# RDS file name for storage Seurat objects
save_Int<- Args$SaveRDSFile

# Load Seurat object containing all samles
MergeData<-readRDS(Args$MergedDataFile)
# Removing cells originating from baseline samples
MergeData<-subset(x = MergeData, Stimulation != "Baseline")
DefaultAssay(MergeData)<-'RNA'

# Copy data to new Seurat object
scrna.list <- MergeData

print("### RNA ###")
scrna.list[["RNA"]] <- split(MergeData[["RNA"]], f = MergeData@meta.data[[Args$SplitArg]])

print("### NORMALIZATION ###")
scrna.list <- Seurat::NormalizeData(scrna.list)

print("### VARIABLE FEATURE SELECTION ###")
scrna.list <- Seurat::FindVariableFeatures(scrna.list)

print("### SCALING ###")
scrna.list <- Seurat::ScaleData(scrna.list, features = rownames(scrna.list))

print("### DIMENSIONALITY REDUCTION ###")
scrna.list <- Seurat::RunPCA(object = scrna.list, assay="RNA", reduction.name = "pca.unintegrated.rna")

scrna.list <- FindNeighbors(scrna.list, dims = 1:20, reduction = "pca.unintegrated.rna")

scrna.list <- FindClusters(scrna.list, resolution = 2, cluster.name = "unintegrated_clusters.rna")

scrna.list <- RunUMAP(scrna.list, dims = 1:20, reduction = "pca.unintegrated.rna", reduction.name = "umap.unintegrated.rna")

print("### CCA INTEGRATION ###")
scrna.list <- IntegrateLayers(
  object = scrna.list, method = CCAIntegration,
  orig.reduction = "pca.unintegrated.rna", new.reduction = "integrated.cca.rna",
  verbose = FALSE
)

print("### RPCA INTEGRATION ###")
scrna.list <- IntegrateLayers(
  object = scrna.list, method = RPCAIntegration,
  orig.reduction = "pca.unintegrated.rna", new.reduction = "integrated.rpca.rna",
  verbose = FALSE
)

print("### HARMONY INTEGRATION ###")
scrna.list <- IntegrateLayers(
  object = scrna.list, method = HarmonyIntegration,
  orig.reduction = "pca.unintegrated.rna",
  verbose = FALSE
)
# Harmony reduction was saved as harmony (no adaptations possible)
# Move and delete
scrna.list@reductions[["harmony.rna"]]<-scrna.list@reductions[["harmony"]]
scrna.list@reductions[["harmony"]]<-NULL

print("### DIMENSIONALITY REDUCTION ON CCA INTEGRATED DATA ###")
scrna.list <- FindNeighbors(scrna.list, reduction = "integrated.cca.rna", dims = 1:30, graph.name = "CCA_Int.rna")

# Store multiple resolutions for ClusTree (under or over clustering?)
for (val_res in seq(0.5,2.5, by=0.5)){
  ResVal<-paste0("cca_clusters_rna_res_", val_res)
  print(paste0("Busy working on CCA clustering with resolution ", val_res))
  scrna.list<-FindClusters(scrna.list,resolution=val_res,verbose = FALSE, cluster.name = ResVal,
	       		   graph.name = "CCA_Int.rna")}

scrna.list <- RunUMAP(scrna.list, reduction = "integrated.cca.rna", dims = 1:30, reduction.name = "umap.cca.rna")

print("### DIMENSIONALITY REDUCTION ON HARMONY INTEGRATED DATA ###")
scrna.list <- FindNeighbors(scrna.list, reduction = "harmony.rna", dims = 1:30, graph.name = "Harmony_Int.rna")

# Store multiple  resolutions for ClusTree (under of over clustering?)
for (val_res in seq(0.5,2.5, by=0.5)){
  ResVal<-paste0("harmony_clusters_rna_res_", val_res)
  print(paste0("Busy working on Harmony clustering with resolution ", val_res))
  scrna.list<-FindClusters(scrna.list,resolution=val_res,verbose = FALSE, cluster.name = ResVal,
               		   graph.name = "Harmony_Int.rna")}

scrna.list <- RunUMAP(scrna.list, reduction = "harmony.rna", dims = 1:30, reduction.name = "umap.harmony.rna")

print("### DIMENSIONALITY REDUCTION ON RPCA INTEGRATED DATA ###")
scrna.list <- FindNeighbors(scrna.list, reduction = "integrated.rpca.rna", dims = 1:30, graph.name = "RPCA_Int.rna")

# Store multiple resolution for ClusTree(under or over clustering?)
for (val_res in seq(0.5,2.5, by=0.5)){
  ResVal<-paste0("rpca_clusters_rna_res_", val_res)
  print(paste0("Busy working on RPCA clustering with resolution ", val_res))
  scrna.list<-FindClusters(scrna.list,resolution=val_res,verbose = FALSE, cluster.name = ResVal,
               		   graph.name = "RPCA_Int.rna")}

scrna.list <- RunUMAP(scrna.list, reduction = "integrated.rpca.rna", dims = 1:30, reduction.name = "umap.rpca.rna")

print("### JOINING LAYERS ###")
scrna.list <- JoinLayers(scrna.list)

print("### ADT ###")
scrna.list[["ADT"]] <- split(MergeData[["ADT"]], f = MergeData@meta.data[[Args$SplitArg]])

DefaultAssay(scrna.list) <- 'ADT'

print("### NORMALIZATION ###")
scrna.list <- Seurat::NormalizeData(scrna.list, normalization.method = "CLR", margin = 2)

print("### VARIABLE FEATURE SELECTION ###")
scrna.list <- Seurat::FindVariableFeatures(scrna.list)

print ("### SCALING ###")
scrna.list <- Seurat::ScaleData(scrna.list, features = rownames(scrna.list))

print("### DIMENSIONALITY REDUCTION ###")
scrna.list <- Seurat::RunPCA(object = scrna.list, assay="ADT", reduction.name="pca.unintegrated.adt")

scrna.list <- FindNeighbors(scrna.list, dims = 1:30, reduction = "pca.unintegrated.adt")

scrna.list <- FindClusters(scrna.list, resolution = 2, cluster.name = "unintegrated_clusters.adt")

scrna.list <- RunUMAP(scrna.list, dims = 1:30, reduction = "pca.unintegrated.adt", reduction.name = "umap.unintegrated.adt")

print("### CCA INTEGRATION ###")
scrna.list <- IntegrateLayers(
  object = scrna.list, method = CCAIntegration,
  orig.reduction = "pca.unintegrated.adt", new.reduction = "integrated.cca.adt",
  verbose = FALSE
)

print("### RPCA INTEGRATION ###")
scrna.list <- IntegrateLayers(
  object = scrna.list, method = RPCAIntegration,
  orig.reduction = "pca.unintegrated.adt", new.reduction = "integrated.rpca.adt",
  verbose = FALSE
)

print("### HARMONY INTEGRATION ###")
scrna.list <- IntegrateLayers(
  object = scrna.list, method = HarmonyIntegration,
  orig.reduction = "pca.unintegrated.adt",
  verbose = FALSE
)
# Harmony reduction was saved as harmony (no adaptations possible)
# Move and delete
scrna.list@reductions[["harmony.adt"]]<-scrna.list@reductions[["harmony"]]
scrna.list@reductions[["harmony"]]<-NULL

print("### DIMENSIONALITY REDUCTION ON RPCA INTEGRATED DATA ###")
scrna.list <- FindNeighbors(scrna.list, reduction = "integrated.rpca.adt", dims = 1:30, graph.name = "RPCA_Int.adt")
for (val_res in seq(0.5,2.5, by=0.5)){
  ResVal<-paste0("rpca_clusters_adt_res_", val_res)
  print(paste0("Busy working on RPCA clustering with resolution ", val_res))
  scrna.list<-FindClusters(scrna.list,resolution=val_res,verbose = FALSE, cluster.name = ResVal,
              		   graph.name = "RPCA_Int.adt")}

#scrna.list <- FindClusters(scrna.list, resolution = 2, cluster.name = "rpca_adt_clusters")

scrna.list <- RunUMAP(scrna.list, reduction = "integrated.rpca.adt", dims = 1:30, reduction.name = "umap.rpca.adt")

print("### DIMENSIONALITY REDUCTION ON CCA INTEGRATED DATA ###")
scrna.list <- FindNeighbors(scrna.list, reduction = "integrated.cca.adt", dims = 1:30, graph.name = "CCA_Int.adt")
for (val_res in seq(0.5,2.5, by=0.5)){
  ResVal<-paste0("cca_clusters_adt_res_", val_res)
  print(paste0("Busy working on CCA clustering with resolution ", val_res))
  scrna.list<-FindClusters(scrna.list,resolution=val_res,verbose = FALSE, cluster.name = ResVal,
              		   graph.name = "CCA_Int.adt")}

#scrna.list <- FindClusters(scrna.list, resolution = 2, cluster.name = "cca_adt_clusters")

scrna.list <- RunUMAP(scrna.list, reduction = "integrated.cca.adt", dims = 1:30, reduction.name = "umap.cca.adt")

print("### DIMENSIONALITY REDUCTION ON HARMONY INTEGRATED DATA ###")
scrna.list <- FindNeighbors(scrna.list, reduction = "harmony.adt", dims = 1:30, graph.name = "Harmony_Int.adt")
for (val_res in seq(0.5,2.5, by=0.5)){
  ResVal<-paste0("harmony_clusters_adt_res_", val_res)
  print(paste0("Busy working on Harmony clustering with resolution ", val_res))
  scrna.list<-FindClusters(scrna.list,resolution=val_res,verbose = FALSE, cluster.name = ResVal,
              		   graph.name = "Harmony_Int.adt")}

#scrna.list <- FindClusters(scrna.list, resolution = 2, cluster.name = "harmony_adt_clusters")

scrna.list <- RunUMAP(scrna.list, reduction = "harmony.adt", dims = 1:30, reduction.name = "umap.harmony.adt")

print("### JOINING LAYERS ###")
scrna.list <- JoinLayers(scrna.list)

saveRDS(scrna.list, file = save_Int)
