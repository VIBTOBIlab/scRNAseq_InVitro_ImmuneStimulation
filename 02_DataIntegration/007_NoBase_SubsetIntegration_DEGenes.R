##### EXTRACT THE INPUT ARGUMENTS #####
library(R.utils)
Args<-commandArgs(trailingOnly = TRUE, asValues = TRUE)
string<-"##################\nInput Arguments:\n##################\n"
for (i in names(Args)){
  string<-paste0(string, i, " : ", Args[i], "\n")}
string<-paste0(string,"#################\n")
cat(string)

# Correct path for libraries
.libPaths(Args$LibFolder)

##### LIBRARY LOADING #####
message("LIBRARY LOADING")
library(Seurat)
library(clustree)
# General data manipulation
library(tidyr)
library(dplyr)
library(tibble)
# Plot generation
library(ggplot2)
library(ggpubr)
library(xlsx)
message("THEME LOADING")
# Themes
source("/kyukon/data/gent/vo/000/gvo00027/TOBI/CVID/Scripts/PlotThemes_sc_CVID.R")
colorsStim<-c(colorsStim, c("CD3/CD28" = colorsStim[["CD3/CD28_beads"]]))
# Functions
source("/kyukon/data/gent/vo/000/gvo00027/TOBI/CVID/Scripts/Functions_sc_CVID.R")

message("LOAD THE INTEGRATED DATA (all cells)")
IntData<-readRDS(file = Args$AllCells_Annotated)

message("EXTRACT CELL SUBSETS (8h and 30h)")
IntData.B.8h<-subset(IntData, TimePoint == "8" & NoDbl_ManAnnot_L1_RNA == "B")
IntData.B.30h<-subset(IntData, TimePoint == "30" & NoDbl_ManAnnot_L1_RNA == "B")
IntData.T_NK.8h<-subset(IntData, TimePoint == "8" & NoDbl_ManAnnot_L1_RNA != "B" & NoDbl_ManAnnot_L1_RNA != "Myeloid")
IntData.T_NK.30h<-subset(IntData, TimePoint == "30" & NoDbl_ManAnnot_L1_RNA != "B" & NoDbl_ManAnnot_L1_RNA != "Myeloid")
rm(IntData)

##### B CELLS 8H #####
message("B CELLS \nINTEGRATION - 8H")
DefaultAssay(IntData.B.8h)<-"RNA"
# Split the data according to sampleID (CEV number) --> only RNA integration
# IntData.B[["RNA"]]<-split(IntData.B[["RNA"]], f = IntData.B$samplename)
IntData.B.8h[["RNA"]]<-split(IntData.B.8h[["RNA"]], f = IntData.B.8h$samplename)
# Normalization
IntData.B.8h<-Seurat::NormalizeData(IntData.B.8h)
# Variable Features
IntData.B.8h<-Seurat::FindVariableFeatures(IntData.B.8h)
# Scaling
IntData.B.8h<-Seurat::ScaleData(IntData.B.8h, features = rownames(IntData.B.8h))
# Dimensionality reduction (unintegrated)
IntData.B.8h<-Seurat::RunPCA(IntData.B.8h, assay = "RNA", reduction.name = "B_pca.unintegrated.rna")
IntData.B.8h<-Seurat::FindNeighbors(IntData.B.8h, dims = 1:20, reduction = "B_pca.unintegrated.rna")
IntData.B.8h<-Seurat::FindClusters(IntData.B.8h, resolution = 2, cluster.name = "B_pca_unintegrated_clusters.rna")
IntData.B.8h<-Seurat::RunUMAP(IntData.B.8h, dims = 1:20, reduction = "B_pca.unintegrated.rna", 
                           reduction.name = "B_umap.unintegrated.rna")
# CCA integration of RNA layer
IntData.B.8h<-IntegrateLayers(IntData.B.8h, method = CCAIntegration,
                           orig.reduction = "B_pca.unintegrated.rna", new.reduction = "B_integrated.cca.rna",
                           verbose = F)
# Joining layers
IntData.B.8h[["RNA"]]<-JoinLayers(IntData.B.8h[["RNA"]])

##### CLUSTERING OF B CELL SUBSET #####
message("CLUSTERING OF B CELLS - 8H")
# Nearest neighbours graph generation based on CCA (only RNA) integration method
IntData.B.8h <- FindNeighbors(IntData.B.8h, 
                           reduction = "B_integrated.cca.rna", 
                           dims = 1:30, graph.name = "B_integrated.cca.rna_graph")
# Clustering using multiple resolution values
for (val_res in seq(0.5,2.5, by=0.5)){
  ResVal<-paste0("Singlets_B_clusters_res_", val_res)
  message(paste0("Busy working on B cell clustering with resolution ", val_res))
  IntData.B.8h<-FindClusters(IntData.B.8h,resolution=val_res,
                          verbose = FALSE, cluster.name = ResVal,
                          graph.name = "B_integrated.cca.rna_graph")}
# Clustree to check the select the appropriate resolution
clustree(IntData.B.8h, 
         # Select only the correct clusters based on the resolution
         prefix = "Singlets_B_clusters_res_",
         # Size of the nodes
         node_size = 6,
         # Edges cross as little a possible
         layout = "sugiyama",
         # Coloring the nodes based on fraction of B cells
         # Per cluster 
         node_colour = "lightblue",
         show_axis = TRUE)+
  scale_color_gradient(low = "lightgrey",
                       high = "#A00D0D")+
  scale_edge_color_continuous(low = "lightgrey",
                              high = "#690F82")+
  labs(edge_alpha = "Proportion moved cells",
       edge_colour = "Counts moved cells",
       colour = "B cell fraction",
       title = "Singlets B cells",
       subtitle = "RNA")+
  ylab("Resolution")+
  theme_clustree
ggsave(filename = paste("NoBase_Integrated",
                        "ClusTree_BCells_8h_CCA.jpeg",
                        sep = "_"),
       path = file.path(Args$PathResults, "ClusTrees"),
       bg = "white", height = 6, width = 15, limitsize = FALSE)
# Resolution of 2 seems okay 
# Better to overcluster, merging clusters is possible if necessary
Res_Method<-c("Resolution" = 2)
Stats_Subset<-tibble(Description = c("Final integrated B cell subset - 8h", "B cell cluster resolution - 8h"),
                     Value = c(dim(IntData.B.8h)[2], Res_Method[["Resolution"]]),
                     Unit = c("Cells", ""))
Res_Method[]<-paste0("Singlets_B_clusters_res_",Res_Method)

IntData.B.8h<-SelectMetaData(object = IntData.B.8h,
                          to_stay = c("B_clusters" = Res_Method[["Resolution"]]),
                          pattern_remove = "Singlets_B_clusters_res_")

IntData.B.8h <- RunUMAP(IntData.B.8h, 
                     reduction = "B_integrated.cca.rna", 
                     dims = 1:30, reduction.name = "B_integrated.umap.rna")

message("SAVE CLUSTER SIZES - 8H")
FracClus<-table(IntData.B.8h$B_clusters)/(dim(IntData.B.8h)[2])
Stats_Subset<-rbind(Stats_Subset,
                    tibble(Description = "Number of RNA CCA clusters in B - 8h",
                           Value = FracClus%>%rownames()%>%as.numeric()%>%max()+1,
                           Unit = "Clusters"))
Stats_Subset<-rbind(Stats_Subset,
                    tibble("Description" = paste0("Fraction B cell cluster ", c(0:Stats_Subset[Stats_Subset["Description"]=="Number of RNA CCA clusters in B - 8h",][["Value"]]), " - 8h"),
                           "Value" = FracClus[c(0:Stats_Subset[Stats_Subset["Description"]=="Number of RNA CCA clusters in B - 8h",]
                                                [["Value"]])%>%as.character()]*100,
                           "Unit" = "%"))

Stats_Subset<-rbind(Stats_Subset,
                    tibble("Description" = paste0("Cells in B cell cluster ", c(0:Stats_Subset[Stats_Subset["Description"]=="Number of RNA CCA clusters in B - 8h",][["Value"]]), " - 8h"),
                           "Value" = table(IntData.B.8h$B_clusters)[c(0:Stats_Subset[Stats_Subset["Description"]=="Number of RNA CCA clusters in B - 8h",]
                                                                 [["Value"]])%>%as.character()],
                           "Unit" = "Cells"))
rm(FracClus)

message("VISUALISATION - 8H")
##### VISUALISATION
ggarrange(plotlist = list(
  # Dimensionality reduction plot with cluster labels (RNA only)
  DimPlot(object = IntData.B.8h,
          group.by = "B_clusters",
          reduction = "B_integrated.umap.rna", 
          label = TRUE)+
    labs(title = "B cell clusters",
         subtitle = "New clustering, only B cells")+
    xlab("UMAP1")+
    ylab("UMAP2")+
    themeG+
    NoLegend(),
  # Dimensionality reduction plot with previous cluster labels (RNA only)
  DimPlot(object = IntData.B.8h,
          group.by = "cca_clusters_rna",
          reduction = "B_integrated.umap.rna", 
          label = TRUE)+
    labs(title = "B cell clusters",
         subtitle = "Old clustering, all cells")+
    xlab("UMAP1")+
    ylab("UMAP2")+
    themeG+
    NoLegend(),
  DimPlot(object = IntData.B.8h,
          group.by = "Stimulation",
          reduction = "B_integrated.umap.rna", 
          label = F)+
    labs(title = "B cell clusters",
         subtitle = "Stimulatory condition")+
    xlab("UMAP1")+
    ylab("UMAP2")+
    themeG+
    scale_color_manual(values = colorsStim)+
    theme(legend.justification=c(1,1.14),
          legend.position.inside=c(1,1.14))+
    guides(color=guide_legend(override.aes = list(size = 3),
                              ncol = 2)),
  DimPlot(object = IntData.B.8h,
          group.by = "TimePoint",
          reduction = "B_integrated.umap.rna", 
          label = F)+
    labs(title = "B cell clusters",
         subtitle = "Incubation time")+
    xlab("UMAP1")+
    ylab("UMAP2")+
    themeG+
    scale_color_manual(values = colorsTime)+
    theme(legend.justification=c(1,1.14),
          legend.position.inside = c(1,1.14))+
    guides(color=guide_legend(override.aes = list(size = 3),
                              ncol = 2))),
  ncol = 2, nrow = 2)

ggsave(filename = paste("NoBase_Integrated",
                        "UMAPs_IntegratedB_8h.jpeg",
                        sep  ="_"),
       path = file.path(Args$PathResults, "BCells"),
       device = "jpeg", bg = "white",
       width = 13, height = 13)

saveRDS(IntData.B.8h,
        file = file.path(Args$PathResults, 
                         "SeuratObjects", 
                         paste0("seuratObj_NoBase_Integrated_BClusters_8h.rds")))

message("DIFFERENTIAL GENE EXPRESSION")
# Identification of differentially expressed genes/proteins
Idents(IntData.B.8h)<-"B_clusters"
# RNA
DiffMarkers_AllG_B.8h <- FindAllMarkers(IntData.B.8h, 
                                     # Minimal gene fraction as condition to test
                                     min.pct = 0.10, 
                                     min.diff.pct = 0.30, 
                                     logfc.threshold = 0.30,
                                     # pvalue threshold to show gene
                                     return.thresh = 0.01, 
                                     only.pos = FALSE, 
                                     assay="RNA",
                                     test.use = "poisson",
                                     group.by = "B_clusters")

DiffMarkers_AllG_B.8h<-DiffMarkers_AllG_B.8h%>%
  mutate(cluster=factor(cluster,
                        # Factorization of the clusters
                        levels = levels(DiffMarkers_AllG_B.8h$cluster)%>%as.numeric()%>%sort(),
                        ordered=TRUE),
         # DE score
         Score = pct.1 / (pct.2 + 0.01) *avg_log2FC)%>%
  arrange(cluster, -Score)

# Save the top 5 genes/proteins in data frame
TopDE_All_B.8h<-data.frame()
for (i in levels(IntData.B.8h$B_clusters)){
  TopDE_All_B.8h<-rbind(TopDE_All_B.8h,
                     data.frame(cluster = i,
                                # Necessary for wider format
                                Feature = c(paste0("gene", c(1:10))),
                                Name = c((DiffMarkers_AllG_B.8h%>%filter(cluster == i))
                                         # Extract genes
                                         [1:10, "gene"])))
}
# To wider format
TopDE_All_B.8h<-TopDE_All_B.8h%>%
  pivot_wider(values_from = Name,
              names_from = Feature)%>%
  mutate(cluster=factor(x = cluster,
                        levels = cluster%>%
                          as.numeric()%>%sort(),
                        ordered = TRUE))%>%
  arrange(cluster)

# Save and clear workspace
# All markers
write.xlsx(x = DiffMarkers_AllG_B.8h, file = file.path(Args$PathResults, "BCells", 
                                                    paste("NoBase_Integrated", "DEGenes_B_8h_RNA.xlsx", 
                                                          sep = "_")),
           sheetName = "DE_Genes", col.names  = T, row.names = F)

# Top 10 markers
write.xlsx(x = as.data.frame(TopDE_All_B.8h),
           file = file.path(Args$PathResults, "BCells",
                            paste("NoBase_Integrated", "Top10_DEGenes_B_8h_RNA.xlsx",
                                  sep = "_")),
           sheetName = "TopDE_Genes", col.names = T, row.names = F, showNA = F)

# Clear workspace
rm(DiffMarkers_AllG_B.8h, TopDE_All_B.8h)

# Clusters based on combination of ADT and RNA (WNN)
# Calculate total number of cells in clusters
# Calculate fractions of each cell type in layer 1
AzimuthProp_B_L2.8h<-left_join(x = CalClusCount(object = IntData.B.8h,
                                             clusters = "B_clusters",
                                             feature = NULL, barplot = F)[[1]],
                            y = CalClusFrac(object = IntData.B.8h, 
                                            clusters = "B_clusters",
                                            feature = "predicted.celltype.l2",
                                            barplot = F)[[1]]%>%
                              mutate(Fraction = Fraction * 100)%>%
                              pivot_wider(names_from = predicted.celltype.l2,
                                          values_from = Fraction),
                            by = "B_clusters")%>%
  select(clusters = B_clusters, 
         `Number of cells` = Counts,
         everything())%>%
  arrange(clusters)

write.xlsx(x = AzimuthProp_B_L2.8h,
           file = file.path(Args$PathResults, "BCells",
                            paste("NoBase_Integrated", "Top10_DEGenes_B_8h_RNA.xlsx",
                                  sep = "_")),
           sheetName = "CellType_L2_Fractions", col.names = T, row.names = F, showNA = F,
           append = T)

rm(AzimuthProp_B_L2.8h, IntData.B.8h)

##### B CELLS 30H #####
message("INTEGRATION - 30H")
DefaultAssay(IntData.B.30h)<-"RNA"
# Split the data according to sampleID (CEV number) --> only RNA integration
# IntData.B[["RNA"]]<-split(IntData.B[["RNA"]], f = IntData.B$samplename)
IntData.B.30h[["RNA"]]<-split(IntData.B.30h[["RNA"]], f = IntData.B.30h$samplename)
# Normalization
IntData.B.30h<-Seurat::NormalizeData(IntData.B.30h)
# Variable Features
IntData.B.30h<-Seurat::FindVariableFeatures(IntData.B.30h)
# Scaling
IntData.B.30h<-Seurat::ScaleData(IntData.B.30h, features = rownames(IntData.B.30h))
# Dimensionality reduction (unintegrated)
IntData.B.30h<-Seurat::RunPCA(IntData.B.30h, assay = "RNA", reduction.name = "B_pca.unintegrated.rna")
IntData.B.30h<-Seurat::FindNeighbors(IntData.B.30h, dims = 1:20, reduction = "B_pca.unintegrated.rna")
IntData.B.30h<-Seurat::FindClusters(IntData.B.30h, resolution = 2, cluster.name = "B_pca_unintegrated_clusters.rna")
IntData.B.30h<-Seurat::RunUMAP(IntData.B.30h, dims = 1:20, reduction = "B_pca.unintegrated.rna", 
                              reduction.name = "B_umap.unintegrated.rna")
# CCA integration of RNA layer
IntData.B.30h<-IntegrateLayers(IntData.B.30h, method = CCAIntegration,
                              orig.reduction = "B_pca.unintegrated.rna", new.reduction = "B_integrated.cca.rna",
                              verbose = F)
# Joining layers
IntData.B.30h[["RNA"]]<-JoinLayers(IntData.B.30h[["RNA"]])

##### CLUSTERING OF B CELL SUBSET #####
message("CLUSTERING OF B CELLS - 30H")
# Nearest neighbours graph generation based on CCA (only RNA) integration method
IntData.B.30h <- FindNeighbors(IntData.B.30h, 
                              reduction = "B_integrated.cca.rna", 
                              dims = 1:30, graph.name = "B_integrated.cca.rna_graph")
# Clustering using multiple resolution values
for (val_res in seq(0.5,2.5, by=0.5)){
  ResVal<-paste0("Singlets_B_clusters_res_", val_res)
  message(paste0("Busy working on B cell clustering with resolution ", val_res))
  IntData.B.30h<-FindClusters(IntData.B.30h,resolution=val_res,
                             verbose = FALSE, cluster.name = ResVal,
                             graph.name = "B_integrated.cca.rna_graph")}
# Clustree to check the select the appropriate resolution
clustree(IntData.B.30h, 
         # Select only the correct clusters based on the resolution
         prefix = "Singlets_B_clusters_res_",
         # Size of the nodes
         node_size = 6,
         # Edges cross as little a possible
         layout = "sugiyama",
         # Coloring the nodes based on fraction of B cells
         # Per cluster 
         node_colour = "lightblue",
         show_axis = TRUE)+
  scale_color_gradient(low = "lightgrey",
                       high = "#A00D0D")+
  scale_edge_color_continuous(low = "lightgrey",
                              high = "#690F82")+
  labs(edge_alpha = "Proportion moved cells",
       edge_colour = "Counts moved cells",
       colour = "B cell fraction",
       title = "Singlets B cells",
       subtitle = "RNA")+
  ylab("Resolution")+
  theme_clustree
ggsave(filename = paste("NoBase_Integrated",
                        "ClusTree_BCells_30h_CCA.jpeg",
                        sep = "_"),
       path = file.path(Args$PathResults, "ClusTrees"),
       bg = "white", height = 6, width = 15, limitsize = FALSE)
# Resolution of 2 seems okay 
# Better to overcluster, merging clusters is possible if necessary
Res_Method<-c("Resolution" = 2)
Stats_Subset<-rbind(Stats_Subset,
                    tibble(Description = c("Final integrated B cell subset - 30h", "B cell cluster resolution - 30h"),
                            Value = c(dim(IntData.B.30h)[2], Res_Method[["Resolution"]]),
                            Unit = c("Cells", "")))
Res_Method[]<-paste0("Singlets_B_clusters_res_",Res_Method)

IntData.B.30h<-SelectMetaData(object = IntData.B.30h,
                             to_stay = c("B_clusters" = Res_Method[["Resolution"]]),
                             pattern_remove = "Singlets_B_clusters_res_")


IntData.B.30h <- RunUMAP(IntData.B.30h, 
                        reduction = "B_integrated.cca.rna", 
                        dims = 1:30, reduction.name = "B_integrated.umap.rna")

message("SAVE CLUSTER SIZES - 30H")
FracClus<-table(IntData.B.30h$B_clusters)/(dim(IntData.B.30h)[2])
Stats_Subset<-rbind(Stats_Subset,
                    tibble(Description = "Number of RNA CCA clusters in B - 30h",
                           Value = FracClus%>%rownames()%>%as.numeric()%>%max()+1,
                           Unit = "Clusters"))
Stats_Subset<-rbind(Stats_Subset,
             tibble("Description" = paste0("Fraction B cell cluster ", c(0:Stats_Subset[Stats_Subset["Description"]=="Number of RNA CCA clusters in B - 30h",][["Value"]]), " - 30h"),
                    "Value" = FracClus[c(0:Stats_Subset[Stats_Subset["Description"]=="Number of RNA CCA clusters in B - 30h",]
                                         [["Value"]])%>%as.character()]*100,
                    "Unit" = "%"))

Stats_Subset<-rbind(Stats_Subset,
             tibble("Description" = paste0("Cells in B cell cluster ", c(0:Stats_Subset[Stats_Subset["Description"]=="Number of RNA CCA clusters in B - 30h",][["Value"]]), " - 30h"),
                    "Value" = table(IntData.B.30h$B_clusters)[c(0:Stats_Subset[Stats_Subset["Description"]=="Number of RNA CCA clusters in B - 30h",]
                                                          [["Value"]])%>%as.character()],
                    "Unit" = "Cells"))
rm(FracClus)

message("VISUALISATION - 30H")
##### VISUALISATION
ggarrange(plotlist = list(
  # Dimensionality reduction plot with cluster labels (RNA only)
  DimPlot(object = IntData.B.30h,
          group.by = "B_clusters",
          reduction = "B_integrated.umap.rna", 
          label = TRUE)+
    labs(title = "B cell clusters",
         subtitle = "New clustering, only B cells")+
    xlab("UMAP1")+
    ylab("UMAP2")+
    themeG+
    NoLegend(),
  # Dimensionality reduction plot with previous cluster labels (RNA only)
  DimPlot(object = IntData.B.30h,
          group.by = "cca_clusters_rna",
          reduction = "B_integrated.umap.rna", 
          label = TRUE)+
    labs(title = "B cell clusters",
         subtitle = "Old clustering, all cells")+
    xlab("UMAP1")+
    ylab("UMAP2")+
    themeG+
    NoLegend(),
  DimPlot(object = IntData.B.30h,
          group.by = "Stimulation",
          reduction = "B_integrated.umap.rna", 
          label = F)+
    labs(title = "B cell clusters",
         subtitle = "Stimulatory condition")+
    xlab("UMAP1")+
    ylab("UMAP2")+
    themeG+
    scale_color_manual(values = colorsStim)+
    theme(legend.justification=c(1,1.14),
          legend.position.inside=c(1,1.14))+
    guides(color=guide_legend(override.aes = list(size = 3),
                              ncol = 2)),
  DimPlot(object = IntData.B.30h,
          group.by = "TimePoint",
          reduction = "B_integrated.umap.rna", 
          label = F)+
    labs(title = "B cell clusters",
         subtitle = "Incubation time")+
    xlab("UMAP1")+
    ylab("UMAP2")+
    themeG+
    scale_color_manual(values = colorsTime)+
    theme(legend.justification=c(1,1.14),
          legend.position.inside=c(1,1.14))+
    guides(color=guide_legend(override.aes = list(size = 3),
                              ncol = 2))),
  ncol = 2, nrow = 2)

ggsave(filename = paste("NoBase_Integrated",
                        "UMAPs_IntegratedB_30h.jpeg",
                        sep  ="_"),
       path = file.path(Args$PathResults, "BCells"),
       device = "jpeg", bg = "white",
       width = 13, height = 13)

saveRDS(IntData.B.30h,
        file = file.path(Args$PathResults, 
                         "SeuratObjects", 
                         paste0("seuratObj_NoBase_Integrated_BClusters_30h.rds")))

message("DIFFERENTIAL GENE EXPRESSION")
# Identification of differentially expressed genes/proteins
Idents(IntData.B.30h)<-"B_clusters"
# RNA
DiffMarkers_AllG_B.30h <- FindAllMarkers(IntData.B.30h, 
                                        # Minimal gene fraction as condition to test
                                        min.pct = 0.10, 
                                        min.diff.pct = 0.30, 
                                        logfc.threshold = 0.30,
                                        # pvalue threshold to show gene
                                        return.thresh = 0.01, 
                                        only.pos = FALSE, 
                                        assay="RNA",
                                        test.use = "poisson",
                                        group.by = "B_clusters")

DiffMarkers_AllG_B.30h<-DiffMarkers_AllG_B.30h%>%
  mutate(cluster=factor(cluster,
                        # Factorization of the clusters
                        levels = levels(DiffMarkers_AllG_B.30h$cluster)%>%as.numeric()%>%sort(),
                        ordered=TRUE),
         # DE score
         Score = pct.1 / (pct.2 + 0.01) *avg_log2FC)%>%
  arrange(cluster, -Score)

# Save the top 5 genes/proteins in data frame
TopDE_All_B.30h<-data.frame()
for (i in levels(IntData.B.30h$B_clusters)){
  TopDE_All_B.30h<-rbind(TopDE_All_B.30h,
                        data.frame(cluster = i,
                                   # Necessary for wider format
                                   Feature = c(paste0("gene", c(1:10))),
                                   Name = c((DiffMarkers_AllG_B.30h%>%filter(cluster == i))
                                            # Extract genes
                                            [1:10, "gene"])))
}
# To wider format
TopDE_All_B.30h<-TopDE_All_B.30h%>%
  pivot_wider(values_from = Name,
              names_from = Feature)%>%
  mutate(cluster=factor(x = cluster,
                        levels = cluster%>%
                          as.numeric()%>%sort(),
                        ordered = TRUE))%>%
  arrange(cluster)

# Save and clear workspace
# All markers
write.xlsx(x = DiffMarkers_AllG_B.30h, file = file.path(Args$PathResults, "BCells", 
                                                       paste("NoBase_Integrated", "DEGenes_B_30h_RNA.xlsx", 
                                                             sep = "_")),
           sheetName = "DE_Genes", col.names  = T, row.names = F)

# Top 10 markers
write.xlsx(x = as.data.frame(TopDE_All_B.30h),
           file = file.path(Args$PathResults, "BCells",
                            paste("NoBase_Integrated", "Top10_DEGenes_B_30h_RNA.xlsx",
                                  sep = "_")),
           sheetName = "TopDE_Genes", col.names = T, row.names = F, showNA = F)

# Clear workspace
rm(DiffMarkers_AllG_B.30h, TopDE_All_B.30h)

# Clusters based on combination of ADT and RNA (WNN)
# Calculate total number of cells in clusters
# Calculate fractions of each cell type in layer 1
AzimuthProp_B_L2.30h<-left_join(x = CalClusCount(object = IntData.B.30h,
                                                clusters = "B_clusters",
                                                feature = NULL, barplot = F)[[1]],
                               y = CalClusFrac(object = IntData.B.30h, 
                                               clusters = "B_clusters",
                                               feature = "predicted.celltype.l2",
                                               barplot = F)[[1]]%>%
                                 mutate(Fraction = Fraction * 100)%>%
                                 pivot_wider(names_from = predicted.celltype.l2,
                                             values_from = Fraction),
                               by = "B_clusters")%>%
  select(clusters = B_clusters, 
         `Number of cells` = Counts,
         everything())%>%
  arrange(clusters)

write.xlsx(x = AzimuthProp_B_L2.30h,
           file = file.path(Args$PathResults, "BCells",
                            paste("NoBase_Integrated", "Top10_DEGenes_B_30h_RNA.xlsx",
                                  sep = "_")),
           sheetName = "CellType_L2_Fractions", col.names = T, row.names = F, showNA = F,
           append = T)

rm(AzimuthProp_B_L2.30h, IntData.B.30h)

##### T AND NK CELLS 8H #####
message("T AND NK CELLS\nINTEGRATION - 8H")
DefaultAssay(IntData.T_NK.8h)<-"RNA"
# Split the data according to sampleID (CEV number) --> only RNA integration
# IntData.T_NK[["RNA"]]<-split(IntData.T_NK[["RNA"]], f = IntData.T_NK$samplename)
IntData.T_NK.8h[["RNA"]]<-split(IntData.T_NK.8h[["RNA"]], f = IntData.T_NK.8h$samplename)
# Normalization
IntData.T_NK.8h<-Seurat::NormalizeData(IntData.T_NK.8h)
# Variable Features
IntData.T_NK.8h<-Seurat::FindVariableFeatures(IntData.T_NK.8h)
# Scaling
IntData.T_NK.8h<-Seurat::ScaleData(IntData.T_NK.8h, features = rownames(IntData.T_NK.8h))
# Dimensionality reduction (unintegrated)
IntData.T_NK.8h<-Seurat::RunPCA(IntData.T_NK.8h, assay = "RNA", reduction.name = "T_NK_pca.unintegrated.rna")
IntData.T_NK.8h<-Seurat::FindNeighbors(IntData.T_NK.8h, dims = 1:20, reduction = "T_NK_pca.unintegrated.rna")
IntData.T_NK.8h<-Seurat::FindClusters(IntData.T_NK.8h, resolution = 2, cluster.name = "T_NK_pca_unintegrated_clusters.rna")
IntData.T_NK.8h<-Seurat::RunUMAP(IntData.T_NK.8h, dims = 1:20, reduction = "T_NK_pca.unintegrated.rna", 
                              reduction.name = "T_NK_umap.unintegrated.rna")
# CCA integration of RNA layer
IntData.T_NK.8h<-IntegrateLayers(IntData.T_NK.8h, method = CCAIntegration,
                              orig.reduction = "T_NK_pca.unintegrated.rna", new.reduction = "T_NK_integrated.cca.rna",
                              verbose = F)
# Joining layers
IntData.T_NK.8h[["RNA"]]<-JoinLayers(IntData.T_NK.8h[["RNA"]])

##### CLUSTERING OF T_NK CELL SUBSET #####
message("CLUSTERING OF T AND NK CELLS - 8H")
# Nearest neighbours graph generation based on CCA (only RNA) integration method
IntData.T_NK.8h <- FindNeighbors(IntData.T_NK.8h, 
                              reduction = "T_NK_integrated.cca.rna", 
                              dims = 1:30, graph.name = "T_NK_integrated.cca.rna_graph")
# Clustering using multiple resolution values
for (val_res in seq(0.5,2.5, by=0.5)){
  ResVal<-paste0("Singlets_T_NK_clusters_res_", val_res)
  message(paste0("Busy working on T and NK cell clustering with resolution ", val_res))
  IntData.T_NK.8h<-FindClusters(IntData.T_NK.8h,resolution=val_res,
                             verbose = FALSE, cluster.name = ResVal,
                             graph.name = "T_NK_integrated.cca.rna_graph")}
# Clustree to check the select the appropriate resolution
clustree(IntData.T_NK.8h, 
         # Select only the correct clusters based on the resolution
         prefix = "Singlets_T_NK_clusters_res_",
         # Size of the nodes
         node_size = 6,
         # Edges cross as little a possible
         layout = "sugiyama",
         # Coloring the nodes based on fraction of T and NK cells
         # Per cluster 
         node_colour = "lightblue",
         show_axis = TRUE)+
  scale_color_gradient(low = "lightgrey",
                       high = "#A00D0D")+
  scale_edge_color_continuous(low = "lightgrey",
                              high = "#690F82")+
  labs(edge_alpha = "Proportion moved cells",
       edge_colour = "Counts moved cells",
       colour = "T and NK cell fraction",
       title = "Singlets T and NK cells",
       subtitle = "RNA")+
  ylab("Resolution")+
  theme_clustree
ggsave(filename = paste("NoBase_Integrated",
                        "ClusTree_T_NKCells_8h_CCA.jpeg",
                        sep = "_"),
       path = file.path(Args$PathResults, "ClusTrees"),
       bg = "white", height = 6, width = 15, limitsize = FALSE)
# Resolution of 2 seems okay 
# Better to overcluster, merging clusters is possible if necessary
Res_Method<-c("Resolution" = 2)
Stats_Subset<-rbind(Stats_Subset,
                    tibble(Description = c("Final integrated T and NK cell subset - 8h", "T and NK cell cluster resolution - 8h"),
                            Value = c(dim(IntData.T_NK.8h)[2], Res_Method[["Resolution"]]),
                            Unit = c("Cells", "")))
Res_Method[]<-paste0("Singlets_T_NK_clusters_res_",Res_Method)

IntData.T_NK.8h<-SelectMetaData(object = IntData.T_NK.8h,
                             to_stay = c("T_NK_clusters" = Res_Method[["Resolution"]]),
                             pattern_remove = "Singlets_T_NK_clusters_res_")


IntData.T_NK.8h <- RunUMAP(IntData.T_NK.8h, 
                        reduction = "T_NK_integrated.cca.rna", 
                        dims = 1:30, reduction.name = "T_NK_integrated.umap.rna")

message("SAVE CLUSTER SIZES - 8H")
FracClus<-table(IntData.T_NK.8h$T_NK_clusters)/(dim(IntData.T_NK.8h)[2])
Stats_Subset<-rbind(Stats_Subset,
                    tibble(Description = "Number of RNA CCA clusters in T/NK - 8h",
                           Value = FracClus%>%rownames()%>%as.numeric()%>%max()+1,
                           Unit = "Clusters"))
Stats_Subset<-rbind(Stats_Subset,
                    tibble(Description = paste0("Fraction T/NK cell cluster ", c(0:Stats_Subset[Stats_Subset["Description"]=="Number of RNA CCA clusters in T/NK - 8h",][["Value"]]), " - 8h"),
                           Value = FracClus[c(0:Stats_Subset[Stats_Subset["Description"]=="Number of RNA CCA clusters in T/NK - 8h",]
                                                [["Value"]])%>%as.character()]*100,
                           Unit = "%"))

Stats_Subset<-rbind(Stats_Subset,
                    tibble(Description = paste0("Cells in B cell cluster ", c(0:Stats_Subset[Stats_Subset["Description"]=="Number of RNA CCA clusters in T/NK - 8h",][["Value"]]), " - 8h"),
                           Value = table(IntData.T_NK.8h$T_NK_clusters)[c(0:Stats_Subset[Stats_Subset["Description"]=="Number of RNA CCA clusters in T/NK - 8h",]
                                                                 [["Value"]])%>%as.character()],
                           Unit = "Cells"))
rm(FracClus)

message("VISUALISATION - 8H")
##### VISUALISATION
ggarrange(plotlist = list(
  # Dimensionality reduction plot with cluster labels (RNA only)
  DimPlot(object = IntData.T_NK.8h,
          group.by = "T_NK_clusters",
          reduction = "T_NK_integrated.umap.rna", 
          label = TRUE)+
    labs(title = "T and NK cell clusters",
         subtitle = "New clustering, only T and NK cells")+
    xlab("UMAP1")+
    ylab("UMAP2")+
    themeG+
    NoLegend(),
  # Dimensionality reduction plot with previous cluster labels (RNA only)
  DimPlot(object = IntData.T_NK.8h,
          group.by = "cca_clusters_rna",
          reduction = "T_NK_integrated.umap.rna", 
          label = TRUE)+
    labs(title = "T and NK cell clusters",
         subtitle = "Old clustering, all cells")+
    xlab("UMAP1")+
    ylab("UMAP2")+
    themeG+
    NoLegend(),
  DimPlot(object = IntData.T_NK.8h,
          group.by = "Stimulation",
          reduction = "T_NK_integrated.umap.rna", 
          label = F)+
    labs(title = "T and NK cell clusters",
         subtitle = "Stimulatory condition")+
    xlab("UMAP1")+
    ylab("UMAP2")+
    themeG+
    scale_color_manual(values = colorsStim)+
    theme(legend.justification=c(1,1.14),
          legend.position.inside=c(1,1.14))+
    guides(color=guide_legend(override.aes = list(size = 3),
                              ncol = 2)),
  DimPlot(object = IntData.T_NK.8h,
          group.by = "TimePoint",
          reduction = "T_NK_integrated.umap.rna", 
          label = F)+
    labs(title = "T and NK cell clusters",
         subtitle = "Incubation time")+
    xlab("UMAP1")+
    ylab("UMAP2")+
    themeG+
    scale_color_manual(values = colorsTime)+
    theme(legend.justification=c(1,1.14),
          legend.position.inside=c(1,1.14))+
    guides(color=guide_legend(override.aes = list(size = 3),
                              ncol = 2))),
  ncol = 2, nrow = 2)

ggsave(filename = paste("NoBase_Integrated",
                        "UMAPs_IntegratedT_NK_8h.jpeg",
                        sep  ="_"),
       path = file.path(Args$PathResults, "T_NKCells"),
       device = "jpeg", bg = "white",
       width = 13, height = 13)

saveRDS(IntData.T_NK.8h,
        file = file.path(Args$PathResults, 
                         "SeuratObjects", 
                         paste0("seuratObj_NoBase_Integrated_T_NKClusters_8h.rds")))

message("DIFFERENTIAL GENE EXPRESSION")
# Identification of differentially expressed genes/proteins
Idents(IntData.T_NK.8h)<-"T_NK_clusters"
# RNA
DiffMarkers_AllG_T_NK.8h <- FindAllMarkers(IntData.T_NK.8h, 
                                        # Minimal gene fraction as condition to test
                                        min.pct = 0.10, 
                                        min.diff.pct = 0.30, 
                                        logfc.threshold = 0.30,
                                        # pvalue threshold to show gene
                                        return.thresh = 0.01, 
                                        only.pos = FALSE, 
                                        assay="RNA",
                                        test.use = "poisson",
                                        group.by = "T_NK_clusters")

DiffMarkers_AllG_T_NK.8h<-DiffMarkers_AllG_T_NK.8h%>%
  mutate(cluster=factor(cluster,
                        # Factorization of the clusters
                        levels = levels(DiffMarkers_AllG_T_NK.8h$cluster)%>%as.numeric()%>%sort(),
                        ordered=TRUE),
         # DE score
         Score = pct.1 / (pct.2 + 0.01) *avg_log2FC)%>%
  arrange(cluster, -Score)

# Save the top 5 genes/proteins in data frame
TopDE_All_T_NK.8h<-data.frame()
for (i in levels(IntData.T_NK.8h$T_NK_clusters)){
  TopDE_All_T_NK.8h<-rbind(TopDE_All_T_NK.8h,
                        data.frame(cluster = i,
                                   # Necessary for wider format
                                   Feature = c(paste0("gene", c(1:10))),
                                   Name = c((DiffMarkers_AllG_T_NK.8h%>%filter(cluster == i))
                                            # Extract genes
                                            [1:10, "gene"])))
}
# To wider format
TopDE_All_T_NK.8h<-TopDE_All_T_NK.8h%>%
  pivot_wider(values_from = Name,
              names_from = Feature)%>%
  mutate(cluster=factor(x = cluster,
                        levels = cluster%>%
                          as.numeric()%>%sort(),
                        ordered = TRUE))%>%
  arrange(cluster)

# Save and clear workspace
# All markers
write.xlsx(x = DiffMarkers_AllG_T_NK.8h, file = file.path(Args$PathResults, "T_NKCells", 
                                                       paste("NoBase_Integrated", "DEGenes_T_NK_8h_RNA.xlsx", 
                                                             sep = "_")),
           sheetName = "DE_Genes", col.names  = T, row.names = F)

# Top 10 markers
write.xlsx(x = as.data.frame(TopDE_All_T_NK.8h),
           file = file.path(Args$PathResults, "T_NKCells",
                            paste("NoBase_Integrated", "Top10_DEGenes_T_NK_8h_RNA.xlsx",
                                  sep = "_")),
           sheetName = "TopDE_Genes", col.names = T, row.names = F, showNA = F)

# Clear workspace
rm(DiffMarkers_AllG_T_NK.8h, TopDE_All_T_NK.8h)

# Clusters based on combination of ADT and RNA (WNN)
# Calculate total number of cells in clusters
# Calculate fractions of each cell type in layer 1
AzimuthProp_T_NK_L2.8h<-left_join(x = CalClusCount(object = IntData.T_NK.8h,
                                                clusters = "T_NK_clusters",
                                                feature = NULL, barplot = F)[[1]],
                               y = CalClusFrac(object = IntData.T_NK.8h, 
                                               clusters = "T_NK_clusters",
                                               feature = "predicted.celltype.l2",
                                               barplot = F)[[1]]%>%
                                 mutate(Fraction = Fraction * 100)%>%
                                 pivot_wider(names_from = predicted.celltype.l2,
                                             values_from = Fraction),
                               by = "T_NK_clusters")%>%
  select(clusters = T_NK_clusters, 
         `Number of cells` = Counts,
         everything())%>%
  arrange(clusters)

write.xlsx(x = AzimuthProp_T_NK_L2.8h,
           file = file.path(Args$PathResults, "T_NKCells",
                            paste("NoBase_Integrated", "Top10_DEGenes_T_NK_8h_RNA.xlsx",
                                  sep = "_")),
           sheetName = "CellType_L2_Fractions", col.names = T, row.names = F, showNA = F,
           append = T)

rm(AzimuthProp_T_NK_L2.8h, IntData.T_NK.8h)

##### T_NK CELLS 30H #####
message("INTEGRATION - 30H")
DefaultAssay(IntData.T_NK.30h)<-"RNA"
# Split the data according to sampleID (CEV number) --> only RNA integration
# IntData.T_NK[["RNA"]]<-split(IntData.T_NK[["RNA"]], f = IntData.T_NK$samplename)
IntData.T_NK.30h[["RNA"]]<-split(IntData.T_NK.30h[["RNA"]], f = IntData.T_NK.30h$samplename)
# Normalization
IntData.T_NK.30h<-Seurat::NormalizeData(IntData.T_NK.30h)
# Variable Features
IntData.T_NK.30h<-Seurat::FindVariableFeatures(IntData.T_NK.30h)
# Scaling
IntData.T_NK.30h<-Seurat::ScaleData(IntData.T_NK.30h, features = rownames(IntData.T_NK.30h))
# Dimensionality reduction (unintegrated)
IntData.T_NK.30h<-Seurat::RunPCA(IntData.T_NK.30h, assay = "RNA", reduction.name = "T_NK_pca.unintegrated.rna")
IntData.T_NK.30h<-Seurat::FindNeighbors(IntData.T_NK.30h, dims = 1:20, reduction = "T_NK_pca.unintegrated.rna")
IntData.T_NK.30h<-Seurat::FindClusters(IntData.T_NK.30h, resolution = 2, cluster.name = "T_NK_pca_unintegrated_clusters.rna")
IntData.T_NK.30h<-Seurat::RunUMAP(IntData.T_NK.30h, dims = 1:20, reduction = "T_NK_pca.unintegrated.rna", 
                               reduction.name = "T_NK_umap.unintegrated.rna")
# CCA integration of RNA layer
IntData.T_NK.30h<-IntegrateLayers(IntData.T_NK.30h, method = CCAIntegration,
                               orig.reduction = "T_NK_pca.unintegrated.rna", new.reduction = "T_NK_integrated.cca.rna",
                               verbose = F)
# Joining layers
IntData.T_NK.30h[["RNA"]]<-JoinLayers(IntData.T_NK.30h[["RNA"]])

##### CLUSTERING OF T_NK CELL SUBSET #####
message("CLUSTERING OF T_NK CELLS - 30H")
# Nearest neighbours graph generation based on CCA (only RNA) integration method
IntData.T_NK.30h <- FindNeighbors(IntData.T_NK.30h, 
                               reduction = "T_NK_integrated.cca.rna", 
                               dims = 1:30, graph.name = "T_NK_integrated.cca.rna_graph")
# Clustering using multiple resolution values
for (val_res in seq(0.5,2.5, by=0.5)){
  ResVal<-paste0("Singlets_T_NK_clusters_res_", val_res)
  message(paste0("Busy working on T and NK cell clustering with resolution ", val_res))
  IntData.T_NK.30h<-FindClusters(IntData.T_NK.30h,resolution=val_res,
                              verbose = FALSE, cluster.name = ResVal,
                              graph.name = "T_NK_integrated.cca.rna_graph")}
# Clustree to check the select the appropriate resolution
clustree(IntData.T_NK.30h, 
         # Select only the correct clusters based on the resolution
         prefix = "Singlets_T_NK_clusters_res_",
         # Size of the nodes
         node_size = 6,
         # Edges cross as little a possible
         layout = "sugiyama",
         # Coloring the nodes based on fraction of T_NK cells
         # Per cluster 
         node_colour = "lightblue",
         show_axis = TRUE)+
  scale_color_gradient(low = "lightgrey",
                       high = "#A00D0D")+
  scale_edge_color_continuous(low = "lightgrey",
                              high = "#690F82")+
  labs(edge_alpha = "Proportion moved cells",
       edge_colour = "Counts moved cells",
       colour = "T and NK cell fraction",
       title = "Singlets T and NK cells",
       subtitle = "RNA")+
  ylab("Resolution")+
  theme_clustree
ggsave(filename = paste("NoBase_Integrated",
                        "ClusTree_T_NKCells_30h_CCA.jpeg",
                        sep = "_"),
       path = file.path(Args$PathResults, "ClusTrees"),
       bg = "white", height = 6, width = 15, limitsize = FALSE)
# Resolution of 2 seems okay 
# Better to overcluster, merging clusters is possible if necessary
Res_Method<-c("Resolution" = 2)
Stats_Subset<-rbind(Stats_Subset,
                    tibble(Description = c("Final integrated T and NK cell subset - 30h", "T and NK cell cluster resolution - 30h"),
                            Value = c(dim(IntData.T_NK.30h)[2], Res_Method[["Resolution"]]),
                            Unit = c("Cells", "")))
Res_Method[]<-paste0("Singlets_T_NK_clusters_res_",Res_Method)

IntData.T_NK.30h<-SelectMetaData(object = IntData.T_NK.30h,
                              to_stay = c("T_NK_clusters" = Res_Method[["Resolution"]]),
                              pattern_remove = "Singlets_T_NK_clusters_res_")


IntData.T_NK.30h <- RunUMAP(IntData.T_NK.30h, 
                         reduction = "T_NK_integrated.cca.rna", 
                         dims = 1:30, reduction.name = "T_NK_integrated.umap.rna")

message("SAVE CLUSTER SIZES - 30H")
FracClus<-table(IntData.T_NK.30h$T_NK_clusters)/(dim(IntData.T_NK.30h)[2])
Stats_Subset<-rbind(Stats_Subset,
                    tibble(Description = "Number of RNA CCA clusters in T/NK - 30h",
                           Value = FracClus%>%rownames()%>%as.numeric()%>%max()+1,
                           Unit = "Clusters"))
Stats_Subset<-rbind(Stats_Subset,
                    tibble(Description = paste0("Fraction T/NK cell cluster ", c(0:Stats_Subset[Stats_Subset["Description"]=="Number of RNA CCA clusters in T/NK - 30h",][["Value"]]), " - 30h"),
                           Value = FracClus[c(0:Stats_Subset[Stats_Subset["Description"]=="Number of RNA CCA clusters in T/NK - 30h",]
                                              [["Value"]])%>%as.character()]*100,
                           Unit = "%"))

Stats_Subset<-rbind(Stats_Subset,
                    tibble(Description = paste0("Cells in B cell cluster ", c(0:Stats_Subset[Stats_Subset["Description"]=="Number of RNA CCA clusters in T/NK - 30h",][["Value"]]), " - 30h"),
                           Value = table(IntData.T_NK.30h$T_NK_clusters)[c(0:Stats_Subset[Stats_Subset["Description"]=="Number of RNA CCA clusters in T/NK - 30h",]
                                                                          [["Value"]])%>%as.character()],
                           Unit = "Cells"))
rm(FracClus)

message("VISUALISATION - 30H")
##### VISUALISATION
ggarrange(plotlist = list(
  # Dimensionality reduction plot with cluster labels (RNA only)
  DimPlot(object = IntData.T_NK.30h,
          group.by = "T_NK_clusters",
          reduction = "T_NK_integrated.umap.rna", 
          label = TRUE)+
    labs(title = "T and NK cell clusters",
         subtitle = "New clustering, only T and NK cells")+
    xlab("UMAP1")+
    ylab("UMAP2")+
    themeG+
    NoLegend(),
  # Dimensionality reduction plot with previous cluster labels (RNA only)
  DimPlot(object = IntData.T_NK.30h,
          group.by = "cca_clusters_rna",
          reduction = "T_NK_integrated.umap.rna", 
          label = TRUE)+
    labs(title = "T and NK cell clusters",
         subtitle = "Old clustering, all cells")+
    xlab("UMAP1")+
    ylab("UMAP2")+
    themeG+
    NoLegend(),
  DimPlot(object = IntData.T_NK.30h,
          group.by = "Stimulation",
          reduction = "T_NK_integrated.umap.rna", 
          label = F)+
    labs(title = "T and NK cell clusters",
         subtitle = "Stimulatory condition")+
    xlab("UMAP1")+
    ylab("UMAP2")+
    themeG+
    scale_color_manual(values = colorsStim)+
    theme(legend.justification=c(1,1.14),
          legend.position.inside=c(1,1.14))+
    guides(color=guide_legend(override.aes = list(size = 3),
                              ncol = 2)),
  DimPlot(object = IntData.T_NK.30h,
          group.by = "TimePoint",
          reduction = "T_NK_integrated.umap.rna", 
          label = F)+
    labs(title = "T and NK cell clusters",
         subtitle = "Incubation time")+
    xlab("UMAP1")+
    ylab("UMAP2")+
    themeG+
    scale_color_manual(values = colorsTime)+
    theme(legend.justification=c(1,1.14),
          legend.position.inside=c(1,1.14))+
    guides(color=guide_legend(override.aes = list(size = 3),
                              ncol = 2))),
  ncol = 2, nrow = 2)

ggsave(filename = paste("NoBase_Integrated",
                        "UMAPs_IntegratedT_NK_30h.jpeg",
                        sep  ="_"),
       path = file.path(Args$PathResults, "T_NKCells"),
       device = "jpeg", bg = "white",
       width = 13, height = 13)

saveRDS(IntData.T_NK.30h,
        file = file.path(Args$PathResults, 
                         "SeuratObjects", 
                         paste0("seuratObj_NoBase_Integrated_T_NKClusters_30h.rds")))

message("DIFFERENTIAL GENE EXPRESSION")
# Identification of differentially expressed genes/proteins
Idents(IntData.T_NK.30h)<-"T_NK_clusters"
# RNA
DiffMarkers_AllG_T_NK.30h <- FindAllMarkers(IntData.T_NK.30h, 
                                         # Minimal gene fraction as condition to test
                                         min.pct = 0.10, 
                                         min.diff.pct = 0.30, 
                                         logfc.threshold = 0.30,
                                         # pvalue threshold to show gene
                                         return.thresh = 0.01, 
                                         only.pos = FALSE, 
                                         assay="RNA",
                                         test.use = "poisson",
                                         group.by = "T_NK_clusters")

DiffMarkers_AllG_T_NK.30h<-DiffMarkers_AllG_T_NK.30h%>%
  mutate(cluster=factor(cluster,
                        # Factorization of the clusters
                        levels = levels(DiffMarkers_AllG_T_NK.30h$cluster)%>%as.numeric()%>%sort(),
                        ordered=TRUE),
         # DE score
         Score = pct.1 / (pct.2 + 0.01) *avg_log2FC)%>%
  arrange(cluster, -Score)

# Save the top 5 genes/proteins in data frame
TopDE_All_T_NK.30h<-data.frame()
for (i in levels(IntData.T_NK.30h$T_NK_clusters)){
  TopDE_All_T_NK.30h<-rbind(TopDE_All_T_NK.30h,
                         data.frame(cluster = i,
                                    # Necessary for wider format
                                    Feature = c(paste0("gene", c(1:10))),
                                    Name = c((DiffMarkers_AllG_T_NK.30h%>%filter(cluster == i))
                                             # Extract genes
                                             [1:10, "gene"])))
}
# To wider format
TopDE_All_T_NK.30h<-TopDE_All_T_NK.30h%>%
  pivot_wider(values_from = Name,
              names_from = Feature)%>%
  mutate(cluster=factor(x = cluster,
                        levels = cluster%>%
                          as.numeric()%>%sort(),
                        ordered = TRUE))%>%
  arrange(cluster)

# Save and clear workspace
# All markers
write.xlsx(x = DiffMarkers_AllG_T_NK.30h, file = file.path(Args$PathResults, "T_NKCells", 
                                                        paste("NoBase_Integrated", "DEGenes_T_NK_30h_RNA.xlsx", 
                                                              sep = "_")),
           sheetName = "DE_Genes", col.names  = T, row.names = F)

# Top 10 markers
write.xlsx(x = as.data.frame(TopDE_All_T_NK.30h),
           file = file.path(Args$PathResults, "T_NKCells",
                            paste("NoBase_Integrated", "Top10_DEGenes_T_NK_30h_RNA.xlsx",
                                  sep = "_")),
           sheetName = "TopDE_Genes", col.names = T, row.names = F, showNA = F)

# Clear workspace
rm(DiffMarkers_AllG_T_NK.30h, TopDE_All_T_NK.30h)

# Clusters based on combination of ADT and RNA (WNN)
# Calculate total number of cells in clusters
# Calculate fractions of each cell type in layer 1
AzimuthProp_T_NK_L2.30h<-left_join(x = CalClusCount(object = IntData.T_NK.30h,
                                                 clusters = "T_NK_clusters",
                                                 feature = NULL, barplot = F)[[1]],
                                y = CalClusFrac(object = IntData.T_NK.30h, 
                                                clusters = "T_NK_clusters",
                                                feature = "predicted.celltype.l2",
                                                barplot = F)[[1]]%>%
                                  mutate(Fraction = Fraction * 100)%>%
                                  pivot_wider(names_from = predicted.celltype.l2,
                                              values_from = Fraction),
                                by = "T_NK_clusters")%>%
  select(clusters = T_NK_clusters, 
         `Number of cells` = Counts,
         everything())%>%
  arrange(clusters)

write.xlsx(x = AzimuthProp_T_NK_L2.30h,
           file = file.path(Args$PathResults, "T_NKCells",
                            paste("NoBase_Integrated", "Top10_DEGenes_T_NK_30h_RNA.xlsx",
                                  sep = "_")),
           sheetName = "CellType_L2_Fractions", col.names = T, row.names = F, showNA = F,
           append = T)

rm(AzimuthProp_T_NK_L2.30h, IntData.T_NK.30h)

message("SAVE SUBSET PROPERTIES")
write.table(x = Stats_Subset, 
            file = file.path(Args$PathResults, 
                             "OverviewFeaturesDataSubSets.txt"), append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

message("FINISHED")
sessionInfo()
