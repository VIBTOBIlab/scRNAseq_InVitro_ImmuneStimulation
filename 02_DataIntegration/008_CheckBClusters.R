# Read in all the necessary datasets and variables using 
# 006_Visualisation_DoubletCheck_ClusterAnalysis_IntegrationV5.Rmd

##### QUALITY CHECK OF CLUSTERS #####
pdf("/kyukon/data/gent/vo/000/gvo00027/TOBI/CVID/Scripts/scCITEAnalysis/Integration/CheckAnnotation_BCellClusters.pdf",
    width = 12, height = 7)

(VlnfeaturePlot_Categories(object = IntData.B.8h, 
                           features = c("percent_mito", "nUMI"),
                           metaData = "B_clusters", numSorting = T)+
  labs(title = "8h")) / 
(VlnfeaturePlot_Categories(object = IntData.B.30h, 
                            features = c("percent_mito", "nUMI"),
                            metaData = "B_clusters", numSorting = T)+
  labs(title = "30h"))

(FeaturePlot(IntData.B.8h, reduction = "B_integrated.umap.rna", 
            features = "percent_mito")+
  scale_color_gradientn(colours = c("white", "pink", "red")) +
  themeG+theme(plot.title = element_blank()) +
  labs(x = "UMAP1", y = "UMAP2", subtitle = "8h")) |
(FeaturePlot(IntData.B.30h, reduction = "B_integrated.umap.rna", 
              features = "percent_mito")+
  scale_color_gradientn(colours = c("white", "pink", "red")) +
  themeG+theme(plot.title = element_blank()) +
  labs(x = "UMAP1", y = "UMAP2", subtitle = "30h")) |
plot_layout(widths = c(1, 1), guides = "collect") |
plot_annotation(title = "% Mitochondrial Genes")

(FeaturePlot(IntData.B.8h, reduction = "B_integrated.umap.rna", 
             features = "nUMI")+
    scale_color_gradientn(colours = c("white", "pink", "red"),
                          limits = c(1000, 3600)) +
    themeG+theme(plot.title = element_blank()) +
    labs(x = "UMAP1", y = "UMAP2", subtitle = "8h")) |
  (FeaturePlot(IntData.B.30h, reduction = "B_integrated.umap.rna", 
               features = "nUMI")+
     scale_color_gradientn(colours = c("white", "pink", "red"),
                           limits = c(1000, 3600)) +
     themeG+theme(plot.title = element_blank()) +
     labs(x = "UMAP1", y = "UMAP2", subtitle = "30h")) |
  plot_layout(widths = c(1, 1), guides = "collect") |
  plot_annotation(title = "Number of UMIs")

##### CHECK T CELL AND INTERMEDIATE GENES #####
FileMarkersB<-"/scratch/gent/vo/000/gvo00027/TOBI/Projects/CVID/OverviewMarkers.xlsx"
B_L2_Feat<-(read.xlsx(FileMarkersB,sheetName = "MarkerGenes")%>%
              select(Gene, BAnnot = Second.level.annotation.B.cell.clusters, BSub = Second.level.annotation.B.cell.clusters.1)%>%
              filter(!is.na(BAnnot)))%>%
  mutate(Plasmablast = case_when(
    grepl(pattern = "plasmablast", x = BSub, ignore.case = T) ~ "X"),
    Naive = case_when(
      grepl(pattern = "naive", x = BSub, ignore.case = T) ~ "X"),
    Memory = case_when(
      grepl(pattern = "memory", x = BSub, ignore.case = T) ~ "X"),
    ABC = case_when(
      grepl(pattern = "ABC", x = BSub, ignore.case = T) ~ "X"))%>%
  select(-BSub, -BAnnot)%>%
  full_join(data.frame(Gene = c("MS4A1", "CD79A", "RALGPS2", "CD79B"),
                   Memory = c("X", "X", "X", NA),
                   Intermediate = rep("X", times = 4)))%>%
  full_join(data.frame(Gene = Sel_Genes%>%
                         filter(CellType == "T")%>%
                         select(Gene),
                       TCell = "X"))%>%
  arrange(Naive, Intermediate, Memory, Plasmablast, ABC, TCell, Gene)

DotPlot_MultipleFeat(object = IntData.B.8h, metaData = "B_clusters",
                     features = B_L2_Feat$Gene, 
                     numSorting = c(0,3,6,16,1,4,13,9,22,21,24,25,  # Naive
                                    2,5,8,18,19,20,7,15,10,11,12,14, # Memory
                                    17, # Other
                                    23))+
  labs(subtitle = "8h")
DotPlot_MultipleFeat(object = IntData.B.30h, metaData = "B_clusters",
                     features = B_L2_Feat$Gene, 
                     numSorting = c(0,4,5,7,1,9,10,6,12,14,15,16,17,21, # Naive
                                    2,3,8,11,18,13,20,22, # Memory
                                    19, # Other
                                    23))+
  labs(subtitle = "30h")

FeaturePlot(IntData.B.8h, features = (B_L2_Feat%>%filter(TCell == "X"))[,"Gene"], 
            reduction = "B_integrated.umap.rna") & 
  labs(x = "UMAP1", y = "UMAP2", subtitle = "8h") & themeG 

FeaturePlot(IntData.B.30h, features = (B_L2_Feat%>%filter(TCell == "X"))[,"Gene"], 
            reduction = "B_integrated.umap.rna") & 
  labs(x = "UMAP1", y = "UMAP2", subtitle = "30h") & themeG 

FeaturePlot(x, features = (B_L2_Feat%>%filter(TCell == "X"))[,"Gene"], 
            reduction = "B_integrated.umap.rna") & 
  labs(x = "UMAP1", y = "UMAP2") & themeG + plot_annotation(title = "8h")

FeaturePlot(object = IntData.B.8h, features = "scDblFinder.score", 
            reduction = "B_integrated.umap.rna")+
  themeG+labs(subtitle = "8h", x = "UMAP1", y = "UMAP2")+theme(plot.title = element_blank())+
  scale_color_gradientn(colours = c("white", "darkred"),
                        limits = c(0,0.65)) |
FeaturePlot(object = IntData.B.30h, features = "scDblFinder.score", 
              reduction = "B_integrated.umap.rna")+
  themeG+labs(subtitle = "30h", x = "UMAP1", y = "UMAP2")+theme(plot.title = element_blank()) +
  scale_color_gradientn(colours = c("white", "darkred"),
                        limits = c(0,0.65)) | 
  plot_layout(guides = "collect")|plot_annotation(title = "scDblFinder score")

(RidgePlot(IntData.B.8h, features = "scDblFinder.score", sort = T,
          cols = alpha(colorRampPalette("#B7C0FF")(26),0.5))+
  theme(legend.position = "none", plot.title = element_blank(), axis.title = element_blank())+
  labs(subtitle = "8h")) |
(RidgePlot(IntData.B.30h, features = "scDblFinder.score", sort = T,
            cols = alpha(colorRampPalette("#B7C0FF")(31),0.5))+
  theme(legend.position = "none", plot.title = element_blank(), axis.title = element_blank())+
  labs(subtitle = "30h")) |
  plot_annotation(title = "scDblFinder score")

dev.off()


##### CLUSTERING #####

### 8h ###
# Remove the clusters high in mitochondrial genes + T marker gene expression
Idents(IntData.B.8h)<-"B_clusters"
IntData.B.8h<-subset(IntData.B.8h, idents = c(0:16,18:21,23:25))

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
IntData.B.8h<-Seurat::RunPCA(IntData.B.8h, assay = "RNA", reduction.name = "B_Sel_pca.unintegrated.rna")
IntData.B.8h<-Seurat::FindNeighbors(IntData.B.8h, dims = 1:20, reduction = "B_Sel_pca.unintegrated.rna")
IntData.B.8h<-Seurat::FindClusters(IntData.B.8h, resolution = 2, cluster.name = "B_Sel_pca_unintegrated_clusters.rna")
IntData.B.8h<-Seurat::RunUMAP(IntData.B.8h, dims = 1:20, reduction = "B_Sel_pca.unintegrated.rna", 
                              reduction.name = "B_Sel_umap.unintegrated.rna")
# CCA integration of RNA layer
IntData.B.8h<-IntegrateLayers(IntData.B.8h, method = CCAIntegration,
                              orig.reduction = "B_Sel_pca.unintegrated.rna", new.reduction = "B_Sel_integrated.cca.rna",
                              verbose = F)
# Joining layers
IntData.B.8h[["RNA"]]<-JoinLayers(IntData.B.8h[["RNA"]])

# UMAP generation
IntData.B.8h <- RunUMAP(IntData.B.8h, 
                        reduction = "B_Sel_integrated.cca.rna", 
                        dims = 1:30, reduction.name = "B_Sel_integrated.umap.rna")

##### CLUSTERING OF B CELL SUBSET #####
message("CLUSTERING OF B CELLS - 8H")
# Nearest neighbours graph generation based on CCA (only RNA) integration method
IntData.B.8h <- FindNeighbors(IntData.B.8h, 
                              reduction = "B_Sel_integrated.cca.rna", 
                              dims = 1:30, graph.name = "B_Sel_integrated.cca.rna_graph")
# Clustering using multiple resolution values
for (val_res in c(0.1, seq(0.5,2.5, by=0.5))){
  ResVal<-paste0("Singlets_B_Sel_clusters_res_", val_res)
  message(paste0("Busy working on B cell clustering with resolution ", val_res))
  IntData.B.8h<-FindClusters(IntData.B.8h,resolution=val_res,
                             verbose = FALSE, cluster.name = ResVal,
                             graph.name = "B_Sel_integrated.cca.rna_graph")}

# Clustree to check the select the appropriate resolution
clustree(IntData.B.8h, 
         # Select only the correct clusters based on the resolution
         prefix = "Singlets_B_Sel_clusters_res_",
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

IntData.B.8h<-SelectMetaData(object = IntData.B.8h,
                             to_stay = c("B_Sel_clusters_0.1" = "Singlets_B_Sel_clusters_res_0.1",
                                         "B_Sel_clusters_1.5" = "Singlets_B_Sel_clusters_res_1.5"),
                             pattern_remove = "Singlets_B_Sel_clusters_res_")

# Visualisation of the clusters 
DimPlot(object = IntData.B.8h, reduction = "B_Sel_integrated.umap.rna",
        group.by = "B_Sel_clusters_0.1", label = T)+
  DimPlot(object = IntData.B.8h, reduction = "B_Sel_integrated.umap.rna",
          group.by = "B_Sel_clusters_1.5", label = T) &
  labs(title = "8h", x = "UMAP1", y = "UMAP2") & themeG &
  theme(plot.subtitle = element_blank(), 
        legend.position = "none")

### 30h ###
# Remove the clusters high in mitochondrial genes + T marker gene expression
Idents(IntData.B.30h)<-"B_clusters"
IntData.B.30h<-subset(IntData.B.30h, idents = c(0,2:18,20:23))

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
IntData.B.30h<-Seurat::RunPCA(IntData.B.30h, assay = "RNA", reduction.name = "B_Sel_pca.unintegrated.rna")
IntData.B.30h<-Seurat::FindNeighbors(IntData.B.30h, dims = 1:20, reduction = "B_Sel_pca.unintegrated.rna")
IntData.B.30h<-Seurat::FindClusters(IntData.B.30h, resolution = 2, cluster.name = "B_Sel_pca_unintegrated_clusters.rna")
IntData.B.30h<-Seurat::RunUMAP(IntData.B.30h, dims = 1:20, reduction = "B_Sel_pca.unintegrated.rna", 
                              reduction.name = "B_Sel_umap.unintegrated.rna")
# CCA integration of RNA layer
IntData.B.30h<-IntegrateLayers(IntData.B.30h, method = CCAIntegration,
                              orig.reduction = "B_Sel_pca.unintegrated.rna", new.reduction = "B_Sel_integrated.cca.rna",
                              verbose = F)
# Joining layers
IntData.B.30h[["RNA"]]<-JoinLayers(IntData.B.30h[["RNA"]])

# UMAP generation
IntData.B.30h <- RunUMAP(IntData.B.30h, 
                        reduction = "B_Sel_integrated.cca.rna", 
                        dims = 1:30, reduction.name = "B_Sel_integrated.umap.rna")

##### CLUSTERING OF B CELL SUBSET #####
message("CLUSTERING OF B CELLS - 30h")
# Nearest neighbours graph generation based on CCA (only RNA) integration method
IntData.B.30h <- FindNeighbors(IntData.B.30h, 
                              reduction = "B_Sel_integrated.cca.rna", 
                              dims = 1:30, graph.name = "B_Sel_integrated.cca.rna_graph")
# Clustering using multiple resolution values
for (val_res in c(0.1, seq(0.5,2.5, by=0.5))){
  ResVal<-paste0("Singlets_B_Sel_clusters_res_", val_res)
  message(paste0("Busy working on B cell clustering with resolution ", val_res))
  IntData.B.30h<-FindClusters(IntData.B.30h,resolution=val_res,
                             verbose = FALSE, cluster.name = ResVal,
                             graph.name = "B_Sel_integrated.cca.rna_graph")}

# Clustree to check the select the appropriate resolution
clustree(IntData.B.30h, 
         # Select only the correct clusters based on the resolution
         prefix = "Singlets_B_Sel_clusters_res_",
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

IntData.B.30h<-SelectMetaData(object = IntData.B.30h,
                             to_stay = c("B_Sel_clusters_0.1" = "Singlets_B_Sel_clusters_res_0.1",
                                         "B_Sel_clusters_2" = "Singlets_B_Sel_clusters_res_2"),
                             pattern_remove = "Singlets_B_Sel_clusters_res_")

# Visualisation of the clusters 
DimPlot(object = IntData.B.30h, reduction = "B_Sel_integrated.umap.rna",
        group.by = "B_Sel_clusters_0.1", label = T)+
DimPlot(object = IntData.B.30h, reduction = "B_Sel_integrated.umap.rna",
          group.by = "B_Sel_clusters_2", label = T) &
  labs(title = "30h", x = "UMAP1", y = "UMAP2") & themeG &
  theme(plot.subtitle = element_blank(), 
        legend.position = "none")

##### DIFFERENTIAL GENE ANALYSIS #####
### 8h ###
Idents(IntData.B.8h)<-"B_Sel_clusters_1.5"
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
                                        group.by = "B_Sel_clusters_1.5")

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
for (i in levels(IntData.B.8h$B_Sel_clusters_1.5)){
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

# Save file
write.xlsx(x = as.data.frame(TopDE_All_B.8h),
           file = file.path(path_Results, "Markers/BCells",
                            paste("Integrated", "Top10_DEGenes_B_8h_RNA.xlsx",
                                  sep = "_")),
           sheetName = "TopDE_Genes", col.names = T, row.names = F, showNA = F)

### 30h ###
Idents(IntData.B.30h)<-"B_Sel_clusters_2"
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
                                        group.by = "B_Sel_clusters_2")

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
for (i in levels(IntData.B.30h$B_Sel_clusters_2)){
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

# Save file
write.xlsx(x = as.data.frame(TopDE_All_B.30h),
           file = file.path(path_Results, "Markers/BCells",
                            paste("Integrated", "Top10_DEGenes_B_30h_RNA.xlsx",
                                  sep = "_")),
           sheetName = "TopDE_Genes", col.names = T, row.names = F, showNA = F)

rm(TopDE_All_B.30h, TopDE_All_B.8h)

saveRDS(IntData.B.8h,
        file = file.path(path_Results, 
                         "SeuratObjects", 
                         paste0("seuratObj_Integrated_BClusters_8h.rds")))

saveRDS(IntData.B.30h,
        file = file.path(path_Results, 
                         "SeuratObjects", 
                         paste0("seuratObj_Integrated_BClusters_30h.rds")))