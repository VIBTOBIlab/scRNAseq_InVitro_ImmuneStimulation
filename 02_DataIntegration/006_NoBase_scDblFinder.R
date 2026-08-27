library(SingleCellExperiment)
library(scDblFinder, lib.loc = "/kyukon/data/gent/vo/000/gvo00027/TOBI/CVID/Rtools/")
library(tidyr)
library(dplyr)
library(tibble)

# Load in the seurat object
message("Read seurat object")
IntData.sce <- readRDS("/kyukon/scratch/gent/vo/000/gvo00027/TOBI/Projects/CVID/Results/scCITE/Exp_Nov2022_HCSamples_CEV017-22/NoBaseline/SeuratObjects/Intermediate_SCE_Pre-scDblFinder.rds")

# Doublet detection with scDblFinder
message("Execute scDblFinder")
IntData.sce<-scDblFinder(IntData.sce, samples="samplename",
                         clusters = "cca_clusters_rna")
message("Save scDblFinder output")
out_DblFinder<-data.frame(CellID = colnames(IntData.sce),
                          scDblFinder.class = IntData.sce@colData@listData[["scDblFinder.class"]],
                          scDblFinder.score = IntData.sce@colData@listData[["scDblFinder.score"]])
write.csv(out_DblFinder,
          file = "/kyukon/scratch/gent/vo/000/gvo00027/TOBI/Projects/CVID/Results/scCITE/Exp_Nov2022_HCSamples_CEV017-22/NoBaseline/scDblFinder_DoubletCells.csv")         