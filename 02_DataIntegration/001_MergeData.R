print("### SAMPLE PATH ###")
# Save path to sample data
samples<-paste0(rep("CEV0", times = 6),
                c(17:22))
PathFiles<-'/scratch/gent/vo/000/gvo00027/TOBI/Projects/CVID/Results/scCITE/Exp_Nov2022_HCSamples_CEV017-22/NF_OutputData'
FileSeurat<-file.path(PathFiles,
                      samples, "QC",
                      paste0("CVID_scCITESeq_", samples),
                      "Robjects", paste0("seuratObj_", samples,'_labeled.rds'))

# Read in data (list as output)
MergeData.list<-lapply(FileSeurat,
                       readRDS)

print("### MERGE DATA ###")
# Merge the datasets
MergeData<-merge(x = MergeData.list[[1]],
                 y = MergeData.list[2:length(FileSeurat)],
                 merge.data = TRUE,
                 add.cell.ids=samples)

print("### SAVE DATA ###")
saveRDS(MergeData, file = "/kyukon/data/gent/vo/000/gvo00027/TOBI/CVID/Results/scCITE/Exp_Nov2022_HCSamples_CEV17-22/IntegrationV5/SeuratObjects/MergedSeuratObj_AllSamples.rds")

