<div style="text-align: justify">

# scRNAseq_InVitro_ImmuneStimulation
 This repository encloses interesting code supporting the manuscript entitled "Stimulation strategy matters: A B-cell focused single-cell transcriptomic study on in vitro stimulated immune cells". This repository contains four folders, each of them going more in depth on specific analysis parts. 

- **01_PreprocessingPipeline**: This folder contains the code used for processing the raw sequencing data files.
- **02_DataIntegration**: This folder contains the code used for integration of the data prior downstream processing.
- **03_MiloR_AbundancyTesting**: This folder contains the code used for the MiloR neighbourhood definition, abundancy testing and refinement of the cell annotation.
- **04_MultiNichNetR**: This folder contains the code used for cell-cell communication inference and visualization using the NicheNetR package.

## Preprocessing of the pipeline.
The code and necessary additional data files for preprocessing of the raw sequencing data files are available in the folder **`01_PreprocessingPipeline`**. An overview of the files and their goal is given in the table below. Furthermore, in each of the scripts, additional information on how the code works and and which parameters are expected is added.

First of all, two files were added in the **`SampleInformation`** folder that contained information on the samples, necessary for downstream analysis.

| File name                             | Information                              |
| :----:                                | :----                                   |
|`20230301_SampleInformation/CVID_scCITE.csv`    |Sequencing indexes linked to the CEV pools, necessary for demultiplexing of the raw sequencing data.|
|`20230301_SampleInformation/SampleAnnotation_CVID.txt`|In depth overview of the sample metadata, similar to what was provided in the Supplementary Table S3 of the manuscript.|

Furthermore, the files in the **`01_PreprocessingPipeline`** are elaborated on in the table below.

| File name                             | Information                              |
| :----:                                | :----                                   |
| `scCITE_NF_CVID.nf` | The NextFlow preprocessing pipeline used for preprocessing of raw sequencing files, starting from demultiplexing up to single-cell preprocessing steps using CellRanger and Freemuxlet.|
| `scCITE_NF_CVID.config` | The configuration file necessary to run the NextFlow preprocessing pipeline. |
| `RunNF.pbs` | This script is used to finally run the NextFlow preprocessing pipeline using the correct configuration file. |
| `005_QC_scCITESeq.Rmd`| 



</div>