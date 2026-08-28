<div style="text-align: justify">

# scRNAseq_InVitro_ImmuneStimulation
 This repository encloses interesting code supporting the manuscript entitled "Stimulation strategy matters: A B-cell focused single-cell transcriptomic study on in vitro stimulated immune cells". This repository contains four folders, each of them going more in depth on specific analysis parts. 

- **`01_PreprocessingPipeline`**: This folder contains the code used for processing the raw sequencing data files.
- **`02_DataIntegration`**: This folder contains the code used for integration of the data prior downstream processing.
- **`03_MiloR_DiffGeneExp`**: This folder contains the code used for the MiloR neighbourhood definition, abundancy testing, refinement of the cell annotation and differential gene expression.
- **`04_MultiNicheNetR`**: This folder contains the code used for cell-cell communication inference and visualization using the NicheNetR package.

## Reference to manuscript and data
The manuscript is accessible via <t style="color:red;">XXX</t>. Data is accessible on EGA via accession ID <t style="color:red;">XXX</t>.

## Preprocessing of the pipeline.
The code and necessary additional data files for preprocessing of the raw sequencing data files are available in the folder **`01_PreprocessingPipeline`**. An overview of the files and their goal is given in the table below. Furthermore, in each of the scripts, additional information on how the code works and and which parameters are expected is added.

First of all, two files were added in the **`SampleInformation`** folder that contained information on the samples, necessary for downstream analysis.

| File name                             | Information                              |
| :----:                                | :---------------------------------------------------------------------------------------|
|`20230301_SampleInformation/CVID_scCITE.csv`    |Sequencing indexes linked to the CEV pools, necessary for demultiplexing of the raw sequencing data.|
|`20230301_SampleInformation/SampleAnnotation_CVID.txt`|In depth overview of the sample metadata, similar to what was provided in the Supplementary Table S3 of the manuscript.|

Furthermore, the files in the **`01_PreprocessingPipeline`** are elaborated on in the table below.

| File name                             | Information                              |
| :----:                                | :---------------------------------------------------------------------------------------|
| `001_scCITE_NF_CVID.nf` | The NextFlow preprocessing pipeline used for preprocessing of raw sequencing files, starting from demultiplexing up to single-cell preprocessing steps using CellRanger and Freemuxlet.|
| `001_scCITE_NF_CVID.config` | The configuration file necessary to run the NextFlow preprocessing pipeline. |
| `001_RunNF.pbs` | This script is used to finally run the NextFlow preprocessing pipeline using the correct configuration file. |
| `002_QC_scCITESeq.Rmd`| This RMarkdown script contains necessary QC checks and formatting of the final dataset into a Seurat object that is used for data integration.| 


## Data integration
Seurat was used to integrate the data obtain in distinct, so called 'CEV' single-cell sample pools. This was performed using the following scripts available in the **`02_DataIntegration`** folder. 

| File name                             | Information                              |
| :----:                                | :---------------------------------------------------------------------------------------|
| `001_MergeData.R` | Merge the sample data sets into one Seurat object.|
| `001_MergeData.sh`| A shell script to run `001_MergeData.R` script easily.|
| `003_IntegrationV5_RNA_ADT.R`| Integration of the distinct samples using Seurat V5. |
| `003_IntegrationV5_RNA_ADT_PatStim.sh` | A shell script to run `003_IntegrationV5_RNA_ADT.R` script easily.|
| `004_Visualisation_ClusterAnalysis_IntegrationV5.Rmd` | An RMarkdown script to explore, visualize and cluster the integrated data.|
| `006_NoBase_scDblFinder.R` | Doublet detection using the scDblFinder tool. |
| `006_NoBase_scDblFinder.R` | A shell script to run `006_NoBase_scDblFinder.R` script easily. |
| `006_NoBase_Visualisation_DoubletCheck_ClusterAnalysis_IntegrationV5.Rmd` | Visualization and removal of the doublets prior to cell clustering and additional dataset visualization. |
| `007_NoBase_SubsetIntegration_DEGenes.R` | This script contains additional integration of the B cell subset finally resulting in the final B cell datasets for each 8 and 30 hours of incubation. |
| `007_NoBase_SubsetIntegration_DEGenes.sh` | A shell script ro run `00_7_NoBase_SubsetIntegration_DEGenes.R` easily. | 
| `008_CheckBClusters.R` | This script explains the additional removal of specific B cell clusters in the final B cell datasets.| 
| `012_BCellFiltering_DEGenes.Rmd` | Final filtering step of the B cell subset as in depth explained in the manuscript. _Note: Differential expression as coded in this script was optimized after MiloR neighbourhood abundance testing and cell annotation refinement_ | 

## MiloR neighbourhood definition and abundance testing followed by differential gene expression
This code was used for the analysis based on neighbourhood definition using MiloR as stated in the manuscript. Furthermore, cell annotation was refined after MiloR analysis and was used for final differential gene expression between conditions. The script mentioned below are available in the **`03_MiloR_DiffGeneExp`** folder.

| File name                             | Information                              |
| :----:                                | :---------------------------------------------------------------------------------------|
| `013_MiloR_DA_FilteredBCells.Rmd` | MiloR analysis. |
| `014_DEGenes_SelectedNhoods.Rmd` | Cell annotation refinement and differential gene expression. |

## Cell-cell communication analysis with MultiNicheNetR 
Cell-cell communication analysis was performed using MulitNicheNetR as stated in the manuscript. Following code files are available in the **`04_MultiNicheNetR`** folder.

| File name                             | Information                              |
| :----:                                | :---------------------------------------------------------------------------------------|
| `003_MultiNicheNetObjectInitialization_LowerDEScenario.Rmd` | Preprocessing needed to obtain the necessary MultiNicheNet object for downstream analysis.|
| `004_MultiNicheNet_LowerDE_CorrelationBasedFiltering.Rmd` | Downstream analysis of the MultiNicheNet results. |
| `AdaptedFunction_MultiNicheNet.R` | Some functions of the MultiNicheNetR package were adapted for improved visualization, which can be found in this R script |

## Manuscript figure generation
In the file `FiguresManuscript.Rmd` all the code used to generate the main and supplementary figures is merged.

</div>