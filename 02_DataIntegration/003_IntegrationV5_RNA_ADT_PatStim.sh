#!/bin/sh

#SBATCH --job-name=IntegrationV5_RNA_ADT_PatStim
#SBATCH --ntasks-per-node=18
#SBATCH --time=7:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tine.dhamers@ugent.be
#SBATCH --chdir=/data/gent/vo/000/gvo00027/TOBI/CVID/Scripts/scCITEAnalysis/Integration
#SBATCH -o 003_IntegrationV5_RNA_ADT_PatStim.out
#SBATCH -e 003_IntegrationV5_RNA_ADT_PatStim.err
#SBATCH --mem=200gb

# Load the R bundle
ml R-bundle-Bioconductor/3.16-foss-2022b-R-4.2.2
# Navigate to the correct working directory
cd /data/gent/vo/000/gvo00027/TOBI/CVID/Scripts/scCITEAnalysis/Integration

# Use the correct inputs for the R script
LibFolder="/data/gent/vo/000/gvo00027/SingleCell10X/seuratintegrate/Rpackages/R_4.2.1-foss-2022a"
SeuratVersion="v5"
SaveRDSFile="/scratch/gent/vo/000/gvo00027/TOBI/Projects/CVID/Results/scCITE/Exp_Nov2022_HCSamples_CEV017-22/NoBaseline/SeuratObjects/seuratObj_CVID_PatStimSplit_AdjDim_integratedV5_RNA_ADT_noBaseline.rds"
MergedDataFile="/kyukon/data/gent/vo/000/gvo00027/TOBI/CVID/Results/scCITE/Exp_Nov2022_HCSamples_CEV17-22/IntegrationV5/SeuratObjects/MergedSeuratObj_AllSamples.rds"
SplitArg="PatStim"

# Execute R script with the necessary arguments
Rscript "/kyukon/data/gent/vo/000/gvo00027/TOBI/CVID/Scripts/scCITEAnalysis/Integration/003_IntegrationV5_RNA_ADT.R" \
	--LibFolder=$LibFolder \
	--SeuratVersion=$SeuratVersion \
	--SaveRDSFile=$SaveRDSFile \
	--MergedDataFile=$MergedDataFile \
	--SplitArg=$SplitArg

<<JobInfo
Used walltime       : 03:14:39
Used CPU time       : 03:13:33
% User (Computation): 99.17%
% System (I/O)      :  0.83%
Mem reserved        : 500G
Max Mem used        : 82.96G (node4109.gallade.os)
Max Disk Write      : 245.76K (node4109.gallade.os)
Max Disk Read       : 19.82M (node4109.gallade.os)
JobInfo

<<JobInfo_CEVSplit
Used walltime       : 01:51:19
Used CPU time       : 01:50:31
% User (Computation): 98.97%
% System (I/O)      :  1.03%
Mem reserved        : 500G
Max Mem used        : 65.87G (node4115.gallade.os)
Max Disk Write      : 245.76K (node4115.gallade.os)
Max Disk Read       : 19.83M (node4115.gallade.os)
JobInfo_CEVSplit

<<JobInfo_PatStimSplit_NoBaseline
Used walltime       : 02:42:33.0
Used CPU time       : 02:41:41.0
% User (Computation): 96.04
% System (I/O)      :  3.96
Mem reserved        : 200G
Max Mem used        : 84.25G (node4116.gallade.os)
Max Disk Write      : 16.76G (node4116.gallade.os)
Max Disk Read       : 5.21G (node4116.gallade.os)
JobInfo_PatStimSplit_NoBaseline
