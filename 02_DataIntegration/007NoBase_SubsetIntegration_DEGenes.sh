#!/bin/sh

#SBATCH --job-name=Integration_Subsets_NoBase
#SBATCH --ntasks-per-node=18
#SBATCH --time=2:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tine.dhamers@ugent.be
#SBATCH --chdir=/data/gent/vo/000/gvo00027/TOBI/CVID/Scripts/scCITEAnalysis/Integration
#SBATCH -o 007NoBase_SubsetIntegration_DEGenes.out
#SBATCH -e 007NoBase_SubsetIntegration_DEGenes.err
#SBATCH --mem=90gb

# Load the R bundle
ml R-bundle-Bioconductor/3.16-foss-2022b-R-4.2.2
# Navigate to the correct working directory
cd /data/gent/vo/000/gvo00027/TOBI/CVID/Scripts/scCITEAnalysis/Integration

# Use the correct inputs for the R script
LibFolder="/data/gent/vo/000/gvo00027/SingleCell10X/seuratintegrate/Rpackages/R_4.2.1-foss-2022a"
AllCells_Annotated="/kyukon/scratch/gent/vo/000/gvo00027/TOBI/Projects/CVID/Results/scCITE/Exp_Nov2022_HCSamples_CEV017-22/NoBaseline/SeuratObjects/seuratObj_PatStim_NoBase_DoubletsRemoved.rds"
PathResults="/kyukon/scratch/gent/vo/000/gvo00027/TOBI/Projects/CVID/Results/scCITE/Exp_Nov2022_HCSamples_CEV017-22/NoBaseline/"

# Execute R script with the necessary arguments
Rscript '/kyukon/data/gent/vo/000/gvo00027/TOBI/CVID/Scripts/scCITEAnalysis/Integration/007NoBase_SubsetIntegration_DEGenes.R' \
       	--LibFolder=$LibFolder \
	--AllCells_Annotated=$AllCells_Annotated \
	--PathResults=$PathResults

<<JobInfo
Used CPU time       : 00:27:45.178
% User (Computation): 97.56
% System (I/O)      :  2.44
Mem reserved        : 90G
Max Mem used        : 45.79G (node4115.gallade.os)
Max Disk Write      : 12.81G (node4115.gallade.os)
Max Disk Read       : 14.55G (node4115.gallade.os)
JobInfo
