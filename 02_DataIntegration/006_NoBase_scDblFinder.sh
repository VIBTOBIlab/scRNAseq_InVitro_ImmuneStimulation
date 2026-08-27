#!/bin/sh

#SBATCH --job-name=scDblFinder_NoBaseline
#SBATCH --ntasks-per-node=2
#SBATCH --time=1:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tine.dhamers@ugent.be
#SBATCH --chdir=/data/gent/vo/000/gvo00027/TOBI/CVID/Scripts/scCITEAnalysis/Integration
#SBATCH -o 006NoBase_scDblFinder.out
#SBATCH -e 006NoBase_scDblFinder.err
#SBATCH --mem=90gb

# Load the R bundle
ml R-bundle-Bioconductor/3.16-foss-2022b-R-4.2.2
# Navigate to the correct working directory
cd /data/gent/vo/000/gvo00027/TOBI/CVID/Scripts/scCITEAnalysis/Integration


Rscript "/kyukon/data/gent/vo/000/gvo00027/TOBI/CVID/Scripts/scCITEAnalysis/Integration/006NoBase_scDblFinder.R"

