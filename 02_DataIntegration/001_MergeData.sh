#!/bin/sh

#SBATCH --job-name=MergeData
#SBATCH --ntasks-per-node=1
#SBATCH --time=0:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tine.dhamers@ugent.be
#SBATCH --chdir=/data/gent/vo/000/gvo00027/TOBI/CVID/Scripts/scCITEAnalysis/Integration
#SBATCH -o MergeData.out
#SBATCH -e MergeData.err
#SBATCH --mem=100gb

ml R-bundle-Bioconductor/3.16-foss-2022b-R-4.2.2
cd /data/gent/vo/000/gvo00027/TOBI/CVID/Scripts/scCITEAnalysis/Integration
Rscript "/kyukon/data/gent/vo/000/gvo00027/TOBI/CVID/Scripts/scCITEAnalysis/Integration/MergeData.R"
