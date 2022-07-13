#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name="CNV_access"
#SBATCH --mail-user=martin.gonzalez@vetsuisse.unibe.ch
#SBATCH --mail-type=all
#SBATCH --time=0:29:59

module load Workspace
module load Anaconda3

eval "$(conda shell.bash hook)"

conda activate cnvkitenv

RESTUM=${WORKSPACE}/resistant_tumors

cd $RESTUM

cnvkit.py access $RESTUM/refgenome_testall_fixed.fa -o $RESTUM/access.FVB-NJ-fa.bed
