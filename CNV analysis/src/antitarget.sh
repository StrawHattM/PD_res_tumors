#!/bin/bash

#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name="CNV_antitarget"
#SBATCH --mail-user=martin.gonzalez@vetsuisse.unibe.ch
#SBATCH --mail-type=all
#SBATCH --time=11:59:59

module load Workspace
module load Anaconda3

eval "$(conda shell.bash hook)"

conda activate cnvkitenv

RESTUM=${WORKSPACE}/resistant_tumors

cd $RESTUM

cnvkit.py antitarget $RESTUM/FixedGene_GCA_001624535.1_FVB_NJ_v1.ensGene.bed \
	-g $RESTUM/access.FVB-NJ-fa.bed \
	-o antitargets.bed
