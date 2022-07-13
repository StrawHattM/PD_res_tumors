#!/bin/bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name="CNV_target"
#SBATCH --mail-user=martin.gonzalez@vetsuisse.unibe.ch
#SBATCH --mail-type=all
#SBATCH --time=0:59:59

module load Workspace
module load Anaconda3

eval "$(conda shell.bash hook)"

conda activate cnvkitenv

RESTUM=${WORKSPACE}/resistant_tumors

cd $RESTUM

cnvkit.py target $RESTUM/FixedGene_GCA_001624535.1_FVB_NJ_v1.ensGene.bed --annotate $RESTUM/gene_symbol_corrected.refflat -o targets.bed
