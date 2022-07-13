#!/bin/bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name="CNV_DLref"
#SBATCH --mail-user=martin.gonzalez@vetsuisse.unibe.ch
#SBATCH --mail-type=all
#SBATCH --time=2:59:59

module load Workspace
module load Anaconda3

eval "$(conda shell.bash hook)"

conda activate cnvkitenv

RESTUM=${WORKSPACE}/resistant_tumors

cd $RESTUM

cnvkit.py reference $RESTUM/batch_trial/DL_targetcoverage.cnn -o $RESTUM/batch_trial/DL_reference.cnn /
		-f $WORKSPACE/FVB-NJ_reference_genome/Mus_musculus_fvbnj.FVB_NJ_v1.dna.toplevel.fa --no-edge
