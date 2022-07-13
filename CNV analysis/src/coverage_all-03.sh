#!/bin/bash

#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name="CNV_DRTcov"
#SBATCH --mail-user=martin.gonzalez@vetsuisse.unibe.ch
#SBATCH --mail-type=all
#SBATCH --time=11:59:59

module load Workspace
module load Anaconda3

eval "$(conda shell.bash hook)"

conda activate cnvkitenv

RESTUM=${WORKSPACE}/resistant_tumors

cd $RESTUM

for BAM in $RESTUM/BAM_files/*.bam
do

NAME=${BAM#*files/}
NAME=${NAME/%_merged*.bam/_targetcoverage.cnn}
NAME=${NAME/%_recal*.bam/_targetcoverage.cnn}

cnvkit.py coverage $BAM $RESTUM/GCA_001624535.1_FVB_NK_v1.ensGene.bed \
	-f $WORKSPACE/FVB-NJ_reference_genome/Mus_musculus_fvbnj.FVB_NJ_v1.dna.toplevel.fa \
	-o $RESTUM/batch_trial/$NAME -p 32
	
done
