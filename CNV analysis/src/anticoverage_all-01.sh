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
NAME=${NAME/%_merged*.bam/.antitargetcoverage.cnn}
NAME=${NAME/%_recal*.bam/.antitargetcoverage.cnn}

cnvkit.py coverage $BAM $RESTUM/antitargets.bed \
	-f $WORKSPACE/FVB-NJ_reference_genome/Mus_musculus_fvbnj.FVB_NJ_v1.dna.toplevel.fa \
	-o $RESTUM/bin_coverages/$NAME -p 32
	
done
