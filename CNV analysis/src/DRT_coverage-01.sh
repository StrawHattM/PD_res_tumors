#!/bin/bash

#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name="CNV_cov"
#SBATCH --mail-user=martin.gonzalez@vetsuisse.unibe.ch
#SBATCH --mail-type=all
#SBATCH --time=23:59:59

module load Workspace
module load vital-it
module load UHTS/Analysis/cnvkit/0.9.6

cd $WORKSPACE/resistant_tumors

for BAM in $WORKSPACE/resistant_tumors/BAM_files/*.bam
do

cnvkit.py coverage $BAM $WORKSPACE/resistant_tumors/GCA_001624535.1_FVB_NK_v1.ensGene.bed \
	-f $WORKSPACE/FVB-NJ_reference_genome/Mus_musculus_fvbnj.FVB_NJ_v1.dna.toplevel.fa \
	-o $WORKSPACE/resistant_tumors/batch_trial/${BAM}_targetcoverage.cnn -p 32
done
