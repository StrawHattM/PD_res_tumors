#!/bin/bash

#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name="CNV_batch_trial"
#SBATCH --mail-user=martin.gonzalez@vetsuisse.unibe.ch
#SBATCH --mail-type=all
#SBATCH --time=23:59:59

module load Workspace
module load vital-it
module load UHTS/Analysis/cnvkit/0.9.6

cd $WORKSPACE/resistant_tumors

cnvkit.py batch $WORKSPACE/resistant_tumors/BAM_files/*_DRT_*.bam \
	--normal $WORKSPACE/resistant_tumors/BAM_files/*_DL_*.bam \
	--targets $WORKSPACE/resistant_tumors/GCA_001624535.1_FVB_NK_v1.ensGene.bed \
	--annotate $WORKSPACE/resistant_tumors/gene_symbol_corrected.refflat \
    	--fasta $WORKSPACE/FVB-NJ_reference_genome/Mus_musculus_fvbnj.FVB_NJ_v1.dna.toplevel.fa \
    	--output-reference DL_Reference.cnn --output-dir $WORKSPACE/resistant_tumors/batch_trial/\
	--processes 32 --method wgs \
    	--diagram --scatter
