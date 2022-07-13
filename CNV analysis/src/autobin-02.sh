#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name="CNV_autobin"
#SBATCH --mail-user=martin.gonzalez@vetsuisse.unibe.ch
#SBATCH --mail-type=all
#SBATCH --time=23:59:59

module load Workspace
module load vital-it
module load UHTS/Analysis/cnvkit/0.9.6

cd $WORKSPACE/resistant_tumors

cnvkit.py autobin $WORKSPACE/resistant_tumors/BAM_files/D10SRO406_DL_merged_fixed.bam -f $WORKSPACE/FVB-NJ_reference_genome/Mus_musculus_fvbnj.FVB_NJ_v1.dna.toplevel.fa -g $WORKSPACE/resistant_tumors/access.FVB-NJ-fa.bed -t $WORKSPACE/resistant_tumors/GCA_001624535.1_FVB_NK_v1.ensGene.bed --annotate $WORKSPACE/resistant_tumors/gene_symbol_corrected.refflat --m wgs --short-names
