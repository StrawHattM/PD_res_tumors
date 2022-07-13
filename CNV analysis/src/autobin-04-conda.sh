#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name="CNV_autobin"
#SBATCH --mail-user=martin.gonzalez@vetsuisse.unibe.ch
#SBATCH --mail-type=all
#SBATCH --time=23:59:59

module load Workspace
module load Anaconda3

RESTUM=${WORKSPACE}/resistant_tumors

cd $RESTUM

eval "$(conda shell.bash hook)"

conda activate cnvkitenv

module list

cnvkit.py autobin $RESTUM/BAM_files/D10SRO406_DL_merged_fixed.bam \
                -f $WORKSPACE/FVB-NJ_reference_genome/Mus_musculus_fvbnj.FVB_NJ_v1.dna.toplevel.fa \
                -g $RESTUM/access.FVB-NJ-fa.bed \
                -t $RESTUM/GCA_001624535.1_FVB_NK_v1.ensGene.bed \
                --annotate $RESTUM/gene_symbol_corrected.refflat \
                --m wgs --short-names
