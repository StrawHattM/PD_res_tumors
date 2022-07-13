#!/bin/bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name="CNV_D1ref"
#SBATCH --mail-user=martin.gonzalez@vetsuisse.unibe.ch
#SBATCH --mail-type=all
#SBATCH --time=3:59:59

module load Workspace
module load Anaconda3

eval "$(conda shell.bash hook)"

conda activate cnvkitenv

RESTUM=${WORKSPACE}/resistant_tumors

cd $RESTUM

cnvkit.py reference $RESTUM/bin_coverages/D1_SPLEEN.{,anti}targetcoverage.cnn -o $RESTUM/D1_fixed_reference.cnn \
		-f $RESTUM/refgenome_testall_fixed.fa --no-edge
