#!/usr/bin/bash

#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name="PD_fix"
#SBATCH --mail-user=martin.gonzalez@vetsuisse.unibe.ch
#SBATCH --mail-type=all
#SBATCH --time=11:59:59

module load Workspace
module load Anaconda3

eval "$(conda shell.bash hook)"

conda activate cnvkitenv

RESTUM=${WORKSPACE}/resistant_tumors

cd $RESTUM

for CNN in $RESTUM/bin_coverages/*_PD.targetcoverage.cnn
do

NAME=${CNN#*coverages/}
NAME=${NAME/.targetcoverage.cnn/.cnr}

cnvkit.py fix $CNN ${CNN/targetcoverage/antitargetcoverage} $RESTUM/D1_fixed_reference.cnn --no-edge -o $RESTUM/copy_number_ratios/${NAME}

done
