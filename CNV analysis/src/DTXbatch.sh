#!/bin/bash

#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name="DTXbatch"
#SBATCH --mail-user=martin.gonzalez@vetsuisse.unibe.ch
#SBATCH --mail-type=all
#SBATCH --time=1:59:59

module load Workspace
module load Anaconda3

eval "$(conda shell.bash hook)"

conda activate cnvkitenv

RESTUM=${WORKSPACE}/resistant_tumors

cd $RESTUM

cnvkit.py reference $RESTUM/bin_coverages/*_DST.{,anti}targetcoverage.cnn -o $RESTUM/DST_reference.cnn \
		-f $RESTUM/refgenome_testall_fixed.fa --no-edge

mkdir $RESTUM/DTX_comp

for CNN in $RESTUM/bin_coverages/*_DRT.targetcoverage.cnn
do

NAME=${CNN#*coverages/}
NAME=${NAME/.targetcoverage.cnn/.cnr}

cnvkit.py fix $CNN ${CNN/targetcoverage/antitargetcoverage} $RESTUM/DST_reference.cnn --no-edge -o $RESTUM/DTX_comp/${NAME}

cnvkit.py segment $RESTUM/DTX_comp/${NAME} -o $RESTUM/DTX_comp/${NAME/.cnr/_toDST-cbs.cns} -m cbs -t 1e-6 -p 128
cnvkit.py segment $RESTUM/DTX_comp/${NAME} -o $RESTUM/DTX_comp/${NAME/.cnr/_toDST-hmm.cns} -m hmm -p 128

done
