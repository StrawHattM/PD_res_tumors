#!/usr/bin/bash

#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name="HMMsegment"
#SBATCH --mail-user=martin.gonzalez@vetsuisse.unibe.ch
#SBATCH --mail-type=all
#SBATCH --time=23:59:59

module load Workspace
module load Anaconda3
module load R

eval "$(conda shell.bash hook)"

conda activate cnvkitenv

RESTUM=${WORKSPACE}/resistant_tumors

cd $RESTUM

for CNR in $RESTUM/copy_number_ratios/*.cnr
do

NAME=${CNR#*ratios/}
NAME=${NAME/.cnr/.cns}

cnvkit.py segment $CNR -o $RESTUM/segment_cns/HMM/$NAME -m hmm -p 64 --drop-low-coverage

done
