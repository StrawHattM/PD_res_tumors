#!/usr/bin/bash

#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name="calls"
#SBATCH --mail-user=martin.gonzalez@vetsuisse.unibe.ch
#SBATCH --mail-type=all
#SBATCH --time=11:59:59

module load Workspace
module load Anaconda3



eval "$(conda shell.bash hook)"

conda activate cnvkitenv

RESTUM=${WORKSPACE}/resistant_tumors

cd $RESTUM

mkdir $RESTUM/called_cns

for FOLDER in CBS DTX_comp HMM
do

mkdir $RESTUM/called_cns/ÂFOLDER

for CNS in $RESTUM/segment_cns/${FOLDER}/*.cns
do

NAME=${CNS#*${FOLDER}/}
NAME=${NAME/.cns/_called.cns}

echo $NAME

cnvkit.py call $CNS -o $RESTUM/called_cns/${FOLDER}/${NAME}

done

done
