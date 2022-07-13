#!/usr/bin/bash

mkdir metrics_reports
mkdir metrics_reports/CBS
mkdir metrics_Reports/HMM

WORKD=/mnt/d/Bibliotecas/Biologia/Experimentos/WGS_Genevia_docetaxel\&RT-resistant_mouse_tumors/matt_analysis

cd $WORKD

conda activate cnvkitenv

for CNR in $WORKD/copy_number_ratios/*.cnr
do

NAME=${CNR#*ratios/}

echo $NAME

cnvkit.py genemetrics ${CNR} -s $WORKD/called_cns/CBS/${NAME/.cnr/_called.cns} -o metrics_reports/CBS/gene_report_${NAME/.cnr/.txt} -m 2 --drop-low-coverage
cnvkit.py genemetrics ${CNR} -s $WORKD/called_cns/HMM/${NAME/.cnr/_called.cns} -o metrics_reports/HMM/gene_report_${NAME/.cnr/.txt} -m 2 --drop-low-coverage

done


