#!/usr/bin/bash

mkdir metrics_reports

for file in annotated_cns/*.cns

do
echo ${file}
NAME=${file##*/}
NAME=${NAME%.cns}
echo $NAME

cnvkit.py genemetrics ${file} -o metrics_reports/gene_report_${NAME}.txt --drop-low-coverage

done

cnvkit.py metrics annotated_cns/*.cns -o metrics_reports/general_metrics.txt
