#!/usr/bin/bash

mkdir annotated_cns

for file in *.cns
do

echo ${file}

cnv_annotate.py ../gene_symbol_corrected.refflat ${file} -o annotated_cns/anno_${file}

done
