#!/usr/bash/bin

for file in *Exclusive.cns

do
echo ${file}
cnvkit.py genemetrics ${file} -o ${file}_genereport_m3.txt -m 3
cnvkit.py genemetrics ${file} -o ${file}_genereport_m2.txt -m 2

done
