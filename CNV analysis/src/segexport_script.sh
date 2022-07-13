#!/usr/bin/bash

mkdir seg_files

for file in *.cns
do

echo ${file}

cnvkit.py export seg ${file} -o seg_files/${file}_segments.seg

done

cnvkit.py export seg *_CIS.call.cns -o CIS_segments.seg
cnvkit.py export seg *_CARBO.call.cns -o CARBO_segments.seg
