#!/usr/bin/bash

mkdir scatters
mkdir diagrams

for file in *.cns
do

echo ${file}

cnvkit.py scatter -s ${file} -o scatters/${file}_scatter.pdf
cnvkit.py diagram -s ${file} -o diagrams/${file}_diagram.pdf

done
