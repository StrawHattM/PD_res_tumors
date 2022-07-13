#!/usr/bin/bash

mkdir heatmaps

cnvkit.py heatmap *.cns -d -o heatmap.pdf

for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 X

do

cnvkit.py heatmap *.cns -d -o heatmaps/heatmap_chr${i}.pdf -c chr${i}

done
