#!/usr/bin/bash

for dir in docetaxel radiotherapy platinum_drugs
do

cd ~/wgs_cnv/${dir}

bash annotate_script.sh
wait
bash genemetrics_script.sh
wait

cd annotated_cns

bash ../../graph_script.sh

done
