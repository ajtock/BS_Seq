#!/bin/bash

#######################################################################################
# Convert bedgraph files to bed files containing coordinates and normalised coverage  #
#######################################################################################

for inFile in "$@"
do
( awk '{print $1, $3, $3, $5}' ${inFile} > ${inFile}.gr.bed
  sed 's/ /\t/g' ${inFile}.gr.bed > ${inFile}.gr.tab.bed
  rm ${inFile}.gr.bed ) &
done
wait

