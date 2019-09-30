#!/bin/bash

#######################################################################################
# Convert bedgraph files to bed files containing coordinates and normalised coverage  #
#######################################################################################


for i in `ls GSM981015_drm12* GSM981025_hen1* GSM981057_kyp* GSM981058_suvh5* GSM981059_suvh6*`
  do
    awk '{print $1, $3, $3, $5}' ${i} > ${i}.gr.bed
    sed 's/ /\t/g' ${i}.gr.bed > ${i}.gr.tab.bed
    rm ${i}.gr.bed
  done

