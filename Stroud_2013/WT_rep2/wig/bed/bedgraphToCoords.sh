#!/bin/bash

#######################################################################################
# Convert bedgraph files to bed files containing coordinates and normalised coverage  #
#######################################################################################


for i in `ls GSM980986_WT_rep2* GSM981003_cmt3* GSM981009_ddm1* GSM981017_drm3* GSM981031_met1*`
  do
    awk '{print $1, $3, $3, $5}' ${i} > ${i}.gr.bed
    sed 's/ /\t/g' ${i}.gr.bed > ${i}.gr.tab.bed
    rm ${i}.gr.bed
  done

