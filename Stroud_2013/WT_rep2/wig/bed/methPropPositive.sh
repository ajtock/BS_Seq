#!/bin/bash

#######################################################################################
# Convert negative to positive methylation proportions                                #
#######################################################################################


for i in `ls GSM981003_cmt3* GSM981009_ddm1* GSM981015_drm12* GSM981017_drm3* GSM981025_hen1* GSM981057_kyp* GSM981031_met1* GSM981058_suvh5* GSM981059_suvh6*`
#GSM980986_WT_rep2` 
  do
    sed 's/-//g' ${i} > ${i}.2
    mv ${i}.2 ${i}
  done

