#!/bin/bash

for i in `ls GSM980986_WT_rep2* GSM981003_cmt3* GSM981009_ddm1* GSM981017_drm3* GSM981031_met1*`
do
echo ${i}
wig2bed < ${i} > ${i}.bed
done

