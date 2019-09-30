#!/bin/bash

for i in `ls GSM981015_drm12* GSM981025_hen1* GSM981057_kyp* GSM981058_suvh5* GSM981059_suvh6*`
do
echo ${i}
wig2bed < ${i} > ${i}.bed
done

