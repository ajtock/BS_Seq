#!/bin/bash

#############################################################
# Convert wig files to TDF format for visualisation in IGV  #
#############################################################

# Usage on node7:
# csmit -m 20G -c 1 "bash ./wigToTDF.sh GSM980986_WT_rep2_CG"

i=$1

igvtools toTDF ${i}.wig ${i}.wig.tdf /projects/ajt200/TAIR10/tair10.chrom.sizes
