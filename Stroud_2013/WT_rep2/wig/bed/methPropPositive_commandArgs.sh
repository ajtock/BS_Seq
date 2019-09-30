#!/bin/bash

#######################################################################################
# Convert negative to positive methylation proportions                                #
#######################################################################################

for inFile in "$@"
do
( sed 's/-//g' ${inFile} > ${inFile}.2
  mv ${inFile}.2 ${inFile} ) &
done
wait

