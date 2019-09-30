#!/bin/bash

for inFile in "$@"
do
( wig2bed < ${inFile}.wig > ${inFile}.bed ) &
done
wait

