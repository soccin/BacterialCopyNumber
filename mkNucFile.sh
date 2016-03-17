#!/bin/bash

BEDTOOLS_GENOME=$1
GENOMEFASTA=$2
GENOMETAG=$(basename $GENOMEFASTA | sed 's/.fa.*//')
echo $GENOMETAG

bedtools makewindows -w 1000 -s 100 -g $BEDTOOLS_GENOME \
    | bedtools nuc -fi $GENOMEFASTA -bed - \
    | gzip -c - > ${GENOMETAG}_1000by100_TileBin.nuc.gz


