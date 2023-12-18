#!/bin/sh

module load ngs/bedtools2/2.28.0
module load ngs/Homer/4.10

cd ./HOMER

samples=`ls *.bed`

for sample in ${samples};do
	FILEBASE=`echo $sample|sed -e "s/.bed//g"`
	annotatePeaks.pl ${FILEBASE}.bed dm6 > ${FILEBASE}_peakAnno.txt 	
done

