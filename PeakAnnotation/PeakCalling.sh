#!/bin/sh

module load ngs/samtools/1.9
module load ngs/Homer/4.9
module load ngs/UCSCutils
module load ngs/bedtools2
module load ngs/MACS2/2.1.2



mkdir ./Peaks_Narrow
mkdir ./Peaks_Narrow_stringent
mkdir ./Peaks_Broad


#### list bigwigs generated from merged replicates of H3K36me1/2/3
SAMPLES=`ls *.bw`



### convert bigwigs to bedgraphs 
for SAMPLE in ${SAMPLES}; do

	BASE=`echo ${SAMPLE}|sed -e 's/.bw//g'`
	bigWigToBedGraph ${BASE}.bw ${BASE}.bedgraph
done

SAMPLES=`ls *.bedgraph`	



for SAMPLE in ${SAMPLES};do
	
	BASE=`echo ${SAMPLE}|sed -e 's/.bedgraph//g'`
	
	### explore different peak calling parameters by MACS2 broad or narrow peaks 
	macs2 bdgpeakcall -l 1000 -i ${BASE}.bedgraph --outdir ./Peaks_Narrow --o-prefix ${BASE}
	macs2 bdgpeakcall -c 3 -l 1000 -i ${BASE}.bedgraph --outdir ./Peaks_Narrow_stringent --o-prefix ${BASE}
	#macs2 bdgpeakcall -c 4 -l 1000 -i ${BASE}.bedgraph --outdir ./Peaks_Narrow_stringent --o-prefix ${BASE}	
	macs2 bdgbroadcall -l 1000 -i ${BASE}.bedgraph --outdir ./Peaks_Broad --o-prefix ${BASE}
done


