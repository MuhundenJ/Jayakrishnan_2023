#!/bin/sh



#######################################################################################################################
################### prepare genomes ###################
module load ngs/bowtie2
module load ngs/samtools
module load ngs/bedtools2
module load ngs/Homer
module load ngs/UCSCutils/3.4.1

#######################################################################################################################
#get filename 

cd ./Homer/


#######################################################################################################################
################### makeUCSCfile ###################


#select inputs normalized to vir or mel

### explore different kinds of normalization with or without inputs  
mkdir ./bedgraphs_Input_Pseudo
mkdir ./bedgraphs_NoInputNorm
mkdir ./bedgraphs_Inputs

gen=("dmelNorm" "dvirNorm")

### change here depending on ChIP conditions 
conditions=("GST" "NSD" "Set2" "Ash1")


for a in ${gen[@]}
do

	for b in ${conditions[@]}
	do
		
		inp_gen_cond=`ls -d *.dir|grep ${a}|grep ${b}|grep "Input"|sed "s/://"`
		samples_gen_cond=`ls -d *.dir|grep ${a}|grep ${b}|grep -v "Input"|sed "s/://"`	

		for sample in ${samples_gen_cond}
		do
			
			echo ${sample} ${inp_gen_cond}

			makeUCSCfile ${sample} -o ./bedgraphs_NoInputNorm/${sample}.bedgraph
        		makeUCSCfile ${sample} -i ${inp_gen_cond} -pseudo 1 -o ./bedgraphs_Input_Pseudo/${sample}.bedgraph

		done
		makeUCSCfile ${inp_gen_cond} -o ./bedgraphs_Inputs/${inp_gen_cond}.bedgraph
	done

done


#### generate bigwigs for visualization

cd ./bedgraphs_Input_Pseudo

for file in *.bedgraph.gz; do
        base=`echo ${file} | sed "s/.bedgraph.gz//"`
        gunzip ${base}.bedgraph.gz
        bedGraphToBigWig ${base}.bedgraph /work/project/becbec_005/ChIP_Seq/Drosophila_melanogaster_UCSC_dm6/dm6.chrom.sizes ${base}.bw

done

cd ../bedgraphs_NoInputNorm

for file in *.bedgraph.gz; do
        base=`echo ${file} | sed "s/.bedgraph.gz//"`
        gunzip ${base}.bedgraph.gz
	bedGraphToBigWig ${base}.bedgraph /work/project/becbec_005/ChIP_Seq/Drosophila_melanogaster_UCSC_dm6/dm6.chrom.sizes ${base}.bw
	
done


cd ../bedgraphs_Inputs

for file in *.bedgraph.gz; do
        base=`echo ${file} | sed "s/.bedgraph.gz//"`
        gunzip ${base}.bedgraph.gz
	bedGraphToBigWig ${base}.bedgraph /work/project/becbec_005/ChIP_Seq/Drosophila_melanogaster_UCSC_dm6/dm6.chrom.sizes ${base}.bw
done

#############################:q############################################################################################
