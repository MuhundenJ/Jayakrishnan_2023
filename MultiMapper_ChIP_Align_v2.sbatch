#!/bin/sh
# 
# bowtie.sbatch
#
#SBATCH -J bowtie_RepetitiveGenome   # A single job name for the array
#SBATCH -p slim16                  # Partition
#SBATCH -n 10                        # 12 cores
#SBATCH -N 1                         # one node ?required
#SBATCH -t 0-10:00                    # Running time of 2 hours
#SBATCH --mem 20000                  # Memory request
#SBATCH -o raw/out/bowtie_%A_%a.out          # Standard output
#SBATCH -e raw/out/bowtie_%A_%a.err          # Standard error





#######################################################################################################################
################### prepare genomes ###################
module load ngs/bowtie2
module load ngs/samtools
module load ngs/bedtools2
module load ngs/Homer
module load ngs/UCSCutils/3.4.1

#######################################################################################################################
#get filename 

cd ./raw/

FILENAME=`ls *_1.fastq.gz.fastqsanger.gz | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`
FILEBASE=`echo ${FILENAME} | sed -e "s/_1.fastq.gz.fastqsanger.gz//g"`

#######################################################################################################################
################### bowtie2 alignment ###################




	# align to dmel and dvir regular genome - need to get read counts for normalization  

bowtie_index="/work/data/genomes/fly/Drosophila_melanogaster/UCSC/dm6/Sequence/Bowtie2Index/genome"
	
BOWTIE_PARAM="-p 12 --end-to-end --very-sensitive --no-unal --no-mixed --no-discordant -I 10 -X 500"	

bowtie2 ${BOWTIE_PARAM} -x $bowtie_index -1 ${FILEBASE}_1.fastq.gz.fastqsanger.gz -2 ${FILEBASE}_2.fastq.gz.fastqsanger.gz > ${FILEBASE}.fullgen_mel.sam 2> ${FILEBASE}.fullgen_mel.stats

samtools view -bS -@ 12 -q 10 ${FILEBASE}.fullgen_mel.sam | samtools sort -@ 12 - | tee ${FILEBASE}.fullgen.dmel.bam | samtools index - ${FILEBASE}.fullgen.dmel.bam.bai

mapped_mel=`samtools view -c -F 4 ${FILEBASE}.fullgen.dmel.bam | bc`

echo ${FILEBASE}.fullgen_dmel $mapped_mel >> MappedReads.txt

rm ${FILEBASE}.fullgen_mel.sam



 	##vir

bowtie_index="/work/project/becbec_005/ChIP_Seq/Drosophila_virilis_UCSC_droVir3/bt2buildvir3/genome"

BOWTIE_PARAM="-p 12 --end-to-end --very-sensitive --no-unal --no-mixed --no-discordant -I 10 -X 500"

bowtie2 ${BOWTIE_PARAM} -x $bowtie_index -1 ${FILEBASE}_1.fastq.gz.fastqsanger.gz -2 ${FILEBASE}_2.fastq.gz.fastqsanger.gz > ${FILEBASE}.fullgen_vir.sam 2> ${FILEBASE}.fullgen_vir.stats

samtools view -bS -@ 12 -q 10 ${FILEBASE}.fullgen_vir.sam | samtools sort -@ 12 - | tee ${FILEBASE}.fullgen.dvir.bam | samtools index - ${FILEBASE}.fullgen.dvir.bam.bai

mapped_vir=`samtools view -c -F 4 ${FILEBASE}.fullgen.dvir.bam | bc`

echo ${FILEBASE}.fullgen_dvir $mapped_vir >> MappedReads.txt

rm ${FILEBASE}.fullgen_vir.sam



        # align to dmel - rep genmel - rep genome -

bowtie_index="/work/project/becbec_005/ChIP_Seq/MultiMapper_Analysis/Scripts/genome/repetitivegenome"

bowtie2 ${BOWTIE_PARAM} -x $bowtie_index -1 ${FILEBASE}_1.fastq.gz.fastqsanger.gz -2 ${FILEBASE}_2.fastq.gz.fastqsanger.gz > ${FILEBASE}.repgen.sam 2> ${FILEBASE}.repgen.stats

samtools view -bS -@ 12 -q 10 ${FILEBASE}.repgen.sam | samtools sort -@ 12 - | tee ${FILEBASE}.repgen.bam | samtools index - ${FILEBASE}.repgen.bam.bai

rm ${FILEBASE}.repgen.sam



#######################################################################################################################


#######################################################################################################################
################### makeTagDirectory ###################



cd ../Homer
mv ../raw/${FILEBASE}.*.bam .
mv ../raw/${FILEBASE}.*.bai .





        #dir_base=`echo $base | sed 's/_[G,A,T,C][G,A,T,C][G,A,T,C][G,A,T,C].*//'`

        # calculate total number of mapped reads

#mapped=`samtools view -c -F 4 ${FILEBASE}.dvir.bam | bc`

        # use 3 kinds of normalization -> normaization to dmel fullgen reads, dvir fullgen reads and dmel repeatgenome total reads

makeTagDirectory ${FILEBASE}.dvir_fullgen_Norm.dir -single -totalReads $mapped_vir ${FILEBASE}.repgen.bam

makeTagDirectory ${FILEBASE}.dmel_fullgen_Norm.dir -single -totalReads $mapped_mel ${FILEBASE}.repgen.bam

makeTagDirectory ${FILEBASE}.dmel_repgen_Norm.dir ${FILEBASE}.repgen.bam




#######################################################################################################################
