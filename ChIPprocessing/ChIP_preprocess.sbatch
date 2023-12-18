#!/bin/sh
# 
# bowtie.sbatch
#
#SBATCH -J bowtie_array   # A single job name for the array
#SBATCH -p slim18                  # Partition
#SBATCH -n 10                        # 12 cores
#SBATCH -N 1                         # one node ?required
#SBATCH -t 0-10:00                    # Running time of 2 hours
#SBATCH --mem 20000                  # Memory request
#SBATCH -o ./processed/out/bowtie_%A_%a.out          # Standard output
#SBATCH -e ./processed/out/bowtie_%A_%a.err          # Standard error





#######################################################################################################################
################### prepare genomes ###################
module load ngs/bowtie2
module load ngs/samtools
module load ngs/bedtools2
module load ngs/Homer
module load ngs/UCSCutils/3.4.1

#######################################################################################################################
#get filename 

cd ./processed/

FILENAME=`ls *_1.fastq.gz.fastqsanger.gz | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`
FILEBASE=`echo ${FILENAME} | sed -e "s/_1.fastq.gz.fastqsanger.gz//g"`

#######################################################################################################################
################### bowtie2 alignment ###################




	# align to dmel

bowtie_index="/work/data/genomes/fly/Drosophila_melanogaster/UCSC/dm6/Sequence/Bowtie2Index/genome"
	
BOWTIE_PARAM="-p 12 --end-to-end --very-sensitive --no-unal --no-mixed --no-discordant -I 10 -X 500"	

bowtie2 ${BOWTIE_PARAM} -x $bowtie_index -1 ${FILEBASE}_1.fastq.gz.fastqsanger.gz -2 ${FILEBASE}_2.fastq.gz.fastqsanger.gz > ${FILEBASE}.dmel.sam 2> ${FILEBASE}.dmel.stats

samtools view -bS -@ 12 -q 10 ${FILEBASE}.dmel.sam | samtools sort -@ 12 - | tee ${FILEBASE}.dmel.bam | samtools index - ${FILEBASE}.dmel.bam.bai

mapped=`samtools view -c -F 4 ${FILEBASE}.dmel.bam | bc` 

echo ${FILEBASE}.dmel $mapped >> MappedReads.txt 

rm ${FILEBASE}.dmel.sam 

        # align to dvir

bowtie_index="/work/project/becbec_005/ChIP_Seq/Drosophila_virilis_UCSC_droVir3/bt2buildvir3/genome"

bowtie2 ${BOWTIE_PARAM} -x $bowtie_index -1 ${FILEBASE}_1.fastq.gz.fastqsanger.gz -2 ${FILEBASE}_2.fastq.gz.fastqsanger.gz > ${FILEBASE}.dvir.sam 2> ${FILEBASE}.dvir.stats

samtools view -bS -@ 12 -q 10 ${FILEBASE}.dvir.sam | samtools sort -@ 12 - | tee ${FILEBASE}.dvir.bam | samtools index - ${FILEBASE}.dvir.bam.bai

mapped=`samtools view -c -F 4 ${FILEBASE}.dvir.bam | bc` 

echo ${FILEBASE}.dvir $mapped >> MappedReads.txt 

rm ${FILEBASE}.dvir.sam





#######################################################################################################################


#######################################################################################################################



cd ../Homer

################################################

mv ../processed/${FILEBASE}.*.bam .
mv ../processed/${FILEBASE}.*.bai .



#######################################################################################################################
################### makeTagDirectory ###################




        #dir_base=`echo $base | sed 's/_[G,A,T,C][G,A,T,C][G,A,T,C][G,A,T,C].*//'`

        # calculate total number of mapped reads

mapped=`samtools view -c -F 4 ${FILEBASE}.dvir.bam | bc`

        # dmel reads scaled according to matched virilis sample  for normalization

makeTagDirectory ${FILEBASE}.dvirNorm.dir -single -totalReads $mapped ${FILEBASE}.dmel.bam

	# dmel reads scaled to total library size for normalization 

makeTagDirectory ${FILEBASE}.dmelNorm.dir ${FILEBASE}.dmel.bam




######################################################################################################################
