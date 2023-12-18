#! /bin/bash
### run bowtie array ###

cd ./processed/
FILES=($(ls -1 *_1.fastq.gz.fastqsanger.gz))

# get size of array
NUMFASTQ=${#FILES[@]}

mkdir -p out
mkdir ../Homer

cd ../

# now submit to SLURM
if [ $NUMFASTQ -ge 0 ]; then
	sbatch --array=1-$NUMFASTQ ChIP_preprocess.sbatch
fi
