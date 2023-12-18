#!/bin/bash
#SBATCH --time=02-10:00
#SBATCH --job-name="Intron_Exon"
#SBATCH --mem=50GB
#SBATCH --output=long-job.out
#SBATCH -p slim16

module load R/3.6.0

### insert appropriate script name for cluster execution

Rscript --vanilla R_cluster.R
