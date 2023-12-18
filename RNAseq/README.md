# Universal Snakemake pipeline for processing paired-end RNA-seq data deposited at SRA

This pipeline will download the data, align the reads agaist any chosen genome and report gene counts as well as TPM for all samples.

Multiple runs (technical replicates) will be merged.


## Setup

* clone this directory. 
* download a RunTable from SRA that comprises the samples of interest and place it into this directory.
* edit the file config.yaml providing the web locations of genome fasta and gff file.
* create a conda environment specified in sra_star_rsem.yaml.

```
conda env create --name sra_star_rsem --file sra_star_rsem.yaml
```

## Running

activate the conda environment

```
conda activate sra_star_rsem
```

### pre-flight checks

```
snakemake -np
snakemake --dag | dot -Tsvg > dag.svg
```

### on local machine
with 16 cores

```
snakemake --cores 16
```
### on the cluster

execution has to be started on master node as internet connection is needed for genome downlaod and SRA access.

```
snakemake -j 999 --cluster-config cluster.json --cluster "sbatch -p {cluster.partition} -n {cluster.n} --mem {cluster.mem}"
````
