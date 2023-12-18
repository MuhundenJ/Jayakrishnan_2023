# conda env create --name sra_star_rsem --file sra_star_rsem.yaml
conda activate sra_star_rsem

#snakemake -np
#snakemake --dag | dot -Tsvg > dag.svg

## on local machine
#snakemake --cores 16

## on the cluster
snakemake -j 999 --cluster-config cluster.json --cluster "sbatch -p {cluster.partition} -n {cluster.n} --mem {cluster.mem}"
