### Script for generation of intermediate files

######### initialize general packages

rm(list=ls())


library(wesanderson)
library(stringr)
library(rtracklayer)
library(GenomicFeatures)
library(RColorBrewer)
library(matrixStats)
library(grid)
library(gridBase)
library(gridExtra)
library(tsTools)
library(BiocParallel)
library(GenomicAlignments)
library(LSD)
library(csaw)
library(pheatmap)
library(Vennerable)
library(writexl)
library(LSD)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(IRanges)
library(ShortRead)
library(ComplexHeatmap)
library(circlize)
library(dendextend)
library(genefilter)
library(tidyverse)
library(zoo)
library(matrixStats)
library(tsTools)
library(Vennerable)
library(HelpersforChIPSeq)
library(ggpmisc)
library(plotly)

source("./functions/functions.R")
source("./functions/chromoMap_custom.R")

`%nin%` = Negate(`%in%`)

######################## load genomes and genebodyfeatures

#GeneBodyFeatures <- readRDS(file="./GeneBodyFeatures.rds")


my_chromosomes_dm6 <- c("chr2L","chr2R","chr3L","chr3R","chrX","chrY","chr4")

my_lengths_dm6 <- seqlengths(keepSeqlevels(TxDb.Dmelanogaster.UCSC.dm6.ensGene, my_chromosomes_dm6, pruning.mode = "coarse"))

### load detailed gene annotaions from geome.gtf -- BDGP6 v104 ENSEMBL
my_allgenes <- makeTxDbFromGFF("./genome.gtf",format="gtf")
my_allgenes <- genes(my_allgenes)
seqlevelsStyle(my_allgenes) <- "UCSC"
my_allgenes <- keepSeqlevels(my_allgenes,  my_chromosomes, pruning.mode="coarse")

my_allgenes <- my_allgenes[grepl("FBgn",my_allgenes$gene_id)] ### filter out genes without Fbgn ID -- these are typically localized transposons etc


#### need to load previously generated bigwig coverages from 3 replicates, make average tables and merge coverages 

coverages_rep1 <- list.files(path = ".rep1/coverage/",pattern = "^coverage.*rda$", full.names = T)
coverages_rep2 <- list.files(path = ".rep2/coverage/",pattern = "^coverage.*rda$", full.names = T)
coverages_rep3 <- list.files(path = ".rep3/coverage/",pattern = "^coverage.*rda$", full.names = T)

### make average tables - load Rdata files with correct prefix

loadRData <- function(filePath){
  #loads an RData file, and returns it
  load(filePath)
  get(ls()[ls() != "filePath"])
}


## load each experiment from separate parent directories
i <- 1

for (i in seq_along(coverages_rep1)){
  
  fileName <- gsub(".rda","",gsub(".*coverage.","coverage.Rep1.",coverages_rep1[i]))
  assign(fileName,loadRData(coverages_rep1[i]),envir = .GlobalEnv)
}

i <- 1

for (i in seq_along(coverages_rep2)){
  
  fileName <- gsub(".rda","",gsub(".*coverage.","coverage.Rep2.",coverages_rep2[i]))
  assign(fileName,loadRData(coverages_rep2[i]),envir = .GlobalEnv)
}

i <- 1

for (i in seq_along(coverages_rep3)){
  
  fileName <- gsub(".rda","",gsub(".*coverage.","coverage.Rep3.",coverages_rep3[i]))
  assign(fileName,loadRData(coverages_rep3[i]),envir = .GlobalEnv)
}

## merge replicates 
my_covs <- ls(pattern = "^coverage")

#### define condiitons and IPs 
conditions <- c("GFP|GST","Ash1","Set2","NSD","Comb") 
chips <- c("Msl3","JASPer","K36Me1","K36Me2","K36Me3","NSD","ePol") 

i <- 1
j <- 1

for (i in seq_along(conditions)){
  my_covs_conds <- my_covs[grep(conditions[i],my_covs)]
  
  for (j in seq_along(chips)){
    my_covs_conds_chip <- my_covs_conds[grep(chips[j],my_covs_conds)]
    my_name <- gsub("coverage.Rep..","",my_covs_conds_chip[1])
    
    nreps <- length(my_covs_conds_chip)
    
    ## average coverages over replicates -- K36me1, NSD and ePol are only 2 replicates so special condition 
    
    if(nreps==2){
      my_cov <- (get(my_covs_conds_chip[1]) + get(my_covs_conds_chip[2]))/2
      assign(paste0("coverage.Merged.",my_name),my_cov)
      export.bw(con=paste0("./merged_bigwig/coverage.Merged.",my_name,".bw"), object=my_cov)
    }else{
      my_cov <- (get(my_covs_conds_chip[1]) + get(my_covs_conds_chip[2]) + get(my_covs_conds_chip[3]))/3
      assign(paste0("coverage.Merged.",my_name),my_cov)
      export.bw(con=paste0("./merged_bigwig/coverage.Merged.",my_name,".bw"), object=my_cov)
    }
    
  }
}


################################################################################################
############ Generate lists of genebody average coverages for all drosophila genes ###############
################################################################################################


############ Use this code to generate genebody average ChIP signal for factors of interest -- For example, KDM4A and KDM2 shown in Fig S8A 
###### can be later combined with cluster assignment for the genes to visualize cluster/supercluster trends

my_covs_v2 <- ls(pattern="^coverage.Rep")
my_covs_v2 <- my_covs_v2[!grepl(paste(chips,collapse = "|"),my_covs_v2)]
my_covs_v2 <-c(my_covs_v2,ls(pattern="^coverage.Merged"))


#### make gene body average matrices 

mat_path <- file.path("./GB_list_avg")

if(!dir.exists(mat_path))
  dir.create(mat_path)

i = 1
j = 1

genesets <- c("my_allgenes")

### parallelize code 

parallel::mclapply(seq_along(my_covs_v2), mc.cores = 4, FUN = function(i){
  for (j in seq_along(genesets)){
    
    my_name <-  paste("mat.GB",gsub("coverage.","", my_covs_v2[i]), sep=".")
    
    my_cov <- get(my_covs_v2[i])
    
    genelist <- get(genesets[j])
    
    ### subset coverages over gene bodies
    
    my_mat <- my_cov[genelist]
    
    names(my_mat) <- genelist$gene_id
    
    ##### calculate mean and log2 transform
    
    my_mat <- lapply(my_mat, function(x){ mean(x)})
    
    my_mat <- as.matrix(log2(unlist(my_mat)+0.001))
    
    assign(my_name, my_mat)
    
    fname <- file.path(mat_path, my_name)
    save(list = my_name, file = paste0(fname, ".rda"))
    print(paste(my_name,"created"))
    
  }
})

#### load saved genebody lists 

my_mat_file_names <- list.files(path = "./GB_list_avg/",pattern = "^mat\\..*rda$")

for(i in seq_along(my_mat_file_names)){
  
  load(file.path("./GB_list_avg/", my_mat_file_names[i]))
  
}


## make average table -- Note that below is a general purpose code for averaging arbitrarily different windows  
## example - 1kb windows around TSS, 500 windows of scaled genebodies etc. -- In our case, we use a singular window representing averaged genebody signal 

##################make average tables for GB
my_sub_range <- 1 #1 is for GB
my_mats <- ls(pattern = "^mat.GB\\.")
site <- "GB"


my_mat_avg <- my_mats
ave.table <- matrix(nrow = nrow(as.matrix(unlist(get(my_mat_avg[1])))), ncol = length(my_mat_avg))
colnames(ave.table) <- gsub("mat.GB.", "",my_mat_avg)
rownames(ave.table) <- rownames(as.matrix(unlist(get(my_mat_avg[1]))))
  
  
################################## Below code just aggregates all lists representing genebody averaged chip signals into a single table
  
j=1
  
for(j in 1:ncol(ave.table)){
    
  my_mat <- as.matrix(unlist(get(my_mat_avg[j])))
    
  stopifnot(
    identical(rownames(my_mat), rownames(as.matrix(unlist(get(my_mat_avg[1]))))) &
      identical(colnames(ave.table)[j], gsub("mat.GB.", "",my_mat_avg)[j])
    )
    
    
    if(site == "TSS"){my_sub_range <- (ncol(my_mat)/2):(ncol(my_mat)) }
    if(site == "TTS"){my_sub_range <- (1:(ncol(my_mat)/2))}
    if(site == "GB"){my_sub_range <- 1}
    
    ave.table[,j] <- rowMeans(as.matrix(my_mat[,my_sub_range]))
}


#### Use ave.table for clustering and other downstream analyses  
