### Script for calculating HMT signal using 2kb windows around signal maxima on genebody- Fig 5A
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

#genebody features contains useful information for all genes including cluster membership etc

GeneBodyFeatures <- readRDS(file="./GeneBodyFeatures.rds")


my_chromosomes <- c("chr2L","chr2R","chr3L","chr3R","chrX","chrY","chr4")

my_lengths <- seqlengths(keepSeqlevels(TxDb.Dmelanogaster.UCSC.dm6.ensGene, my_chromosomes_dm6, pruning.mode = "coarse"))

### load detailed gene annotaions from geome.gtf -- BDGP6 v104 ENSEMBL
my_allgenes <- makeTxDbFromGFF("./genome.gtf",format="gtf")
my_allgenes <- genes(my_allgenes)
seqlevelsStyle(my_allgenes) <- "UCSC"
my_allgenes <- keepSeqlevels(my_allgenes,  my_chromosomes, pruning.mode="coarse")

my_allgenes <- my_allgenes[grepl("FBgn",my_allgenes$gene_id)] ### filter out genes without Fbgn ID -- these are typically localized transposons etc


#### code below edited to calulcate mean over 2kb window around the max signal point 

#### load HMT bigwig merged coverages prior to this -- Note that Ash1 profile is from Huang et al., 2017
my_covs_HMT <- c("coverage.ePol_merged.bw","coverage.GSM2443786_Ash1-in-WT.rep1_dm6.bw","coverage.NSD_merged.bw")

mat_path <- file.path("./Coverage_HMT/GB_list_Mean_aroundMax")

if(!dir.exists(mat_path))
  dir.create(mat_path)

i = 1
j = 1

genesets <- c("my_allgenes")

parallel::mclapply(seq_along(my_covs_HMT), mc.cores = 4, FUN = function(i){
  for (j in seq_along(genesets)){
    
    my_name <-  paste("mat.GB_meanMax_2kb",gsub("coverage.","", my_covs_HMT[i]), sep=".")
    
    my_cov <- get(my_covs_HMT[i])
    
    genelist <- get(genesets[j])
    
    my_mat <- my_cov[genelist]
    
    names(my_mat) <- genelist$gene_id
    
    ## above steps extracts per bp coverage scores for each gene
    ## this step below retreives which bp position within the gene has highest score (its an integer lying between 1:width of gene)
    my_mat_max <- data.frame(which.max(my_mat)) %>% rownames_to_column(var="gene_id")
    
    ## use this as a new 'center' to create a new granges object containg the coordinates of the new center (represented as a new start and end with +-1000bp around center)
    my_ranges_centerwindow <- data.frame(genelist) %>% left_join(my_mat_max,by="gene_id") %>% mutate(center=start-1+which.max.my_mat.,start=center-1000,end=center+1000) %>%
                            makeGRangesFromDataFrame(keep.extra.columns = T)
    
    seqlengths(my_ranges_centerwindow) <- my_lengths[c(1,2,3,4,7,5,6)]  ### force appropriate chromosome seqlengths onto new ranges 
    
    ######this allows you to trim ranges which lie outside seqlengths
    
    my_ranges_centerwindow <- trim(my_ranges_centerwindow)
  
    my_mat_centerwindow <- my_cov[my_ranges_centerwindow]
    
    names(my_mat_centerwindow) <- my_ranges_centerwindow$gene_id
    
    ### average signal over entire window and log2 transform as usual 
    
    my_mat_centerwindow <- lapply(my_mat_centerwindow, function(x){ mean(x)})  
    
    my_mat_centerwindow <- as.matrix(log2(unlist(my_mat_centerwindow)+0.001))
    
    assign(my_name, my_mat_centerwindow)
    
    fname <- file.path(mat_path, my_name)
    save(list = my_name, file = paste0(fname, ".rda"))
    print(paste(my_name,"created"))
    
  }
})

#### These values can be now grouped by clusters/superclusters and plotted using dplyr! 

my_mat_file_names <- list.files(path = "./Coverage_HMT/GB_list_Mean_aroundMax//",pattern = "^mat\\..*rda$")

for(i in seq_along(my_mat_file_names)){
  load(file.path("./Coverage_HMT/GB_list_Mean_aroundMax//", my_mat_file_names[i]))
  }


## make average table -- .ie. aggregate different samples into single ave.table object

##################make average tables for GB
my_sub_range <- 1 #1 is for GB
my_mats <- ls(pattern = "^mat.GB") 
site <- "GB"

my_mat_avg <- my_mats
max.table_HMT <- matrix(nrow = nrow(as.matrix(unlist(get(my_mat_avg[1])))), ncol = length(my_mat_avg))
colnames(max.table_HMT) <- gsub("mat.GB_meanMax_2kb.", "",my_mat_avg)
rownames(max.table_HMT) <- rownames(as.matrix(unlist(get(my_mat_avg[1]))))

################################## below is general purpose code for averaging arbitrarily different windows   

j=1

for(j in 1:ncol(max.table_HMT)){
  my_mat <- as.matrix(unlist(get(my_mat_avg[j])))
  stopifnot(
    identical(rownames(my_mat), rownames(as.matrix(unlist(get(my_mat_avg[1]))))) &
      identical(colnames(max.table_HMT)[j], gsub("mat.GB_meanMax_2kb.", "",my_mat_avg)[j])
  )
  
  
  if(site == "TSS"){my_sub_range <- (ncol(my_mat)/2):(ncol(my_mat)) }
  if(site == "TTS"){my_sub_range <- (1:(ncol(my_mat)/2))}
  if(site == "GB"){my_sub_range <- 1}
  
  max.table_HMT[,j] <- rowMeans(as.matrix(my_mat[,my_sub_range]))
}

#######

### define subsets of genes -- use all K36me genes--- Plot !

### select superclusters as well as 'background' set of randomly samples 3000 genes

dens_plot_list <- GeneBodyFeatures %>% filter(K36Me_bound_new_v2_Clusters_merged %in% c("1","2_Eu","2_Het","3")) %>% select(gene_id,K36Me_bound_new_v2_Clusters_merged)
dens_plot_list <- rbind(dens_plot_list,setNames(sample_n(data.frame(GeneBodyFeatures[GeneBodyFeatures$K36Me_bound_new_v2_Clusters_merged=="0",c("gene_id","K36Me_bound_new_v2_Clusters_merged")]),size=3000),colnames(dens_plot_list)))

### modify avetable to dataframe
max.table_HMT <- as.data.frame(max.table_HMT) %>% rownames_to_column(var="gene_id")
dens_plot_list <- dens_plot_list %>% left_join(max.table_HMT) %>% pivot_longer(cols=c(3,4,5),names_to="HMT_ChIP",values_to="Mean_ChIP_coverage")


pdf(file="HMTChIP_ClusterEnrichment_meanMax_around2kb.pdf",width=7,height=5)

i <- unique(dens_plot_list$HMT_ChIP)[1]

### limits are different for Ash1 profile as it was generated by Huang et al.,2017

for (i in unique(dens_plot_list$HMT_ChIP)){
  if (grepl("Ash1",i)){
    lim=c(2,8)
  } else {
    lim=c(-0.5,1)
  }
  print(dens_plot_list %>% filter(HDM_ChIP==i) %>% ggplot(aes(x=K36Me_bound_new_v2_Clusters,y=Mean_ChIP_coverage,fill=K36Me_bound_new_v2_Clusters)) + geom_violin(trim=T) + geom_boxplot(width=0.1)+theme_minimal() + ylim(lim[1],lim[2]) +labs(col="Cluster")+ ggtitle(i))
}
dev.off()
