

####### Plot chromoMaps for Control ChIP profiles used in Fig 1B




rm(list=ls())

##### initialize general packages

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
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(tsTools)
library(Vennerable)
library(HelpersforChIPSeq)

source("./functions/functions.R")
source("./functions/chromoMap_custom.R")

######################## load genomes and genebodyfeatures

#GeneBodyFeatures <- readRDS(file="./GeneBodyFeatures.rds")


my_chromosomes_dm6 <- c("chr2L","chr2R","chr3L","chr3R","chrX","chrY","chr4")

my_lengths_dm6 <- seqlengths(keepSeqlevels(TxDb.Dmelanogaster.UCSC.dm6.ensGene, my_chromosomes_dm6, pruning.mode = "coarse"))

### load detailed gene annotaions from geome.gtf -- BDGP6 v104 ENSEMBL
my_allgenes_dm6 <- makeTxDbFromGFF("./genome.gtf",format="gtf")
my_allgenes_dm6 <- genes(my_allgenes_dm6)
seqlevelsStyle(my_allgenes_dm6) <- "UCSC"
my_allgenes_dm6 <- keepSeqlevels(my_allgenes_dm6,  my_chromosomes_dm6, pruning.mode="coarse")

#### load additional packages

library("chromoMap")
library("htmlwidgets")
chr_file <- read.delim("Fly/Further_Analysis/dm6.chrom.sizes",header=F)

### generate annotation file as described by package along with 'pseudo'-chromocenters

chr_file <- chr_file %>% filter(V1 %in% my_chromosomes_dm6) 
chr_file <- add_column(chr_file,d=rep(1,nrow(chr_file)),.after = 1) %>% mutate(V4=c(1,28110227,1,23542271,23513712,NA,1)) 
chr_file <- chr_file[order(chr_file$V1),]
write_delim(chr_file,file="Fly/Further_Analysis/chr_file_anno.txt",quote="none",col_names = FALSE,delim="\t")


### annotation files -> read in selected bedgraphs as dataframes
###select interested ChIPs for visualization 
chromoMaps_list <- c("Ash1","ePol","H3K36Me1","JASPer","Msl3","K36Me2","K36Me3","H3K9me2","NSD")

chromoMaps_files <- list.files(path="./BigWigs_CorrMap/", full.names=T)

i <- 7


#### plot and save chromoMaps
#### NOTE -- you can add parameter plots="bar" to chromoMap to obtain barplots (coverage) of ChIP signal overlayed on top of chromoMap, like shown in Fig 3

for (i in seq_along(chromoMaps_list)){
  
  anno_chip <- read.delim(chromoMaps_files[grep(chromoMaps_list[i],chromoMaps_files)],header=F)
  anno_chip <- add_column(anno_chip, name=paste0(anno_chip$V1,":",anno_chip$V2,"-",anno_chip$V3),.after=0) %>% filter(V1 %in% my_chromosomes_dm6)
  write_delim(anno_chip, file=paste0("./chromoMaps_2023/",chromoMaps_list[i],"_anno.text"),quote="none",col_names=F,delim="\t")
  
  
  x <- chromoMap("./chromoMaps_2023/chr_file_anno.txt",paste0("Fly/Further_Analysis/chromoMaps_2023/",chromoMaps_list[i],"_anno.text"),data_based_color_map=T,data_type="numeric",data_colors = list(c("white","orange","red")), title=chromoMaps_list[i], legend=T,title_font_size = 15, lg_y=350)
 
  saveWidget(x,paste0("./chromoMaps_2023/",chromoMaps_list[i],".html"),selfcontained = F)
  
  }
