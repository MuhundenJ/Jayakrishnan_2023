

####### Plot chromatin state enrichment of chromatin marks and proteins




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


###### Chrom State Enrichment Maps 

### read in MODENCODE 9-state chromatin annotation 
Chrom_State <- read.table("9STATE_S2_NArepl_IGV.bed")[-c(5,6,7,8)] %>% makeGRangesFromDataFrame(seqnames.field = "V1",start.field = "V2",end.field = "V3",keep.extra.columns = T)

## load files (self generated as well as MODEncode) to be plotted on enrichment map. 
#Note that this loop takes really long, so recommended to resize bedgraphs to lower resolution bins before proceeding
corr_files <- list.files(path="./BigWigs_CorrMap/",pattern=".*_resize.bedgraph")
corr_files_path <- list.files(path="./BigWigs_CorrMap/",pattern=".*_resize.bedgraph", full.names = T)


i <- 1
j <- 1

#initialize empty Dataframe

ME_ChromState_enrichment <- data.frame(matrix(nrow=9,ncol=length(corr_files)))
for (i in 1:9){
  my_region <- Chrom_State[Chrom_State$V4==i]
  for (j in seq_along(corr_files)){
    my_name <- gsub("_resize.bedgraph","",corr_files[j])
    my_bedgraph <- rtracklayer::import(corr_files_path[j])
    seqlevelsStyle(my_bedgraph) <- "UCSC"
    my_bedgraph <- keepSeqlevels(my_bedgraph, my_chromosomes_dm6, pruning.mode = "coarse")
    seqlengths(my_bedgraph) <- my_lengths_dm6
    my_cov <- coverage(my_bedgraph, weight = "score")
    ME_ChromState_enrichment[i,j] <- mean(unlist(my_cov[my_region]))
    colnames(ME_ChromState_enrichment)[j] <- my_name
  }
}
saveRDS(ME_ChromState_enrichment,file="./ME_ChromState_ChIPEnrichment.rda")
ME_ChromState_enrichment <- readRDS("./ME_ChromState_ChIPEnrichment.rda")

###scale to improve visuliazation 
scale(ME_ChromState_enrichment)


colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
colorRamp2(seq(-1.75,2.15,length=11),rev(brewer.pal(n = 11, name ="RdYlBu")))

### plot subsetted heatmaps of the plot using ComplexHeatmap package

ME_ChromState_HM1 <- ComplexHeatmap::Heatmap(scale(ME_ChromState_enrichment[,c(9,17,18)]),cluster_rows = F, cluster_columns = F, name="log2(ChIP/Input) Z-scores",col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100), width= 0.5)
ME_ChromState_HM2 <- ComplexHeatmap::Heatmap(scale(ME_ChromState_enrichment[,c(1,4,23)]),cluster_rows = F, cluster_columns = F, name="log2(ChIP/Input) Z-scores",col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),width=0.5)
ME_ChromState_HM3 <- ComplexHeatmap::Heatmap(scale(ME_ChromState_enrichment[,-c(1,4,23,9,17,18)]),cluster_rows = F, cluster_columns = T, name="log2(ChIP/Input) Z-scores",col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),width =2.5)
ME_ChromState_HM <- ME_ChromState_HM1 + ME_ChromState_HM2 + ME_ChromState_HM3
pdf(file="ME_ChromState_enrichment.pdf",height=5,width=10)
draw(ME_ChromState_HM)
dev.off()


