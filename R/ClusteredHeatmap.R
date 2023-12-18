### Script for clustered heatmap shown in Figure 4A
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


my_chromosomes <- c("chr2L","chr2R","chr3L","chr3R","chrX","chrY","chr4")

my_lengths <- seqlengths(keepSeqlevels(TxDb.Dmelanogaster.UCSC.dm6.ensGene, my_chromosomes_dm6, pruning.mode = "coarse"))

### load detailed gene annotaions from geome.gtf -- BDGP6 v104 ENSEMBL
my_allgenes <- makeTxDbFromGFF("./genome.gtf",format="gtf")
my_allgenes <- genes(my_allgenes)
seqlevelsStyle(my_allgenes) <- "UCSC"
my_allgenes <- keepSeqlevels(my_allgenes,  my_chromosomes, pruning.mode="coarse")

my_allgenes <- my_allgenes[grepl("FBgn",my_allgenes$gene_id)] ### filter out genes without Fbgn ID -- these are typically localized transposons etc


#### LOAD ave.table OBJECT GENERATED FROM GENERATEINTERMEDIATEFILES.R


#### clustered heatmaps 

ave.table_clustering <- ave.table[rownames(ave.table) %in% GeneBodyFeatures[GeneBodyFeatures$K36Me_bound_new_v2==1,]$gene_id,grepl("K36me|K36Me",colnames(ave.table))][,c(4,8,6,2,3,7,5,1,11,14,13,10,9)]
ave.table_clustering_scaled <- cbind(t(scale(t(ave.table_clustering[,1:4]))), t(scale(t(ave.table_clustering[,5:8]))),t(scale(t(ave.table_clustering[,9:13]))))

NA_genes <- ave.table_clustering_scaled[rowSums(is.na(ave.table_clustering_scaled)) > 0, ]  
#NA_genes_unscaled <- ave.table_clustering[rowSums(is.na(ave.table_clustering_scaled)) > 0, ] ### around 90 genes give NAs - exclude them for clustering

#ave.table_clustering_scaled <- ave.table_clustering_scaled[rownames(ave.table_clustering_scaled) %nin% rownames(NA_genes),]
ave.table_clustering <- ave.table_clustering[rownames(ave.table_clustering) %nin% rownames(NA_genes),]

## initialize new matrix which will be corrected for clustering - This row order can be used to plot the original uncorrected raw signal matrix

ave.table_clustering_robdist <- ave.table_clustering
col_fun_scaled <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
#col_fun <- colorRamp2(c(-3, 0, 8), c("blue", "yellow", "red"))
col_fun_unscaled <- colorRamp2(seq(-2,2.75,length=11),rev(brewer.pal(n = 11, name ="RdYlBu")))


####### load cluster dendrogram object generated from K=12, Complete linkage and Minkowski distance HM -> generated on the computational cluster using HeatMap_HPC.R

HM_dend_12k <- readRDS("HM_Dend_12k.rds")  

cluster_list <- cutree(HM_dend_12k,k=12) ## define k clusters 


pdf("HM_K36Me_final.pdf",height=15,width=7)
p <-Heatmap(ave.table_clustering, cluster_columns= F,row_labels = rep("",nrow(ave.table_clustering)),split=cluster_list,col = col_fun_unscaled,column_title="12k",row_split=factor(cluster_list, levels = c(3,8,10,7,1,5,12,11,2,6,4,9)),cluster_row_slices = F, width=20)
p
dev.off()


row_order_final <- unlist(row_order(p),use.names = F)

#### Group clusters into Superclusters 

cluster_df <- as.data.frame(cluster_list) %>% mutate(Merged_cluster=as.factor(case_when(cluster_list %in% c(3,8,10,7)~1,
                                                                              cluster_list %in% c(1,5,12,11)~2,
                                                                              cluster_list %in% c(2,4,6)~3,
                                                                              cluster_list %in% c(9)~4)))

#### annotate colors for superclusters 
p_clustermerge <- Heatmap(cluster_df[,"Merged_cluster",drop=F],cluster_columns= F,row_order=row_order_final,row_labels = rep("",nrow(cluster_df)),width=1,col=c("red","green","#1E88E5","gray"))

p_1 <- p + p_clustermerge

pdf("HM_K36Me_final_clusters.pdf",width=7,height=20)
draw(p_1)
dev.off()


########## Make heatmaps for JASPer and MSL3 using same order - shown in Fig 6A

### generate ave.table for readers 

ave.table_clustering_readers <- ave.table[rownames(ave.table) %in% rownames(ave.table_clustering),grepl("JASPer|Msl3",colnames(ave.table))]


### Select appropriate columns and cluster using same order 
p_cluster_jasper <- Heatmap(ave.table_clustering_readers[,c(3,7,5,1)],cluster_columns= F,row_order=row_order_final,row_labels = rep("",nrow(ave.table_clustering_readers)),row_split=factor(cluster_list, levels = c(3,8,10,7,1,5,12,11,2,6,4,9)),width=15,col=colorRamp2(seq(-2,2.75,length=11),rev(brewer.pal(n = 11, name ="RdYlBu"))))
p_cluster_msl3 <- Heatmap(ave.table_clustering_readers[,c(4,8,6,2)],cluster_columns= F,row_order=row_order_final,row_labels = rep("",nrow(ave.table_clustering_readers)),row_split=factor(cluster_list, levels = c(3,8,10,7,1,5,12,11,2,6,4,9)),width=15,col=colorRamp2(seq(-2,4,length=11),rev(brewer.pal(n = 11, name ="RdYlBu"))))
p_reader <- p_cluster_jasper + p_cluster_msl3

pdf("HM_Readers_final_clusters.pdf",width=7,height=20)
draw(p_reader)
dev.off()
