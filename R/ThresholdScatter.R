### Script to generate scatter plots comparing the relationship between reader binding changes (Z-score) with absolute K36me2/3 scores upon RNAi
### Also formally defines threshold 
### As shown in Figures S7B,C,D

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

all_genes <- my_allgenes ### alternate name compatible with other scripts 

##### load ave.table objects for readers and K36me1/2/3 generated for clustering heatmaps and other downstream analyses in GenerateIntermediateFiles.R
#### this object contains genebody averaged coverages for all chips in all conditions 

#### generate temp dataframe that contains JASPer changes (Z-scores) along with aboslute K36me2/3 scores for Set2 RNAi condition only -- This is limited to Supercluster 1 genes

temp <- data.frame(Set2_K36Me2=ave.table_clustering[,c("Merged.Set2_K36Me2")],Set2_K36Me3=ave.table_clustering[,c("Merged.Set2_K36Me3")]) %>% rownames_to_column("gene_id") %>% right_join(ave.table_readers_zscore %>% select(gene_id,Set2_JASPer,K36Me_bound_new_v2_Clusters_merged), by="gene_id") %>% filter(K36Me_bound_new_v2_Clusters_merged=="1") %>% select(-gene_id,-K36Me_bound_new_v2_Clusters_merged)

##### What is the nature of relationship between JASPer changes and K36me2/3 signal in Set2 RNAi condition? Linear or higher order? 

######### Test polynomial regression for few degrees with K-fold cross validations 

temp_shuffled <- temp[sample(nrow(temp)),]  ##randomly shuffle to get rid of sampling bias 

k <- 10   ## K for k-fold cross validation

degree <- 5 ### number of maximum degrees to fit 

folds <- cut(seq(1,nrow(temp_shuffled)),breaks=k,labels=FALSE)   ###create k folds representing testing sequence

mse = matrix(data=NA,nrow=k,ncol=degree)     ## create empty matrix to store mse

### perform k fold cross validation

for(i in 1:k){
  
  #define training and testing data
  testIndexes <- which(folds==i,arr.ind=TRUE)
  testData <- temp_shuffled[testIndexes, ]
  trainData <- temp_shuffled[-testIndexes, ]
  
  #use k-fold cv to evaluate models
  for (j in 1:degree){
    fit.train = lm(Set2_JASPer ~ polym(Set2_K36Me2,Set2_K36Me3,degree=j),data=trainData)
    fit.test = predict(fit.train, newdata=testData)
    mse[i,j] = mean((fit.test-testData$Set2_JASPer)^2) 
  }
}

colMeans(mse) ### 2 degree model has the best fit - but only slightly better RMSE - Not worth the additional complexity !

####### Stick to simple linear regression 

pdf(file="./JASPerZ_K36me23_Set2KD_SCI.pdf",width=5,height=3)
### plot JASPer Z vs Me3, color by Me2. Fit a smooth lm line based on our previous interpretations
ggplot(temp,aes(x=Set2_K36Me3,y=Set2_JASPer,col=Set2_K36Me2)) + geom_point() + scale_color_gradientn(colours=viridis(n=100),limits=c(-1,2))+ geom_smooth(method='lm')  + theme_classic() + geom_hline(yintercept=0,lty="dotted",size=1) + geom_vline(xintercept=1.35,lty="dotted",size=1, color="red") + xlab("log2(Set2 KD Me3/Inp)") + ylab("JASPer Set2 KD Z-score")
dev.off()

#### These results suggest Me3 must fall below 1.35 ChIP/Inp to start losing JASPer in Set2 RNAi 

### Does the same threshold apply in JASPer at SCIIA genes upon NSD knockdown? 

temp <- data.frame(NSD_K36Me2=ave.table_clustering[,c("Merged.NSD_K36Me2")],NSD_K36Me3=ave.table_clustering[,c("Merged.NSD_K36Me3")]) %>% rownames_to_column("gene_id") %>% right_join(ave.table_readers_zscore %>% select(gene_id,NSD_JASPer,K36Me_bound_new_v2_Clusters_merged), by="gene_id") %>% filter(K36Me_bound_new_v2_Clusters_merged=="2_Eu") %>% select(-gene_id,-K36Me_bound_new_v2_Clusters_merged)

pdf(file="./JASPerZ_K36me23_NSDKD_SCIIA.pdf",width=5,height=3)
ggplot(temp,aes(x=NSD_K36Me3,y=NSD_JASPer, col=NSD_K36Me2)) + geom_point() + scale_color_gradientn(colours=viridis(n=100),limits=c(-1,2))+ geom_smooth(method='lm')  + theme_classic() + geom_hline(yintercept=0,lty="dotted",size=1) + geom_vline(xintercept=1.35,lty="dotted",size=1, color="red") + xlab("log2(NSD KD Me3/Inp)") + ylab("JASPer NSD KD Z-score")
dev.off()

### SAME THRESHOLD WORKS! 

#### WHAT ABOUT MSL3 AT SET2-DEP GENES (RESTRICTED TO CHR X) ? 

temp <- data.frame(Set2_K36Me2=ave.table_clustering[,c("Merged.Set2_K36Me2")],Set2_K36Me3=ave.table_clustering[,c("Merged.Set2_K36Me3")]) %>% rownames_to_column("gene_id") %>% right_join(ave.table_readers_zscore %>% select(gene_id,Set2_Msl3,K36Me_bound_new_v2_Clusters_merged), by="gene_id") %>% filter(K36Me_bound_new_v2_Clusters_merged=="1" &gene_id %in% all_genes[seqnames(all_genes)=="chrX"]$gene_id) %>% select(-gene_id,-K36Me_bound_new_v2_Clusters_merged)
pdf(file="./Msl3Z_K36me23_Set2KD_SCI_X.pdf",width=5,height=3)
ggplot(temp,aes(x=Set2_K36Me3,y=Set2_Msl3, col=Set2_K36Me2)) + geom_point() + scale_color_gradientn(colours=viridis(n=100),limits=c(-1,2))+ geom_smooth(method='lm')  + theme_classic() + geom_hline(yintercept=0,lty="dotted",size=1) + geom_vline(xintercept=1.35,lty="dotted",size=1, color="red") + xlab("log2(Set2 KD Me3/Inp)") + ylab("Msl3 Set2 KD Z-score")
dev.off()
