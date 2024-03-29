#### Exploratory analyses for Figure S9 ### 

### Note -- For expression analysis Fig S9F, see Revision_FigS7 script 

rm(list=ls())

setwd("~/Desktop/mount/ChIP_Seq/240229_Revision")

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
library(ggforce)
source("./functions/functions.R")

`%nin%` = Negate(`%in%`)

###################### Load necessary annotations


my_chromosomes <- c("chr2L","chr2R","chr3L","chr3R","chrX","chrY","chr4")


GeneBodyFeatures <- readRDS("../../GeneBodyFeatures_v2/GeneBodyFeatures_expanded_v20.rds")  ## load most recent GeneBodyFeatures file

my_genes <- genes(TxDb.Dmelanogaster.UCSC.dm6.ensGene)

my_genes <- keepSeqlevels(my_genes,my_chromosomes,pruning.mode = "coarse")

my_lengths <- seqlengths(keepSeqlevels(TxDb.Dmelanogaster.UCSC.dm6.ensGene, my_chromosomes))

### Load blacklisted regions - These are regions with no uniquely mapping reads in MNase input samples (however, upon generating input normlaized
### coverages. it is given a pseudocount which can sometimes affect z-score calculations - especially if we are focussing on heterochromatic regions)

BlackList <- read.delim("../../BlackList/Resized_NoNorm_bedgraphs/BlackList_NoInpReads_1kb.bed",header=F,col.names = c("chr","start","end")) %>% makeGRangesFromDataFrame()


### Load peaks for K27 or K9 methylation - Note that these peaks were called for each RNAi condition and then merged

K27_regions_bed <- read.delim("./peaks/K27me_regions.bed",header=F) %>% makeGRangesFromDataFrame(start.field = "V2",end.field = "V3",seqnames="V1",keep.extra.columns = T)
K9_regions_bed <- read.delim("./peaks/K9me_regions.bed",header=F) %>% makeGRangesFromDataFrame(start.field = "V2",end.field = "V3",seqnames="V1",keep.extra.columns = T)


### Read in relevant 5-kb resized bedgraph files for z-score calculation

my_file_directory <- list.files("./bigwig_resize/",pattern="*_resize_5kb.bedgraph",full.names = T)
my_file_names <-list.files("./bigwig_resize/",pattern="*_resize_5kb.bedgraph")


my_conditions <- c("Ash1","Set2","NSD")

my_chips <- c("K27","K9")


i <- 1

for (i in seq_along(my_chips)){
  
  my_files_chip <- my_file_directory[grepl(my_chips[i],my_file_directory)]
  
  #### tailor filenames 
  my_filenames_chip <- gsub(paste0("_",my_chips[i],".*"),"",my_file_names[grepl(my_chips[i],my_file_names)])
  
  #### read all files using vroom, convert to df and filter useless bins 
  
  raw_scores <- lapply(my_files_chip, function(x) {read.table(file = x, header = F, sep ="\t")})
  
  names(raw_scores) <- my_filenames_chip
  
  raw_scores <- map_df(raw_scores,data.frame,.id="name") %>% pivot_wider(names_from="name",values_from="V4") %>% makeGRangesFromDataFrame(seqnames.field ="V1",start.field = "V2",end.field = "V3",keep.extra.columns = T)  ### this already does a outerjoin based on coordinates !
  
  if (my_chips[i]=="K9"){
    my_regions <- K9_regions_bed
  }else{
    my_regions <- K27_regions_bed
  }
  
  ### standard filtering out of no input read regions. Followed by filtering only regions that overlap a ChIP peak
  
  raw_scores <- raw_scores[!overlapsAny(raw_scores,BlackList)]
  raw_scores_filt <- as.data.frame(raw_scores[overlapsAny(raw_scores,my_regions,minoverlap = 1500)])
  raw_scores_negativeregions <- as.data.frame(raw_scores[!overlapsAny(raw_scores,my_regions,minoverlap = 1500)])
  
  raw_scores_filt <- raw_scores_filt %>% replace(is.na(.),0)
  ### Calculate z-scores for appropriate ChIP-RNAi combination as described by Chaouch et al., 2022 using genebody averaged coverage instead of RPKPM 
  
  
  if (my_chips[i]=="K27"){
    z_scores_df <- raw_scores_filt %>% transmute(chr=seqnames,start=start,end=end,
                                                 avg=rowMeans(select(raw_scores_filt,matches("GST|Ash1|Set2|NSD"))),
                                                 sdev=rowSds(select(raw_scores_filt,matches("GST|Ash1|Set2|NSD"))),
                                                 GST_raw = GST,
                                                 Ash1=(Ash1-GST)/sqrt(Ash1+GST),
                                                 Set2=(Set2-GST)/sqrt(Set2+GST),
                                                 NSD=(NSD-GST)/sqrt(NSD+GST))
    
    ### may filter low GST signal or avg signal
    
  }else{
    z_scores_df <- raw_scores_filt %>% transmute(chr=seqnames,start=start,end=end,
                                                 GST_raw = GST,
                                                 NSD=(NSD-GST)/sqrt(NSD+GST),
                                                 avg=rowMeans(select(raw_scores_filt,matches("GST|Ash1|Set2|NSD"))),
                                                 sdev=rowSds(select(raw_scores_filt,matches("GST|Ash1|Set2|NSD"))))
    
    
    
  }
  
  assign(paste0(my_chips[i],"_z_scores"),z_scores_df)
  assign(paste0(my_chips[i],"_raw_scores_filt"),raw_scores_filt)
  assign(paste0(my_chips[i],"_raw_scores_negativeregions"),raw_scores_negativeregions)
  
}


####################### Preliminary Analysis for K27me3 changes with only 1 replicate !

### classify bins as unchanging or changing upon individual RNAi. Note that only one replicate is available ! - So can't use statistical
### tools like csaw to identify differentially changing regions. Instead, use arbitrary cutoffs based on viewing on IGV 

### In a single replicate, there can be fluctuations in nucleosomal signal sometimes within 5kb domains - For example, a few nucleosomes (approx 1 kb)
### might show increased K27me3 signal within a large 20-50kb domains - Its unclear whether its a reproducible effect from only replicate. 
### Better to manually adjust thresholds so you dont capture these events as changing (underestimation is better than overestimation !)

### to identify denovo formed domains (.ie increasing), add an additional threshold filtering for very low signal in Control cells (basically capture 
### regions that go from 0 signal in control to appreciable signal in RNAi)

## similarly, for lost domains (.ie. decreasing), the regions should start with substantial signal in Control and go to low signal upon RNAi

classification_vec <- as.data.frame(matrix(data = 0, nrow = nrow(K27_z_scores), ncol = 7))

colnames(classification_vec) <- c("Ash1_decrease","Ash1_increase","Set2_decrease","Set2_increase","NSD_decrease","NSD_increase","unchanged")

i <- 1

for (i in 1:nrow(K27_z_scores)){
  
  ### z-scores cutoffs were decided based on inspection of score quantiles, trial and error and IGV viewing 
  
  classification_vec[i,]$Ash1_decrease <- ifelse((K27_z_scores[i,]$Ash1<=-0.32 & K27_z_scores[i,]$GST_raw>1.25),1,0)
  classification_vec[i,]$Ash1_increase <- ifelse((K27_z_scores[i,]$Ash1>=0.275 & K27_z_scores[i,]$GST_raw<1.25),1,0)
  classification_vec[i,]$Set2_decrease <- ifelse((K27_z_scores[i,]$Set2<=-0.32 & K27_z_scores[i,]$GST_raw>1.25),1,0)
  classification_vec[i,]$Set2_increase <- ifelse((K27_z_scores[i,]$Set2>=0.275 & K27_z_scores[i,]$GST_raw<1.25),1,0)
  classification_vec[i,]$NSD_decrease <- ifelse((K27_z_scores[i,]$NSD<=-0.32 & K27_z_scores[i,]$GST_raw>1.25),1,0)
  classification_vec[i,]$NSD_increase <- ifelse((K27_z_scores[i,]$NSD>=0.275 & K27_z_scores[i,]$GST_raw<1.25),1,0)
  
}

## regions not passing threshold in any of the above conditions are classified as unchanging 
classification_vec$unchanged <- sapply(1:nrow(classification_vec), function(i) ifelse(any(classification_vec[i,1:6]==1),0,1))

K27_z_scores_v4 <- cbind(K27_z_scores,classification_vec)


########### look at counts as table 


K27_z_scores_v4 %>% filter((.$unchanged==0|(.$unchanged==1 & .$avg>1))) %>% select(10:16) %>% pivot_longer(cols=1:7,names_to = "Region_Type",values_to = "Change") %>% filter(Change==1) %>% select(Region_Type) %>% table()
K27_z_scores_v4 %>% select(10:16) %>% pivot_longer(cols=1:7,names_to = "Region_Type",values_to = "Change") %>% filter(Change==1) %>% select(Region_Type) %>% table()           


##### Are there overlapping windows between diff 'Up' categories ? 

K27_z_scores_v4$name <- paste0(K27_z_scores_v4$chr,":",K27_z_scores_v4$start,"-",K27_z_scores_v4$end)


pdf(paste0("./venn.HMTKD_K27me3_updomains.pdf"), width = 6, height = 6)

my_overlaps_list <- list(K27_z_scores_v4[K27_z_scores_v4$Ash1_increase==1,"name"],K27_z_scores_v4[K27_z_scores_v4$Set2_increase==1,"name"],K27_z_scores_v4[K27_z_scores_v4$NSD_increase==1,"name"])

names(my_overlaps_list) <- c("Ash1-KD","Set2-KD","NSD-KD")
my_overlaps_venn <- Venn(my_overlaps_list)

my_overlaps_venn_plot <- compute.Venn(my_overlaps_venn,doWeights = T)
SetLabels <- VennGetSetLabels(my_overlaps_venn_plot)
SetFaceLabels <- VennGetFaceLabels(my_overlaps_venn_plot)

my_overlaps_venn_plot <- Vennerable:::VennSetFaceLabels(my_overlaps_venn_plot,SetFaceLabels)
gp <- VennThemes(my_overlaps_venn_plot)

plot(my_overlaps_venn_plot, show = list(Faces = FALSE), gp = gp)

dev.off()

################# NOT MUCH OVERLAP !

################### Classify bins as genic SC0-III or intergenic


### in the bottom loop, we attempt to do two things: 
### 1) Bins-based analysis: For each 5kb bin, count number of SC0-III genes it contains and store it in gene_SC_K27 object. If no gene, then classify as intergenic 

### 2) Genes-based analysis: For each 5kb bin, find the gene_id of the gene it intersects and store it in the genenames_SC_K27 object. This is necessary because
### often long genes span multiple bins. Consider the example where we identify 10 bins which gain K27me3 upon Set2 RNAi and 10 bins that gain K27me3 in an Ash1 RNAi
### A proportion analysis will conclude that the similar proportion of the genome contains Set2 and Ash1 depenedent K27me3. However, imagine that 10 small genes 
### are contained with the 10 set2-dependent bins but only 1 superlong gene in contained with the 10 Ash-dependent bins. A gene based proportion analyses
### will conclude that Set2 dependent K27me3 gaining genes are 10 times more than Ash1-dep genes ! -- Both bins and genes informatio are useful for 
### biological interpretation


## loop through each window, find which genes it intersects, take the appropriate gene(s)' supercluster ID and append to K27_z_scores df 
## Initilize new df 
gene_SC_K27 <- as.data.frame(matrix(data=0,nrow=nrow(K27_z_scores_v4),ncol=7)) ### this is to store counts
colnames(gene_SC_K27) <- c("Intergenic","SC_0","SC_1","SC_2_Eu","SC_2_Het","SC_3","SC_4")

genenames_SC_K27 <- as.data.frame(matrix(data=0,nrow=nrow(K27_z_scores_v4),ncol=6)) ### this is to store names to genes - this is to later collapse bins with same genes into one to get gene counts 
### only 6 columns as there is no intergenic genes !
colnames(genenames_SC_K27) <- c("SC_0","SC_1","SC_2_Eu","SC_2_Het","SC_3","SC_4")

i <- 1

for (i in 1:nrow(K27_z_scores_v4)){
  my_range <- K27_z_scores_v4[i,1:3] %>% makeGRangesFromDataFrame()
  
  ## for a gene to be part of the window, atleast 500bp should overlap - note that I cant increase the number too much as dmel genes are kinda small, especially houskeeping genes
  ## This we way underestimate intergenic windows -- Play around with this number. Alternatively, can also try using 2.5kb bins in the beginning 
  my_overlap_genes <- my_genes[overlapsAny(my_genes,my_range,minoverlap = 500)]$gene_id
  
  ## if no overlap with genes, then integenic
  if (is_empty(my_overlap_genes)){
    gene_SC_K27[i,]$Intergenic <- 1
  } else {
    j <-1
    ## 5kb bin may overlap multiple genes. Loop through each gene to get its supercluster annotation
    ## note that they are quite likely to be the same, as supercluster annotations were decided based on K36me ChIP patterns which are frequently 
    ## same for neighbouring genes 
    for (j in 1:length(my_overlap_genes)){
      my_gene <- my_overlap_genes[j]
      ### Pull out Supercluster annotation associated with the gene_id and classify the window as containing that particular cluster gene
      if (is_empty(GeneBodyFeatures[GeneBodyFeatures$gene_id == my_gene,]$K36Me_bound_new_v2_Clusters_merged)) next
      
      ### this increments the count depending on how many genes of a particular cluster intersect with the bin
      gene_SC_K27[i,paste0("SC_",GeneBodyFeatures[GeneBodyFeatures$gene_id == my_gene,]$K36Me_bound_new_v2_Clusters_merged)] <- gene_SC_K27[i,paste0("SC_",GeneBodyFeatures[GeneBodyFeatures$gene_id == my_gene,]$K36Me_bound_new_v2_Clusters_merged)]+1
      
      ## this is to store gene_ids which will later collapse to find unique genes
      genenames_SC_K27[i,paste0("SC_",GeneBodyFeatures[GeneBodyFeatures$gene_id == my_gene,]$K36Me_bound_new_v2_Clusters_merged)] <- paste(genenames_SC_K27[i,paste0("SC_",GeneBodyFeatures[GeneBodyFeatures$gene_id == my_gene,]$K36Me_bound_new_v2_Clusters_merged)],my_gene,sep = ",")
      
    }
  }
}

### Note,SC_4 (around 200) are super long genes which have a tiny speck of K36me somewhere. Almost all SC_4 bins are unchanging polycomb -> Force these counts (names) to SC_0 

gene_SC_K27 <- gene_SC_K27 %>% mutate(SC_0=SC_0+SC_4) %>% select(-c("SC_4"))
genenames_SC_K27 <- genenames_SC_K27 %>% mutate(SC_0=paste(SC_0,SC_4,sep=",")) %>% select(-c("SC_4"))

################ BINS BASED ANALYSES ##################

#### Note that bins based analyses is crucial as many K27me3 domains are intergenic (see Main Fig 1C -- you can see substantial K27me3 signal in 
#### Chromatin state 9 which are usually gene poor ) - a gene based analyses will miss all of these non-genic polycomb domains

#### background distribution of genome (used as reference in the barplot in Fig S9D) 

background_dist <- GeneBodyFeatures %>% group_by(K36Me_bound_new_v2_Clusters_merged) %>% select(width,K36Me_bound_new_v2_Clusters_merged) %>% summarize(n=sum(width,na.rm = T))

## estimate intergenic regions - this way you get rid of overlapping genes and stuff 
intergenic <- gaps(IRanges::reduce(my_genes,ignore.strand=T))
intergenic <- intergenic[strand(intergenic)=="*"]

background_dist[7,2] <- sum(width(intergenic))
background_dist[7,1] <- "Intergenic"

### force SC_4 gene widths to SC_0
background_dist[1,2] <- background_dist[1,2] + background_dist[6,2]

background_dist <- background_dist %>% filter(K36Me_bound_new_v2_Clusters_merged != 4) %>% transmute(Region_Type="Genome",Annotation=ifelse(!grepl("Intergenic",K36Me_bound_new_v2_Clusters_merged),paste0("SC_",K36Me_bound_new_v2_Clusters_merged),"Intergenic"),counts=n/sum(n))

### VISUALIZE CLASSIFICATION OF BINS 
plot_df <- cbind(K27_z_scores_v4,gene_SC_K27) %>% pivot_longer(cols=10:16,names_to = "Region_Type",values_to = "Change") %>% filter(Change==1) %>% pivot_longer(11:16,names_to = "Annotation",values_to = "Index") %>% filter(Index!=0) %>% group_by(Region_Type,Annotation) %>% summarize(counts=sum(Index)) %>% 
  group_by(Region_Type) %>% mutate(counts=counts/sum(counts)) 

plot_df <- rbind(plot_df,background_dist)


plot_df$Region_Type <- factor(plot_df$Region_Type, levels = c("Genome","unchanged","Ash1_decrease","Ash1_increase","Set2_decrease","Set2_increase","NSD_decrease","NSD_increase"))

pdf(file="RegionType_GenesIntergenes.pdf", width = 10, height=5)

plot_df %>% ggplot(aes(x=Region_Type,y=counts, fill=Annotation)) + geom_bar(stat="identity") + theme_bw()

dev.off()

#### ################ GENES BASED ANALYSES - Not part of supplementary figures ! ##################
 
### Background proprotion of # of genes

background_dist_genes <- GeneBodyFeatures %>% select(K36Me_bound_new_v2_Clusters_merged) %>% table() %>% as.data.frame() %>% setNames(c("K36Me_bound_new_v2_Clusters_merged","n"))

## force SC4 genes to SC0

background_dist_genes[1,2] <- background_dist_genes[1,2] + background_dist_genes[6,2]

background_dist_genes <- background_dist_genes %>% filter(K36Me_bound_new_v2_Clusters_merged != 4) %>% transmute(Region_Type="Genome",Annotation=paste0("SC_",K36Me_bound_new_v2_Clusters_merged),counts=n/sum(n))


### reduce lists 
summarize_genenames <- function(v) {
  Reduce(f=paste0, x = v)
}

plot_df <- cbind(K27_z_scores_v4,genenames_SC_K27) %>% pivot_longer(cols=10:16,names_to = "Region_Type",values_to = "Change") %>% filter(Change==1) %>% pivot_longer(11:15,names_to = "Annotation",values_to = "Genes") %>% filter(Genes %nin% c("0","0,0")) %>% mutate(Genes=paste0(Genes,",")) %>% group_by(Region_Type,Annotation) %>% summarize(names=summarize_genenames(Genes))

#### you have a long string containing gene names that overlap K27me changing bins. Extract unique gene names from this

plot_df$counts <- unlist(lapply(plot_df$names, function(x){
  length(sapply(unique(unlist(str_split(x,pattern = ","))), function(y) {y[y!="0"]})) - 2    ###### filter 0s, first and last characters 
}))

plot_df_v2 <- plot_df %>% select(-names) %>% group_by(Region_Type) %>% mutate(counts=counts/sum(counts)) %>% rbind(background_dist_genes)

plot_df_v2$Region_Type <- factor(plot_df_v2$Region_Type, levels = c("Genome","unchanged","Ash1_decrease","Ash1_increase","Set2_decrease","Set2_increase","NSD_decrease","NSD_increase"))

pdf(file="GenesType_Genes.pdf", width = 10, height=5)
plot_df_v2 %>% ggplot(aes(x=Region_Type,y=counts, fill=Annotation)) + geom_bar(stat="identity") + theme_bw()
dev.off()

################ Genes based analysis also gives the same conclusions as bins based analysis. ####################  

#### SAVE THIS AS VERSION 1 OF UP DOMAINS -- .ie. genes are classified as gaining K27me3 domains if they contain a bin that gains K27me3
##### These numbers correspond to the v1 shown in Fig S9E table
GeneBodyFeatures <- GeneBodyFeatures %>% mutate(K27_up_Ash1_v1 = ifelse(gene_id %in% unlist(unique(str_split(paste(plot_df[plot_df$Region_Type=="Ash1_increase",]$names,collapse=","),pattern = ","))),1,0),
                                                K27_up_NSD_v1 = ifelse(gene_id %in% unlist(unique(str_split(paste(plot_df[plot_df$Region_Type=="NSD_increase",]$names,collapse=","),pattern = ","))),1,0),
                                                K27_up_Set2_v1 = ifelse(gene_id %in% unlist(unique(str_split(paste(plot_df[plot_df$Region_Type=="Set2_increase",]$names,collapse=","),pattern = ","))),1,0))


#### Compare this to an alternate gene-centric approach, similar to whats used by DePierre et al 2023. In brief, calculate average K27me3 over all 
### genes, calculate genic z-score and pick genes with highest z-score. This approach will find genes with uniform increase in k27me3 over 
### entire gene body (more robust) but may miss genes which show increases only in small region within the gene (maybe just exons or introns etc)


#### load average K27me3 table (similar to whats shown in GENERATEGENEBODYFEATURES.R) - Each row represents a gene, each column is average genebody signal
#### for different RNAi conditions 

ave.table <- readRDS("ave.table_K27me3.rds")

## ave.table contains all Drosophila genes. Filter which overlaps K27me signal
ave.table_filt <- ave.table[rownames(ave.table) %in% my_genes[overlapsAny(my_genes,K27_regions_bed,minoverlap = 500)]$gene_id,]

### classify based on arbitrary thresholds ! 

ave.table_filt <- as.data.frame(ave.table_filt) %>% rownames_to_column("gene_id") %>% transmute(gene_id=gene_id,
                                                                                                K27_GST=KD30.GST_K27me3,
                                                                                                K27_Ash1_z=(KD30.Ash1_K27me3-KD30.GST_K27me3)/sqrt(KD30.Ash1_K27me3+KD30.GST_K27me3),
                                                                                                K27_NSD_z=(KD30.NSD_K27me3-KD30.GST_K27me3)/sqrt(KD30.NSD_K27me3+KD30.GST_K27me3),
                                                                                                K27_Set2_z=(KD30.Set2_K27me3-KD30.GST_K27me3)/sqrt(KD30.Set2_K27me3+KD30.GST_K27me3),
                                                                                                K27_Ash1_up=ifelse((K27_GST<1.25)&(K27_Ash1_z>=0.25),1,0),
                                                                                                K27_NSD_up=ifelse((K27_GST<1.25)&(K27_NSD_z>=0.25),1,0),
                                                                                                K27_Set2_up=ifelse((K27_GST<1.250)&(K27_Set2_z>=0.25),1,0))

### SAVE THESE AS UP DOMAINS V2 - These numbers are rperesented in Fig S9E table as v2
GeneBodyFeatures <- GeneBodyFeatures %>% mutate(K27_up_Ash1_v2 = ifelse(gene_id %in% ave.table_filt[ave.table_filt$K27_Ash1_up==1,]$gene_id,1,0),
                                                K27_up_NSD_v2 = ifelse(gene_id %in% ave.table_filt[ave.table_filt$K27_NSD_up==1,]$gene_id,1,0),
                                                K27_up_Set2_v2 = ifelse(gene_id %in% ave.table_filt[ave.table_filt$K27_Set2_up==1,]$gene_id,1,0))

### similar number of genes -- fairly similar annotations in v1 and v2, but differences do exist. More robust classification will be possible 
### only with more replicates !



