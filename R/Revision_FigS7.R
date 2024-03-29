#### Exploratory analyses for Figure S7 ### 

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


install.packages("ggbeeswarm")
library(ggbeeswarm)
###################### Load necessary annotations


my_chromosomes <- c("chr2L","chr2R","chr3L","chr3R","chrX","chrY","chr4")


GeneBodyFeatures <- readRDS("../../GeneBodyFeatures_v2/GeneBodyFeatures_expanded_v20.rds")  ## load most recent GeneBodyFeatures file

my_genes <- genes(TxDb.Dmelanogaster.UCSC.dm6.ensGene)

my_genes <- keepSeqlevels(my_genes,my_chromosomes,pruning.mode = "coarse")

my_lengths <- seqlengths(keepSeqlevels(TxDb.Dmelanogaster.UCSC.dm6.ensGene, my_chromosomes))


################### Does control average genic K36me1/2/3 signal correlate with S2 cell transcription ? - Fig S7A


pdf("Correlation_TPMvsK36me_K36meboundgenes.pdf",width=5,height=5)
PearsonR_me3 <- round(cor(GeneBodyFeatures[GeneBodyFeatures$K36Me_bound_new_v2==1,]$Average_Exp,GeneBodyFeatures[GeneBodyFeatures$K36Me_bound_new_v2==1,]$GFP_K36Me3,use = "complete.obs"),digits = 2) 
PearsonR_me2 <- round(cor(GeneBodyFeatures[GeneBodyFeatures$K36Me_bound_new_v2==1,]$Average_Exp,GeneBodyFeatures[GeneBodyFeatures$K36Me_bound_new_v2==1,]$GFP_K36Me2,use = "complete.obs"),digits = 2)
PearsonR_me1 <- round(cor(GeneBodyFeatures[GeneBodyFeatures$K36Me_bound_new_v2==1,]$Average_Exp,GeneBodyFeatures[GeneBodyFeatures$K36Me_bound_new_v2==1,]$GFP_K36Me1,use = "complete.obs"),digits = 2)

GeneBodyFeatures %>% filter(K36Me_bound_new_v2==1) %>% select(Average_Exp,GFP_K36Me3) %>% mutate(Average_Exp=log2(Average_Exp+1))%>% ggplot(aes(x=GFP_K36Me3,y=Average_Exp)) + ggpointdensity::geom_pointdensity(size=0.75,alpha=0.8) + viridis::scale_color_viridis() + theme_bw() + annotate(label=paste0("PearsonR:",PearsonR_me3),x=2,y=15,"text")
GeneBodyFeatures %>% filter(K36Me_bound_new_v2==1) %>% select(Average_Exp,GFP_K36Me2) %>% mutate(Average_Exp=log2(Average_Exp+1))%>% ggplot(aes(x=GFP_K36Me2,y=Average_Exp)) + ggpointdensity::geom_pointdensity(size=0.75,alpha=0.8) + viridis::scale_color_viridis() + theme_bw() + annotate(label=paste0("PearsonR:",PearsonR_me2),x=2, y=15,"text")
GeneBodyFeatures %>% filter(K36Me_bound_new_v2==1) %>% select(Average_Exp,GFP_K36Me1) %>% mutate(Average_Exp=log2(Average_Exp+1))%>% ggplot(aes(x=GFP_K36Me1,y=Average_Exp)) + ggpointdensity::geom_pointdensity(size=0.75,alpha=0.8) + viridis::scale_color_viridis() + theme_bw() + annotate(label=paste0("PearsonR:",PearsonR_me1),x=2,y=15,"text")

dev.off()



############## Exploratory analyses of Huang et al and DePierre et al RNAseq datasets. Define filepaths 


RNAseq_filepaths <- list.files("RNAseq/",full.names = T)

######## ONLY 1 REPLICATE IS AVAILABLE ON GEO FOR BOTH HUANG AND DEPIERRE DATASETS !

#### THIS SEVERELY LIMITS OUR ABILITY TO ANALYSE THE DATA IN A STATISTICAL FRAMEWORK USING TOOLS LIKE DESEQ2 (FOR DIFFERENTIAL EXPRESSION ANALYSES)

#### How to compare data? 
##### TPM is good for intra-sample comparisons but not INTERSAMPLE comparisons (see https://www.biostars.org/p/9555918/#9556305 ; https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html)

### DESEQ2's median of ratios (or edgeR's TMM) is better as it accounts for sequencing depth AND compositonal biases 
### However, both of these tools do not allow analyses without replicates (DESEQ2 fails at setupDDS step!)
### Either perform manual median of ratios normlaization as illustrated by toy example here in above hbc link, OR,
### ALTERNATIVELY, you can use Variance Stabilizing Transformation from DESEQ2 dds object with a design matrix of ~1 as suggested my Mike Love (DESEQ2 founder)
### (https://support.bioconductor.org/p/9142214/) -- This transfomration calculates size factors to make the variance indepednent of the mean (homoskedastic) and does log transformation
### Subtraction of appropriate columns results in LFC !

### READ in raw ReadsPerGene objects -- all of them appear to be unstranded 
library(HelpersforDESeq2)
library(DESeq2)

read_RPG <- function(x) {
  base_name <- basename(x) %>% gsub(".ReadsPerGene.out.tab","", x=.)
  read.table(file=x, col.names = c("ID",base_name,"NULL","NULL"),colClasses = c("factor","integer","NULL","NULL"))
}

raw_RPG <- RNAseq_filepaths %>% map_dfc(~ read_RPG(.)) %>% select(-paste("ID...",seq(3,9,by=2), sep="")) 
raw_RPG <- raw_RPG[-(1:4),]

colnames(raw_RPG)[1] <- "ID"

sample_table <- read.delim("RNAseq/SampleTable.txt",header = T,stringsAsFactors = F)


rownames(raw_RPG) <- raw_RPG$ID
raw_RPG <- raw_RPG %>% select(-ID)
raw_RPG$filt <- rowSums(raw_RPG>1)
raw_RPG_filt <- raw_RPG %>% filter(filt>=4) %>% select(-filt)  ### FILTER SAMPLES THAT HAVE ATLEAST 2 READS IN ATLEAST 4 SAMPLES 


#### Setup DDS object (design ~1 as recommended in above link by DESEq2 creator)

dds <- DESeqDataSetFromMatrix(countData = raw_RPG_filt,colData = sample_table, design = ~1)
dds_norm <- assay(vst(dds))   ## apply VST transformation, NOTE: Does blind=F change output? 

### dds_norm <- vst(as.matrix(raw_RPG_filt))  -- THIS ALSO GIVES SAME RESULTS 

### calulcate log2 Fold Change 

lfc_df <- as.data.frame(dds_norm) %>% rownames_to_column(var="gene_id") %>% transmute(gene_id=gene_id, 
                                                                                      lfcRNA_Ash1=GSM2443797-GSM2443796,
                                                                                      lfcRNA_NSD=GSM4411852-GSM4411851,
                                                                                      lfcRNA_Set2=GSM4411853-GSM4411851)

GeneBodyFeatures <- GeneBodyFeatures %>% left_join(lfc_df, by="gene_id")  #save


### Plot expression changes for each supercluster as boxplot for Fig S7B. Force SC4 genes to SC0

### plot 
pdf("ExpressionLFC_clusters.pdf",width=7,height=5)
GeneBodyFeatures %>% select(K36Me_bound_new_v2_Clusters_merged,lfcRNA_Ash1,lfcRNA_NSD,lfcRNA_Set2) %>% 
  mutate(K36Me_bound_new_v2_Clusters_merged = ifelse(K36Me_bound_new_v2_Clusters_merged %in% c("4"),"0",K36Me_bound_new_v2_Clusters_merged)) %>%
  pivot_longer(cols=2:4,values_to = "Expression_Change",names_to = "RNAi") %>% 
  ggplot(aes(x=RNAi,y=Expression_Change,fill=RNAi)) + geom_boxplot() + facet_wrap(vars(K36Me_bound_new_v2_Clusters_merged)) + 
  scale_x_discrete(labels=NULL) + ylim(c(-1.5,1.5))+ theme_bw()
dev.off()

#### Plot expression changes for selected superclusters with cluster information overlayed on top as beeswarm plots - Fig S7C

pdf("ExpressionLFC_subclusters.pdf",height=5)
for (i in 1:3){
  print(GeneBodyFeatures %>% select(K36Me_bound_new_v2_Clusters_merged,lfcRNA_Ash1,lfcRNA_NSD,lfcRNA_Set2,K36Me_bound_new_v2_Clusters) %>% 
          filter(K36Me_bound_new_v2_Clusters_merged %in% c("1","2_Eu","3")) %>%
          pivot_longer(cols=2:4,values_to = "Expression_Change",names_to = "RNAi") %>% 
          ggplot(aes(x=RNAi,y=Expression_Change,col=factor(K36Me_bound_new_v2_Clusters))) + geom_beeswarm(corral="gutter",alpha=0.9,size=0.75,cex=0.5) + scale_color_manual(values=c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black", "gold1", "skyblue2","brown","gray70","darkturquoise")) + theme_bw() +theme(legend.title = element_blank()) + coord_cartesian(ylim=c(-1.5,1.5)) + facet_wrap_paginate(~K36Me_bound_new_v2_Clusters_merged,nrow = 1,ncol = 1, page = i))
}
dev.off()

####### Can part of the expression changes be explained by K27me3 gain seen in Fig S9 -- Experiment selection of genes either from 
#### method v1 or v2 (described in detail in script Revision_FigS9)

plot_df_v2 <- GeneBodyFeatures %>% select(gene_name,gene_id,K36Me_bound_new_v2_Clusters_merged,K27_up_Ash1_v2,lfcRNA_Ash1,K27_up_NSD_v2,lfcRNA_NSD,K27_up_Set2_v2,lfcRNA_Set2) %>% 
  mutate(K36Me_bound_new_v2_Clusters_merged = ifelse(K36Me_bound_new_v2_Clusters_merged %in% c("4"),"0",K36Me_bound_new_v2_Clusters_merged))

### helper function for appropriate plotting of gene number on top of boxplots shown in Fig S9F

give.n <- function(x){
  return(c(y = median(x) + 1, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}


pdf("ExpressionChange_K27me3gainCorr_UpdomainsFromGenicZ.pdf",width=5,height=5)
for (i in unique(plot_df_v2$K36Me_bound_new_v2_Clusters_merged)){
  for (j in my_conditions){
    plot_df_ij <- plot_df_v2 %>% filter(K36Me_bound_new_v2_Clusters_merged == i) %>% select(matches(j)) %>% mutate(RNAi=j) 
    plot_df_ij[,paste0("K27_up_",j,"_v2")] <- as.factor(plot_df_ij[,paste0("K27_up_",j,"_v2")] %>% pull(paste0("K27_up_",j,"_v2")))
    # p <- plot_df_ij %>% ggplot(aes_string(x="RNAi",y=paste0("lfcRNA_",j),col=paste0("K27_up_",j), alpha=paste0("K27_up_",j))) + geom_beeswarm() + scale_color_manual(values=c("#606060","black")) + ylim(c(-2,2)) + theme_bw()  + ggtitle(paste0("Cluster ",i,":",j))
    p <- plot_df_ij %>% ggplot(aes_string(x="RNAi",y=paste0("lfcRNA_",j),fill=paste0("K27_up_",j,"_v2"))) + geom_boxplot(alpha=0.8) +scale_fill_manual(values=c("gray","#CC0000"))+stat_summary(fun.data = give.n, geom = "text", fun.y = mean,position = position_dodge(width = 0.75)) + theme_bw() + ggtitle(paste0("Cluster ",i,":",j)) + coord_cartesian(ylim=c(-1.5,1.5))
    
    print(p)
  }
}
dev.off()


## Above observations suggest that changes in K36me1/2/3 (defined by cluster grouping) correlate poorly to gene expression changes
## What about Me2:Me3 ratios?

###################### Do top 10% upregulated, downregulate and unchanged lfc genes show difference in me2/3 ratios ? - Fig S7D 

#### Perform this analyses for Supercluster 1 genes

quantile(GeneBodyFeatures[GeneBodyFeatures$K36Me_bound_new_v2_Clusters_merged=="1",]$lfcRNA_Set2,na.rm=T,type=1, probs=seq(0,1,0.05))
#quantile(GeneBodyFeatures[GeneBodyFeatures$K36Me_bound_new_v2_Clusters_merged=="2_Eu",]$lfcRNA_Set2,na.rm=T,type=1, probs=seq(0,1,0.05))


pdf("Me23RatioChange_GSTSet2_Top10perUpDownGenes_SC12Eu.pdf",width=5,height=4)

### cutoffs were decided based on quantiles
temp <- GeneBodyFeatures %>% filter(K36Me_bound_new_v2_Clusters_merged=="1") %>% mutate(lfc_Set2_cutoffs=case_when((lfcRNA_Set2< -0.21) ~ "-1", #downregulated
                                                                                                                   ((lfcRNA_Set2< 0.03) & (lfcRNA_Set2> -0.01)) ~ "0", #unchanged
                                                                                                                   (lfcRNA_Set2> 0.24) ~ "1", ## upregulated
                                                                                                                   .default = "R")) %>% filter(lfc_Set2_cutoffs %in% c("-1","0","1")) %>% select(GST_Me23_ratio,Set2_Me23_ratio,gene_id,lfc_Set2_cutoffs) 

### add median FAKE genes for each group - these represent the red lines in the plot 
temp <- rbind(temp, c(median(temp[temp$lfc_Set2_cutoffs=="0",]$GST_Me23_ratio, na.rm=T),median(temp[temp$lfc_Set2_cutoffs=="0",]$Set2_Me23_ratio, na.rm=T),"Fake_unchanged","0"))
temp <- rbind(temp, c(median(temp[temp$lfc_Set2_cutoffs=="1",]$GST_Me23_ratio, na.rm=T),median(temp[temp$lfc_Set2_cutoffs=="1",]$Set2_Me23_ratio, na.rm=T),"Fake_up","1"))
temp <- rbind(temp, c(median(temp[temp$lfc_Set2_cutoffs=="-1",]$GST_Me23_ratio, na.rm=T),median(temp[temp$lfc_Set2_cutoffs=="-1",]$Set2_Me23_ratio, na.rm=T),"Fake_down","-1"))


temp$med_line <- factor(ifelse(grepl("Fake",temp$gene_id),1,0))


print(temp %>% pivot_longer(cols=1:2,names_to = "RNAi",values_to = "Me23_ratio") %>% ggplot(aes(x=factor(RNAi),y=as.numeric(Me23_ratio),group=gene_id,col=med_line,alpha=med_line,size=med_line)) + geom_line() + ylim(c(-3,3)) + scale_color_manual(values=c("#404040","red")) + scale_size_manual(values=c(0.1,1)) + scale_alpha_manual(values=c(0.25,1)) + facet_wrap(~lfc_Set2_cutoffs) + theme_bw() + ggtitle("SC_1"))

dev.off()

### Modify code above if you want to check for other superclusters ! - Not much differnece in conclusions !