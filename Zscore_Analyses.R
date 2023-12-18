### Script for following z-score related analyses 
#### Fig 3 - z-score paired scatter
#### Fig 4C - z-score for Huang K36me2 in Ash1 RNAi 
#### Fig 6 and S7 - z-score densities for selected clusters for K36me1/2/3 and readers
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

########################### Z -score plots to represent global changes
### Prior to this, convert all merged bigwig coverages to bedgraph and resize to 5kb windows (containing average signal within those windows)

my_file_directory <- list.files("Z_scores_temp/",pattern="*_resize.bedgraph",full.names = T)
my_file_names <-list.files("Z_scores_temp/",pattern="*_resize.bedgraph")

#### load a bed file containing merged peak sets of K36me1/2/3 across RNAi conditions. This will be used to filter out bins that have no K36me1/2/3
#### .ie. contains atleast one of K36me1 or 2 or 3
K36me_regions_bed <- read.delim("./H3K36me_regions.bed",header=F) %>% makeGRangesFromDataFrame(start.field = "V2",end.field = "V3",seqnames="V1",keep.extra.columns = T)


my_chromosomes <- c("chr2L","chr2R","chr3L","chr3R","chrX","chrY","chr4")

my_conditions <- c("Ash1","Set2","NSD","Comb")

my_chips <- c("K36Me1","K36Me2","K36Me3")


i <- 1

for (i in seq_along(my_chips)){

  my_files_chip <- my_file_directory[grepl(my_chips[i],my_file_directory)]
  
  #### tailor filenames 
  my_filenames_chip <- gsub("coverage.Merged.","",gsub(paste0("_",my_chips[i],".*"),"",my_file_names[grepl(my_chips[i],my_file_names)]))
  
  #### read all files using vroom, convert to df and filter useless bins 
  
  raw_scores <- lapply(my_files_chip, function(x) {read.table(file = x, header = F, sep ="\t")})
  
  names(raw_scores) <- my_filenames_chip
  
  raw_scores <- map_df(raw_scores,data.frame,.id="name") %>% pivot_wider(names_from="name",values_from="V4") %>% makeGRangesFromDataFrame(seqnames.field ="V1",start.field = "V2",end.field = "V3",keep.extra.columns = T)  ### this already does a outerjoin based on coordinates !

  raw_scores_filt <- as.data.frame(raw_scores[overlapsAny(raw_scores,K36me_regions_bed)])
  
  
  ### Calculate z-scores for appropriate ChIP-RNAi combination as described by Chaouch et al., 2022 using genebody averaged coverage instead of RPKPM 
  
  ### Note that Huang dataset cant be directly compared with our dataset due to difference in protocols 
  ### Calculate Huang Ash1 z-score using the respective matched control -
  
  if (my_chips[i]=="K36Me2"){
    z_scores_df <- raw_scores_filt %>% transmute(chr=seqnames,start=start,end=end,Set2=(Set2-GFP)/sqrt(Set2+GFP),
                                                                    Ash1=(Huang_Ash1-Huang_GST)/sqrt(Huang_Ash1+Huang_GST),
                                                                    NSD=(NSD-GFP)/sqrt(NSD+GFP),
                                                                    Comb=(Comb-GFP)/sqrt(Comb+GFP))
    
    }else if(my_chips[i]=="K36Me3"){
    z_scores_df <- raw_scores_filt %>% transmute(chr=seqnames,start=start,end=end,Set2=(Set2-GFP)/sqrt(Set2+GFP),
                                                                    NSD=(NSD-GFP)/sqrt(NSD+GFP),
                                                                    Comb=(Comb-GFP)/sqrt(Comb+GFP))
    

    }else{
    z_scores_df <- raw_scores_filt %>% transmute(chr=seqnames,start=start,end=end,Set2=(Set2-GFP)/sqrt(Set2+GFP),
                                                                    Ash1=(Ash1-GFP)/sqrt(Ash1+GFP),
                                                                    NSD=(NSD-GFP)/sqrt(NSD+GFP),
                                                                    Comb=(Comb-GFP)/sqrt(Comb+GFP))
    
  
    
  }
  
  assign(paste0(my_chips[i],"_z_scores"),z_scores_df)
  
  }

######### Paired scatter plots 
library(GGally)
library(ggpointdensity)
library(viridis)

i <- 3


for(i in seq_along(my_chips)){
  
  pdf(paste0("pointDens_Z_score.",my_chips[i],".pdf"), width = 7.5, height = 7.5)
  
  ### exclude irreleveant columns that arent necessary for the scatter
  
  if(my_chips[i]=="K36Me2"){
  
  ave.table_pointDens <- get(paste0(my_chips[i],"_z_scores")) %>% select(-c("chr","start","end","Ash1")) %>% drop_na()
  
  } else {
  ave.table_pointDens <- get(paste0(my_chips[i],"_z_scores")) %>% select(-c("chr","start","end")) %>% drop_na()
    
  }
  
  p_ggpairs <- ggpairs(ave.table_pointDens,diag=list(continuous = "blankDiag"),lower=list(continuous = wrap("points",alpha = 0.1,size=0.2)))+ xlim(-1.5,1.5)+ ylim(-1.5,1.5)
  
  ### THIS STEP IS NECESSARY AS OTHERWISE GEOM_POINTDENSITY REPLACES AES FROM GGPAIRS !
  add_to_ggmatrix(p_ggpairs, geom_pointdensity(alpha=0.1,adjust=0.5),location="lower") + scale_color_viridis() +theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + geom_vline(xintercept=0,linetype="dotted") + geom_hline(yintercept=0,linetype="dotted") + geom_abline(slope=1,intercept=0,linetype="dotted")
  
  dev.off()
}


#### Make a zscore plot for GB from Huang dataset -- As shown in Fig 4C

my_file_directory_Huang <- list.files("Z_scores_temp/",pattern="^Huang.*_resize.bedgraph",full.names = T)
my_file_names_Huang <-list.files("Z_scores_temp/",pattern="^Huang.*_resize.bedgraph")


i <- 1

for (i in seq_along(my_file_directory_Huang)){

  my_filename <- gsub("_resize.bedgraph","",my_file_names_Huang[i])
  
  my_cov_bedgraph <- import.bedGraph(my_file_directory_Huang[i])
  seqlevelsStyle(my_cov_bedgraph) <- "UCSC"
  
  my_cov_bedgraph <- keepSeqlevels(my_cov_bedgraph, my_chromosomes, pruning.mode = "coarse")
  seqlengths(my_cov_bedgraph) <- my_lengths
  
  my_cov_bedgraph <- coverage(my_cov_bedgraph, weight = "score")
  
  my_genes <- all_genes[GeneBodyFeatures[GeneBodyFeatures$K36Me_bound_new_v2_Clusters_merged %in% c("1","2_Eu","2_Het","3"),]$gene_id]
  
  my_mat <- my_cov_bedgraph[my_genes]
  
  names(my_mat) <- my_genes$gene_id
  
  my_mat <- do.call(rbind,lapply(my_mat, function(x){ mean(x)}))
  
  colnames(my_mat) <- my_filename
  
  assign(paste0("GB_Zscores_",my_filename),my_mat)
  
}

GB_Zscores_Huang <- as.data.frame(cbind(GB_Zscores_Huang_GST_K36Me2,GB_Zscores_Huang_Ash1_K36Me2)) %>% transmute(z_score=(Huang_Ash1_K36Me2-Huang_GST_K36Me2)/sqrt(Huang_Ash1_K36Me2+Huang_GST_K36Me2)) %>% rownames_to_column(var="gene_id") %>% left_join(GeneBodyFeatures[,c("gene_id","K36Me_bound_new_v2_Clusters_merged")],by="gene_id")
GB_Zscores_Huang_v2 <- GB_Zscores_Huang

GB_Zscores_Huang_v2$K36Me_bound_new_v2_Clusters_merged <- fct_collapse(GB_Zscores_Huang_v2$K36Me_bound_new_v2_Clusters_merged,`2`=c("2_Eu","2_Het"))

pdf("Huang_Ash1_Zscores_clusters_v2.pdf",height=5,width=7.5)

GB_Zscores_Huang_v2 %>% ggplot(aes(x=z_score,col=K36Me_bound_new_v2_Clusters_merged))+geom_density() + geom_density(size=1) + scale_colour_manual(values=c("red","green","#1E88E5")) + scale_x_continuous(breaks=seq(-4,4,1),limits = c(-4,4)) + geom_vline(xintercept = 0,linetype="dotted", size=1) + ggtitle("Huang Ash1 RNAi H3K36me2 Z-scores") + theme_classic() + labs(col="Cluster")

dev.off()




######### Code for Z-score desntiy plots with genebodies as shown in Fig 6 and Supplemantary Fig7

##### load ave.table objects for readers and K36me1/2/3 generated for clustering heatmaps and other downstream analyses in GenerateIntermediateFiles.R
#### this object contains genebody averaged coverages for all chips in all conditions 

#### calculate z-scores. Note that ave.table consists of log2 transformed values which cant be used for z-score calulcation (negative values throw error)
#### convert back to ChIP/Inp simple values

### last step is filtering appropriate clusters of interest -- Note that cluster numbers shown here do not correspond to the one shown on Final Figures 
ave.table_readers_zscore_v2 <- 2^(as.data.frame(ave.table_clustering_readers)) %>% rownames_to_column(var="gene_id") %>% transmute(gene_id=gene_id,
                                                                                                                            Set2_JASPer=(Merged.Set2_JASPer-Merged.GFP_JASPer)/sqrt(Merged.Set2_JASPer+Merged.GFP_JASPer),
                                                                                                                            NSD_JASPer=(Merged.NSD_JASPer-Merged.GFP_JASPer)/sqrt(Merged.NSD_JASPer+Merged.GFP_JASPer),
                                                                                                                            Comb_JASPer=(Merged.Comb_JASPer-Merged.GFP_JASPer)/sqrt(Merged.Comb_JASPer+Merged.GFP_JASPer),
                                                                                                                            Set2_Msl3=(Merged.Set2_Msl3-Merged.GFP_Msl3)/sqrt(Merged.Set2_Msl3+Merged.GFP_Msl3),
                                                                                                                            NSD_Msl3=(Merged.NSD_Msl3-Merged.GFP_Msl3)/sqrt(Merged.NSD_Msl3+Merged.GFP_Msl3),
                                                                                                                            Comb_Msl3=(Merged.Comb_Msl3-Merged.GFP_Msl3)/sqrt(Merged.Comb_Msl3+Merged.GFP_Msl3)) %>%   
                                                                                                                              left_join(GeneBodyFeatures[,c("gene_id","K36Me_bound_new_v2_Clusters")],by="gene_id") %>% filter(K36Me_bound_new_v2_Clusters %in% c("3","8","10","5","12"))


#### plot z-score density plots for JASPer (whole genome, autosomes only) and Msl3 (ChrX only)
pdf("Readers_Zscores_clusters.pdf",height=5,width=7.5)

ave.table_readers_zscore_v2 %>% filter(gene_id %in% all_genes[seqnames(all_genes)!="chrX"]$gene_id)%>% select(matches("JASPer|K36me")) %>% pivot_longer(cols=contains("JASPer"),names_to = "RNAi_Condition", values_to = "Z_Score") %>% ggplot(aes(x=Z_Score,col=RNAi_Condition)) + geom_density(size=1) + scale_colour_manual(values=my_colours[-1]) + scale_x_continuous(breaks=seq(-1,1,0.5),limits = c(-1,1)) + geom_vline(xintercept = 0,linetype="dotted", size=1) +facet_wrap(~K36Me_bound_new_v2_Clusters) + ggtitle("HMT RNAi JASPer Z-scores (ChrA,n=8539)") + theme_classic() + labs(col="RNAi_Condition")
ave.table_readers_zscore_v2 %>% select(matches("JASPer|K36me")) %>% pivot_longer(cols=contains("JASPer"),names_to = "RNAi_Condition", values_to = "Z_Score") %>% ggplot(aes(x=Z_Score,col=RNAi_Condition)) + geom_density(size=1) + scale_colour_manual(values=my_colours[-1]) + scale_x_continuous(breaks=seq(-1,1,0.5),limits = c(-1,1)) + geom_vline(xintercept = 0,linetype="dotted", size=1) +facet_wrap(~K36Me_bound_new_v2_Clusters) + ggtitle("HMT RNAi JASPer Z-scores (n=10477)") + theme_classic() + labs(col="RNAi_Condition")
ave.table_readers_zscore_v2 %>% filter(gene_id %in% all_genes[seqnames(all_genes)=="chrX"]$gene_id) %>% select(matches("Msl3|K36me")) %>% pivot_longer(cols=contains("Msl3"),names_to = "RNAi_Condition", values_to = "Z_Score") %>% ggplot(aes(x=Z_Score,col=RNAi_Condition)) + geom_density(size=1) + scale_colour_manual(values=my_colours[-1]) + scale_x_continuous(breaks=seq(-1.5,1.5,0.5),limits = c(-1.5,1.5)) + geom_vline(xintercept = 0,linetype="dotted", size=1) +facet_wrap(~K36Me_bound_new_v2_Clusters) + ggtitle("HMT RNAi Msl3 Z-scores (ChrX only, n=1663)") + theme_classic() + labs(col="RNAi_Condition")

dev.off()


### make density z score plots for K36me1/2/3 as well


ave.table_K36me_Zscore_v2 <- 2^(as.data.frame(ave.table_clustering)) %>% rownames_to_column(var="gene_id") %>% transmute(gene_id=gene_id,
                                                                                                                                Set2_K36Me3=(Merged.Set2_K36Me3-Merged.GFP_K36Me3)/sqrt(Merged.Set2_K36Me3+Merged.GFP_K36Me3),
                                                                                                                                NSD_K36Me3=(Merged.NSD_K36Me3-Merged.GFP_K36Me3)/sqrt(Merged.NSD_K36Me3+Merged.GFP_K36Me3),
                                                                                                                                Comb_K36Me3=(Merged.Comb_K36Me3-Merged.GFP_K36Me3)/sqrt(Merged.Comb_K36Me3+Merged.GFP_K36Me3),
                                                                                                                                Set2_K36Me2=(Merged.Set2_K36Me2-Merged.GFP_K36Me2)/sqrt(Merged.Set2_K36Me2+Merged.GFP_K36Me2),
                                                                                                                                NSD_K36Me2=(Merged.NSD_K36Me2-Merged.GFP_K36Me2)/sqrt(Merged.NSD_K36Me2+Merged.GFP_K36Me2),
                                                                                                                                Comb_K36Me2=(Merged.Comb_K36Me2-Merged.GFP_K36Me2)/sqrt(Merged.Comb_K36Me2+Merged.GFP_K36Me2),
                                                                                                                                Set2_K36Me1=(Rep3.Set2_K36Me1-Rep3.GFP_K36Me1)/sqrt(Rep3.Set2_K36Me1+Rep3.GFP_K36Me1),
                                                                                                                                NSD_K36Me1=(Rep3.NSD_K36Me1-Rep3.GFP_K36Me1)/sqrt(Rep3.NSD_K36Me1+Rep3.GFP_K36Me1),
                                                                                                                                Comb_K36Me1=(Rep3.Comb_K36Me1-Rep3.GFP_K36Me1)/sqrt(Rep3.Comb_K36Me1+Rep3.GFP_K36Me1),
                                                                                                                                Ash1_K36Me1=(Rep3.Ash1_K36me1-Rep3.GFP_K36Me1)/sqrt(Rep3.Ash1_K36me1+Rep3.GFP_K36Me1)) %>%   
                                                                                                                                    left_join(GeneBodyFeatures[,c("gene_id","K36Me_bound_new_v2_Clusters")],by="gene_id") %>% filter(K36Me_bound_new_v2_Clusters %in% c("3","8","10","5","12"))



pdf("K36me_Zscores_clusters_v2.pdf",height=5,width=7.5)

ave.table_K36me_Zscore_v2 %>% select(matches("clusters|K36me3")) %>% pivot_longer(cols=contains("K36Me3"),names_to = "RNAi_Condition", values_to = "Z_Score") %>% ggplot(aes(x=Z_Score,col=RNAi_Condition)) + geom_density(size=1) + scale_colour_manual(values=my_colours[-1]) + scale_x_continuous(breaks=seq(-1,1,0.5),limits = c(-1,1)) + geom_vline(xintercept = 0,linetype="dotted", size=1) +facet_wrap(~K36Me_bound_new_v2_Clusters) + ggtitle("HMT RNAi K36me3 Z-scores (n=10477)") + theme_classic() + labs(col="RNAi_Condition")
ave.table_K36me_Zscore_v2 %>% select(matches("clusters|K36me2")) %>% pivot_longer(cols=contains("K36Me2"),names_to = "RNAi_Condition", values_to = "Z_Score") %>% ggplot(aes(x=Z_Score,col=RNAi_Condition)) + geom_density(size=1) + scale_colour_manual(values=my_colours[-1]) + scale_x_continuous(breaks=seq(-1,1,0.5),limits = c(-1,1)) + geom_vline(xintercept = 0,linetype="dotted", size=1) +facet_wrap(~K36Me_bound_new_v2_Clusters) + ggtitle("HMT RNAi K36me2 Z-scores (n=10477)") + theme_classic() + labs(col="RNAi_Condition")
ave.table_K36me_ZscoreP_v2 %>% select(matches("clusters|K36me1")) %>% pivot_longer(cols=contains("K36Me1"),names_to = "RNAi_Condition", values_to = "Z_Score") %>% ggplot(aes(x=Z_Score,col=RNAi_Condition)) + geom_density(size=1) + scale_colour_manual(values=my_colours) + scale_x_continuous(breaks=seq(-1,1,0.5),limits = c(-1,1)) + geom_vline(xintercept = 0,linetype="dotted", size=1) +facet_wrap(~K36Me_bound_new_v2_Clusters) + ggtitle("HMT RNAi K36me1 Z-scores (n=10477)") + theme_classic() + labs(col="RNAi_Condition")


ave.table_K36me_Zscore_v2 %>% filter(gene_id %in% all_genes[seqnames(all_genes)=="chrX"]$gene_id)%>% select(matches("clusters|K36me3"))%>% pivot_longer(cols=contains("K36Me3"),names_to = "RNAi_Condition", values_to = "Z_Score") %>% ggplot(aes(x=Z_Score,col=RNAi_Condition)) + geom_density(size=1) + scale_colour_manual(values=my_colours[-1]) + scale_x_continuous(breaks=seq(-1,1,0.5),limits = c(-1,1)) + geom_vline(xintercept = 0,linetype="dotted", size=1) +facet_wrap(~K36Me_bound_new_v2_Clusters) + ggtitle("HMT RNAi K36me3 ChrX Z-scores (n=1663)") + theme_classic() + labs(col="RNAi_Condition")
ave.table_K36me_Zscore_v2 %>% filter(gene_id %in% all_genes[seqnames(all_genes)=="chrX"]$gene_id)%>% select(matches("clusters|K36me2"))%>% pivot_longer(cols=contains("K36Me2"),names_to = "RNAi_Condition", values_to = "Z_Score") %>% ggplot(aes(x=Z_Score,col=RNAi_Condition)) + geom_density(size=1) + scale_colour_manual(values=my_colours[-1]) + scale_x_continuous(breaks=seq(-1,1,0.5),limits = c(-1,1)) + geom_vline(xintercept = 0,linetype="dotted", size=1) +facet_wrap(~K36Me_bound_new_v2_Clusters) + ggtitle("HMT RNAi K36me2 ChrX Z-scores (n=1633)") + theme_classic() + labs(col="RNAi_Condition")
ave.table_K36me_Zscore_v2 %>% filter(gene_id %in% all_genes[seqnames(all_genes)=="chrX"]$gene_id)%>% select(matches("clusters|K36me1"))%>% pivot_longer(cols=contains("K36Me1"),names_to = "RNAi_Condition", values_to = "Z_Score") %>% ggplot(aes(x=Z_Score,col=RNAi_Condition)) + geom_density(size=1) + scale_colour_manual(values=my_colours) + scale_x_continuous(breaks=seq(-1,1,0.5),limits = c(-1,1)) + geom_vline(xintercept = 0,linetype="dotted", size=1) +facet_wrap(~K36Me_bound_new_v2_Clusters) + ggtitle("HMT RNAi K36me1 ChrX Z-scores (n=1633)") + theme_classic() + labs(col="RNAi_Condition")


dev.off()

