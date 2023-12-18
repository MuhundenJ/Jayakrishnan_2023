### Script for analyzing flyAtlas data - Fig 5C
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

########## load raw flyatlas data -- txt file of expression arrays 

### Last 150 rows are messed up and do not have probes associated to them nor Expression values. Another additional row also has special characters instead of numeric values
### Cleanup after manual inspection and keep only columns containing relevant information 

### Note that geneID is absent, rather only Oligos (often multiple spanning a single gene) expression data is present in this table
FlyAtlas_raw <- read.table("FlyAtlas/FlyAtlas.txt",header = T, fill=T)
FlyAtlas_raw <- FlyAtlas_raw[-c(12659,18771:nrow(FlyAtlas_raw)),] %>% select(matches("Oligo|Mean|mean|Present|ratio"))

### load cleaned up excel data matching Oligo IDs to FBgn IDs 
###Here there are two columns that have flybase IDs but sometimes missing in one of them. Use both columns to make consensus set of FlybaseID

FlyAtlas_Anno <- read_excel("FlyAtlas/Drosophila_2.na32.annot.xlsx",col_names=c("probe_ID","gene_id"),skip=1,sheet = 2) %>% transmute(Oligo=gsub(",.*","",.$probe_ID),
                                                                                                                                      gene_id_1=gsub(" /SEG","",str_extract(.$probe_ID,"FBgn.*/SEG")),
                                                                                                                                      gene_id_2=gsub(" //.*","",str_extract(.$gene_id,"FBgn.*//")),
                                                                                                                                      gene_id=ifelse(is.na(gene_id_2),gene_id_1,gene_id_2)) %>% select(Oligo,gene_id)
                                                                                                                                
FlyAtlas_Anno <- FlyAtlas_Anno %>% filter(!is.na(gene_id)) %>% left_join(FlyAtlas_raw,by="Oligo") ## 10% of Flybase IDs have non unique probe ID - can deal later

##########

FlyAtlas_Anno <- FlyAtlas_Anno[,1:(ncol(FlyAtlas_Anno)-5)]    ### skip S2 and whole fly data -- keep only 25 tissues 

### Calculation expression variance (weighted by mean) across 25 tissues 

FlyAtlas_Anno <- FlyAtlas_Anno %>% mutate(Exp_sd=select(.,matches("Mean|mean"))%>%rowSds(na.rm=T),Exp_avg=select(.,matches("Mean|mean"))%>%rowMeans(na.rm=T), z=Exp_sd/Exp_avg) ## variance in expression

##### visualize distribution of expression Z and select appropriate cutoff 
pdf("FlyAtlas/Z_score.pdf",width=7,height=5)
FlyAtlas_Anno %>% ggplot(aes(x=z)) + geom_histogram(aes(y=..density..),alpha=0.3,color="gray50") + geom_density() + geom_vline(xintercept=1.2,linetype="dotted",color="red",size=1.0)
dev.off()

##### cutoff of z=1.2 chosen. However, genes with very low expression in ALL 25 tissues will have low Z 
##### However, it makes sense that these genes are probably uniquely expressed in some cell type not sampled by FlyAtlas 25 tissues .ie. DV and not HK
##### so manually classify these genes (defined by N_TissueExpressed<=3 ; this value was selected by trial and error ) as Developmental

FlyAtlas_Anno <- FlyAtlas_Anno %>% mutate(HK_DV_FlyAtlas=ifelse(z>=1.2|N_TissueExpressed<=3,"Developmental","Housekeeping"))

saveRDS(FlyAtlas_Anno,file="FlyAtlas/FlyAtlas_table.rds")

### Note that these are Probe_IDs and not genes -> About 1.1k genes have more than one probe -> Need to make a reduced list - HOW TO DEAL WITH MIXED 
### GENES THAT HAVE PROBES IN BOTH HOUSEKEEPING AND DEVELOPMENTAL ? 

my_filter <- FlyAtlas_Anno %>% group_by(gene_id) %>% summarize(n=n()) %>% filter(n>1) %>% arrange(desc(n))

final_list <- FlyAtlas_Anno %>% select(gene_id,HK_DV_FlyAtlas) %>% filter(gene_id %nin% my_filter$gene_id)

## this snippet loops through non unique probes and classifies them as HK or DV based on majority classification. If tie then HK

i <- 1

for (i in seq_along(my_filter$gene_id)){

  my_gene <- my_filter$gene_id[i]
  
  my_counts <- FlyAtlas_Anno[FlyAtlas_Anno$gene_id == my_gene ,c("gene_id","HK_DV_FlyAtlas")] %>% select(HK_DV_FlyAtlas) %>% table() 

  if (length(my_counts)>1){
  
    final_list<- rbind(final_list,data.frame(gene_id=my_gene,HK_DV_FlyAtlas=ifelse(my_counts["Housekeeping"]>=my_counts["Developmental"],"Housekeeping","Developmental")))
  
  } else {
    
    final_list<- rbind(final_list,data.frame(gene_id=my_gene,HK_DV_FlyAtlas=names(my_counts)))
  }
    
}


GeneBodyFeatures <- GeneBodyFeatures %>% left_join(final_list,by="gene_id")

##### around 8.2k developmental genes and 6.7k housekeeping genes -- plot distirbution in relation to superclusters
 
pdf("FlyAtlas/HK_DV_Plots.pdf",width=7,height=5)

print(GeneBodyFeatures %>% filter(K36Me_bound_new_v2_Clusters_merged %in% c("0","1","2_Eu","2_Het","3") & !is.na(HK_DV_FlyAtlas)) %>% select(K36Me_bound_new_v2_Clusters_merged,HK_DV_FlyAtlas) %>%
  group_by(K36Me_bound_new_v2_Clusters_merged,HK_DV_FlyAtlas) %>% summarize(n=n()) %>% ggplot(aes(x=K36Me_bound_new_v2_Clusters_merged,y=n,fill=HK_DV_FlyAtlas)) + geom_col(position="fill",color="black")+ scale_fill_brewer(palette = "Pastel1") + coord_flip())

print(GeneBodyFeatures %>% filter(K36Me_bound_new_v2_Clusters_merged %in% c("1","2_Eu","2_Het","3") & !is.na(HK_DV_FlyAtlas)) %>% select(K36Me_bound_new_v2_Clusters_merged,HK_DV_FlyAtlas) %>%
        group_by(K36Me_bound_new_v2_Clusters_merged,HK_DV_FlyAtlas) %>% summarize(n=n()) %>% ggplot(aes(x=K36Me_bound_new_v2_Clusters_merged,y=n,fill=HK_DV_FlyAtlas)) + geom_col(position="fill",color="black")+ scale_fill_brewer(palette = "Pastel1") + coord_flip())

dev.off()

FlyAtlas_Anno <- FlyAtlas_raw %>% select(Oligo,FlyMean) %>% right_join(FlyAtlas_Anno,by="Oligo")
