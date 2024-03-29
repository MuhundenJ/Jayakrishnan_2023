###### Script to analyzse a) Intronic Proportion of genes used in Fig 5 and b) K36me1/2/3 signal across different RNAi conditions separately for Exons and introns as shown in FigS10
#### run on computational cluster!

rm(list = ls())
try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))



library(ComplexHeatmap)
library(dendextend)
library(dplyr)
library(RColorBrewer)
library(circlize)
library(purrr)
library(GenomicRanges)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(rtracklayer)

`%nin%` = Negate(`%in%`)


#### load genebodyfeatures file containing useful information like cluster memberships

GeneBodyFeatures <- readRDS("GeneBodyFeatures.rds")


#### load relevant annotation files 

### load detailed gene annotaions from geome.gtf -- BDGP6 v104 ENSEMBL
my_allgenes <- makeTxDbFromGFF("./genome.gtf",format="gtf")
my_allgenes <- genes(my_allgenes)
seqlevelsStyle(my_allgenes) <- "UCSC"
my_allgenes <- keepSeqlevels(my_allgenes,  my_chromosomes, pruning.mode="coarse")

my_allgenes <- my_allgenes[grepl("FBgn",my_allgenes$gene_id)] ### filter out genes without Fbgn ID -- these are typically localized transposons etc

my_chromosomes <- c("chr2L","chr2R","chr3L","chr3R","chrX","chrY","chr4")

my_lengths <- seqlengths(keepSeqlevels(TxDb.Dmelanogaster.UCSC.dm6.ensGene, my_chromosomes, pruning.mode = "coarse"))

df_allgenes <- data.frame(chr=seqnames(all_genes),
                          start=start(all_genes)-1,
                          end=end(all_genes),
                          strand = strand(all_genes),
                          width = width(all_genes),
                          row.names = all_genes$gene_id,
                          gene_id=all_genes$gene_id)



##### Identify introns and exons 

### get intron and exon info from my_allgenes object using intronicParts function

my_allgenes_introns <- intronicParts(my_allgenes)
df_allgenes <- data.frame(chr=seqnames(all_genes),
                          start=start(all_genes)-1,
                          end=end(all_genes),
                          strand = strand(all_genes),
                          width = width(all_genes),
                          row.names = all_genes$gene_id,
                          gene_id=all_genes$gene_id)seqlevelsStyle(my_allgenes_introns) <- "UCSC"
my_allgenes_introns <- as.data.frame(keepSeqlevels(my_allgenes_introns,my_chromosomes, pruning.mode = "coarse"))

my_allgenes_exons <- exonicParts(my_allgenes)
seqlevelsStyle(my_allgenes_exons) <- "UCSC"
my_allgenes_exons <- as.data.frame(keepSeqlevels(my_allgenes_exons,my_chromosomes, pruning.mode = "coarse"))


########### calculate intronic widths 

#### Note - About 3% of entries in the intron dataset are shred between multiple genes -- Treat them as shared for both
#### For each intron, make a column containing all the FBgn IDs (gene.id) associated as a comma seperated string
#### For calculating intronic widths for a given gene, one could grep the specific gene id through previously defined index column, pull out the 
#### corresponding introns and sum them ! 

intron_geneid <- data.frame(matrix(nrow=length(my_allgenes_introns$gene_id),ncol=1))
colnames(intron_geneid)<-"GeneIDs"

### generate list of all geneids for each intronic column

i <- 1

for (i in 1:nrow(intron_geneid)){
  intron_geneid[i,1] <- paste(my_allgenes_introns$gene_id[i],collapse=";") 
}


df_allgenes_introns <- data.frame(chr=seqnames(my_allgenes_introns),
                                  start=start(my_allgenes_introns),
                                  end=end(my_allgenes_introns),
                                  strand=strand(my_allgenes_introns),
                                  width=width(my_allgenes_introns),
                                  gene_id=intron_geneid)


#in this below loop, we select a geneid from all_genes, search for this ID through the gene_IDs columns of introns to find all introns belonging
## to this geneID. Then sum up the widths of those introns and add it the the df_allgenes

i <- 1
for (i in 1:nrow(df_allgenes)){
  df_allgenes$intron_widths[i]<-sum(df_allgenes_introns$width[grep(df_allgenes$gene_id[i],df_allgenes_introns$GeneIDs)])
  
  #### note that the below step extracts only list of possible introns - however they may contain alternatively spliced exons 
  ### ie. exonics parts and intronic parts for a given gene are not mutually exclusive. However these regions are very small (few hundred bp max) relative to gene length
  df_allgenes$intron_widths[i]<-df_allgenes_introns[grep(df_allgenes$gene_id[i],df_allgenes_introns$GeneIDs),]
  
}

df_allgenes <- df_allgenes %>% mutate(Intron_Prop = intron_widths/width)   #### USE THIS FOR VISULIAZATION !




#############################################################################################################################
#############################################################################################################################

### load K36me1/2/3 coverages across all RNAi -- obtained from GenerateIntermediateFiles.R

coverages_filepath <- list.files(path = "Coverages_Final/",pattern = "^coverage.", full.names = T)
coverages_filenames <- list.files(path = "Coverages_Final/",pattern = "^coverage.", full.names =F)

for (i in seq_along(coverages_filepath)){
  filename <- gsub(".rds","",coverages_filenames[i])
  file <- readRDS(coverages_filepath[i])
  assign(filename,file,envir=.GlobalEnv)
}


my_covs <- ls(pattern="^coverage\\.")

print(my_covs)

my_chips <- c("K36Me1","K36Me2","K36Me3")

################# for each chip, loop through RNAi conditions to calculate per gene separate exon and intron ChIP coverage 

########### Of note, we explored some edge cases to ensure our observations were robust. For instance, some alternatively spliced exonic regions appear in
########## both introns and exon defintions. However these are rare and these regions are really small so do not affect whole exon or whole intron averages

########### Further, eliminating nested genes (where an exon of Gene A may intersect with intron of Gene B) also didn't cause much difference, 
########## likely because they are distributed fairly evenly across superclusters

i<- 1

for (i in seq_along(my_chips)){
  my_covs_chip <- my_covs[grep(my_chips[i],my_covs)]
  my_conditions <- c("GFP","Ash1","Set2","NSD","Comb")
  
  int_ex_table <- matrix(as.numeric(0),nrow=length(GeneBodyFeatures$gene_id),ncol=1+length(my_conditions)*2)
  
  int_ex_table[,1] <- GeneBodyFeatures$gene_id
  
  #int_ex_table <- data.frame(gene_id=GeneBodyFeatures$gene_id)
  
  rownames(int_ex_table) <- GeneBodyFeatures$gene_id
  
  colnames(int_ex_table) <- c("gene_id",apply(expand.grid(my_conditions,c("Exon","Intron")),1,paste,collapse="_"))

  int_ex_table <- as.data.frame(int_ex_table)
  
  int_ex_table[,-1] <- sapply(int_ex_table[,-1],as.numeric)
  
  j <- 1

  for (j in seq_along(my_conditions)){

    my_covs_chip_cond <- my_covs_chip[grep(my_conditions[j],my_covs_chip)]

    if (!isEmpty(my_covs_chip_cond)){
        
      
      
       k <-1

     
      for (k in seq_along(int_ex_table$gene_id)){
        geneid <- int_ex_table$gene_id[k]

        ### filter features below 500 bp 

        my_intron_ranges <- my_allgenes_introns[grepl(geneid,my_allgenes_introns$gene_id),] %>%
        dplyr::select(-tx_id,-tx_name) %>% dplyr::filter(as.numeric(width)>500) %>% makeGRangesFromDataFrame(keep.extra.columns = T)
        my_exon_ranges <- my_allgenes_exons[grepl(geneid,my_allgenes_exons$gene_id),] %>%
        dplyr::select(-tx_id,-tx_name,-exon_id,-exon_name,-exon_rank)%>% dplyr::filter(width>500) %>% makeGRangesFromDataFrame(keep.extra.columns = T)



        int_ex_table[k,paste0(my_conditions[j],"_Intron")]<-as.numeric(mean(unlist(get(my_covs_chip_cond)[my_intron_ranges])))
        int_ex_table[k,paste0(my_conditions[j],"_Exon")]<-mean(as.numeric(unlist(get(my_covs_chip_cond)[my_exon_ranges])))

        print("done") 

       }
    
        
      
    }

     
    
  }
  assign(paste0("Int_Ex_table_",my_chips[i]),int_ex_table)
  
  ### K36me2/3 do not have Ash1 RNAi condition - so eliminate those empty columns 
  
  if (my_chips[i]!="K36Me1"){
  int_ex_table <- int_ex_table %>% dplyr::select(-Ash1_Intron,-Ash1_Exon)
  saveRDS(int_ex_table,paste0("Int_Ex_table_",my_chips[i],".rds"))
  } else{
  saveRDS(int_ex_table,paste0("Int_Ex_table_",my_chips[i],".rds"))
  }
}


###### plot 


pdf("Intron_Exon/Me1_filt_KDclusters.pdf",width=10,height=7.5)
readRDS("Int_Ex_table_K36Me1.rds") %>% left_join(GeneBodyFeatures[,c("gene_id","K36Me_bound_new_v2_Clusters_merged")], by="gene_id") %>% filter(K36Me_bound_new_v2_Clusters_merged %in% c("1","2_Eu","2_Het","3")) %>% pivot_longer(-c(K36Me_bound_new_v2_Clusters_merged,gene_id),names_to="name",values_to = "Mean_coverage_Me1") %>%
  mutate(Intron_Exon=ifelse(grepl("Exon",.$name),"Exon","Intron"),Condition=gsub("_Exon|_Intron","",.$name)) %>% ggplot(aes(x=factor(Condition,levels=c("GFP","Set2","NSD","Comb","Ash1")),y=Mean_coverage_Me1,fill=as.factor(Intron_Exon))) + geom_split_violin(width=1.5,draw_quantiles = c(0.5),alpha=0.8) + ylim(c(0,5)) + facet_wrap(~K36Me_bound_new_v2_Clusters_merged) + theme_minimal() +scale_fill_manual(values=c('#999999','#E69F00'))
dev.off()


pdf("Intron_Exon/Me2_filt_KDclusters.pdf",width=10,height=7.5)
readRDS("Int_Ex_table_K36Me2.rds") %>% left_join(GeneBodyFeatures[,c("gene_id","K36Me_bound_new_v2_Clusters_merged")], by="gene_id") %>% filter(K36Me_bound_new_v2_Clusters_merged %in% c("1","2_Eu","2_Het","3")) %>% pivot_longer(-c(K36Me_bound_new_v2_Clusters_merged,gene_id),names_to="name",values_to = "Mean_coverage_Me2") %>%
  mutate(Intron_Exon=ifelse(grepl("Exon",.$name),"Exon","Intron"),Condition=gsub("_Exon|_Intron","",.$name)) %>% ggplot(aes(x=factor(Condition,levels=c("GFP","Set2","NSD","Comb")),y=Mean_coverage_Me2,fill=as.factor(Intron_Exon))) + geom_split_violin(width=1.5,draw_quantiles = c(0.5),alpha=0.8) + ylim(c(0,5)) + facet_wrap(~K36Me_bound_new_v2_Clusters_merged) + theme_minimal() +scale_fill_manual(values=c('#999999','#E69F00'))
dev.off()

pdf("Intron_Exon/Me3_filt_KDclusters.pdf",width=10,height=7.5)
readRDS("Int_Ex_table_K36Me3.rds") %>% left_join(GeneBodyFeatures[,c("gene_id","K36Me_bound_new_v2_Clusters_merged")], by="gene_id") %>% filter(K36Me_bound_new_v2_Clusters_merged %in% c("1","2_Eu","2_Het","3")) %>% pivot_longer(-c(K36Me_bound_new_v2_Clusters_merged,gene_id),names_to="name",values_to = "Mean_coverage_Me3") %>%
  mutate(Intron_Exon=ifelse(grepl("Exon",.$name),"Exon","Intron"),Condition=gsub("_Exon|_Intron","",.$name)) %>% ggplot(aes(x=factor(Condition,levels=c("GFP","Set2","NSD","Comb")),y=Mean_coverage_Me3,fill=as.factor(Intron_Exon))) + geom_split_violin(width=1.5,draw_quantiles = c(0.5),alpha=0.8) + ylim(c(0,5)) + facet_wrap(~K36Me_bound_new_v2_Clusters_merged) + theme_minimal() +scale_fill_manual(values=c('#999999','#E69F00'))
dev.off()

