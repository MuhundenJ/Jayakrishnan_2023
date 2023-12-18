### Script for visualizing K36me1/2/3 coverage over genebodies as shown in Fig4E
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

##### make cumulative plots for each cluster to distribution over bodies 

### Note -- This plot uses inner unscaled windows of 500bp as well as central genebody of 1000 bp 
### .ie. gene must be of minimal length 2000bp -- Filter accordingly 

### counts indicate that not many genes are eliminated, especially for superclusters II and III where genes tended to be longer

GeneBodyFeatures %>% mutate(small=ifelse(width<2000,1,0)) %>% group_by(K36Me_bound_new_v2_Clusters_merged) %>% tally(small)

GeneBodyFeatures %>% mutate(outlier=ifelse(width>100000 | width<2000,1,0)) %>% group_by(K36Me_bound_new_v2_Clusters_merged) %>% tally(outlier)


i <- 1

### Nested genes may alter plot profile -- identify and verify 

GeneBodyFeatures$nested_unstranded <- NA

for (i in 1:length(all_genes)){
  query <- all_genes[i]
  subject <- all_genes[-i]
  GeneBodyFeatures[GeneBodyFeatures$gene_id==query$gene_id,"nested_unstranded"]<- ifelse(overlapsAny(query,subject,type="within",ignore.strand=T),1,0)
}

saveRDS(GeneBodyFeatures,file="GeneBodyFeatures.rds")

#GeneBodyFeatures %>% filter(!is.na(nested_unstranded)) %>% group_by(K36Me_bound_new_v2_Clusters_merged) %>% select(K36Me_bound_new_v2_Clusters_merged,nested_unstranded) %>% table()

## not too many nested genes, and distributed fairly evenly in all clusters 

rm(list=ls(pattern="^mat.GB"))

### make coverage plots 

mat_path <- file.path("./matrices_genebody_innerwindow")

if(!dir.exists(mat_path))
  dir.create(mat_path)

#coveragegenebodyscaled doesnt work with granges. Create dataframes for coordinates

df_all_genes_filt <- data.frame(chr=seqnames(my_allgenes),
                          start=start(my_allgenes)-1,
                          end=end(my_allgenes),
                          strand = strand(my_allgenes),
                          width = width(my_allgenes),
                          gene_id=my_allgenes$gene_id,
                          row.names = my_allgenes$gene_id) %>% filter(width>=2000 & width<=100000)

df_all_genes_filt <- GeneBodyFeatures %>% select(K36Me_bound_new_v2_Clusters_merged,gene_id,Average_Exp) %>% right_join(df_all_genes_filt,by="gene_id")

row.names(df_all_genes_filt) <- df_all_genes_filt$gene_id

margin_outer = 500
margin_inner = 500
genebody_scaler = 1000

##### select only K36me1/2/3 replicate merged coverages for control conditions

my_covs_CM <- my_covs_v2[grep("GFP_K36",my_covs_v2)]

my_clusters <- unique(df_all_genes_filt$K36Me_bound_new_v2_Clusters_merged)

i <- 1
j <- 1

parallel::mclapply(seq_along(my_covs_CM), mc.cores = 2, FUN = function(i){
  
  for (j in seq_along(my_clusters)){

    genesets <- df_all_genes_filt %>% filter(K36Me_bound_new_v2_Clusters_merged==my_clusters[j])

    genesets <- data.frame(genesets)
    rownames(genesets) <- genesets$gene_id
    
    my_name <-  paste("mat.GENEBODY_Cluster",my_clusters[j], gsub("coverage.","", my_covs_CM[i]), sep=".")

    my_cov <- get(my_covs_CM[i])

    my_mat <- HelpersforChIPSeq::coverageGeneBodyScaled(my_coverage = my_cov,
                                                        my_coordinates = genesets,
                                                        margin_outer = margin_outer,
                                                        margin_inner = margin_inner,
                                                        genebody_scaler = genebody_scaler)


    assign(my_name, my_mat)
    fname <- file.path(mat_path,my_name)
    save(list = my_name, file = paste0(fname, ".rda", sep=""))
  }
})

#load mat.genebody files
my_mat_file_names <- list.files(path = "matrices_genebody_innerwindow/",pattern = "^mat.GENEBODY")

for(i in seq_along(my_mat_file_names)){
  load(file.path("matrices_genebody_innerwindow//", my_mat_file_names[i]))
}


########### log2 transform

i <- 1

for (i in seq_along(my_mat_file_names)){
  
  my_mat_name <- gsub(".rda$","",my_mat_file_names[i])
  my_mat <- get(my_mat_name)
  
  if (class(my_mat)=="list"){
    my_mat <- as.matrix(log2(unlist(my_mat)+0.001))
  } else {
    my_mat <- log2(my_mat+0.001)
  }
  
  assign(my_mat_name,my_mat)
  
}

########plot genebody cumulative plot

my_mats <- ls(pattern="mat.GENEBODY")
my_chips <- unique(gsub("^.*\\_","", my_mats))

j <- 1

### plot for superclusters as well as background cluster 0
my_clusters_v2 <- my_clusters


for (j in seq_along(my_clusters_v2)){
  my_mat_geneset <- my_mats[grep(paste0("Cluster.",my_clusters_v2[j]),my_mats)]

  my_mat_chip <- my_mat_geneset

  pdf(paste0(my_clusters_v2[j],".GBScaledCompsite_Clusters.pdf"), height = 5, width = 5)
    
  par(mfrow =c(1,1), mar = c(4,4,2,2), oma = c(1,1,1,1), mgp = c(2,1,0))
    
  HelpersforChIPSeq::plotComposite(my_sample_mats = my_mat_chip, 
                                     my_sub_range = 1:ncol(get(my_mat_chip[1])), 
                                     ylims =  c(-2,2), 
                                     my_binning = 1, 
                                     my_colors_composite = c("black", "red", "gray"), 
                                     my_title = paste(my_clusters_v2[j]),
                                     add_range = NULL, 
                                     smoother = 31, 
                                     line_lwd = 2, 
                                     add_axis = FALSE)  
    
    
  mtext(text = "log2 ChIP/Input", side = 2, line = 2.5)
    
    #  add axis (extra function for gene body scaled matrix)
  axisGeneBodyScaled(margin_outer = margin_outer,
                       margin_inner = margin_inner,
                       genebody_scaler = genebody_scaler,
                       add_line = TRUE)
    
  plotLegendtoPosition(conditions = factor(my_chips, levels = unique(my_chips)),
                         legend_colors = c("black", "red", "gray"), 
                         horizontal = F,
                         position = c(0.65, 0.85))

  dev.off() 
      
}

