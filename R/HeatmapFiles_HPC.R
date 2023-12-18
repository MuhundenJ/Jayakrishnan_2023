
###### Script to run on cluster to test optimal heatmap clustering algorithm and generating intermediate files

rm(list = ls())
try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))

#### preinstall packages onto computational cluster 

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

# load genebody features table

GeneBodyFeatures <- readRDS("GeneBodyFeatures.rds")
# 

# load ave.table object containing genebody averaged signals generated in GenerateIntermediateFiles.R

ave.table <- readRDS("ave.table.rds")


# # #### clustered heatmaps 

### filter only genes bound by atleast one of K36me1/2/3 (as defined by K36Me_bound_new_v2 column of GeneBodyFeatures) for clustering -- Rearrange 13 columns in specified order

ave.table_clustering <- ave.table[rownames(ave.table) %in% GeneBodyFeatures[GeneBodyFeatures$K36Me_bound_new_v2==1,]$gene_id,grepl("K36me|K36Me",colnames(ave.table))][,c(4,8,6,2,3,7,5,1,11,14,13,10,9)]

### colors for ComplexHeatmap

col_fun_unscaled <- colorRamp2(seq(-2,3,length=11),rev(brewer.pal(n = 11, name ="RdYlBu")))


#### In pilot tests, various distance and linkage methods (shown below) as well as number of subclusters were explored to pick optimal
#### clustering algorithm based on manual inspection -- The pilot test code was skipped for conciseness

# dist_methods <- c("euclidean","maximum","manhattan","canberra","minkowski","pearson","spearman","kendall")
#link_methods <- c("ward","complete","average","mcquitty","median","centroid") # drop single

#### AFTER TESTING IT WAS DETERMINED THAT OPTIMAL CONDITIONS ARE MINKOWSKI COMPLETE WITH ROBUST DISTANCE - 12 CLUSTERS

p_comp_mink <- Heatmap(ave.table_clustering, cluster_columns= F, cluster_rows = T,clustering_method_rows = "complete",clustering_distance_rows = "minkowski",row_labels = rep("",nrow(ave.table_clustering)),col = col_fun_unscaled)
 
saveRDS(p_comp_mink,"p_complete_mink.rds")

HM_dend <- row_dend(p_comp_mink)   #extract dendgram object
saveRDS(HM_dend,"HM_dend.rds")

# HM_dend <- readRDS("HM_dend_noK.rds")
# str(HM_dend)

##### 12 clusters is a good starting point to merge clusters 

#### color 12 branches for visualization of dendrogram tree

HM_12k <- color_branches(HM_dend,k=12)

pdf(file=paste0("Heatmap_12k.pdf"), height = 15, width=7)
print(Heatmap(ave.table_clustering, cluster_columns= F, cluster_rows = HM_12k,row_labels = rep("",nrow(ave.table_clustering)),col = col_fun_unscaled,column_title="12k",row_split=12))
dev.off()

##### Below code allows us to extract the cluster identity for all genes used in clustering 

HM_12k_cutree <- cutree(HM_dend,k=12)  

saveRDS(HM_12k,"HM_Dend_12k.rds")

HM <- Heatmap(ave.table_clustering, cluster_columns= F, cluster_rows = HM_12k,row_labels = rep("",nrow(ave.table_clustering)),col = col_fun_unscaled,column_title="12k",row_split=12)
cluster_list <- row_order(HM)

saveRDS(cluster_list, file="Cluster_list.rda")

clu_df <- lapply(names(cluster_list), function(i){
  out <- data.frame(GeneID = rownames(ave.table_clustering)[cluster_list[[i]]],
                    Cluster = i,
                   stringsAsFactors = FALSE)
   return(out)
}) %>%
   do.call(rbind, .)

print(clu_df)

#### final list of cluster identities 

saveRDS(clu_df,file="Cluster_list.rda")

