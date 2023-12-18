################## Analysis of mouse cell line datasets from Weinberg et al 2019 and Sun et al 2023 - Fig S6


rm(list=ls())

setwd("/Users/ra36doj/Desktop/mount/ChIP_Seq/RelevantDatasets_HMTKO/Weinberg_2019/Weinberg2019_Analysis/")



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
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
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

source("./functions/functions.R")


################################  Peaks were calleed at 1kb min length for Me2 and Me3 across all KO conditions using MACS2 and 
######################### merged using BEDTOOLs.. However you still see spurious peaks which are tiny- 1.5kb length
##### manually filter it out and check 

PeakRegions <- read.table("./Peaks_Narrow/K36Me_mergedPeaks.bed")
PeakRegion_1.5 <- PeakRegions %>% filter(V3-V2>1500)
PeakRegion_2 <- PeakRegions %>% filter(V3-V2>2000)

write.table(PeakRegion_1.5,"./Peaks_Narrow/K36Me_mergedPeaks_1500bp.bed", quote=FALSE,row.names = F,col.names=F)
write.table(PeakRegion_2,"./Peaks_Narrow/K36Me_mergedPeaks_2000bp.bed", quote=FALSE,row.names = F,col.names=F)

############ Around 30% for both combined Peaksets - Similar number of genes under peaks as well

my_chromosomes_mm10 <- c(paste("chr",1:19,sep=""),"chrX","chrY")

my_lengths <- seqlengths(keepSeqlevels(TxDb.Mmusculus.UCSC.mm10.knownGene, my_chromosomes_mm10, pruning.mode = "coarse"))

my_allgenes_mm10 <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)

### use this definition 

bound_genes_2 <- my_allgenes_mm10[overlapsAny(my_allgenes_mm10,makeGRangesFromDataFrame(PeakRegion_2,seqnames.field = "V1",start.field = "V2",end.field = "V3"))]


Bound_genes_index <- ifelse(my_allgenes_mm10$gene_id %in% bound_genes_2$gene_id, 1,0)
names(Bound_genes_index) <- my_allgenes_mm10$gene_id

####################### generate coverages using Bigwig files, make average tables (containing genebody averaged ChiP signal) for K36me bound genes 
#################as performed in GenerateIntermediateFiles.R - Load those intermediate files 


################################ heatmaps - first cluster using ChIP data for K36me3 -- Use same row order to observe K36me2 as well as K36me3 CutNRun  

#### Load average table
my_avetable <- ave.my_mats

#reorder columns at this step if necessary --- First use only ChIP data 
V <- my_avetable[,grep("ChIP",colnames(my_avetable))]

#### remove 0 variance rows --- Often genes with no signal that erroneously gets flagged as being under K36me peaks (because its overlapping only a few bp etc)

var_v <- apply(V, 1, var)
var_v <- var_v[var_v != 0]
selected_avetable <- V[names(var_v),]

selectedGeneNames <- names(var_v[order(var_v, decreasing = T)][1:50])
row.labs <- ifelse(rownames(selected_avetable) %in% selectedGeneNames, rownames(selected_avetable), "")

selected_avetable_Me3 <- selected_avetable[,grep(".K36me3",colnames(selected_avetable))]  ### Filter only ChIP K36me3 columns from Weinberg et al

############# IMPORTANT : Note that scale function is applied to this specific dataset -- .ie. we aren't visualizing raw ChIP signal but rows are transformed
########## to z-scores unlike Drosophila heatmap shown in Fig 4a-- This is because unlike MNase ChIP seq, sonication based ChIP produces 3-5x lower ChIP signal
######## for NSD dependent genes relative to SetD2 dependent genes - hence to focus on patterns upon Knockout rather than abosulte signal, scaling was used
######## But to reiterate, NSD dependent genes are still classified as K36me bound from peak calling

yy_1 <- Heatmap(t(scale(t(selected_avetable_Me3))), cluster_columns= FALSE, width=8,use_raster=F,row_labels = row.labs, col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100))

Me3_cluster <- row_order(yy_1) ## store row order 


### Plot Me2 using using the clustering order derived from ChIP Me3

selected_avetable_Me2 <- selected_avetable[,grep("ChIP.K36me2",colnames(selected_avetable))]

yy_2 <- Heatmap(t(scale(t(selected_avetable_Me2))), cluster_columns= FALSE, width=8,use_raster=F,row_order=Me3_cluster,row_labels = row.labs, col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100))


### Plot Sun et al 2023 CutNRun data using the clustering order derived from ChIP Me3 -- visualize

selected_avetable_CnR <- my_avetable[,grep("CUTNRUN",colnames(my_avetable))]

## Make sure rows are ordered according to original ave.table from which clustering orders were generated
selected_avetable_CnR <- selected_avetable_CnR[rownames(selected_avetable),]   

yy_3 <- Heatmap(t(scale(t(selected_avetable_CnR))), cluster_columns= FALSE, width=8,use_raster=F,row_order=Me3_cluster,row_labels = row.labs, col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100))

########### Do genes in these clusters contain particular signatures of Intronic Proportion, transcription and gene lengths? Perform same analyses as what was done for Drosophila ############

#### Using Expression data from Weinberg et al 2019 to generate expression quantiles; Genes with no detected reads are given quantile value of 0 ############
#### while rest of the 'expressed' genes are split into 4 equal quantiles (low, medium, high, very high) of expression  ################
#### Calculate intronic proportions for Mouse genome based on pipeline used in IntronExon_Analyses.R script  ############

#### These are stored in the mm10.Features.my_mats object ! ####################


### Weinberg et al have 3 replicates for expression. Calculate quantiles for each replicate and select modal value across all replicates ############

Modal_value <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

mm10.Features.my_mats[["Exp_Quantile"]] <- rep(NA,nrow(mm10.Features.my_mats))

########### Final expression quantile is modal value across 3 replicates################

for (i in 1:nrow(mm10.Features.my_mats)){
  mm10.Features.my_mats[i,]$Exp_Quantile <- Modal_value(c(mm10.Features.my_mats[i,]$Quantile_r1,mm10.Features.my_mats[i,]$Quantile_r2,mm10.Features.my_mats[i,]$Quantile_r3))
}

## visualize 

colors_features = setNames(c("gray","blue1","steelblue1","#FFFF00","red"),c("0","1","2","3","4"))

############################################################################################################################
##### Filter and reorder the Tx Quantile dataframe to the same order used in original selected_avetable for ChIP heatmap ################

mm10.Features.my_mats_HM <- mm10.Features.my_mats[(rownames(mm10.Features.my_mats) %in% rownames(selected_avetable)),]

rownames(mm10.Features.my_mats_HM) <- mm10.Features.my_mats_HM$gene_id

mm10.Features.my_mats_HM <- mm10.Features.my_mats_HM[rownames(selected_avetable),]   



############################################################################
######## Plot individual heatmaps using same order used as before ########

tx_HM <- ComplexHeatmap::Heatmap(mm10.Features.my_mats_HM$Exp_Quantile,cluster_rows = F,name="Expression Quantile",cluster_columns = F, row_order = Me3_cluster, col=colors_features, width=2)

genLen_HM <- ComplexHeatmap::Heatmap(mm10.Features.my_mats_HM$GeneLengthQuantile,cluster_rows = F,name="Gene Length Quantile",cluster_columns = F, row_order = Me3_cluster, col=colors_features, width=2)

intronProp_HM <- ComplexHeatmap::Heatmap(mm10.Features.my_mats_HM$IntronPropQuantile,cluster_rows = F,name="Intron Prop Quantile",cluster_columns = F, row_order = Me3_cluster, col=colors_features, width=2)


heatmaps_list <- yy_1 + yy_2 + yy_3 + tx_HM + genLen_HM + intronProp_HM

pdf(file="./Me23_scaled_GenicFeatures.pdf", height = 15, width=7)

draw(heatmaps_list, row_title="Me3/2 bound genes filtered (n=16368)", column_title="Me3/2 correlation with Genic Features")

dev.off()

##########
