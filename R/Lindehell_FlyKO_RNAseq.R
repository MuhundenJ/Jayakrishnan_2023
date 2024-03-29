######### #Script to analyze RNAseq data from Lindehell 2021 et al --- Venn diagrams in Fig S3

library(tidyverse)
library(devtools)
library(purrr)
library(HelpersforDESeq2)
library("DESeq2")

setwd("./Lindehell_2021_SciAdv/Lindehell_2021_Analysis")

###############get directories of STAR aligner output (see RNAseq snakemake file) and import only the unstranded ReadsPerGene for all samples  ########################

dir <- dir(path="../", pattern="_out")

i <- 1

test_file_path <- vector(mode="list", length=length(dir))

for (i in seq_along(dir)){

test_file_path[i] <- list.files(path = paste("../",dir[i],"/",sep=""), pattern="ReadsPerGene.out.tab",full.names=TRUE)

}

read_RPG <- function(x) {
  base_name <- basename(x) %>% gsub(".ReadsPerGene.out.tab","", x=.)
  read.table(file=x, col.names = c("ID",base_name,"NULL","NULL"),colClasses = c("factor","integer","NULL","NULL"))
}

raw_RPG <- test_file_path %>% map_dfc(~ read_RPG(.)) %>% select(-paste("ID...",seq(3,55,by=2), sep=""))

colnames(raw_RPG)[1] = "ID"

write.table(raw_RPG, file = "raw_RPG.txt", quote =FALSE, row.names = FALSE)

####### setup DESeq2

sample_table <- read.delim(file="SampleTable_v2.txt", header = TRUE, stringsAsFactors = FALSE)
sample_table$Condition <- as.factor(sample_table$Condition)
sample_table$Batch <- as.factor(sample_table$Batch)

raw_RPG_genes <- raw_RPG[-1:-4,]

######## make compatible format for HelpersforDeseq2 - named matrix 

rownames(raw_RPG_genes) <- raw_RPG_genes$ID
raw_RPG_genes_processed <- raw_RPG_genes %>% select(-ID)

#change colnames with Sample ID from sample table 
names(raw_RPG_genes_processed)[match(sample_table$Index, names(raw_RPG_genes_processed))] <- as.character(sample_table$`ID`)

########## setup DDS function ################################################################################

#note that here we don't know the true batch of the experiments 
### also n_samples_for_filtering filters out genes which have less than min_reads in 25% of samples) 

if(identical(colnames(raw_RPG_genes_processed), sample_table$ID)){
  
  dds <- setupDDS(CountTableName = "raw_RPG_genes_processed", 
                  SampleTableName = "sample_table",
                  SampleIdName = "ID", 
                  ConditionName = "Condition", 
                  BatchName = NULL,
                  n_samples_for_filtering = ncol(raw_RPG_genes_processed)*0.75,
                  min_number_of_reads = 1)
}

save(dds, file="dds_Uncorrected.rda")
############### explore counts################################################################################

my_conditions <- colData(dds)$Sample

################ check removing multimapper reads !

library(sva)

## log2 normalization

log2_counts_uncor <- log2(counts(dds, normalized = TRUE)+1)


###### Note that a modified getResults function is used which is edited from original version in HelpersFromDESEQ2 package to subset contrasts (see Functions.R)
##### Original function assumes a different format for the GTF file -- As well as uses an outdated shrinkage method 
#### Note that AshR shrinkage is used (instead of 'normal') as now 'normal' is deprecated based on description by creators of DESEQ2 

#### Return results table for interested contrasts -- KOs vs WT

contrast_list <- list(c("Set21","OregonR"),
                      c("Nsd46","OregonR"),
                      c("Ash1","OregonR"))
                      

for(i in seq_along(contrast_list)){
  
  res <- getResults_v2_shrink(dds = dds,
                              contrast = contrast_list[i],
                              lfc_cutoff = 0,
                              shrink = TRUE,
                              annotation = "gtf",
                              anno_symbol = "gene_name",
                              anno_id = "gene_id")
  
  
  res_name <- paste0("res_ashrShrink_", as.character(map(contrast_list[i],1)),"-",as.character(map(contrast_list[i],2)))
  res$symbol <- rownames(res)
  assign(res_name, res)
  
  write.table(res,file = paste0(res_name, ".txt"), 
              quote = F, sep = "\t", row.names = T, col.names = NA) 
  
}

##### Explore results -- Here only Venn Plotting is shown

i = 1

res_names <- ls(pattern = "^res_ashrShrink")   ### Load results 

### convert results to DataFrame and keep only significant ones -- Atleast 25 % down

library(Vennerable)
library(grid)
library(gridExtra)

## keep only FBgn IDs --- Sigifnicntly Downregulated genes 

Ash1_down <- as.data.frame(`res_ashrShrink_Ash1-OregonR`) %>% filter(padj<0.05 & log2FoldChange < -0.41 & grepl("FBgn",symbol)) %>% select(symbol)  ## decrease by 25% 
Set2_down <- as.data.frame(`res_ashrShrink_Set21-OregonR`) %>% filter(padj<0.05 & log2FoldChange < -0.41 & grepl("FBgn",symbol)) %>% select(symbol)
NSD_down <- as.data.frame(`res_ashrShrink_Nsd46-OregonR`) %>% filter(padj<0.05 & log2FoldChange < -0.41 & grepl("FBgn",symbol)) %>% select(symbol)

pdf(paste0("./venn.downregulated_KO.pdf"), width = 6, height = 6)


my_overlaps_list <- list(Ash1_down$symbol,Set2_down$symbol,NSD_down$symbol)
names(my_overlaps_list) <- c("Ash1-KO","Set2-KO","NSD-KO")
my_overlaps_venn <- Venn(my_overlaps_list)

my_overlaps_venn_plot <- compute.Venn(my_overlaps_venn,doWeights = T)
SetLabels <- VennGetSetLabels(my_overlaps_venn_plot)
SetFaceLabels <- VennGetFaceLabels(my_overlaps_venn_plot)

my_overlaps_venn_plot <- Vennerable:::VennSetFaceLabels(my_overlaps_venn_plot,SetFaceLabels)
gp <- VennThemes(my_overlaps_venn_plot)

plot(my_overlaps_venn_plot, show = list(Faces = FALSE), gp = gp)

dev.off()


## keep only FBgn IDs --- Significantly Upregulated genes 

Ash1_up <- as.data.frame(`res_ashrShrink_Ash1-OregonR`) %>% filter(padj<0.05 & log2FoldChange > 0.32 & grepl("FBgn",symbol)) %>% select(symbol)  ### increase by 25% 
Set2_up <- as.data.frame(`res_ashrShrink_Set21-OregonR`) %>% filter(padj<0.05 & log2FoldChange > 0.32 & grepl("FBgn",symbol)) %>% select(symbol)
NSD_up <- as.data.frame(`res_ashrShrink_Nsd46-OregonR`) %>% filter(padj<0.05 & log2FoldChange > 0.32 & grepl("FBgn",symbol)) %>% select(symbol)

pdf(paste0("./venn.upregulated_KO.pdf"), width = 6, height = 6)

my_overlaps_list <- list(Ash1_up$symbol,Set2_up$symbol,NSD_up$symbol)
names(my_overlaps_list) <- c("Ash1-KO","Set2-KO","NSD-KO")
my_overlaps_venn <- Venn(my_overlaps_list)

my_overlaps_venn_plot <- compute.Venn(my_overlaps_venn,doWeights = T)
SetLabels <- VennGetSetLabels(my_overlaps_venn_plot)
SetFaceLabels <- VennGetFaceLabels(my_overlaps_venn_plot)

my_overlaps_venn_plot <- Vennerable:::VennSetFaceLabels(my_overlaps_venn_plot,SetFaceLabels)
gp <- VennThemes(my_overlaps_venn_plot)

plot(my_overlaps_venn_plot, show = list(Faces = FALSE), gp = gp)

dev.off()

