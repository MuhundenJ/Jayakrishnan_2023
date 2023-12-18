#### Script to visualize ChIP signal at transposons 


rm(list=ls())

setwd("/Users/ra36doj/Desktop/mount/ChIP_Seq/MultiMapper_Analysis/MultiMapper_Analysis")

library(tidyverse)
library(GenomicAlignments)
library(Rsamtools)
library(rtracklayer)
library(seqinr)
library(ape) 
library(seqinr)
library(stringr)
library(RColorBrewer)
library(matrixStats)
library(IRanges)
library(ShortRead)
library(zoo)
library(RColorBrewer)


`%!in%` = Negate(`%in%`)

## read in virilis normalized bigwig files -- Of note, we verified that virilis scaling didn't change any major patterns 
### we also explored normalizing to total number of reads aligned to either a) repetitive genome only or b) repetitive+non repetitive dm6 genomes. These also didn't make differences 

#################### load BW files - dVir normalized ############################################################

bigwig_dir <-"./bedgraphs_Input_Pseudo/"

bigwig_files <- list.files(bigwig_dir, "dvir_fullgen_Norm.dir.bw")
bigwig_files_path <- file.path(bigwig_dir, bigwig_files)

cov_path <- file.path("./coverage")

if(!dir.exists(cov_path))
  dir.create(cov_path)

i <- 1

parallel::mclapply(seq_along(bigwig_files), mc.cores = 2, FUN = function(i){
  
  
  my_name <- paste("coverage.",stringr::str_extract(bigwig_files[i], "[^_]*_[^_]*"), sep="")
  
  my_bigwig <- import.bw(bigwig_files_path[i])
  
  ### filter unnecessary elements later 

  my_cov <- coverage(my_bigwig, weight = "score")
  
  assign(my_name, my_cov)
  
  fname <- file.path(cov_path, my_name)
  save(list = my_name, file = paste0(fname,".rda"))
  
})


my_coverage_file_names <- list.files(path = "coverage/",pattern = "^coverage.*rda$")

for(i in seq_along(my_coverage_file_names)){
  
  load(file.path("coverage/", my_coverage_file_names[i]))
  
}

my_covs <- ls(pattern = "^coverage\\.")


### some elements are missing from certain coverages - check

TP_unique_list <- c()
TP_drop <- c()

for (i in seq_along(my_covs)){
  
  TP_unique_list <- unique(c(TP_unique_list,seqlevels(get(my_covs[i]))))
  
}

for (i in seq_along(my_covs)){
  
  TP_drop<- unique(c(TP_drop,TP_unique_list[TP_unique_list %!in% seqlevels(get(my_covs[i]))]))
  
}

TP_final <- TP_unique_list[TP_unique_list %!in% TP_drop]


### looks like 4 elements are missing - these are the custom Gypsi, Hoppel-1 Dvir//Uvir, TARTdvir that we included -- Exclude for further analyses
###### make average list 
i=1

my_TP_table <- data.frame(TP_id = TP_final)

for (i in seq_along(my_covs)){
  
    my_name <-  paste("mat.TP",gsub("coverage.","", my_covs[i]), sep=".")
    
    my_cov <- get(my_covs[i])
    
    my_cov <- my_cov[seqlevels(my_cov) %in% TP_final]

    my_mat <- as.data.frame(mean(my_cov))
    
    colnames(my_mat) <-gsub("coverage.","",my_covs[i])
    
    my_mat$TP_id <- rownames(my_mat)
    
    my_TP_table <<- my_TP_table %>% left_join(my_mat, by="TP_id")
    
}

rownames(my_TP_table) <- gsub("gb.*\\|","",my_TP_table$TP_id)

my_TP_table <- my_TP_table %>% select(!TP_id) 

my_TP_table <- my_TP_table[,c(4,8,6,2,3,7,5,1)]

saveRDS(my_TP_table, file = paste0("Transposon_AveCov", ".rda"))

####### Heatmaps

library(ComplexHeatmap)
library(circlize)



me3_HM <- Heatmap(my_TP_table[,1:4], cluster_columns= FALSE, cluster_rows = TRUE, name="K36Me3 occupancy",row_labels = rownames(my_TP_table),width = 10,col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100), row_names_gp = gpar(fontsize = 8))

me3_cluster <- row_order(me3_HM)

me2_HM <- Heatmap(my_TP_table[,5:8], cluster_columns= FALSE, cluster_rows = FALSE,row_order = me3_cluster, name="K36Me2 occupancy",row_labels = rownames(my_TP_table),width = 10,col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100), row_names_gp = gpar(fontsize = 8))

my_HM <- me3_HM + me2_HM

pdf(file="Transposon_Me23Coverage.pdf", height = 15, width=10)

draw(my_HM,row_title="Transposons (n=139)", column_title="Me2/3 distribution at Transposons")

dev.off()