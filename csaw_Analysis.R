#### Scripts for implementation of CSAW analytical framework to identify differential binding regions as shown in Fig 3

#### reference (https://bioconductor.org/books/release/csawBook/) 

######### initialize packages

rm(list=ls())

setwd("/Users/ra36doj/Desktop/mount/ChIP_Seq/csaw")

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
library(edgeR)
library(csaw)
source("./csaw/functions/functions.R")
source("./csaw/functions/chromoMap_custom.R")


####################### load BAM Files for all individual replicates

my_ChIP <- c("K36Me3","K36Me2","K36Me1|K36me1")
my_control <- "Input"


bam.files <- list.files("./DiffNorm/", pattern = ".bam$")
bam.files.path <-list.files("./DiffNorm/", pattern = ".bam$",full.names = T)

##clean names
my_samples <- gsub("_[G,A,T,C].*","",bam.files)

`%nin%` = Negate(`%in%`)


########################### Counting

param <- readParam(minq=0, pe="both", max.frag = 500)  #### define params for counting - no quality filtering (as already filtered during samtools), pe=both readmates

### 2 different counting approaches - small windows and large bins. Large bins will be used for normalization later 

win.data <- windowCounts(bam.files.path, param=param, width=250, ext=150) 
bin.data <- windowCounts(bam.files.path, bin=TRUE, param=param, width=10000)   ### this large bin calculation will be used later for normlization


saveRDS(win.data, file = paste("./csaw/win.data.rds", sep="."))
saveRDS(bin.data, file = paste("./csaw/bin.data.rds", sep="."))

win.data

metadata(win.data)

colData(win.data)

rowRanges(win.data)

head(assay(win.data))

####### split data
chip <- "K36Me1|K36me1"

for (chip in my_ChIP){

  chip.win.data <- win.data[,(grepl(chip, bam.files))]
  inp.win.data <- win.data[,grep("Input", bam.files)]
  
  chip.bin.data <- bin.data[,(grepl(chip, bam.files))]
  inp.bin.data <- bin.data[,grep("Input", bam.files)]
  
  
  ### filter based on enrichment over Input and Prior minimum read counts
  filter.stat <- filterWindowsControl(data = chip.win.data, 
                                      background = inp.win.data,
                                      prior.count = 5, 
                                      scale.info = scaleControlFilter(chip.bin.data, inp.bin.data))
  
  min.fc <- 2  ## atleast 2 fold more reads in sample
  keep <- filter.stat$filter > log2(min.fc)
  
  chip.filt.data <- chip.win.data[keep,]

  saveRDS(chip.filt.data, paste0("./csaw/",chip,".filt.data.rds"))
  
  summary(keep) ### Number of bins retained should roughly match genome coverage of called peaks - around 20%
 
  ################ Normalization ###########
  
  #### CHOICE OF NORMALIZATION DEPENDS STRONGLY ON BIOLOGICAL PROBLEM 
  ### USE WHOLE GENOME TMM NORMALIZATION AS SUGGESTED IN CSAW BOOK FOR CBP KO OR K27ME3 KO EXAMPLES
  
  # calculate normFactors using whole genome bins 
  chip.win.data <- normFactors(chip.bin.data,se.out=chip.win.data)
  normfacs <- chip.win.data$norm.factors
  bin.ab <- scaledAverage(chip.bin.data)
  adjc <- calculateCPM(chip.bin.data,use.norm.factors=FALSE)
  
  
  ######### statistical modelling 
  library(edgeR)
  y <- asDGEList(chip.filt.data)
  summary(y)
  
  
  
  if (chip=="K36Me1|K36me1"){
    filenames_design <- colData(chip.win.data)$bam.files 
    
    genotype <- factor(case_when(grepl("GFP|GST",filenames_design)~"Control",
                                 grepl("Ash1",filenames_design)~"Ash1_KD",
                                 grepl("Set2",filenames_design)~"Set2_KD",
                                 grepl("NSD",filenames_design)~"NSD_KD",
                                 grepl("Comb",filenames_design)~"DKD"))
    
    ##for K36me1, only 2 replicates are available
    ## adjust here if replicates are labelled with some other unique ID
    batch <- factor(case_when(grepl("rep1",filenames_design)~"1",
                              grepl("rep2",filenames_design)~"2"))
                            
    #### batch variable important especially for K36me1 as GST from 2nd replicate is quite noisy !
    
    design <- model.matrix(~0+genotype+batch)
    colnames(design)[1:5] <- levels(genotype) 
    
  
  
  }else{
    filenames_design <- colData(chip.win.data)$bam.files 
    
    genotype <- factor(case_when(grepl("GFP",filenames_design)~"Control",
                                 grepl("Set2",filenames_design)~"Set2_KD",
                                 grepl("NSD",filenames_design)~"NSD_KD",
                                 grepl("Comb",filenames_design)~"DKD"))
    
    batch <- factor(case_when(grepl("rep1",filenames_design)~"1",
                              grepl("rep2",filenames_design)~"2",
                              grepl("rep3",filenames_design)~"3"))
    
    
    design <- model.matrix(~0+genotype+batch)
    #design <- model.matrix(~0+genotype)
    colnames(design)[1:4] <- levels(genotype) 
    
  }
  
  design
  
  y <- estimateDisp(y, design)  
  summary(y$trended.dispersion)
  
  #### diagnostic plots
  pdf(paste0("BCVplot_",chip,".pdf"))
  plotBCV(y)
  dev.off()

  fit <- glmQLFit(y, design, robust=TRUE)
  summary(fit$df.prior)

  pdf(paste0("QLDisp_",chip,".pdf"))
  plotQLDisp(fit)
  dev.off()

  pdf(paste0("MDS_",chip,".pdf"))
  plotMDS(cpm(y, log=TRUE), top=10000, labels=genotype,
          col=c("black","red", "blue","orange","brown")[as.integer(genotype)])
  dev.off()
  
  ##### Diff bind testing 
  
  #### Need to run individual contrasts separately according to BioConductor. Running it on total design returns only 1 pvalue signifying variance
  
  if (chip=="K36Me1|K36me1"){
    my_conditions <- c("Set2_KD","Ash1_KD","NSD_KD","DKD")
  } else {
    my_conditions <- c("Set2_KD","NSD_KD","DKD")
  }
  
  ### loop through condition and set contrast as 'KD-Control' as input for makeContrasts
  
  condition <- "Set2_KD"
  for (condition in my_conditions){
    contrast_cond <- makeContrasts(paste0(condition,"-Control"), levels=design)
    res_cond <- glmQLFTest(fit, contrast=contrast_cond)
    
    ### merge small 150bp windows to max of 3kb 
    ### 3kb was chosen based on trial and error. Using too large values (10s of kb) causes chaining of domains, especially because Drosophila genome is quite small
    ### for instance, 2 Set2 dependent regions on either side of a small NSD region will chain to form a large Set2 domain - This leads to artificial 
    ### overlap between conditions 
    
    merged_cond<- mergeResults(chip.filt.data, res_cond$table, tol=100, 
                               merge.args=list(max.width=3000))
    
    
    out.ranges_cond <- merged_cond$regions
    summary_cond <- merged_cond$combined 
    mcols(out.ranges_cond) <-data.frame(summary_cond)
    
    #### filter for only significant windows and output in BED9 format for IGV visualization
    
    out.ranges_cond <- data.frame(out.ranges_cond) %>% filter(FDR<=0.05) %>% transmute(chrom=seqnames,
                                                                                        chromStart=start,
                                                                                        chromEnd=end,
                                                                                        name=paste0(seqnames,":",start,"-",end),
                                                                                        score=rep.logFC,
                                                                                        strand="*",
                                                                                        thickStart=start,
                                                                                        thickEnd=end,
                                                                                        itemRGB=case_when(direction=="down"~"0,0,255",
                                                                                                          direction=="up"~"255,0,0",
                                                                                                          direction=="mixed"~"0,0,0"))
    
    write_delim(out.ranges_cond,paste0("./csaw/DiffRegions_",chip,"_",condition,".bed"),col_names=F, quote="none",delim="\t")
    
    }
}

##### generate chromoMaps of differential regions 


library("chromoMap")
library("htmlwidgets")
## use Chr annotation file from before (Control_chromoMaps.R)

my_beds <- list.files("./csaw/",pattern="*.bed$",full.names=T)  ### list bed files generated by CSAW 

i <- 1

for (i in seq_along(my_beds)){

  my_bed <- read_delim(my_beds[i],col_select = c(4,1,2,3,5),col_names = F) %>% write_delim(paste0("./csaw/chromoMaps/",gsub(".bed","",gsub(".*_tol100_","",my_beds[i])),".txt"),quote=NULL,delim="\t",col_names = F)
  
  ### caution! As these are high resolution chromoMaps (~2kb), the chromomaps will be really long on html browser. 
  
  x <- chromoMap_custom("./csaw/chromoMaps/chr_file_anno_noCentro.txt",paste0("./csaw/chromoMaps/",gsub(".bed","",gsub(".*_tol100_","",my_beds[i])),".txt"),scale_lims = scale_lims,data_based_color_map=T,data_type="numeric",data_colors = list(c("blue","white","red")), title=gsub(".*_tol100_","",my_beds[i]),title_font_size=15, lg_y=350,legend = T,chr_color = c("white"),win.summary.display = T,interactivity = F, n_win.factor = 10,chr_width = 60,ch_gap = 30)
  
  saveWidget(x,paste0("./csaw/chromoMaps/",gsub(".*_tol100_","",my_beds[i]),".html"),selfcontained = F)
  
  }


### annotation chromomaps for genes and K9me2 peaks

my_genes <- as.data.frame(genes(TxDb.Dmelanogaster.UCSC.dm6.ensGene))[,c("gene_id","seqnames","start","end")]  %>% write_delim("./csaw/chromoMaps/gene_anno.txt",col_names = F,quote=NULL,delim="\t")

my_K9me2 <- read_delim("./csaw/K9me2.narrowPeak",col_names = F)[,c(4,1,2,3)]%>% write_delim("./csaw/chromoMaps/K9me2_anno.txt",col_names = F,quote=NULL,delim="\t")


x <- chromoMap("./csaw/chromoMaps/chr_file_anno.txt","./csaw/chromoMaps/gene_anno.txt",title_font_size=15, lg_y=350,legend = T,anno_col = c("black"),chr_color = c("white"),win.summary.display = T,interactivity = F, n_win.factor = 10,chr_width = 60,ch_gap = 30)
saveWidget(x,"./csaw/chromoMaps/gene_anno.html",selfcontained = F)


x <- chromoMap("./csaw/chromoMaps/chr_file_anno.txt","./csaw/chromoMaps/K9me2_anno.txt",title_font_size=15, lg_y=350,legend = T,anno_col = c("black"),chr_color = c("white"),win.summary.display = T,interactivity = F, n_win.factor = 10,chr_width = 60,ch_gap = 30)
saveWidget(x,"./csaw/chromoMaps/K9me2_anno.html",selfcontained = F)