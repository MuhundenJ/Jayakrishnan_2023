



############################################################################################################################
############################################################################################################################
############################################################################################################################


################################################    make Count Table   ##################################################### 



makeCountTable <- function(count_files, count_file_path, stranded = FALSE){
      
      if(stranded){
            cidx <- 4
      } else {
            cidx <- 2
      }
      
      for(i in seq_along(count_files)){
            
            tmp <- read.table(count_file_path[i])
            
            if(i == 1){
                  my_counts <- tmp[,cidx] 
            } else {
                  my_counts <- cbind(my_counts, tmp[,cidx])
            }
      }
      
      rownames(my_counts) <- tmp[,1]
      colnames(my_counts) <- gsub("_[G,A,T,C].*","", count_files)
      
      return(my_counts)
}




############################################################################################################################
############################################################################################################################
############################################################################################################################


#####################################################      TPM      ######################################################## 





countToTpm <- function(counts, effLen, scaler=1e6)
{
      rate <- log(counts) - log(effLen)
      denom <- log(sum(exp(rate)))
      exp(rate - denom + log(scaler))
}





############################################################################################################################
############################################################################################################################
############################################################################################################################


#####################################################      PCA      ######################################################## 






plottingPCA <- function(my_data, 
                        xcomp = 1,
                        ycomp = 2,
                        color_palette,
                        conditions,
                        quantiles = c(0,1),
                        show_labels = TRUE,
                        point_size = 1.5){
      
      rv <- rowVars(my_data)

      selection <- (rv >  quantile(rv, quantiles[1])  & rv < quantile(rv, quantiles[2]))
      
      
      pca <- prcomp(t(my_data[selection, ]), scale. = TRUE)
      
      percentVar <- round(pca$sdev^2/sum(pca$sdev^2)*100,1)[1:10]

      

      plot(pca$x[, xcomp], pca$x[, ycomp]*-1, 
           col = color_palette[conditions], 
           pch=16, cex = point_size,
           xlab = paste("PC",xcomp," (", percentVar[xcomp], "%)", sep=""),
           ylab = paste("PC",ycomp," (", percentVar[ycomp], "%)", sep=""),
           #xlim= my_limits, ylim=my_limits,
           main="PCA")
      

      if(show_labels){
            text(pca$x[, xcomp], pca$x[, ycomp]*-1, labels = rownames(pca$x), 
                 adj = -0.5, col = "gray32", cex=0.5)      
      }
      
      
      
}



############################################################################################################################
############################################################################################################################
############################################################################################################################

###################################################       dot plots      ################################################### 



plotDots <- function(my_data, 
                     my_title,
                     color_palette,
                     color_groups,
                     conditions,
                     point_size = 1.5,
                     ylims = c(0, 15),
                     x_label = NULL){
      

      
      grouped_data = data.frame(exp = as.numeric(my_data),
                                group =  as.numeric(conditions))
      
      plot(grouped_data$group + runif(length(grouped_data$group), -0.05, 0.05),
           grouped_data$exp,
           pch=19, cex = point_size, ylim =  ylims,
           xaxt = "n", xlab = "", ylab = "",
           main= my_title, col =  color_palette[color_groups])
      
      if(!(is.null(x_label))){
            axis(side = 1, at = seq_along(unique(grouped_data$group)), labels = x_label)
            
      } else {
            axis(side = 1, labels = FALSE)
            
      }
      
  
}











############################################################################################################################
############################################################################################################################
############################################################################################################################

###################################################       Heatmaps       ################################################### 


# plotHeatmap <- function(my_mat,
#                         my_row_order,
#                         my_col_order,
#                         min_value = 0,
#                         max_value = 15,
#                         my_title,
#                         my_color_palette,
#                         show_xaxis = FALSE,
#                         useRaster = TRUE,
#                         cex.main = 1.25){
#       
#       my_mat[my_mat < min_value] <- min_value
#       my_mat[my_mat > max_value] <- max_value
#       
#       image(t(my_mat[my_row_order, my_col_order]), 
#             main = my_title, cex.main = cex.main,
#             col = my_color_palette, 
#             breaks =  seq(min_value, max_value, length.out = 101),
#             axes=FALSE, useRaster = useRaster)
#       
#       if(show_xaxis){
#             axis(side = 1, at = c(0,0.5,1), lwd = 0, lwd.ticks = 1, las=1,tck = -0.15, cex.axis = 0.8,
#                  labels = c( paste(round(as.integer(colnames(my_mat))[1] / 1000), "kb", sep=""), "",
#                              paste("+",round(as.integer(colnames(my_mat))[ncol(my_mat)] / 1000), "kb", sep="")
#                  ))
#       }
#       
#     
# }

plotHeatmap_HCSeq <- function(my_sample_mats,
                        my_sample_names = "",
                        my_site_name = "0",
                        font_size = 0.75,
                        my_colors = colorRampPalette(brewer.pal(9, "Blues"))(100),
                        #min_value = 1,
                        #max_value = 7,
                        my_binning = 10,
                        smoother = 1){
  
  hidx <- seq(0, 1, length.out = length(my_sample_mats)+1)
  
  y_data <- colMeans(get(my_sample_mats[1]), na.rm = TRUE)[1:ncol(get(my_sample_mats[1]))] 
  
  min_value <- min(y_data, na.rm=TRUE)*0.9
  max_value <- max(y_data, na.rm=TRUE)*1.1
  
  for(i in seq_along(my_sample_mats)){
    
    my_sample <- my_sample_mats[i]
    
    my_mat <- get(my_sample)
    x_range <- ncol(my_mat)
    
    ###########################################
    
    my_mat[my_mat < min_value] <- min_value
    my_mat[my_mat > max_value] <- max_value
    
    my_breaks <- seq(min_value, max_value, length.out = 101)
    
    ###########################################
    
    par(mar=c(3,1,4,1), cex = font_size)
    par(fig=c(hidx[i],hidx[i+1],0.075,1), new=TRUE)
    
    if(smoother > 1){
      my_mat <- t(apply(my_mat, 1, function(x){
        zoo::rollmean(c(rep(x[1], (smoother-1)/2), x, rep(x[x_range], (smoother-1)/2)),
                      smoother)}))
    }
    
    
    image(t(my_mat)[,nrow(my_mat):1],
          col = my_colors,
          breaks =  my_breaks,
          axes=FALSE,
          useRaster = TRUE)
    
    
    ###########################################
    
    my_title <- my_sample_names[i]
    
    title(main = my_title, line = 1)
    
    axis(side = 1, at = c(0,0.5,1), labels = c("", gsub(".*\\.","",my_site_name), ""))
    axis(side = 1, at = c(0.1,0.9),
         labels = c(paste("-",round((x_range/2)*my_binning/1000), "kb", sep=""),
                    paste("+",round((x_range/2)*my_binning/1000), "kb", sep="")),
         col.ticks = NA, col = NA)
    
    ###########################################
    
    par(fig=c(hidx[i],hidx[i+1],0,0.12), new=TRUE, mar=c(3,1,2,1))
    
    image(matrix(seq_along(my_colors)), col=my_colors, axes=FALSE, useRaster = TRUE)
    axis(side = 1, at = c(0,0.5,1), labels = round(my_breaks[c(1,51,101)],1))
    
  }
   
}

                                                               





plotHeatmapKey <- function(my_mat,
                        min_value = 0,
                        max_value = 15,
                        my_title,
                        my_color_palette,
                        show_yaxis = FALSE,
                        show_xaxis = FALSE,
                        useRaster = TRUE,
                        cex.main = 1.25,
                        cex.axis = 0.6){
      
      my_mat[my_mat < min_value] <- min_value
      my_mat[my_mat > max_value] <- max_value
      
      image((my_mat), 
            main = my_title, cex.main = cex.main,
            col = my_color_palette, 
            breaks =  seq(min_value, max_value, length.out = 101),
            axes=FALSE, useRaster = useRaster)
      
      if(show_yaxis){
            axis(side = 4, at = c(0,1), lwd = 0, lwd.ticks = 1, las=1, cex.axis = cex.axis, tck = -0.25,
                 labels = round(c(min_value,max_value)))
      }
      
      if(show_xaxis){
            axis(side = 1, at = c(0,1), lwd = 0, lwd.ticks = 1, las=1, cex.axis = cex.axis, tck = -0.25,
                 labels = round(c(min_value,max_value)))
      }
      
      
}



############################################################################################################################
############################################################################################################################
############################################################################################################################

###################################################      GO analysis      ################################################### 





createGOTable <- function(allGenes, 
                        shown_terms = 10){
      
      
      tgd <- new( "topGOdata", 
                  ontology="BP", 
                  allGenes = allGenes, 
                  nodeSize=5,
                  annot=annFUN.org, 
                  mapping="org.Dm.eg.db", 
                  ID = "ensembl" )
      
      resultTopGO <- runTest(tgd, algorithm = "classic", statistic = "Fisher" )
      
      
      
      bp <- GenTable( tgd,
                      Fisher.classic = resultTopGO,
                      orderBy = "Fisher.classic" , topNodes = 200)
      
      
      head(bp, shown_terms)
}



plotTable <- function(my_table,
                      xpos = 0.6,
                      ypos = 1){
      
      mytheme <- gridExtra::ttheme_default(
            core = list(fg_params=list(cex = 0.5)),
            colhead = list(fg_params=list(cex = 0.5)),
            rowhead = list(fg_params=list(cex = 0.5)))
      
      myt <- tableGrob(my_table, theme = mytheme, rows = NULL)
      
      myt$widths <-  unit(c(16,42.5,12.5,12.5,12,15), "mm")
      
      sample_vp <- viewport(x = xpos, y = ypos, width = 0.25, height = 0.25,just = c("left", "top"))
      pushViewport(sample_vp)
      grid.draw(myt)
      popViewport()
      
}



############################################################################################################################
############################################################################################################################
############################################################################################################################


################################################       plot Composite     ################################################## 




plotComposite <- function(my_sample_mats, 
                          my_sub_range = 1:400,
                          #ylims = c(1,5),
                          my_binning = 10,
                          my_colors_composite,
                          my_title = "",
                          site_label = 0,
                          add_line = FALSE,
                          line_lwd = 2,
                          smoother = 11){
      
      
      x_range <- ncol(get(my_sample_mats[1])[,my_sub_range])
      y_data <- colMeans(get(my_sample_mats[1]), na.rm = TRUE)[my_sub_range] 
      
      ylims = c(min(y_data,na.rm=TRUE)*0.9,max(y_data,na.rm=TRUE)*1.1)
      
      plot(1:x_range, y_data, xaxt = "n", 
           main = my_title, xlab = "", ylab = "",
           type="n", ylim = ylims)
      
      #abline(h=1)
      
      
      axis(side = 1, at = seq(1, x_range , length.out = 3),
           labels =  c( paste("-",round((x_range/2)*my_binning/1000), "kb", sep="") ,
                        site_label,
                        paste("+",round((x_range/2)*my_binning/1000), "kb", sep=""))
      )
      
      
      for(i in seq_along(my_sample_mats)){
            
            x_range <- ncol(get(my_sample_mats[i])[,my_sub_range])
            y_data <- colMeans(get(my_sample_mats[i]), na.rm = TRUE)[my_sub_range] 
            
            my_yerror <- apply(get(my_sample_mats[i])[,my_sub_range], 2, function(x){ qt(0.975, df = length(x)-1)*sd(x, na.rm = TRUE)/sqrt(length(x)) })
            
            y_data <- zoo::rollmean(y_data, smoother)
            my_yerror <- zoo::rollmean(my_yerror, smoother)
            
            
            xx <- c(((smoother-1)/2+1):(x_range-((smoother-1)/2)), (x_range-((smoother-1)/2)):((smoother-1)/2+1))
            yy <- c(y_data-my_yerror, rev(y_data+my_yerror))
            
            #polygon(xx, yy, border = NA, col = paste(my_colors_composite[i], "55",sep=""))
          
            
            lines(((smoother-1)/2+1):(x_range-((smoother-1)/2)), y_data, col = my_colors_composite[i], lwd=line_lwd)
            

      }
      
      if(add_line){
            abline(h = c(min(y_data) , max(y_data)), lty=2)
      }
      
}



############################################################################################################################
############################################################################################################################
############################################################################################################################

################################################      plot Coverage      ################################################### 




plotCoverage <- function(my_coverage, 
                         my_region = c("chr2R", 7292189,  7436822),
                         my_bin_size = 100,
                         my_color = rgb(220, 161, 55, maxColorValue = 255),
                         my_ylims = c(1.1,15),
                         show_yaxis = TRUE,
                         my_title = "Track",
                         my_peaks = NULL,
                         scale = FALSE, 
                         my_lwd = 0.1){
      
      par(mar = c(0,3,2,0))
      
      my_signal <- as.numeric(my_coverage[[my_region[1]]])[my_region[2]:my_region[3]]
      
      if(scale){
            my_signal <- scale_min_max(my_signal, lower = 0.001, upper = 0.999)
      }
      
      my_cov_subset <- data.frame(signal = my_signal[1 : (length(my_signal) -  (length(my_signal) %% my_bin_size)  )],
                                  binner = rep(1:floor((length(my_signal)/my_bin_size)), each = my_bin_size))
      
      my_cov_subset <- aggregate(my_cov_subset$signal, by = list(my_cov_subset$binner), FUN = mean, na.rm = TRUE )
      my_cov_subset <- my_cov_subset$x
      
      
      plot(my_cov_subset, type="n", ylim = c(my_ylims[1],my_ylims[2]), xaxt = "n", yaxt = "n", ylab = "", main = "", bty = "n")
      title(main = my_title, line = 1, col.main = my_color, cex.main = 1)
      
      
      if(show_yaxis){
            axis(side = 2, at = c(round(my_ylims[1]), mean(round(my_ylims)), round(my_ylims[2])),
                 labels = c(round(my_ylims[1]), "", round(my_ylims[2],1)))  
      }
      
      
      xx <- c(rev(seq_along(my_cov_subset)), (seq_along(my_cov_subset)))
      
      yy <- c(rep(0, length(my_cov_subset)), my_cov_subset)
      
      
      polygon(xx, yy, col = my_color, border = my_color, lwd =1.5)
      
      
      if(!(is.null(my_peaks))){
            
            
            my_region_gr <- makeGRangesFromDataFrame(data.frame(chr = my_region[1], start = my_region[2], end = my_region[3]))
            
            my_peaks_subset <- subsetByOverlaps(my_peaks, my_region_gr)
            
            if(length(my_peaks_subset) != 0){
                  
                  my_starts <- (start(my_peaks_subset) - as.numeric(my_region[2])) / my_bin_size
                  my_ends   <- (end(my_peaks_subset) - as.numeric(my_region[2])) / my_bin_size
                  
                  par(xpd = NA)
                  rect(xleft = my_starts, xright = my_ends, 
                       ytop = my_ylims[1] + ((my_ylims[2] - my_ylims[1])*1.35), 
                       ybottom = my_ylims[1] + ((my_ylims[2] - my_ylims[1])*1.15), 
                       border = my_color, col = my_color, lwd = my_lwd)    
                  par(xpd = FALSE)
                  
            }
            
      }
      

      
      
}





plotAnnotation <- function(my_genes = my_genes, 
                           my_exons = my_exons,
                           my_region = c("chr2R", 7292189,  7436822),
                           x_scale = 10^5,
                           x_round = 1,
                           my_peaks = NULL,
                           my_ylims =  c(-1,1)){
      
      par(mar = c(0,3,1,0))
      
      my_region_gr <- makeGRangesFromDataFrame(data.frame(chr = my_region[1], start = my_region[2], end = my_region[3]))
      
      my_span <- as.numeric(my_region[3])-as.numeric(my_region[2])
      
      
      plot(1:my_span, rep(0, my_span), ylim = my_ylims,
           type="n", xaxt = "n", yaxt = "n", ylab = "", xlab = "" , bty = "n")
      
      axis(side = 1, at =  seq(1, my_span+1, x_scale), line = -1,
           labels =  round( seq(as.numeric(my_region[2]), as.numeric(my_region[3]), x_scale)/10^6, x_round))
      
      mtext(paste0("chr", my_region[1], " [Mb]"), side = 1, line = 1.5, outer = FALSE, cex = 1.33)
      
      my_genes_subset <- subsetByOverlaps(my_genes, my_region_gr)
      my_starts <- start(my_genes_subset) - as.numeric(my_region[2])
      my_ends <- end(my_genes_subset) - as.numeric(my_region[2])
      
      # text(rowMeans(cbind(my_starts[!(is.na(my_genes_subset$symbol))], my_ends[!(is.na(my_genes_subset$symbol))]))[seq(1, length(my_starts),2)], 
      #      c(-0.1,-0.5), 
      #      labels = my_genes_subset$symbol[!(is.na(my_genes_subset$symbol))][seq(1, length(my_starts),2)], 
      #      col = "darkgrey", cex = 1.25)
      
      
      rect(xleft =  my_starts[as.logical(strand(my_genes_subset) == "+")], 
           xright = my_ends[as.logical(strand(my_genes_subset) == "+")], 
           ytop =  0.75, ybottom = 0.75, border = "darkgrey", lwd=3, col = "darkgrey")
      
      rect(xleft =  my_starts[as.logical(strand(my_genes_subset) == "-")], 
           xright = my_ends[as.logical(strand(my_genes_subset) == "-")], 
           ytop =  0.25, ybottom = 0.25, border = "darkgrey", lwd=3, col = "darkgrey")
      
      my_exons_subset <- subsetByOverlaps(my_exons, my_region_gr)
      my_starts <- start(my_exons_subset) - as.numeric(my_region[2])
      my_ends <- end(my_exons_subset) - as.numeric(my_region[2])                 
      
      rect(xleft =  my_starts[as.logical(strand(my_exons_subset) == "+")], 
           xright = my_ends[as.logical(strand(my_exons_subset) == "+")], 
           ytop =  0.55, ybottom = 1, border = "darkgrey", lwd=1, col = "darkgrey")
      
      rect(xleft =  my_starts[as.logical(strand(my_exons_subset) == "-")], 
           xright = my_ends[as.logical(strand(my_exons_subset) == "-")], 
           ytop =  0, ybottom = 0.45, border = "darkgrey", lwd=1, col = "darkgrey")
 
      
      if(!(is.null(my_peaks))){
            
            
            my_region_gr <- makeGRangesFromDataFrame(data.frame(chr = my_region[1], start = my_region[2], end = my_region[3]))
            
            my_peaks_subset <- subsetByOverlaps(my_peaks, my_region_gr)
            
            if(length(my_peaks_subset) != 0){
                  
                  my_starts <- (start(my_peaks_subset) - as.numeric(my_region[2])) 
                  my_ends   <- (end(my_peaks_subset) - as.numeric(my_region[2])) 
                  
                  par(xpd = NA)
                  rect(xleft = my_starts, xright = my_ends, 
                       ytop = my_ylims[1] + ((my_ylims[2] - my_ylims[1])*1.35), 
                       ybottom = my_ylims[1] + ((my_ylims[2] - my_ylims[1])*1.15), 
                       border = "grey33", col = "grey33", lwd=3)    
                  par(xpd = FALSE)
                  
            }
            
      }
      
}





############################################################################################################################
############################################################################################################################
############################################################################################################################


keepBSgenomeSequences <- function(genome, seqnames)
{
      stopifnot(all(seqnames %in% seqnames(genome)))
      genome@user_seqnames <- setNames(seqnames, seqnames)
      genome@seqinfo <- genome@seqinfo[seqnames]
      genome
}


############################################################################################################################
############################################################################################################################
############################################################################################################################



scale_min_max <- function(x, lower=0, upper=1){
      y <-(x-quantile(x, lower))/(quantile(x, upper)-quantile(x, lower))
      return(y)
}





############################################################################################################################
############################################################################################################################
############################################################################################################################


writeRangesToFasta <- function(my_genome, my_ranges, my_file_path){
      
      seqlevelsStyle(my_ranges) <- "UCSC"
      seqlevelsStyle(my_genome) <- "UCSC"
      
      my_seqs <- getSeq(my_genome, my_ranges)
      names(my_seqs) <- paste(seqnames(my_ranges), start(my_ranges), end(my_ranges), sep="_")
      
      writeXStringSet(my_seqs, paste0(my_file_path, ".fasta"))
      
      export.bed(my_ranges, paste0(my_file_path, ".bed"))
}




############################################################################################################################
############################################################################################################################
############################################################################################################################



plotHeatScatter <- function(my_data,
                            my_subset,
                            my_title,
                            my_stages,
                            my_stage_names,
                            my_legend,
                            my_colors,
                            k_means,
                            my_limits = c(-1,3))
      {
      
      
      
      my_xdata <- rowMeans(my_data[,grep(my_stages[1], colnames(my_data))])[my_subset]
      my_ydata <- rowMeans(my_data[,grep(my_stages[2], colnames(my_data))])[my_subset]
      
      plot(my_xdata, my_ydata,  
            main = "", pch =19, cex = 0.25,
            col = my_colors[k_means],
            xlab = my_stage_names[1], 
            ylab = my_stage_names[2], 
            xlim = my_limits, ylim = my_limits)
      title(my_title)
      
      legend("topleft", legend = "", title = my_legend, bty="n")
      
      abline(coef = c(0,1)) 
      
}











############################################################################################################################
############################################################################################################################
############################################################################################################################














saveVenn3 <- function(my_gene_list = list(Consitutive = names(my_con_genes), 
                                          H4K16ac = names(K16_clust1_genes_HASdist),
                                          H3K36me3 = names(K36_clust1_genes_HASdist)),
                      my_file_name = "clust1_con",
                      face_xadjust = c(-4,+5,0,0,-1,+1,0,0),
                      face_yadjust = c(0,-10,0,0,0,0,0,0),
                      lab_xadjust = c(0,5,-5),
                      lab_yadjust = c(0,0,0),
                      my_color = rep(my_grey_palette[1],3),
                      Faces = FALSE
                      
){
      
      
      
      
      
      my_overlaps_venn <-  Venn(my_gene_list)
      my_overlaps_venn_plot <- compute.Venn(my_overlaps_venn,doWeights=TRUE)
      
      SetLabels <- VennGetSetLabels(my_overlaps_venn_plot)
      
      
      for(i in seq_along(SetLabels[,"x"] )){
            SetLabels[i,"x"] <- SetLabels[i,"x"] + lab_xadjust[i]
      }
      
      for(i in seq_along(SetLabels[,"y"] )){
            SetLabels[i,"y"] <- SetLabels[i,"y"] + lab_yadjust[i]
      }
      
      
      
      my_overlaps_venn_plot <- VennSetSetLabels(my_overlaps_venn_plot,SetLabels)
      
      SetFaceLabels <- VennGetFaceLabels(my_overlaps_venn_plot)
      
      
      
      for(i in seq_along(SetFaceLabels$x)){
            SetFaceLabels$x[i] <- SetFaceLabels$x[i] + face_xadjust[i]
      }
      
      
      for(i in seq_along(SetFaceLabels$y)){
            SetFaceLabels$y[i] <- SetFaceLabels$y[i] + face_yadjust[i]
      }
      
      
      my_overlaps_venn_plot <- Vennerable:::VennSetFaceLabels(my_overlaps_venn_plot,SetFaceLabels)
      
      
      gp <- VennThemes(my_overlaps_venn_plot)
      
      
      gp$Set$Set1$col <- my_color[1]
      gp$Set$Set2$col <- my_color[2]
      gp$Set$Set3$col <- my_color[3]
      
      gp$SetText$Set1$col <- my_color[1]
      gp$SetText$Set2$col <- my_color[2]
      gp$SetText$Set3$col <- my_color[3]
      
      gp$Set$Set1$lwd <- 6
      gp$Set$Set2$lwd <- 6
      gp$Set$Set3$lwd <- 6
      
      gp$SetText$Set1$fontsize <- 30
      gp$SetText$Set2$fontsize <- 30
      gp$SetText$Set3$fontsize <- 30
      
      for(i in seq_along(gp$FaceText)){
            gp$FaceText[[i]]$fontsize <- 30
      }
      
      png(paste0("../data_folder/venn.",my_file_name,".png"), height = 8, width = 8, units = "in", res = 300)
      
      grid.newpage()
      plot(my_overlaps_venn_plot, show = list(Faces = Faces, Universe=FALSE), gp = gp)
      
      dev.off()    
      
      
      
}





############################################################################################################################
############################################################################################################################
############################################################################################################################


my_line <- function(x,y,...){
      points(x,y,...)
      abline(coef = c(0,1),col="black",lty=1)
}



my_hv_line <- function(x,y,...){
      points(x,y,...)
      abline(h=0,v=0,col="black",lty=2)
}

my_dhv_line <- function(x,y,...){
      points(x,y,...)
      abline(coef = c(0,1),h=0,v=0,col="black",lty=2)
}


add_cor_to_plot <- function(x,y,...){
      text(0,0, labels = paste("r =", round(cor(x,y, method = "spearman"),2)), cex = 1.5)
}

add_cor_to_plot_FC <- function(x,y,...){
  text(1.5,1.5, labels = paste("r =", round(cor(x,y, method = "spearman"),2)), cex = 1.5)
}



############################################################################################################################
############################################################################################################################
############################################################################################################################



callback = function(hc, mat){
      sv = svd(t(mat))$v[,1]
      dend = rev(reorder(as.dendrogram(hc), wts = sv))
      as.hclust(dend)
}




############################################################################################################################
############################################################################################################################
############################################################################################################################
###### Modified function for HelpersforDeSEQ2 package #########

getResults_v2_shrink <- function(dds,
                                 contrast,
                                 result_name = NULL,
                                 lfc_cutoff = 0,
                                 shrink = TRUE,
                                 annotation = "gtf",
                                 anno_id = "gene_id",
                                 anno_symbol = "gene_name")
{
  
  if(is.null(result_name)){
    
    res <- DESeq2::results(dds,
                           contrast = c("Sample", as.character(map(contrast,1)), as.character(map(contrast,2))),
                           lfcThreshold = lfc_cutoff,
                           independentFiltering = FALSE)
    
    if(shrink){
      res <- DESeq2::lfcShrink(dds,  contrast = c("Sample", as.character(map(contrast,1)), as.character(map(contrast,2))), type="ashr", res = res)
    }
    
  } else {
    
    res <- DESeq2::results(dds,
                           name = result_name,
                           lfcThreshold = lfc_cutoff,
                           independentFiltering = FALSE)
    
  }
  
  
  
  anno <- get(annotation)
  
  res[anno_symbol] <- mcols(anno)[anno_symbol][,1][match(rownames(res), mcols(anno)[anno_id][,1])]
  res$chr <- as.character(seqnames(anno))[match(rownames(res), mcols(anno)[anno_id][,1])]
  
  res$padj[is.na(res$padj)] <- 1
  res$log2FoldChange[is.na(res$log2FoldChange)] <- 0
  
  res <- res[order(res$pvalue),]
  
  return(res)
  
}

