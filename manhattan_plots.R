library(ggplot2)

#############
# FUNCTIONS #
#############

# Multiple plot function from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_%28ggplot2%29/
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

add_sig_or_nonsig_column <- function(sig_SNP_table, all_SNP_table){
  all_chrpos <- apply(all_SNP_table[1:2], 1, paste, collapse="_")
  sig_chrpos <- apply(sig_SNP_table, 1, paste, collapse="_")
  all_SNP_table$sig_or_nonsig <- "nonsig"
  all_SNP_table$sig_or_nonsig[all_chrpos%in%sig_chrpos] <- "sig"
  return(all_SNP_table[,c(1,2,6,which(colnames(all_SNP_table)=="sig_or_nonsig"))])
}
  
add_cumulative_pos_column <- function(all_SNP_table){

  chrNum=5
  chrs <- c("2L", "2R", "3L", "3R", "X")
  all_SNP_table$cumulative_pos <- all_SNP_table[,2]
  new_col_location <- (which(colnames(all_SNP_table)=="cumulative_pos"))
  for (i in 1:chrNum){
    tempchr <- chrs[i]
    ndx <- which(all_SNP_table[, 1]==tempchr)
    lstMrk <- max(all_SNP_table[ndx, new_col_location])
    #print(lstMrk)
    if (i < chrNum) ndx2 <- which(all_SNP_table[, 1]==chrs[i+1])
    if (i < chrNum) all_SNP_table[ndx2, new_col_location] <- all_SNP_table[ndx2, new_col_location] + lstMrk
  }
  
  bpMidVec <- vector(length=chrNum)
  
  for (i in 1:chrNum){ndx <- which(all_SNP_table[, 1]==chrs[i])
                      posSub <- all_SNP_table[ndx, new_col_location]
                      bpMidVec[i] <- ((max(posSub) - min(posSub))/2) + min(posSub)
  }
  print(cbind(chrs, bpMidVec))
  return(all_SNP_table)
}

colScale <- scale_fill_manual(name = "type",values = exon_utr_colours)

man_plot <- function(i, SNP_table, this_fet_threshold){
  marks <- c(11759252,36148367,62847735,
             92915027,120681945)
  chrs <- c("2L", "2R", "3L", "3R", "X")
  sig_subset <- subset(SNP_table, SNP_table$sig_or_nonsig=="sig")
  SNP_table$chr <- ordered(as.factor(SNP_table$chr), levels=c("2L", "2R", "3L", "3R", "X"))
  p <- ggplot() +
    geom_point(data=SNP_table, aes(x=cumulative_pos, y=-log10(pval), 
                                   colour=factor(chr)), 
                                   alpha=1/3, size=0.2) + 
    scale_size(range=c(0.3, 0.3)) +
    scale_color_manual(values=c('grey', 'black', 'grey', 'black', 'grey')) +
    scale_x_continuous(labels=chrs, breaks=marks) +
    scale_y_continuous(limits=c(0, 40)) +
    geom_point(data=sig_subset, aes(x=cumulative_pos, y=-log10(pval)), colour="red", alpha=0.5, size=0.2) +
    #scale_size(range=c(0.3, 0.3)) +
    #theme_bw(base_size=15)  +
    theme(legend.position='none')  +
    theme(panel.background=element_blank(), 
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.border=element_blank()) + 
    theme(axis.line=element_line(color="black"), plot.margin=unit(c(0.25,0,0.25,0.25),units="cm")) + 
    geom_hline(y=-log10(this_fet_threshold), linetype=2, col='grey', lwd=0.2) +
    xlab('') + 
    ylab(expression(paste('-log10(', italic('P'),')', sep=""), cex=4)) +
    annotate("text", x=-100, y=30, label=LETTERS[which(Sample_code==i)], cex=4) 
}


###########
# MAIN #
########

setwd("~/Documents/Drosophila\ Selection\ Experiment/snp_and_gene_lists/")

sig_dir <- "~/Documents/Drosophila\ Selection\ Experiment/snp_and_gene_lists/"
sync_dir <- "~/Documents/Drosophila\ Selection\ Experiment/pileup_and_sync_files/"

Sample_code <- c("C1", "C2", "C3", "C4", "C5",
                 "D1", "D2", "D3", "D4", "D5")

fet_threshold <- c(3.46E-13, 3.55E-13,  3.56E-13, 3.58E-13, 3.58E-13,
                   3.55E-13, 3.69E-13, 3.46E-13, 3.35E-13, 3.57E-13)
names(fet_threshold) <- Sample_code


for(i in Sample_code){
  temp_sample <- i
  temp_fet_threshold <- fet_threshold[i]
  temp_all_SNPs_name <- paste(temp_sample, "_all_SNPs", sep="")
  temp_sig_SNPs_name <- paste(temp_sample, "_sig_SNPs", sep="")
  temp_table_name <- paste(temp_sample, "_SNP_table", sep="")
  temp_all_SNPs <- read.table(paste(sync_dir, "MB_", temp_sample, 
                                    "_noN_reduced_150818.mpileup.sync", sep=""), 
                              sep="\t", header=FALSE, stringsAsFactors=FALSE)
  temp_sig_SNPs <- read.table(paste(sig_dir, temp_sig_SNPs_name, ".txt", 
                                              sep=""), sep="\t", 
                                        header=FALSE, stringsAsFactors=FALSE)
  print(paste("Data read in for", i))
  temp_table <- add_sig_or_nonsig_column(temp_sig_SNPs, temp_all_SNPs)
  print(head(temp_table))
  rm(temp_all_SNPs)
  rm(temp_sig_SNPs)
  colnames(temp_table)[1:3] <- c("chr", "pos", "pval")
  temp_table <- temp_table[-which(temp_table$chr=="dmel_mitochondrion_genome"),]  
  temp_table <- add_cumulative_pos_column(temp_table)
  print(paste("Formatting complete for", i))
  assign(temp_table_name, temp_table)
}



for(i in Sample_code){
  temp_sample <- i
  temp_fet_threshold <- fet_threshold[i]
  temp_table <- get(paste(temp_sample, "_SNP_table", sep=""))
  temp_plot_name <- paste(temp_sample, "_manhattan", sep="")
  temp_plot <- man_plot(i=i, SNP_table=temp_table, this_fet_threshold=temp_fet_threshold)
  assign(temp_plot_name, temp_plot)
  print(paste("Plot created for", i))
}
  
jpeg(file="Manhattan Plots for Selection Experiment.jpg", width=21, height=15,
     units="cm", res=1200)
multiplot(plotlist=list(C1_manhattan, C2_manhattan, C3_manhattan,
                        C4_manhattan, C5_manhattan, D1_manhattan,
                        D2_manhattan, D3_manhattan, D4_manhattan, D5_manhattan),
          layout=matrix(1:10, nrow=5, byrow=FALSE))
dev.off()



