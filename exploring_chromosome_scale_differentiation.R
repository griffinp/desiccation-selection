# Plotting differentiation across chromosomes in more detail

#############
# FUNCTIONS #
#############

add_sig_or_nonsig_column <- function(sig_SNP_table, all_SNP_table){
  all_chrpos <- apply(all_SNP_table[1:2], 1, paste, collapse="_")
  sig_chrpos <- apply(sig_SNP_table, 1, paste, collapse="_")
  all_SNP_table$sig_or_nonsig <- "nonsig"
  all_SNP_table$sig_or_nonsig[all_chrpos%in%sig_chrpos] <- "sig"
  return(all_SNP_table[,c(1,2,6,which(colnames(all_SNP_table)=="sig_or_nonsig"))])
}

add_cumulative_pos_column <- function(all_SNP_table, include_X=TRUE){
  if (include_X==TRUE){ chrNum <- 5 }
  else {chrNum <- 4 }
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
  bpMidVec[i] <- ((max(posSub, na.rm=TRUE) - min(posSub, na.rm=TRUE))/2) + min(posSub, na.rm=TRUE)
  }
  print(cbind(chrs, bpMidVec))
  return(all_SNP_table)
}



###########
# MAIN #
########


################################

setwd("~/Documents/Documents_from_Hoffmann_Computer/Documents/Drosophila\ Selection\ Experiment/snp_and_gene_lists/")

sig_dir <- "~/Documents/Documents_from_Hoffmann_Computer/Documents/Drosophila\ Selection\ Experiment/snp_and_gene_lists/"
sync_dir <- "~/Documents/Documents_from_Hoffmann_Computer/Documents/Drosophila\ Selection\ Experiment/pileup_and_sync_files/"
freq_dir <- "~/Documents/Documents_from_Hoffmann_Computer/Documents/Drosophila\ Selection\ Experiment/allele_frequency_difference_testing/"

Sample_code <- c("C1", "C2", "C3", "C4", "C5",
                 "D1", "D2", "D3", "D4", "D5")

fet_threshold <- c(3.46E-13, 3.55E-13,  3.56E-13, 3.58E-13, 3.58E-13,
                   3.55E-13, 3.69E-13, 3.46E-13, 3.35E-13, 3.57E-13)
names(fet_threshold) <- Sample_code



# for(i in Sample_code){
#   # Don't need to run this every time (it's quite slow...)
#   # See below for importing files instead
#   temp_sample <- i
#   temp_fet_threshold <- fet_threshold[i]
#   temp_all_SNPs_name <- paste(temp_sample, "_all_SNPs", sep="")
#   temp_sig_SNPs_name <- paste(temp_sample, "_sig_SNPs", sep="")
#   temp_table_name <- paste(temp_sample, "_SNP_table", sep="")
#   temp_all_SNPs <- read.table(paste(sync_dir, "MB_", temp_sample, 
#                                     "_noN_reduced_150818.mpileup.sync", sep=""), 
#                               sep="\t", header=FALSE, stringsAsFactors=FALSE)
#   temp_sig_SNPs <- read.table(paste(sig_dir, temp_sig_SNPs_name, ".txt", 
#                                     sep=""), sep="\t", 
#                               header=FALSE, stringsAsFactors=FALSE)
#   print(paste("Data read in for", i))
#   temp_table <- add_sig_or_nonsig_column(temp_sig_SNPs, temp_all_SNPs)
#   print(head(temp_table))
#   rm(temp_all_SNPs)
#   rm(temp_sig_SNPs)
#   colnames(temp_table)[1:3] <- c("chr", "pos", "pval")
#   temp_table <- temp_table[-which(temp_table$chr=="dmel_mitochondrion_genome"),]  
#   temp_table <- add_cumulative_pos_column(temp_table)
#   for(tempchr in c("2L", "2R", "3L", "3R", "X")){
#     tempsubset <- temp_table[which(temp_table$chr==tempchr),]
#     print(paste("Merging freq info for", temp_sample, tempchr))
#     freqfile <- read.table(paste(freq_dir, "MB_", tempchr, "_all_startfreq.txt", sep=""),
#                            sep="\t", header=TRUE, stringsAsFactors=FALSE)
#     names(freqfile) <- c("chr", "pos", "start_freq")
#     tempsub_inc_freq <- merge(tempsubset, freqfile, by=c("chr", "pos"))
#     assign(paste("subset_", tempchr, sep=""), tempsub_inc_freq)
#     rm(freqfile)
#   }
#   table_output <- rbind(subset_2L, subset_2R, subset_3L, subset_3R, subset_X)
#   print(paste("Formatting complete for", i))
#   assign(temp_table_name, table_output)
# }
# 
# for(each_table in paste(Sample_code, "_SNP_table", sep="")){
#   temptab <- get(each_table)
#   tempsavename <- paste(sig_dir, each_table, "_with_startfreq.txt", sep="")
#   write.table(temptab, file=tempsavename, sep="\t", col.names=TRUE,
#               row.names=FALSE, quote=FALSE)
# }

for(sample in Sample_code){
  tempsavename <- paste(sig_dir, sample, "_SNP_table_with_startfreq.txt", sep="")
  temptabname <- paste(sample, "_SNP_table", sep="")
  temptab <- read.table(file=tempsavename, header=TRUE, stringsAsFactors=FALSE,
                        sep="\t")
  assign(temptabname, temptab)
}

binsize <- 100000
#binsize <- 1000000
min_SNPs_in_bin <- 10
quantile_threshold <- 0.8
#sampletab <- get(paste(sample, "_binned_SNP_table", sep=""))



for(i in Sample_code){
  temp_table <- get(paste(i, "_SNP_table", sep=""))
  bin_results <- data.frame(chr=character(),
                            binstart=numeric(),
                            binend=numeric(),
                            binmid=numeric(),
                            meanp=numeric(),
                            number_sig_snps=numeric(),
                            total_snps=numeric(), 
                            mean_start_freq=numeric(), 
                            number_startfreq_below_0.1=numeric(),
                            number_startfreq_below_0.05=numeric(),
                            startfreq_ratio=numeric(),
                            stringsAsFactors = FALSE)
  results_name <- paste(i, "_binned_SNP_table", sep="")
  #results_name <- paste(i, "_broad_binned_SNP_table", sep="")
  existing_bin_no <- 0
  
  for(tempchr in c("2L", "2R", "3L", "3R", "X")){
    #for each chr, calculate bin starts
    #for each bin in a chr, 
    #test whether at least a min. no. SNPs in bin
    #calculate av. p-val and number sig
    tempsubset <- subset(temp_table, temp_table$chr==tempchr)
    minbase <- min(tempsubset$pos)
    maxbase <- max(tempsubset$pos)
    binstarts <- seq(minbase, maxbase, by=binsize)
    binends <- seq(minbase+binsize-1, maxbase+binsize-1, by=binsize)
    bincoords <- data.frame(binstarts, binends)
    
    for(k in 1:nrow(bincoords)){
      bin_number <- k+existing_bin_no
      if(k%%(nrow(bincoords)/10)==0){print(paste("Processing bin #", k, "for", tempchr, "in", i))}
      binstart <- bincoords[k, "binstarts"]
      binend <- bincoords[k, "binends"]
      binmid <- binstart+(binsize/2)
      bincontents <- subset(tempsubset, tempsubset$pos>=binstart&tempsubset$pos<binend)
      meanp <- mean(bincontents$pval)
      meanstartfreq <- mean(bincontents$start_freq)
      number_below_0.1 <- length(which(bincontents$start_freq<0.1))
      number_below_0.05 <- length(which(bincontents$start_freq<0.05))
      number_sig <- length(which(bincontents$sig_or_nonsig=="sig"))
      # calculate the ratio of mean start freq in sig loci / mean start freq in nonsig loci
      meanstart_sig <- mean(bincontents[which(bincontents$sig_or_nonsig=="sig"),]$start_freq)
      meanstart_nonsig <- mean(bincontents[which(bincontents$sig_or_nonsig=="nonsig"),]$start_freq)
      startfreq_ratio <- meanstart_sig/meanstart_nonsig
      total_number <- nrow(bincontents)
      bin_results[bin_number,1] <- as.character(tempchr)
      if(total_number >= min_SNPs_in_bin){
        bin_results[bin_number,2:11] <- c(binstart, binend, binmid, meanp, number_sig, total_number, meanstartfreq, number_below_0.1, 
                                          number_below_0.05, startfreq_ratio)
      } else {
        bin_results[bin_number,2:11] <- c(binstart, binend, binmid, NA, NA, total_number, NA, NA, NA, NA)
      }
    }
    existing_bin_no <- existing_bin_no+nrow(bincoords)
  }
  prop_sig <- bin_results$number_sig_snps/bin_results$total_snps
  prop_sig_nonzero <- prop_sig[which(prop_sig>0)]
  quant <- quantile(prop_sig_nonzero, probs=quantile_threshold)
  high_sig <- prop_sig>quant
  bin_results$high_sig <- high_sig
  assign(results_name, bin_results)
}


# Loop #2!
no_bins_to_connect <- 10
connect_width <- no_bins_to_connect-1

for(i in Sample_code){
  temp_sample_snp_table <- get(paste(i, "_SNP_table", sep=""))
  temp_bin_sample_table <- get(paste(i, "_binned_SNP_table", sep=""))
  results_name <- paste(i, "_region_binned_SNP_table", sep="")
  existing_bin_no <- 0
  
  for(tempchr in c("2L", "2R", "3L", "3R", "X")){
    print(paste("Working on", tempchr, "for", i))
    regioncoords_name <- paste(i, "_", tempchr, "_region_coords", sep="")
    temp_snp_table <- temp_sample_snp_table[which(temp_sample_snp_table$chr==tempchr),]
    temp_bin_table <- temp_bin_sample_table[which(temp_bin_sample_table$chr==tempchr),]
    temp_bin_results <- temp_bin_table$high_sig
    region_no <- rep(NA, times=length(temp_bin_results))
    region_sig <- rep(NA, times=length(temp_bin_results))

    for(m in 1:length(temp_bin_results)){
      temp_result <- temp_bin_results[m]
      region_sig[m] <- temp_result
      #print(region_sig[m])
      if(m>connect_width){
        if((region_sig[m]|is.na(region_sig[m]))&TRUE%in%region_sig[(m-connect_width):(m-1)]){
          lasttrue <- which(region_sig[(m-connect_width):(m-1)])
          region_sig[(m-no_bins_to_connect+lasttrue[1]):m] <- TRUE
        }
      }
    }
    region_no[1] <- 1
    for(n in 2:length(region_sig)){
      current <- region_sig[n]
      previous <- region_sig[n-1]
      if(current==previous|is.na(current)|is.na(previous)){
        region_no[n] <- region_no[n-1]
      } else{
        region_no[n] <- region_no[n-1]+1
      }
    }
    temp_bin_table$region_no <- region_no
    print(paste(max(region_no), "regions identified"))
    temp_bin_table$region_sig <- region_sig
    region_details <- as.data.frame(xtabs(region_sig~region_no))
    region_coords <- data.frame(chr=character(),
                                region_number=numeric(),
                                region_start=numeric(),
                                region_end=numeric(),
                                region_sig_status=logical(),
                                prop_below_0.1=numeric(),
                                prop_sig_below_0.1=numeric(), stringsAsFactors=FALSE)
    for(o in region_details$region_no){
      tempsub <- temp_bin_table[which(as.character(temp_bin_table$region_no)==o),]
      tempstart <- min(tempsub$binstart)
      tempend <- max(tempsub$binend)
      tempsub2 <- temp_snp_table[which(temp_snp_table$pos>=tempstart&temp_snp_table$pos<tempend),]
      prop_below_0.1 <- length(which(tempsub2$start_freq<0.1&is.na(tempsub2$start_freq)==FALSE))/nrow(tempsub2)
      if(length(which(tempsub2$sig_or_nonsig=="sig"))==0){
        prop_sig_below_0.1<-0
      } else {
        prop_sig_below_0.1 <- length(which(tempsub2$start_freq<0.1&tempsub2$sig_or_nonsig=="sig"&is.na(tempsub2$start_freq)==FALSE))/length(which(tempsub2$sig_or_nonsig=="sig"))
      }
      tempstatus <- tempsub[1, "region_sig"]
      region_coords[o,] <- c(as.character(tempchr), o, tempstart, tempend, tempstatus, prop_below_0.1, prop_sig_below_0.1)
    }
    assign(regioncoords_name, region_coords)
    
  }
  tables_to_combine <- paste(i, "_", c("2L", "2R", "3L", "3R", "X"), "_region_coords", sep="")
  allchrs <- rbind(get(tables_to_combine[1]), get(tables_to_combine[2]), get(tables_to_combine[3]),
                   get(tables_to_combine[4]), get(tables_to_combine[5]), stringsAsFactors=FALSE)
  # Modelling at the whole-genome level (too few regions to do it per chromosome)
  glm1 <- glm(as.factor(region_sig_status)~as.numeric(prop_below_0.1), family=binomial, data=allchrs)
  glm2 <- glm(as.factor(region_sig_status)~as.numeric(prop_below_0.1)+as.numeric(prop_sig_below_0.1), family=binomial, data=allchrs)
  anova_test <- anova(glm1, glm2, test="Chisq")
  print(anova_test)
  if(nrow(allchrs)<5){print("Too few regions to test")}
  if(nrow(allchrs)>=5&anova_test[5][2,1] < 0.05){ print(summary(glm2))}
}




pdf(file="window_freq_plotting_3.pdf", width=25, height=25)
par(mfcol=c(10,5))
for(tempchr in c("2L", "2R", "3L", "3R", "X")){
  for(sample in Sample_code){
    tab <- get(paste(sample, "_binned_SNP_table", sep=""))
    subtab <- tab[which(tab$chr==tempchr),]
    print(paste(nrow(subtab), "windows examined"))
    prop_lowstart <- subtab$number_startfreq_below_0.1/subtab$total_snps
    prop_sig <- subtab$number_sig_snps/subtab$total_snps
    #Linear models 
    model1 <- lm(prop_sig~prop_lowstart)
    rsq_mod1 <- round(summary(model1)$r.squared, 3)
    pval_mod1 <- signif(summary(model1)$coefficients[2,4], 2)
    subtab2 <- subtab[which(subtab$number_sig_snps>0),]
    print(paste(nrow(subtab2), "windows excluding non-sig"))
    prop_lowstart2 <- subtab2$number_startfreq_below_0.1/subtab2$total_snps
    prop_sig2 <- subtab2$number_sig_snps/subtab2$total_snps
    model2 <- lm(prop_sig2~prop_lowstart2)
    rsq_mod2 <- round(summary(model2)$r.squared, 3)
    pval_mod2 <- signif(summary(model2)$coefficients[2,4], 2)
    #Predicting for plotting
    predvalues <- seq(0.1, 0.45, by=0.01)
    mod1pred <- predict(model1, newdata=list(prop_lowstart=predvalues))
    mod2pred <- predict(model2, newdata=list(prop_lowstart2=predvalues))
    #plotting
    plot(x=prop_lowstart, 
         y=prop_sig, type="p", 
         pch=19, cex=0.85, col=rgb(100, 100, 100, 100, maxColorValue = 255),
         xlab="prop. snps with startfreq<0.1", ylab="prop. sig snps in window",
         xlim=c(0, 0.7), ylim=c(-0.004, max(prop_sig, na.rm=TRUE)+0.01), main=paste(sample, "chr", tempchr))
    lines(x=predvalues, y=mod1pred, col="pink", lty=2)
    if(pval_mod1<0.05){
      text(x=max(predvalues)+0.05, y=mod1pred[length(mod1pred)], col="pink", 
           labels=paste("R2 =", rsq_mod1, "P =", pval_mod1), pos=4)
    }
    lines(x=predvalues, y=mod2pred, col="red")
    if(pval_mod2<0.05){
      text(x=max(predvalues)+0.05, y=mod2pred[length(mod2pred)], col="red", 
           labels=paste("R2 =", rsq_mod2, "P =", pval_mod2), pos=4)
    }
  }
}
dev.off()

sig_freq_effect <- matrix(nrow=5, ncol=10, data=rep("nonsig", times=50))
colnames(sig_freq_effect)<-Sample_code
rownames(sig_freq_effect)=c("2L", "2R", "3L", "3R", "X")
sig_freq_effect["3R", "C1"] <- "sig"
sig_freq_effect["3L", "C2"] <- "sig"
sig_freq_effect["2R", "D1"] <- "sig"
sig_freq_effect["3L", "D2"] <- "sig"
sig_freq_effect["X", "D3"] <- "sig"
sig_freq_effect["3L", "D3"] <- "sig"
sig_freq_effect["2L", "D4"] <- "sig"
sig_freq_effect["3R", "D4"] <- "sig"
sig_freq_effect["X", "D4"] <- "sig"
sig_freq_effect["3R", "D5"] <- "sig"
sig_freq_effect["2R", "D5"] <- "sig"

pdf(file="window_sig_plotting.pdf", width=25, height=25)
par(mfcol=c(10,5))
for(tempchr in c("2L", "2R", "3L", "3R", "X")){
  for(sample in Sample_code){
    tab <- get(paste(sample, "_binned_SNP_table", sep=""))
    subtab <- tab[which(tab$chr==tempchr),]
    plot(subtab$binmid, subtab$number_sig_snps/subtab$total_snps, type="l", 
         xlab="Position (bp)", ylab="",
         ylim=c(-0.02, 0.16), main=paste(sample, "chr", tempchr),
         col=rgb(0, 0, 0, 100, maxColorValue = 255))
    mtext(side=2, line=2, text="prop. SNPs sig. diff.", cex=0.5)
    polygon_tab <- get(paste(sample, "_", tempchr, "_region_coords", sep=""))
    sub_polygon <- polygon_tab[which(polygon_tab$region_sig_status==TRUE),]
    segments(x0=as.numeric(sub_polygon$region_start), x1=as.numeric(sub_polygon$region_end), 
             y0=rep(-0.01, length=nrow(sub_polygon)), y1=rep(-0.01, length=nrow(sub_polygon)),
             lwd=2)
    tab2 <- get(paste(sample, "_broad_binned_SNP_table", sep=""))
    subtab2 <- tab2[which(tab2$chr==tempchr),]
    lines(subtab2$binmid, subtab2$number_sig_snps/subtab2$total_snps, type="l", 
          col="black")
    if(sig_freq_effect[tempchr, sample]=="sig"){
      text(x=100, y=0.15, labels="*", cex=3, col="green")
    }
    par(new=TRUE)
    plot(subtab$binmid, subtab$number_startfreq_below_0.1/subtab$total_snps,
          col=rgb(0, 0, 255, 100, maxColorValue = 255), 
         ylim=c(0, 0.4), axes=FALSE, type="l", ylab="", xlab="")
    lines(subtab2$binmid, subtab2$number_startfreq_below_0.1/subtab2$total_snps,
          col=rgb(0, 0, 255, 255, maxColorValue=255))
    axis(side = 4, col=rgb(0, 0, 255, 255, maxColorValue=255), 
         col.axis=rgb(0, 0, 255, 255, maxColorValue=255))
    mtext(side = 4, line = 2, "prop. SNPs with startfreq<0.1", 
          col=rgb(0, 0, 255, 100, maxColorValue=255), cex = 0.5)
  }
}
dev.off()

# model highly-differentiated windows

for(tempchr in c("2L", "2R", "3L", "3R", "X")){
  for(sample in Sample_code){
    tab <- get(paste(sample, "_binned_SNP_table", sep=""))
    subtab <- tab[which(tab$chr==tempchr),]
    prop_lowstart <- subtab$number_startfreq_below_0.1/subtab$total_snps
    prop_sig <- subtab$number_sig_snps/subtab$total_snps
    #Linear models 
    model1 <- lm(prop_sig~prop_lowstart)
    rsq_mod1 <- round(summary(model1)$r.squared, 3)
    pval_mod1 <- signif(summary(model1)$coefficients[2,4], 2)
    subtab2 <- subtab[which(subtab$number_sig_snps>0),]
    prop_lowstart2 <- subtab2$number_startfreq_below_0.1/subtab2$total_snps
    prop_sig2 <- subtab2$number_sig_snps/subtab2$total_snps
    model2 <- lm(prop_sig2~prop_lowstart2)
    rsq_mod2 <- round(summary(model2)$r.squared, 3)
    pval_mod2 <- signif(summary(model2)$coefficients[2,4], 2)
    #Predicting for plotting
    predvalues <- seq(0.1, 0.45, by=0.01)
    mod1pred <- predict(model1, newdata=list(prop_lowstart=predvalues))
    mod2pred <- predict(model2, newdata=list(prop_lowstart2=predvalues))
    #plotting
    plot(x=prop_lowstart, 
         y=prop_sig, type="p", 
         pch=19, cex=0.85, col=rgb(100, 100, 100, 100, maxColorValue = 255),
         xlab="prop. snps with startfreq<0.1", ylab="prop. sig snps in window",
         xlim=c(0, 0.7), ylim=c(0, 0.16), main=paste(sample, "chr", tempchr))
    lines(x=predvalues, y=mod1pred, col="pink", lty=2)
    if(pval_mod1<0.05){
      text(x=max(predvalues)+0.05, y=mod1pred[length(mod1pred)], col="pink", 
           labels=paste("R2 =", rsq_mod1, "P =", pval_mod1), pos=4)
    }
    lines(x=predvalues, y=mod2pred, col="red")
    if(pval_mod2<0.05){
      text(x=max(predvalues)+0.05, y=mod2pred[length(mod2pred)], col="red", 
           labels=paste("R2 =", rsq_mod2, "P =", pval_mod2), pos=4)
    }
  }
}
dev.off()






# pdf(file="window_freq_plotting_4.pdf", width=25, height=50)
# par(mfcol=c(10,5))
# for(tempchr in c("2L", "2R", "3L", "3R", "X")){
#   for(sample in Sample_code){
#     tab <- get(paste(sample, "_binned_SNP_table", sep=""))
#     subtab <- tab[which(tab$chr==tempchr),]
#     prop_lowstart <- subtab$number_startfreq_below_0.05/subtab$total_snps
#     prop_sig <- subtab$number_sig_snps/subtab$total_snps
#     #Linear models 
#     model1 <- lm(prop_sig~prop_lowstart)
#     rsq_mod1 <- round(summary(model1)$r.squared, 3)
#     pval_mod1 <- round(summary(model1)$coefficients[2,4], 3)
#     subtab2 <- subtab[which(subtab$number_sig_snps>0),]
#     prop_lowstart2 <- subtab2$number_startfreq_below_0.05/subtab2$total_snps
#     prop_sig2 <- subtab2$number_sig_snps/subtab2$total_snps
#     model2 <- lm(prop_sig2~prop_lowstart2)
#     rsq_mod2 <- round(summary(model2)$r.squared, 3)
#     pval_mod2 <- round(summary(model2)$coefficients[2,4], 3)
#     #Predicting for plotting
#     predvalues <- seq(0.1, 0.45, by=0.01)
#     mod1pred <- predict(model1, newdata=list(prop_lowstart=predvalues))
#     mod2pred <- predict(model2, newdata=list(prop_lowstart2=predvalues))
#     #plotting
#     plot(x=prop_lowstart, 
#          y=prop_sig, type="p", 
#          pch=19, cex=0.85, col=rgb(100, 100, 100, 100, maxColorValue = 255),
#          xlab=paste(sample, tempchr, "prop. snps with startfreq<0.05"), ylab="prop. sig snps in window",
#          xlim=c(0, 0.7), ylim=c(0, 0.13))
#     lines(x=predvalues, y=mod1pred, col="green", lty=2)
#     if(pval_mod1<0.05){
#       text(x=max(predvalues)+0.05, y=mod1pred[length(mod1pred)], col="green", 
#            labels=paste("R2 =", rsq_mod1, "P =", pval_mod1), pos=4)
#     }
#     lines(x=predvalues, y=mod2pred, col="red")
#     if(pval_mod2<0.05){
#       text(x=max(predvalues)+0.05, y=mod2pred[length(mod2pred)], col="red", 
#            labels=paste("R2 =", rsq_mod2, "P =", pval_mod2), pos=4)
#     }
#   }
# }
# dev.off()


