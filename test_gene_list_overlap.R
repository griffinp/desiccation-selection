
#############
# FUNCTIONS #
#############

get_gene_lengths <- function(gene_list, gene_list_for_conversion){
  gene_list <- as.data.frame(gene_list)
  colnames(gene_list) <- "gene_names"
  merging <- merge(x=gene_list, y=gene_list_for_conversion, 
                   by.x="gene_names", by.y="name",
                   all.x=TRUE, all.y=FALSE, sort=FALSE)
  gene_lengths <- as.numeric(merging$total_length)
  return(gene_lengths)
}

test_if_gene_in_list <- function(gene_list, gene){
  if(gene%in%gene_list){
    return(TRUE)
  } else {return(FALSE)}
}

########
# MAIN #
########

Sample_code <- c("C1", "C2", "C3", "C4", "C5",
                 "D1", "D2", "D3", "D4", "D5")

# Get object containing simulated gene lists and overlap lists
# (this is the version saved *after* parsing 'exon' gene names)

exons_parsed <- readRDS(file="/Users/pgriffin/Documents/Drosophila Selection Experiment/gene_list_overlap_testing/Gene_lists_exons_parsed_including_overlap_from_resampled_SNPs.rds")

# Extract overlap genes called in each combination of D replicates

setwd("/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/")

# Obtain full list of sig C genes

C_sig_genes <- read.csv('Putative_lab_adaptation_gene_list.txt', header=FALSE,
                        stringsAsFactors=FALSE)[,1]

# Now get real gene lists (again, before removing C genes, and after exon parsing)
for(i in Sample_code[1:10]){
  if(i %in% c('D1', 'D2', 'D3', 'D4', 'D5')){
    temp_input_name <- paste('/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/',
                             i, '_sig_withC_gene_list.txt', sep="")
  } else{
    temp_input_name <- paste('/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/',
                             i, '_sig_gene_list.txt', sep="")
  }
  temp_file <- read.csv(temp_input_name, header=FALSE, stringsAsFactors=FALSE)[,1]
  temp_object_name <- paste(i, '_sig_gene_names', sep="")
  assign(temp_object_name, temp_file)
}

D1_noC_gene_names <- setdiff(D1_sig_gene_names, C_sig_genes)
D2_noC_gene_names <- setdiff(D2_sig_gene_names, C_sig_genes)
D3_noC_gene_names <- setdiff(D3_sig_gene_names, C_sig_genes)
D4_noC_gene_names <- setdiff(D4_sig_gene_names, C_sig_genes)
D5_noC_gene_names <- setdiff(D5_sig_gene_names, C_sig_genes)

D1_D2<-intersect(D1_noC_gene_names, D2_noC_gene_names)
D1_D3<-intersect(D1_noC_gene_names, D3_noC_gene_names)
D1_D4<-intersect(D1_noC_gene_names, D4_noC_gene_names)
D1_D5<-intersect(D1_noC_gene_names, D5_noC_gene_names)
D2_D3<-intersect(D2_noC_gene_names, D3_noC_gene_names)
D2_D4<-intersect(D2_noC_gene_names, D4_noC_gene_names)
D2_D5<-intersect(D2_noC_gene_names, D5_noC_gene_names)
D3_D4<-intersect(D3_noC_gene_names, D4_noC_gene_names)
D3_D5<-intersect(D3_noC_gene_names, D5_noC_gene_names)
D4_D5<-intersect(D4_noC_gene_names, D5_noC_gene_names)
D1_D2_D3<-intersect(intersect(D1_noC_gene_names, D2_noC_gene_names), D3_noC_gene_names)
D1_D2_D4<-intersect(intersect(D1_noC_gene_names, D2_noC_gene_names), D4_noC_gene_names)
D1_D2_D5<-intersect(intersect(D1_noC_gene_names, D2_noC_gene_names), D5_noC_gene_names)
D1_D3_D4<-intersect(intersect(D1_noC_gene_names, D3_noC_gene_names), D4_noC_gene_names)
D1_D3_D5<-intersect(intersect(D1_noC_gene_names, D3_noC_gene_names), D5_noC_gene_names)
D1_D4_D5<-intersect(intersect(D1_noC_gene_names, D4_noC_gene_names), D5_noC_gene_names)
D2_D3_D4<-intersect(intersect(D2_noC_gene_names, D3_noC_gene_names), D4_noC_gene_names)
D2_D3_D5<-intersect(intersect(D2_noC_gene_names, D3_noC_gene_names), D5_noC_gene_names)
D2_D4_D5<-intersect(intersect(D2_noC_gene_names, D4_noC_gene_names), D5_noC_gene_names)
D3_D4_D5<-intersect(intersect(D3_noC_gene_names, D4_noC_gene_names), D5_noC_gene_names)
D1_D2_D3_D4<-intersect(intersect(D1_noC_gene_names, D2_noC_gene_names), 
                       intersect(D3_noC_gene_names, D4_noC_gene_names))
D1_D2_D3_D5<-intersect(intersect(D1_noC_gene_names, D2_noC_gene_names),
                       intersect(D3_noC_gene_names, D5_noC_gene_names))
D1_D2_D4_D5<-intersect(intersect(D1_noC_gene_names, D2_noC_gene_names),
                       intersect(D4_noC_gene_names, D5_noC_gene_names))
D1_D3_D4_D5<-intersect(intersect(D1_noC_gene_names, D3_noC_gene_names),
                       intersect(D4_noC_gene_names, D5_noC_gene_names))
D2_D3_D4_D5<-intersect(intersect(D2_noC_gene_names, D3_noC_gene_names),
                       intersect(D4_noC_gene_names, D5_noC_gene_names))
D1_D2_D3_D4_D5<-intersect(intersect(D1_noC_gene_names, D2_noC_gene_names), 
                          intersect(intersect(D3_noC_gene_names, D4_noC_gene_names), D5_noC_gene_names))


###### FIGURE S7 ######

setwd("/Users/pgriffin/Documents/Drosophila Selection Experiment/gene_list_overlap_testing/")

pdf("Histograms of simulated gene number overlap, resampling with position.pdf", width=28, height=16)
par(mfcol=c(4,7))

combn2 <- combn(D_names, 2)
combn2a <- apply(combn2, 2, paste, collapse="_")
combn3 <- combn(D_names, 3)
combn3a <- apply(combn3, 2, paste, collapse="_")
combn4 <- combn(D_names, 4)
combn4a <- apply(combn4, 2, paste, collapse="_")
combn5 <- "D1_D2_D3_D4_D5"
allcombns <- c(combn2a, combn3a, combn4a, combn5)
for(i in 1:length(allcombns)){
  temp_combn <- allcombns[i]
  real_overlap_length <- length(get(temp_combn))
  temp_all_list <- lapply(exons_parsed, "[[", i+16)
  temp_dist <- sapply(temp_all_list, length)
  
  multiplier <- hist(temp_dist, plot=FALSE)$counts / hist(temp_dist, plot=FALSE)$density
  mydensity <- density(temp_dist)
  mydensity$y <- mydensity$y * multiplier[1] 
  
  hist(temp_dist, main="", ylim=c(0, 380), xlim=c(0, max(temp_dist)*1.5),
       xlab=paste("No. genes in", temp_combn, "overlap", sep=" "))
  lines(mydensity)
  lines(x=c(real_overlap_length, real_overlap_length),
        y=c(0, 225), col="red")
  myx <- seq(min(temp_dist), max(temp_dist), length = 100)
  normal <- dnorm(x=myx, mean = mean(temp_dist), sd = sd(temp_dist))
  lines(myx, normal * multiplier[1], col = "blue", lwd = 1)
  normless <- pnorm(q = real_overlap_length, mean=mean(temp_dist), sd=sd(temp_dist))
  #pval <- (1-normless)*2
  text(x=real_overlap_length, y=275, labels=signif(normless, 3))
  text(x=2, y=300, labels=LETTERS[i], cex=2)
}
dev.off() 

#######################
# Investigating gene  #
# length and its      #
# relationship with   #
# frequency of        #
# sampling            #
#######################


gene_list_for_conversion <- read.table(file="All_gene_list_for_name_conversion.txt",
                                        sep="\t", stringsAsFactors=FALSE, header=TRUE)



# test_l <- get_gene_lengths(rearrange_iterations[[1]][[1]], gene_list_for_conversion)
# test_l2 <- get_gene_lengths(D1_gene_names_exons_parsed, gene_list_for_conversion)
# 
# gene_length_breaks <- seq(0, 401000, by=1000)
# 
# genes_in_temp_list <- gene_list_for_conversion[gene_list_for_conversion$name%in%D1_gene_names_exons_parsed,]
# plot(genes_in_temp_list$total_length, y=rep(1, times=nrow(genes_in_temp_list)))

#find count of resampled lists that each gene appears in


presence_count <- matrix(NA, nrow=nrow(gene_list_for_conversion), 
                         ncol=68)
for(i in 1:nrow(gene_list_for_conversion)){
  temp_gene <- gene_list_for_conversion$name[i]
  if(i%%100==0){
    print(paste("Processing gene #", i, sep=" "))
  }
  presabs_matrix <- matrix(NA, nrow=nIter, ncol=68)
  for(j in 1:68){
    temp_list <- lapply(exons_parsed, "[[", j)
    temp_gene_present <- sapply(X=temp_list, FUN=test_if_gene_in_list, gene=temp_gene)
    presabs_matrix[,j] <- temp_gene_present
  }
  colnames(presabs_matrix) <- names(exons_parsed[[1]])
  number_pres <- colSums(presabs_matrix)
  presence_count[i,] <- number_pres
}



rownames(presence_count) <- gene_list_for_conversion$name
colnames(presence_count) <- names(rearrange_iterations[[1]])
saveRDS(presence_count, file="Gene_presence_count_in_simulated_gene_lists_150826.rds")
presence_count <- readRDS(file="Gene_presence_count_in_simulated_gene_lists_150826.rds")

pres_and_length <- cbind(presence_count, gene_list_for_conversion$total_length)

#Remove genes that have length 0 and/all that never get sampled
# (probably due to masking or lack of SNPs in region), 
# which equates to 735 genes

#which(pres_and_length[,69]<1)
#which(rowSums(pres_and_length[,1:68])<1)
exclude_condition <- which(pres_and_length[,69]<1 | rowSums(pres_and_length[,1:68])<1)

pres_and_length <- pres_and_length[-(exclude_condition),]








name_vector <- c("C1_sig_gene_names","C2_sig_gene_names",
                 "C3_sig_gene_names","C4_sig_gene_names",
                 "C5_sig_gene_names","D1_sig_gene_names",
                 "D2_sig_gene_names","D3_sig_gene_names",
                 "D4_sig_gene_names","D5_sig_gene_names",
                 "C_sig_genes", "D1_noC_gene_names",
                 "D2_noC_gene_names", "D3_noC_gene_names",
                 "D4_noC_gene_names", "D5_noC_gene_names"
)
####### FIGURE S9 ########

jpeg(file="Gene sampling frequency vs length 150718.jpeg", width=15, height=18,
     units="in", res=300)
par(mfcol=c(7,6))
for(i in 1:42){
  temp_plot_label <- colnames(pres_and_length)[i]
  
  if(i<=16){
    temp_real_list <- get(name_vector[i])
  } else {temp_real_list <- get(temp_plot_label)}
  plot(pres_and_length[,i]~pres_and_length[,69], 
       col=rgb(100, 100, 100, 100, maxColorValue=255), pch=19, cex=0.8,
       xlab="Gene length (bp)", ylab=paste("Times resampled in", temp_plot_label))
  plot_in_red <- which(rownames(pres_and_length)%in%temp_real_list)
  points(pres_and_length[plot_in_red,i]~pres_and_length[plot_in_red,69],
         col=rgb(255, 0, 0, 100, maxColorValue=255), pch=19, cex=0.2)
  if(length(plot_in_red)>50){
    usespan<-0.1} else {usespan<-0.4}
  lwfull <- loess(pres_and_length[,i]~pres_and_length[,69], span=usespan)
  ordering <- order(pres_and_length[,69])
  lines(pres_and_length[ordering,69],lwfull$fitted[ordering],col="grey",lwd=1)
  lwobs <- loess(pres_and_length[plot_in_red,i]~pres_and_length[plot_in_red,69],
                 span=usespan)
  ordering2 <- order(pres_and_length[plot_in_red,69])
  lines(pres_and_length[plot_in_red,69][ordering2],lwobs$fitted[ordering2],col="darkred",lwd=1)
  
}
dev.off()


jpeg(file="Gene sampling frequency density 150718.jpeg", width=15, height=18,
     units="in", res=300)
par(mfcol=c(7,6))
for(i in 1:42){
  temp_plot_label <- colnames(pres_and_length)[i]  
  if(i<=16){
    temp_real_list <- get(name_vector[i])
  } else {temp_real_list <- get(temp_plot_label)}
  #plot(pres_and_length[,i]~pres_and_length[,69], 
  #     col=rgb(100, 100, 100, 100, maxColorValue=255), pch=19, cex=0.8,
  #     xlab="Gene length (bp)", ylab=paste("Times resampled in", temp_plot_label))
  plot(density(pres_and_length[,i]), 
       xlab=paste("No. times resampled in", temp_plot_label),
       ylab="Density", main=""
  )
  plot_in_red <- which(rownames(pres_and_length)%in%temp_real_list)
  lines(density(pres_and_length[plot_in_red,i]), col="red")
}
dev.off()

jpeg(file="Gene sampling length density.jpeg", width=15, height=18,
     units="in", res=300)
par(mfcol=c(7,6))
for(i in 1:42){
  temp_plot_label <- colnames(pres_and_length)[i]  
  if(i<=16){
    temp_real_list <- get(name_vector[i])
  } else {temp_real_list <- get(temp_plot_label)}
  #plot(pres_and_length[,i]~pres_and_length[,69], 
  #     col=rgb(100, 100, 100, 100, maxColorValue=255), pch=19, cex=0.8,
  #     xlab="Gene length (bp)", ylab=paste("Times resampled in", temp_plot_label))
  plot(density(pres_and_length[,69]), 
       xlab=paste("Gene length", temp_plot_label),
       ylab="Density", main="", xlim=c(0, 40000), ylim=c(0, 0.0001)
  )
  plot_in_red <- which(rownames(pres_and_length)%in%temp_real_list)
  lines(density(pres_and_length[plot_in_red,69]), col="red")
}
dev.off()


jpeg(file="Gene sampling length density scaled.jpeg", width=15, height=18,
     units="in", res=300)
par(mfcol=c(7,6))
for(i in 1:42){
  temp_plot_label <- colnames(pres_and_length)[i]  
  if(i<=16){
    temp_real_list <- get(name_vector[i])
  } else {temp_real_list <- get(temp_plot_label)}
  #plot(pres_and_length[,i]~pres_and_length[,69], 
  #     col=rgb(100, 100, 100, 100, maxColorValue=255), pch=19, cex=0.8,
  #     xlab="Gene length (bp)", ylab=paste("Times resampled in", temp_plot_label))
  plot(density(pres_and_length[,69]*pres_and_length[,i]/1000), 
       xlab=paste("Gene length", temp_plot_label),
       ylab="Density", main="", xlim=c(0, 40000), ylim=c(0, 0.0001)
  )
  plot_in_red <- which(rownames(pres_and_length)%in%temp_real_list)
  lines(density(pres_and_length[plot_in_red,69]), col="red")
}
dev.off()



hist(pres_and_length[,6], breaks=100, ylim=c(0, 100))
hist(pres_and_length[which(gene_list_for_conversion$name%in%D1_gene_names_exons_parsed),6], add=TRUE, col="red", breaks=100)




