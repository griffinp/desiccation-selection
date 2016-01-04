
#############
# FUNCTIONS #
#############


test_if_snp_in_list <- function(snp_list, snp){
  if(snp%in%snp_list){
    return(TRUE)
  } else {return(FALSE)}
}

# test_if_snp_in_D_or_C_lists <- function(snp_name, gene_lists){
#   present <- unlist(lapply(gene_lists, test_if_gene_in_list, gene=gene_name))
#   #output <- cbind(names(gene_lists), present)
#   return(present)
# }




remove_C_SNPs_from_D_list <- function(input_snp_list, C_SNP_table){
  temp_D_SNP_table <- t(as.data.frame(strsplit(input_snp_list, split=":")))
  colnames(temp_D_SNP_table) <- c("chr", "pos")
  for(i in 1:nrow(C_SNP_table)){
    temp_C_chr <- C_SNP_table[i,1]
    temp_C_up1000 <- C_SNP_table[i,4]
    temp_C_down1000 <- C_SNP_table[i,5]
    cond_vector <- temp_D_SNP_table[,1]==temp_C_chr & temp_D_SNP_table[,2]>=temp_C_up1000 & temp_D_SNP_table[,2]<=temp_C_down1000
    #print(summary(cond_vector))
    temp_D_SNP_table <- temp_D_SNP_table[cond_vector==FALSE,]
  }
  #new_object_name <- paste("noC_", input_snp_list, sep="")
  #write.table(temp_D_SNP_table, file=new_file_name, sep="\t", col.names=FALSE,
  #            row.names=FALSE, quote=FALSE)
  processed_snp_list <- paste(temp_D_SNP_table[,1], as.character(as.numeric(temp_D_SNP_table[,2])), sep=":")
  return(processed_snp_list)
}


########
# MAIN #
########

Sample_code <- c("C1", "C2", "C3", "C4", "C5",
                 "D1", "D2", "D3", "D4", "D5")

# Get object containing simulated snp lists and overlap lists
setwd("/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/")

snp_simulation <- readRDS(file="SNP_and_eQTL_lists_from_resampled_SNP_table.rds")
snp_simulation <- lapply(snp_simulation, "[[", 1)

# Get object containing coordinates of C genes (plus 1000 bp up and downstream)

setwd('/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists')
allC_coords <- read.table(file="Putative_lab_adaptation_gene_coords_plus_1000_bp.txt", 
                          header=TRUE, stringsAsFactors=FALSE)


# Extract overlap snps called in each combination of D replicates

setwd("/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/")


# Now get real snp lists (again, before removing SNPs close to C genes)
for(i in Sample_code[1:10]){
  if(i %in% c('D1', 'D2', 'D3', 'D4', 'D5')){
    temp_input_name <- paste('/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/',
                             i, '_sig_SNPs.txt', sep="")
  } else{
    temp_input_name <- paste('/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/',
                             i, '_sig_SNPs.txt', sep="")
  }
  temp_file <- read.csv(temp_input_name, header=FALSE, stringsAsFactors=FALSE)[,1]
  temp_object_name <- paste(i, '_sig_SNPs', sep="")
  assign(temp_object_name, temp_file)
}

# get real overlaps

C1_C2<-intersect(C1_sig_SNPs, C2_sig_SNPs)
C1_C3<-intersect(C1_sig_SNPs, C3_sig_SNPs)
C1_C4<-intersect(C1_sig_SNPs, C4_sig_SNPs)
C1_C5<-intersect(C1_sig_SNPs, C5_sig_SNPs)
C2_C3<-intersect(C2_sig_SNPs, C3_sig_SNPs)
C2_C4<-intersect(C2_sig_SNPs, C4_sig_SNPs)
C2_C5<-intersect(C2_sig_SNPs, C5_sig_SNPs)
C3_C4<-intersect(C3_sig_SNPs, C4_sig_SNPs)
C3_C5<-intersect(C3_sig_SNPs, C5_sig_SNPs)
C4_C5<-intersect(C4_sig_SNPs, C5_sig_SNPs)
C1_C2_C3<-intersect(intersect(C1_sig_SNPs, C2_sig_SNPs), C3_sig_SNPs)
C1_C2_C4<-intersect(intersect(C1_sig_SNPs, C2_sig_SNPs), C4_sig_SNPs)
C1_C2_C5<-intersect(intersect(C1_sig_SNPs, C2_sig_SNPs), C5_sig_SNPs)
C1_C3_C4<-intersect(intersect(C1_sig_SNPs, C3_sig_SNPs), C4_sig_SNPs)
C1_C3_C5<-intersect(intersect(C1_sig_SNPs, C3_sig_SNPs), C5_sig_SNPs)
C1_C4_C5<-intersect(intersect(C1_sig_SNPs, C4_sig_SNPs), C5_sig_SNPs)
C2_C3_C4<-intersect(intersect(C2_sig_SNPs, C3_sig_SNPs), C4_sig_SNPs)
C2_C3_C5<-intersect(intersect(C2_sig_SNPs, C3_sig_SNPs), C5_sig_SNPs)
C2_C4_C5<-intersect(intersect(C2_sig_SNPs, C4_sig_SNPs), C5_sig_SNPs)
C3_C4_C5<-intersect(intersect(C3_sig_SNPs, C4_sig_SNPs), C5_sig_SNPs)
C1_C2_C3_C4<-intersect(intersect(C1_sig_SNPs, C2_sig_SNPs), 
                       intersect(C3_sig_SNPs, C4_sig_SNPs))
C1_C2_C3_C5<-intersect(intersect(C1_sig_SNPs, C2_sig_SNPs),
                       intersect(C3_sig_SNPs, C5_sig_SNPs))
C1_C2_C4_C5<-intersect(intersect(C1_sig_SNPs, C2_sig_SNPs),
                       intersect(C4_sig_SNPs, C5_sig_SNPs))
C1_C3_C4_C5<-intersect(intersect(C1_sig_SNPs, C3_sig_SNPs),
                       intersect(C4_sig_SNPs, C5_sig_SNPs))
C2_C3_C4_C5<-intersect(intersect(C2_sig_SNPs, C3_sig_SNPs),
                       intersect(C4_sig_SNPs, C5_sig_SNPs))
C1_C2_C3_C4_C5<-intersect(intersect(C1_sig_SNPs, C2_sig_SNPs), 
                          intersect(intersect(C3_sig_SNPs, C4_sig_SNPs), C5_sig_SNPs))


D1_D2<-intersect(D1_sig_SNPs, D2_sig_SNPs)
D1_D3<-intersect(D1_sig_SNPs, D3_sig_SNPs)
D1_D4<-intersect(D1_sig_SNPs, D4_sig_SNPs)
D1_D5<-intersect(D1_sig_SNPs, D5_sig_SNPs)
D2_D3<-intersect(D2_sig_SNPs, D3_sig_SNPs)
D2_D4<-intersect(D2_sig_SNPs, D4_sig_SNPs)
D2_D5<-intersect(D2_sig_SNPs, D5_sig_SNPs)
D3_D4<-intersect(D3_sig_SNPs, D4_sig_SNPs)
D3_D5<-intersect(D3_sig_SNPs, D5_sig_SNPs)
D4_D5<-intersect(D4_sig_SNPs, D5_sig_SNPs)
D1_D2_D3<-intersect(intersect(D1_sig_SNPs, D2_sig_SNPs), D3_sig_SNPs)
D1_D2_D4<-intersect(intersect(D1_sig_SNPs, D2_sig_SNPs), D4_sig_SNPs)
D1_D2_D5<-intersect(intersect(D1_sig_SNPs, D2_sig_SNPs), D5_sig_SNPs)
D1_D3_D4<-intersect(intersect(D1_sig_SNPs, D3_sig_SNPs), D4_sig_SNPs)
D1_D3_D5<-intersect(intersect(D1_sig_SNPs, D3_sig_SNPs), D5_sig_SNPs)
D1_D4_D5<-intersect(intersect(D1_sig_SNPs, D4_sig_SNPs), D5_sig_SNPs)
D2_D3_D4<-intersect(intersect(D2_sig_SNPs, D3_sig_SNPs), D4_sig_SNPs)
D2_D3_D5<-intersect(intersect(D2_sig_SNPs, D3_sig_SNPs), D5_sig_SNPs)
D2_D4_D5<-intersect(intersect(D2_sig_SNPs, D4_sig_SNPs), D5_sig_SNPs)
D3_D4_D5<-intersect(intersect(D3_sig_SNPs, D4_sig_SNPs), D5_sig_SNPs)
D1_D2_D3_D4<-intersect(intersect(D1_sig_SNPs, D2_sig_SNPs), 
                       intersect(D3_sig_SNPs, D4_sig_SNPs))
D1_D2_D3_D5<-intersect(intersect(D1_sig_SNPs, D2_sig_SNPs),
                       intersect(D3_sig_SNPs, D5_sig_SNPs))
D1_D2_D4_D5<-intersect(intersect(D1_sig_SNPs, D2_sig_SNPs),
                       intersect(D4_sig_SNPs, D5_sig_SNPs))
D1_D3_D4_D5<-intersect(intersect(D1_sig_SNPs, D3_sig_SNPs),
                       intersect(D4_sig_SNPs, D5_sig_SNPs))
D2_D3_D4_D5<-intersect(intersect(D2_sig_SNPs, D3_sig_SNPs),
                       intersect(D4_sig_SNPs, D5_sig_SNPs))
D1_D2_D3_D4_D5<-intersect(intersect(D1_sig_SNPs, D2_sig_SNPs), 
                          intersect(intersect(D3_sig_SNPs, D4_sig_SNPs), D5_sig_SNPs))


snp_overlap_lengths <- c(length(D1_D2)/length(D1), length(D1_D2)/length(D2),
                         length(D1_D3)/length(D1), length(D1_D3)/length(D3),
                         length(D1_D4)/length(D1), length(D1_D4)/length(D4),
                         length(D1_D5)/length(D1), length(D1_D5)/length(D5),
                         length(D2_D3)/length(D2), length(D2_D3)/length(D3),
                         length(D2_D4)/length(D2), length(D2_D4)/length(D4),
                         length(D2_D5)/length(D2), length(D2_D5)/length(D5),
                         length(D3_D4)/length(D3), length(D3_D4)/length(D4),
                         length(D3_D5)/length(D3), length(D3_D5)/length(D5),
                         length(D4_D5)/length(D4), length(D4_D5)/length(D5))

D1_unique_snps <- length(setdiff(D1, union(union(D2, D3), union(D4, D5))))/length(D1)
D2_unique_snps <- length(setdiff(D2, union(union(D1, D3), union(D4, D5))))/length(D2)
D3_unique_snps <- length(setdiff(D3, union(union(D1, D2), union(D4, D5))))/length(D3)
D4_unique_snps <- length(setdiff(D4, union(union(D2, D3), union(D1, D5))))/length(D4)
D5_unique_snps <- length(setdiff(D5, union(union(D2, D3), union(D4, D1))))/length(D5)

C1_unique_snps <- length(setdiff(C1, union(union(C2, C3), union(C4, C5))))/length(C1)
C2_unique_snps <- length(setdiff(C2, union(union(C1, C3), union(C4, C5))))/length(C2)
C3_unique_snps <- length(setdiff(C3, union(union(C1, C2), union(C4, C5))))/length(C3)
C4_unique_snps <- length(setdiff(C4, union(union(C2, C3), union(C1, C5))))/length(C4)
C5_unique_snps <- length(setdiff(C5, union(union(C2, C3), union(C4, C1))))/length(C5)

unique_comparison <- data.frame(number_unique=c(D1_unique_snps, D2_unique_snps, D3_unique_snps,
                                                D4_unique_snps, D5_unique_snps, C1_unique_snps,
                                                C2_unique_snps, C3_unique_snps, C4_unique_snps,
                                                C5_unique_snps), 
                                C_or_D=c(rep("D", times=5), rep("C", times=5)))

t.test(number_unique~C_or_D, data=unique_comparison)


# now rearrange simulated SNP lists so as to create overlap lists

rearrange_snp_iterations <- list()
for(j in 1:nIter) {
  if(j%%50==0){
    print(paste("Rearranging iteration #", j, sep=" "))
  }
  temp_name <- paste("Iter", j, sep="_")
  rearrange_snp_iterations[[temp_name]] <- lapply(snp_simulation, "[[", j)
  
  names(rearrange_snp_iterations[[j]]) <- Sample_code
  #temp_allC <- 
  #rearrange_snp_iterations[[j]][["allC"]] <- union(union(union(rearrange_snp_iterations[[j]][[1]], 
  #                                                           rearrange_snp_iterations[[j]][[2]]), 
  #                                                     rearrange_snp_iterations[[j]][[3]]), 
  #                                               union(rearrange_snp_iterations[[j]][[4]], 
  #                                                     rearrange_snp_iterations[[j]][[5]]))
  D_names <- c("D1", "D2", "D3", "D4", "D5")
  #   for(i in 1:5){
  #     temp_D <- D_names[i]
  #     #temp_noC_name <- paste(temp_D, "noC", sep="_")
  #     rearrange_snp_iterations[[j]][[temp_D]] <- rearrange_snp_iterations[[j]][[temp_D]]
  #   }
  combn2 <- combn(D_names, 2)
  combn2a <- apply(combn2, 2, paste, collapse="_")
  for(i in 1:ncol(combn2)){
    temp_Da <- rearrange_snp_iterations[[j]][[combn2[1,i]]]
    temp_Db <- rearrange_snp_iterations[[j]][[combn2[2,i]]]
    temp_combn_name <- combn2a[i]
    rearrange_snp_iterations[[j]][[temp_combn_name]] <- intersect(temp_Da, temp_Db)
  }
  combn3 <- combn(D_names, 3)
  combn3a <- apply(combn3, 2, paste, collapse="_")
  for(i in 1:ncol(combn3)){
    temp_Da <- rearrange_snp_iterations[[j]][[combn3[1,i]]]
    temp_Db <- rearrange_snp_iterations[[j]][[combn3[2,i]]]
    temp_Dc <- rearrange_snp_iterations[[j]][[combn3[3,i]]]
    temp_combn_name <- combn3a[i]
    rearrange_snp_iterations[[j]][[temp_combn_name]] <- intersect(temp_Da, intersect(temp_Db, temp_Dc))
  } 
  combn4 <- combn(D_names, 4)
  combn4a <- apply(combn4, 2, paste, collapse="_")
  for(i in 1:ncol(combn4)){
    temp_Da <- rearrange_snp_iterations[[j]][[combn4[1,i]]]
    temp_Db <- rearrange_snp_iterations[[j]][[combn4[2,i]]]
    temp_Dc <- rearrange_snp_iterations[[j]][[combn4[3,i]]]
    temp_Dd <- rearrange_snp_iterations[[j]][[combn4[4,i]]]
    temp_combn_name <- combn4a[i]
    rearrange_snp_iterations[[j]][[temp_combn_name]] <- intersect(intersect(temp_Da, temp_Db), intersect(temp_Dc, temp_Dd))
  } 
  combn5 <- "D1_D2_D3_D4_D5"
  rearrange_snp_iterations[[j]][[combn5]] <- intersect(intersect(rearrange_snp_iterations[[j]][["D1"]], 
                                                                 intersect(rearrange_snp_iterations[[j]][["D2"]], 
                                                                           rearrange_snp_iterations[[j]][["D3"]])), 
                                                       intersect(rearrange_snp_iterations[[j]][["D4"]], 
                                                                 rearrange_snp_iterations[[j]][["D5"]]))
  C_names <- c("C1", "C2", "C3", "C4", "C5")
  #   for(i in 1:5){
  #     temp_C <- C_names[i]
  #   }
  combn2 <- combn(C_names, 2)
  combn2a <- apply(combn2, 2, paste, collapse="_")
  for(i in 1:ncol(combn2)){
    temp_Ca <- rearrange_snp_iterations[[j]][[combn2[1,i]]]
    temp_Cb <- rearrange_snp_iterations[[j]][[combn2[2,i]]]
    temp_combn_name <- combn2a[i]
    rearrange_snp_iterations[[j]][[temp_combn_name]] <- intersect(temp_Ca, temp_Cb)
  }
  combn3 <- combn(C_names, 3)
  combn3a <- apply(combn3, 2, paste, collapse="_")
  for(i in 1:ncol(combn3)){
    temp_Ca <- rearrange_snp_iterations[[j]][[combn3[1,i]]]
    temp_Cb <- rearrange_snp_iterations[[j]][[combn3[2,i]]]
    temp_Cc <- rearrange_snp_iterations[[j]][[combn3[3,i]]]
    temp_combn_name <- combn3a[i]
    rearrange_snp_iterations[[j]][[temp_combn_name]] <- intersect(temp_Ca, intersect(temp_Cb, temp_Cc))
  } 
  combn4 <- combn(C_names, 4)
  combn4a <- apply(combn4, 2, paste, collapse="_")
  for(i in 1:ncol(combn4)){
    temp_Ca <- rearrange_snp_iterations[[j]][[combn4[1,i]]]
    temp_Cb <- rearrange_snp_iterations[[j]][[combn4[2,i]]]
    temp_Cc <- rearrange_snp_iterations[[j]][[combn4[3,i]]]
    temp_Cd <- rearrange_snp_iterations[[j]][[combn4[4,i]]]
    temp_combn_name <- combn4a[i]
    rearrange_snp_iterations[[j]][[temp_combn_name]] <- intersect(intersect(temp_Ca, temp_Cb), intersect(temp_Cc, temp_Cd))
  } 
  combn5 <- "C1_C2_C3_C4_C5"
  rearrange_snp_iterations[[j]][[combn5]] <- intersect(intersect(rearrange_snp_iterations[[j]][["C1"]], 
                                                                 intersect(rearrange_snp_iterations[[j]][["C2"]], 
                                                                           rearrange_snp_iterations[[j]][["C3"]])), 
                                                       intersect(rearrange_snp_iterations[[j]][["C4"]], 
                                                                 rearrange_snp_iterations[[j]][["C5"]]))
  
  
  
}











###### FIGURE S7 ######

setwd("/Users/pgriffin/Documents/Drosophila Selection Experiment/gene_list_overlap_testing/")

pdf("Histograms of simulated snp number overlap, resampling with position.pdf", width=14, height=16)
par(mfcol=c(4,3))

combn2 <- combn(D_names, 2)
combn2a <- apply(combn2, 2, paste, collapse="_")
#combn3 <- combn(D_names, 3)
#combn3a <- apply(combn3, 2, paste, collapse="_")
#combn4 <- combn(D_names, 4)
#combn4a <- apply(combn4, 2, paste, collapse="_")
#combn5 <- "D1_D2_D3_D4_D5"
allcombns <- c(combn2a)
for(i in 1:length(allcombns)){
  temp_combn <- allcombns[i]
  real_overlap_length <- length(get(temp_combn))
  temp_all_list <- lapply(rearrange_snp_iterations, "[[", i+10)
  temp_dist <- sapply(temp_all_list, length)
  
  multiplier <- hist(temp_dist, plot=FALSE)$counts / hist(temp_dist, plot=FALSE)$density
  mydensity <- density(temp_dist)
  mydensity$y <- mydensity$y * multiplier[1] 
  
  hist(temp_dist, main="", ylim=c(0, 600), 
       xlim=c(0, real_overlap_length+10),
       #xlim=c(0, 700),
       xlab=paste("No. SNPs in", temp_combn, "overlap", sep=" "))
  lines(mydensity)
  lines(x=c(real_overlap_length, real_overlap_length),
        y=c(0, 225), col="red")
  myx <- seq(min(temp_dist), max(temp_dist), length = 100)
  normal <- dnorm(x=myx, mean = mean(temp_dist), sd = sd(temp_dist))
  lines(myx, normal * multiplier[1], col = "blue", lwd = 1)
  normless <- pnorm(q = real_overlap_length, mean=mean(temp_dist), sd=sd(temp_dist))
  #pval <- (1-normless)*2
  text(x=real_overlap_length, y=275, labels=signif(normless, 3))
  text(x=2, y=500, labels=LETTERS[i], cex=2)
}
dev.off() 

######## Figure showing histograms of simulated overlap for C replicates ########

setwd("/Users/pgriffin/Documents/Drosophila Selection Experiment/gene_list_overlap_testing/")

pdf("Histograms of simulated snp number overlap for C replicates, resampling with position.pdf", width=14, height=16)
par(mfcol=c(4,3))
C_names <- Sample_code[1:5]
combn2 <- combn(C_names, 2)
combn2a <- apply(combn2, 2, paste, collapse="_")
combn3 <- combn(C_names, 3)
combn3a <- apply(combn3, 2, paste, collapse="_")
combn4 <- combn(C_names, 4)
combn4a <- apply(combn4, 2, paste, collapse="_")
combn5 <- "C1_C2_C3_C4_C5"
allcombns <- c(combn2a)
for(i in 1:length(allcombns)){
  temp_combn <- allcombns[i]
  real_overlap_length <- length(get(temp_combn))
  temp_all_list <- lapply(rearrange_snp_iterations, "[[", i+36)
  temp_dist <- sapply(temp_all_list, length)
  
  multiplier <- hist(temp_dist, plot=FALSE)$counts / hist(temp_dist, plot=FALSE)$density
  mydensity <- density(temp_dist)
  mydensity$y <- mydensity$y * multiplier[1] 
  
  hist(temp_dist, main="", ylim=c(0, 1000), xlim=c(0, 13),
       xlab=paste("No. SNPs in", temp_combn, "overlap", sep=" "))
  #lines(mydensity)
  lines(x=c(real_overlap_length, real_overlap_length),
        y=c(0, 700), col="red")
  myx <- seq(min(temp_dist), max(temp_dist), length = 100)
  normal <- dnorm(x=myx, mean = mean(temp_dist), sd = sd(temp_dist))
  #lines(myx, normal * multiplier[1], col = "blue", lwd = 1)
  normless <- pnorm(q = real_overlap_length, mean=mean(temp_dist), sd=sd(temp_dist))
  #pval <- (1-normless)*2
  #text(x=real_overlap_length, y=275, labels=signif(normless, 3))
  text(x=2, y=800, labels=LETTERS[i], cex=2)
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
                                       sep="\t", stringsAsFactors=FALSE, header=TRUE,
                                       quote="\"")



# test_l <- get_gene_lengths(rearrange_iterations[[1]][[1]], gene_list_for_conversion)
# test_l2 <- get_gene_lengths(D1_gene_names_exons_parsed, gene_list_for_conversion)
# 
# gene_length_breaks <- seq(0, 401000, by=1000)
# 
# genes_in_temp_list <- gene_list_for_conversion[gene_list_for_conversion$name%in%D1_gene_names_exons_parsed,]
# plot(genes_in_temp_list$total_length, y=rep(1, times=nrow(genes_in_temp_list)))

#find count of resampled lists that each gene appears in

nIter=1000

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
colnames(presence_count) <- names(exons_parsed[[1]])
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

jpeg(file="Gene sampling frequency vs length 151011.jpeg", width=15, height=18,
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


jpeg(file="Gene sampling frequency density 151011.jpeg", width=15, height=18,
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

jpeg(file="Gene sampling length density 151011.jpeg", width=15, height=18,
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


jpeg(file="Gene sampling length density scaled 151011.jpeg", width=15, height=18,
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




