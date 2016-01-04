setwd("~/Documents/Drosophila Selection Experiment/snp_and_gene_lists")

Sample_code <- c("C1", "C2", "C3", "C4", "C5",
                 "D1", "D2", "D3", "D4", "D5")

eqtl_table <- read.table("eQTLs_Huang_et_al_2015.txt", sep="\t",
                         stringsAsFactors=FALSE, header=TRUE)

SNPs_only <- subset(eqtl_table, eqtl_table$SNP_or_DEL=="SNP")
#all_eQTLs <- unique(SNPs_only$v6_coord)

zfor(i in Sample_code){
  if(which(Sample_code==i) > 5){
    #temp_sig_SNPs_file_name <- paste("noC_", i, "_sig_SNPs.txt", sep="")
    # USE FULL SNP LISTS, NOT ONES WITH LAB-ADAPTATION SNPS REMOVED
    temp_sig_SNPs_file_name <- paste(i, "_sig_SNPs.txt", sep="")
  } else {
    temp_sig_SNPs_file_name <- paste(i, "_sig_SNPs.txt", sep="")
  }
  temp_sig_eQTLs_df_name <- paste(i, "_sig_eQTLs", sep="")
  temp_sig_eQTL_gene_list <- paste(i, "_sig_eQTL_genes", sep="")
  temp_sig_SNPs <- read.csv(temp_sig_SNPs_file_name, header=FALSE, sep="\t",
                          stringsAsFactors=FALSE, strip.white=TRUE)
  
  temp_sig_SNPs[,2] <- as.character(as.numeric(temp_sig_SNPs[,2]))
  
  temp_sig <- apply(temp_sig_SNPs, 1, paste0, collapse=":")
  length(which(temp_sig%in%SNPs_only$v6_coord))
  
  temp_sig_eQTLs <- SNPs_only[which(SNPs_only$v6_coord%in%temp_sig),]
  assign(temp_sig_eQTLs_df_name, temp_sig_eQTLs)
  temp_sig_eQTL_genes_pre <- levels(as.factor(temp_sig_eQTLs$Gene.Symbol))
  temp_sig_eQTL_genes <- subset(temp_sig_eQTL_genes_pre, temp_sig_eQTL_genes_pre!="-")
  assign(temp_sig_eQTL_gene_list, temp_sig_eQTL_genes)
  
  temp_important <- temp_sig_eQTLs[,c("v6_coord", "Locus_Type", "Sex",
                                      "FlyBase.ID", "Gene.Symbol")]
  write.table(temp_important, file=paste(i, "_eQTL_withC_sig_SNPs.txt", sep=""),
              sep="\t", row.names=FALSE,
              col.names=TRUE, quote=FALSE)
#   temp_all_SNPs_file_name <- paste(i, "_all_SNPs.txt", sep="")
#   temp_all_SNPs <- read.table(temp_all_SNPs_file_name, header=FALSE, sep="\t",
#                               stringsAsFactors=FALSE, strip.white=TRUE)
#   temp_all_SNPs[,2] <- as.character(as.numeric(temp_all_SNPs[,2]))
#   temp_all <- apply(temp_all_SNPs, 1, paste0, collapse=":")
#   
}

#############
# QUESTIONS #
#############

# Do more SNPs hit known eQTLs than expected by chance? 
# (use simulations? Alternatively could build into snpEff somehow?)

# From my SNP simulations, I only saved gene lists previously, not SNP lists
# Therefore I need to redo them!

#############
# FUNCTIONS #
#############

extract_SNPs_from_table <- function(snpsift_table){
  temp_snp_vector <- paste(snpsift_table[,1], snpsift_table[,2], sep=":")
  return(temp_snp_vector)
}

########
# MAIN #
########

#Resample simulated "significant" SNPs from the file of all SNPs

Sample_code <- c("C1", "C2", "C3", "C4", "C5", 
                 "D1", "D2", "D3", "D4", "D5")

setwd("~/Documents/Drosophila Selection Experiment/snp_and_gene_lists")

# Now look at the results
colname_vector <- c("chr", "pos", paste(rep("genename", times=88), 1:88, sep=""))

C1_sig_ex <- read.table("MB_C1_sig_genes_extracted.txt", fill=TRUE, header=FALSE, 
                        colClasses="character", col.names=colname_vector)
C2_sig_ex <- read.table("MB_C2_sig_genes_extracted.txt", fill=TRUE, header=FALSE, 
                        colClasses="character", col.names=colname_vector)
C3_sig_ex <- read.table("MB_C3_sig_genes_extracted.txt", fill=TRUE, header=FALSE, 
                        colClasses="character", col.names=colname_vector)
C4_sig_ex <- read.table("MB_C4_sig_genes_extracted.txt", fill=TRUE, header=FALSE, 
                        colClasses="character", col.names=colname_vector)
C5_sig_ex <- read.table("MB_C5_sig_genes_extracted.txt", fill=TRUE, header=FALSE, 
                        colClasses="character", col.names=colname_vector)
D1_sig_ex <- read.table("MB_D1_sig_genes_extracted.txt", fill=TRUE, header=FALSE, 
                        colClasses="character", col.names=colname_vector)
D2_sig_ex <- read.table("MB_D2_sig_genes_extracted.txt", fill=TRUE, header=FALSE, 
                        colClasses="character", col.names=colname_vector)
D3_sig_ex <- read.table("MB_D3_sig_genes_extracted.txt", fill=TRUE, header=FALSE, 
                        colClasses="character", col.names=colname_vector)
D4_sig_ex <- read.table("MB_D4_sig_genes_extracted.txt", fill=TRUE, header=FALSE, 
                        colClasses="character", col.names=colname_vector)
D5_sig_ex <- read.table("MB_D5_sig_genes_extracted.txt", fill=TRUE, header=FALSE, 
                        colClasses="character", col.names=colname_vector)

C1_sig_SNPs <- paste(C1_sig_ex[,1], as.character(as.numeric(C1_sig_ex[,2])), sep=":")
C2_sig_SNPs <- paste(C2_sig_ex[,1], as.character(as.numeric(C2_sig_ex[,2])), sep=":")
C3_sig_SNPs <- paste(C3_sig_ex[,1], as.character(as.numeric(C3_sig_ex[,2])), sep=":")
C4_sig_SNPs <- paste(C4_sig_ex[,1], as.character(as.numeric(C4_sig_ex[,2])), sep=":")
C5_sig_SNPs <- paste(C5_sig_ex[,1], as.character(as.numeric(C5_sig_ex[,2])), sep=":")
D1_sig_SNPs <- paste(D1_sig_ex[,1], as.character(as.numeric(D1_sig_ex[,2])), sep=":")
D2_sig_SNPs <- paste(D2_sig_ex[,1], as.character(as.numeric(D2_sig_ex[,2])), sep=":")
D3_sig_SNPs <- paste(D3_sig_ex[,1], as.character(as.numeric(D3_sig_ex[,2])), sep=":")
D4_sig_SNPs <- paste(D4_sig_ex[,1], as.character(as.numeric(D4_sig_ex[,2])), sep=":")
D5_sig_SNPs <- paste(D5_sig_ex[,1], as.character(as.numeric(D5_sig_ex[,2])), sep=":")

resample_names <- c("C1", "C2", "C3", "C4", "C5","D1", "D2", "D3", "D4", "D5"
)
nIter=1000

### Resampling SNPs where position
### is taken into account. 

for (i in 1:length(resample_names)) {
  temp_name <- resample_names[i]
  temp_filename <- paste("/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/MB_",
                         temp_name,
                         "_all_genes_extracted.txt", sep="")
  temp_full <- read.table(temp_filename, fill=TRUE, header=FALSE, 
                          colClasses="character", col.names=colname_vector)
  temp_full <- temp_full[temp_full$chr%in%c("2L", "2R", "3L", "3R", "X", "dmel_mitochondrion_genome"),]
  # SORT BY CHR AND POS! 
  temp_full <- temp_full[with(temp_full, order(temp_full$chr, as.integer(temp_full$pos))),]
  temp_full$chrpos <- paste(temp_full$chr, temp_full$pos, sep="_")
  temp_sig <- get(paste(temp_name, "_sig_ex", sep=""))
  temp_sig$chrpos <- paste(temp_sig$chr, temp_sig$pos, sep="_")
  #temp_length <- nrow(temp_sig)
  find_sig_rows <- which(temp_full$chrpos %in% temp_sig$chrpos)
  resample_snp_table=as.list(rep(NA, length=nIter))
  resample_eQTL_table=as.list(rep(NA, length=nIter))
  resample_eQTL_gene_table=as.list(rep(NA, length=nIter))
  for (j in 1:nIter) {
    
    temp_resample_index_shift <- sample(1:nrow(temp_full), size=1)
    temp_resample_index <- find_sig_rows + temp_resample_index_shift
    # 'wrap around' indices that would extend off the table
    temp_resample_index_wrapped <- temp_resample_index[temp_resample_index > nrow(temp_full)]
    temp_resample_index_to_use <- c(temp_resample_index_wrapped - nrow(temp_full),
                                    temp_resample_index[temp_resample_index < nrow(temp_full)])
    #print(cbind(find_sig_rows, temp_resample_index, temp_resample_index_to_use))
    temp_resample <- temp_full[temp_resample_index_to_use,]
    temp_snp_list <- extract_SNPs_from_table(temp_resample)
    temp_eQTL_hits <- temp_snp_list[which(temp_snp_list%in%SNPs_only$v6_coord)]
    temp_eQTL_genes_pre <- SNPs_only[which(SNPs_only$v6_coord%in%temp_snp_list),2]
    temp_eQTL_genes <- unique(subset(temp_eQTL_genes_pre, temp_eQTL_genes_pre!="-"))
    resample_snp_table[[j]] <- temp_snp_list
    resample_eQTL_table[[j]] <- temp_eQTL_hits
    resample_eQTL_gene_table[[j]] <- temp_eQTL_genes
    if(j%%50==0){
      print(paste(temp_name, " Resampling iteration #", j, sep=""))
      print(paste(nrow(temp_sig), "sig. SNPs;", 
                  length(find_sig_rows), "sig. SNPs in full data;",
                  temp_resample_index_shift, "is the index shift;",
                  nrow(temp_resample), "rows in resampled table;",
                  "Equates to", length(temp_eQTL_hits), "eQTLs hitting",
                  length(temp_eQTL_genes), "genes",
                  sep=" "))
    }
  }
  assign(paste(temp_name, "resampled_SNPs_with_position", sep="_"), resample_snp_table)
  assign(paste(temp_name, "resampled_SNPs_hitting_eQTLs", sep="_"), resample_eQTL_table)
  assign(paste(temp_name, "resampled_eQTL_genes_with_position", sep="_"), resample_eQTL_gene_table)
  remove(temp_full)
}

full_resample_snp_table <- list(C1=list(C1_resampled_SNPs_with_position, C1_resampled_SNPs_hitting_eQTLs, 
                                     C1_resampled_eQTL_genes_with_position),
                                     C2=list(C2_resampled_SNPs_with_position, C2_resampled_SNPs_hitting_eQTLs, 
                                     C2_resampled_eQTL_genes_with_position),
                                     C3=list(C3_resampled_SNPs_with_position, C3_resampled_SNPs_hitting_eQTLs, 
                                     C3_resampled_eQTL_genes_with_position),
                                     C4=list(C4_resampled_SNPs_with_position, C4_resampled_SNPs_hitting_eQTLs, 
                                     C4_resampled_eQTL_genes_with_position),
                                     C5=list(C5_resampled_SNPs_with_position, C5_resampled_SNPs_hitting_eQTLs, 
                                     C5_resampled_eQTL_genes_with_position),
                                     D1=list(D1_resampled_SNPs_with_position, D1_resampled_SNPs_hitting_eQTLs, 
                                     D1_resampled_eQTL_genes_with_position),
                                     D2=list(D2_resampled_SNPs_with_position, D2_resampled_SNPs_hitting_eQTLs, 
                                     D2_resampled_eQTL_genes_with_position),
                                     D3=list(D3_resampled_SNPs_with_position, D3_resampled_SNPs_hitting_eQTLs, 
                                     D3_resampled_eQTL_genes_with_position),
                                     D4=list(D4_resampled_SNPs_with_position, D4_resampled_SNPs_hitting_eQTLs, 
                                     D4_resampled_eQTL_genes_with_position),
                                     D5=list(D5_resampled_SNPs_with_position, D5_resampled_SNPs_hitting_eQTLs, 
                                     D5_resampled_eQTL_genes_with_position))

setwd("/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/")
saveRDS(full_resample_snp_table, file="SNP_and_eQTL_lists_from_resampled_SNP_table.rds")

#setwd("/Users/pgriffin/Documents/Drosophila Selection Experiment/gene_list_overlap_testing/")

pdf("Histograms of eQTL SNP hits resampling with position.pdf", width=8, height=12)

par(mfcol=c(5,2))
for(i in 1:10){
  temp_sample <- Sample_code[i]
  temp_observed_list <- get(paste(temp_sample, "_sig_eQTLs", sep=""))
  temp_all_list <- full_resample_snp_table[[i]]
  names(temp_all_list) <- c("number_resampled_snps", "number_resampled_snps_hitting_eQTLs", "number_genes_hit_by_resampled_eQTL_snps")
  temp_snps_hitting_eQTL_list <- temp_all_list[[2]]
  temp_all_length <- sapply(temp_snps_hitting_eQTL_list, length)

  temp_dist <- temp_all_length
  
  multiplier <- hist(temp_dist, plot=FALSE)$counts / hist(temp_dist, plot=FALSE)$density
  mydensity <- density(temp_dist)
  mydensity$y <- mydensity$y * multiplier[1] 
  #offset <- sd(temp_dist)*2.5
  hist(temp_dist, 
       xlim=c(0, mean(temp_dist)*4), ylim=c(0, 400), main="",
       xlab=paste("No.", temp_sample, "SNPs hitting eQTLs or veQTLs", sep=" "))
  lines(mydensity)
  lines(y=c(0, 250), 
        x=c(nrow(temp_observed_list), nrow(temp_observed_list)),
        col="red")
  myx <- seq(min(temp_dist), max(temp_dist), length = 100)
  normal <- dnorm(x=myx, mean = mean(temp_dist), sd = sd(temp_dist))
  lines(myx, normal * multiplier[1], col = "blue", lwd = 1)
  normless <- pnorm(q = nrow(temp_observed_list), mean=mean(temp_dist), sd=sd(temp_dist))
  #pval <- (1-normless)*2
  text(x=nrow(temp_observed_list), y=300, labels=signif(normless, 3))
  text(x=0, y=300, labels=LETTERS[i], cex=2)
  
}
dev.off()


# ANSWER: Yes, it seems that eQTL SNPs are significantly overrepresented in our candidate lists...
# for all D replicates and some of the C replicates too (reps C3, C5)
# *REMEMBER* this is using the D SNP lists without excluding any putative lab-adaptation SNPs

##############
# QUESTION 2 #
##############

# Do the same eQTL SNPs get hit more than expected by chance across replicates?
# related: Do eQTL SNPs *for the same genes* get hit more than expected by chance across replicates?

# OVERLAP AMONG OBSERVED EQTL SNPS

C1_C2<-intersect(C1_sig_eQTLs$v6_coord, C2_sig_eQTLs$v6_coord)
C1_C3<-intersect(C1_sig_eQTLs$v6_coord, C3_sig_eQTLs$v6_coord)
C1_C4<-intersect(C1_sig_eQTLs$v6_coord, C4_sig_eQTLs$v6_coord)
C1_C5<-intersect(C1_sig_eQTLs$v6_coord, C5_sig_eQTLs$v6_coord)
C2_C3<-intersect(C2_sig_eQTLs$v6_coord, C3_sig_eQTLs$v6_coord)
C2_C4<-intersect(C2_sig_eQTLs$v6_coord, C4_sig_eQTLs$v6_coord)
C2_C5<-intersect(C2_sig_eQTLs$v6_coord, C5_sig_eQTLs$v6_coord)
C3_C4<-intersect(C3_sig_eQTLs$v6_coord, C4_sig_eQTLs$v6_coord)
C3_C5<-intersect(C3_sig_eQTLs$v6_coord, C5_sig_eQTLs$v6_coord)
C4_C5<-intersect(C4_sig_eQTLs$v6_coord, C5_sig_eQTLs$v6_coord)
C1_C2_C3<-intersect(intersect(C1_sig_eQTLs$v6_coord, C2_sig_eQTLs$v6_coord), C3_sig_eQTLs$v6_coord)
C1_C2_C4<-intersect(intersect(C1_sig_eQTLs$v6_coord, C2_sig_eQTLs$v6_coord), C4_sig_eQTLs$v6_coord)
C1_C2_C5<-intersect(intersect(C1_sig_eQTLs$v6_coord, C2_sig_eQTLs$v6_coord), C5_sig_eQTLs$v6_coord)
C1_C3_C4<-intersect(intersect(C1_sig_eQTLs$v6_coord, C3_sig_eQTLs$v6_coord), C4_sig_eQTLs$v6_coord)
C1_C3_C5<-intersect(intersect(C1_sig_eQTLs$v6_coord, C3_sig_eQTLs$v6_coord), C5_sig_eQTLs$v6_coord)
C1_C4_C5<-intersect(intersect(C1_sig_eQTLs$v6_coord, C4_sig_eQTLs$v6_coord), C5_sig_eQTLs$v6_coord)
C2_C3_C4<-intersect(intersect(C2_sig_eQTLs$v6_coord, C3_sig_eQTLs$v6_coord), C4_sig_eQTLs$v6_coord)
C2_C3_C5<-intersect(intersect(C2_sig_eQTLs$v6_coord, C3_sig_eQTLs$v6_coord), C5_sig_eQTLs$v6_coord)
C2_C4_C5<-intersect(intersect(C2_sig_eQTLs$v6_coord, C4_sig_eQTLs$v6_coord), C5_sig_eQTLs$v6_coord)
C3_C4_C5<-intersect(intersect(C3_sig_eQTLs$v6_coord, C4_sig_eQTLs$v6_coord), C5_sig_eQTLs$v6_coord)
C1_C2_C3_C4<-intersect(intersect(C1_sig_eQTLs$v6_coord, C2_sig_eQTLs$v6_coord), 
                       intersect(C3_sig_eQTLs$v6_coord, C4_sig_eQTLs$v6_coord))
C1_C2_C3_C5<-intersect(intersect(C1_sig_eQTLs$v6_coord, C2_sig_eQTLs$v6_coord),
                       intersect(C3_sig_eQTLs$v6_coord, C5_sig_eQTLs$v6_coord))
C1_C2_C4_C5<-intersect(intersect(C1_sig_eQTLs$v6_coord, C2_sig_eQTLs$v6_coord),
                       intersect(C4_sig_eQTLs$v6_coord, C5_sig_eQTLs$v6_coord))
C1_C3_C4_C5<-intersect(intersect(C1_sig_eQTLs$v6_coord, C3_sig_eQTLs$v6_coord),
                       intersect(C4_sig_eQTLs$v6_coord, C5_sig_eQTLs$v6_coord))
C2_C3_C4_C5<-intersect(intersect(C2_sig_eQTLs$v6_coord, C3_sig_eQTLs$v6_coord),
                       intersect(C4_sig_eQTLs$v6_coord, C5_sig_eQTLs$v6_coord))
C1_C2_C3_C4_C5<-intersect(intersect(C1_sig_eQTLs$v6_coord, C2_sig_eQTLs$v6_coord), 
                          intersect(intersect(C3_sig_eQTLs$v6_coord, C4_sig_eQTLs$v6_coord), C5_sig_eQTLs$v6_coord))


D1_D2<-intersect(D1_sig_eQTLs$v6_coord, D2_sig_eQTLs$v6_coord)
D1_D3<-intersect(D1_sig_eQTLs$v6_coord, D3_sig_eQTLs$v6_coord)
D1_D4<-intersect(D1_sig_eQTLs$v6_coord, D4_sig_eQTLs$v6_coord)
D1_D5<-intersect(D1_sig_eQTLs$v6_coord, D5_sig_eQTLs$v6_coord)
D2_D3<-intersect(D2_sig_eQTLs$v6_coord, D3_sig_eQTLs$v6_coord)
D2_D4<-intersect(D2_sig_eQTLs$v6_coord, D4_sig_eQTLs$v6_coord)
D2_D5<-intersect(D2_sig_eQTLs$v6_coord, D5_sig_eQTLs$v6_coord)
D3_D4<-intersect(D3_sig_eQTLs$v6_coord, D4_sig_eQTLs$v6_coord)
D3_D5<-intersect(D3_sig_eQTLs$v6_coord, D5_sig_eQTLs$v6_coord)
D4_D5<-intersect(D4_sig_eQTLs$v6_coord, D5_sig_eQTLs$v6_coord)
D1_D2_D3<-intersect(intersect(D1_sig_eQTLs$v6_coord, D2_sig_eQTLs$v6_coord), D3_sig_eQTLs$v6_coord)
D1_D2_D4<-intersect(intersect(D1_sig_eQTLs$v6_coord, D2_sig_eQTLs$v6_coord), D4_sig_eQTLs$v6_coord)
D1_D2_D5<-intersect(intersect(D1_sig_eQTLs$v6_coord, D2_sig_eQTLs$v6_coord), D5_sig_eQTLs$v6_coord)
D1_D3_D4<-intersect(intersect(D1_sig_eQTLs$v6_coord, D3_sig_eQTLs$v6_coord), D4_sig_eQTLs$v6_coord)
D1_D3_D5<-intersect(intersect(D1_sig_eQTLs$v6_coord, D3_sig_eQTLs$v6_coord), D5_sig_eQTLs$v6_coord)
D1_D4_D5<-intersect(intersect(D1_sig_eQTLs$v6_coord, D4_sig_eQTLs$v6_coord), D5_sig_eQTLs$v6_coord)
D2_D3_D4<-intersect(intersect(D2_sig_eQTLs$v6_coord, D3_sig_eQTLs$v6_coord), D4_sig_eQTLs$v6_coord)
D2_D3_D5<-intersect(intersect(D2_sig_eQTLs$v6_coord, D3_sig_eQTLs$v6_coord), D5_sig_eQTLs$v6_coord)
D2_D4_D5<-intersect(intersect(D2_sig_eQTLs$v6_coord, D4_sig_eQTLs$v6_coord), D5_sig_eQTLs$v6_coord)
D3_D4_D5<-intersect(intersect(D3_sig_eQTLs$v6_coord, D4_sig_eQTLs$v6_coord), D5_sig_eQTLs$v6_coord)
D1_D2_D3_D4<-intersect(intersect(D1_sig_eQTLs$v6_coord, D2_sig_eQTLs$v6_coord), 
                       intersect(D3_sig_eQTLs$v6_coord, D4_sig_eQTLs$v6_coord))
D1_D2_D3_D5<-intersect(intersect(D1_sig_eQTLs$v6_coord, D2_sig_eQTLs$v6_coord),
                       intersect(D3_sig_eQTLs$v6_coord, D5_sig_eQTLs$v6_coord))
D1_D2_D4_D5<-intersect(intersect(D1_sig_eQTLs$v6_coord, D2_sig_eQTLs$v6_coord),
                       intersect(D4_sig_eQTLs$v6_coord, D5_sig_eQTLs$v6_coord))
D1_D3_D4_D5<-intersect(intersect(D1_sig_eQTLs$v6_coord, D3_sig_eQTLs$v6_coord),
                       intersect(D4_sig_eQTLs$v6_coord, D5_sig_eQTLs$v6_coord))
D2_D3_D4_D5<-intersect(intersect(D2_sig_eQTLs$v6_coord, D3_sig_eQTLs$v6_coord),
                       intersect(D4_sig_eQTLs$v6_coord, D5_sig_eQTLs$v6_coord))
D1_D2_D3_D4_D5<-intersect(intersect(D1_sig_eQTLs$v6_coord, D2_sig_eQTLs$v6_coord), 
                          intersect(intersect(D3_sig_eQTLs$v6_coord, D4_sig_eQTLs$v6_coord), D5_sig_eQTLs$v6_coord))


rearrange_snp_iterations <- list()
for(j in 1:nIter) {
  if(j%%50==0){
    print(paste("Rearranging iteration #", j, sep=" "))
  }
  temp_name <- paste("Iter", j, sep="_")
  rearrange_snp_iterations[[temp_name]] <- list(C1_resampled_SNPs_hitting_eQTLs[[j]],
                                            C2_resampled_SNPs_hitting_eQTLs[[j]],
                                            C3_resampled_SNPs_hitting_eQTLs[[j]],
                                            C4_resampled_SNPs_hitting_eQTLs[[j]],
                                            C5_resampled_SNPs_hitting_eQTLs[[j]],
                                            D1_resampled_SNPs_hitting_eQTLs[[j]],
                                            D2_resampled_SNPs_hitting_eQTLs[[j]],
                                            D3_resampled_SNPs_hitting_eQTLs[[j]],
                                            D4_resampled_SNPs_hitting_eQTLs[[j]],
                                            D5_resampled_SNPs_hitting_eQTLs[[j]])
  
  names(rearrange_snp_iterations[[j]]) <- resample_names
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

saveRDS(rearrange_snp_iterations, file="eQTL_hit_lists_including_overlap_from_resampled_SNPs.rds")

################
# OVERLAP AMONG GENES HIT BY EQTL LOCI
################

C1_C2<-intersect(C1_sig_eQTL_genes, C2_sig_eQTL_genes)
C1_C3<-intersect(C1_sig_eQTL_genes, C3_sig_eQTL_genes)
C1_C4<-intersect(C1_sig_eQTL_genes, C4_sig_eQTL_genes)
C1_C5<-intersect(C1_sig_eQTL_genes, C5_sig_eQTL_genes)
C2_C3<-intersect(C2_sig_eQTL_genes, C3_sig_eQTL_genes)
C2_C4<-intersect(C2_sig_eQTL_genes, C4_sig_eQTL_genes)
C2_C5<-intersect(C2_sig_eQTL_genes, C5_sig_eQTL_genes)
C3_C4<-intersect(C3_sig_eQTL_genes, C4_sig_eQTL_genes)
C3_C5<-intersect(C3_sig_eQTL_genes, C5_sig_eQTL_genes)
C4_C5<-intersect(C4_sig_eQTL_genes, C5_sig_eQTL_genes)
C1_C2_C3<-intersect(intersect(C1_sig_eQTL_genes, C2_sig_eQTL_genes), C3_sig_eQTL_genes)
C1_C2_C4<-intersect(intersect(C1_sig_eQTL_genes, C2_sig_eQTL_genes), C4_sig_eQTL_genes)
C1_C2_C5<-intersect(intersect(C1_sig_eQTL_genes, C2_sig_eQTL_genes), C5_sig_eQTL_genes)
C1_C3_C4<-intersect(intersect(C1_sig_eQTL_genes, C3_sig_eQTL_genes), C4_sig_eQTL_genes)
C1_C3_C5<-intersect(intersect(C1_sig_eQTL_genes, C3_sig_eQTL_genes), C5_sig_eQTL_genes)
C1_C4_C5<-intersect(intersect(C1_sig_eQTL_genes, C4_sig_eQTL_genes), C5_sig_eQTL_genes)
C2_C3_C4<-intersect(intersect(C2_sig_eQTL_genes, C3_sig_eQTL_genes), C4_sig_eQTL_genes)
C2_C3_C5<-intersect(intersect(C2_sig_eQTL_genes, C3_sig_eQTL_genes), C5_sig_eQTL_genes)
C2_C4_C5<-intersect(intersect(C2_sig_eQTL_genes, C4_sig_eQTL_genes), C5_sig_eQTL_genes)
C3_C4_C5<-intersect(intersect(C3_sig_eQTL_genes, C4_sig_eQTL_genes), C5_sig_eQTL_genes)
C1_C2_C3_C4<-intersect(intersect(C1_sig_eQTL_genes, C2_sig_eQTL_genes), 
                       intersect(C3_sig_eQTL_genes, C4_sig_eQTL_genes))
C1_C2_C3_C5<-intersect(intersect(C1_sig_eQTL_genes, C2_sig_eQTL_genes),
                       intersect(C3_sig_eQTL_genes, C5_sig_eQTL_genes))
C1_C2_C4_C5<-intersect(intersect(C1_sig_eQTL_genes, C2_sig_eQTL_genes),
                       intersect(C4_sig_eQTL_genes, C5_sig_eQTL_genes))
C1_C3_C4_C5<-intersect(intersect(C1_sig_eQTL_genes, C3_sig_eQTL_genes),
                       intersect(C4_sig_eQTL_genes, C5_sig_eQTL_genes))
C2_C3_C4_C5<-intersect(intersect(C2_sig_eQTL_genes, C3_sig_eQTL_genes),
                       intersect(C4_sig_eQTL_genes, C5_sig_eQTL_genes))
C1_C2_C3_C4_C5<-intersect(intersect(C1_sig_eQTL_genes, C2_sig_eQTL_genes), 
                          intersect(intersect(C3_sig_eQTL_genes, C4_sig_eQTL_genes), C5_sig_eQTL_genes))


D1_D2<-intersect(D1_sig_eQTL_genes, D2_sig_eQTL_genes)
D1_D3<-intersect(D1_sig_eQTL_genes, D3_sig_eQTL_genes)
D1_D4<-intersect(D1_sig_eQTL_genes, D4_sig_eQTL_genes)
D1_D5<-intersect(D1_sig_eQTL_genes, D5_sig_eQTL_genes)
D2_D3<-intersect(D2_sig_eQTL_genes, D3_sig_eQTL_genes)
D2_D4<-intersect(D2_sig_eQTL_genes, D4_sig_eQTL_genes)
D2_D5<-intersect(D2_sig_eQTL_genes, D5_sig_eQTL_genes)
D3_D4<-intersect(D3_sig_eQTL_genes, D4_sig_eQTL_genes)
D3_D5<-intersect(D3_sig_eQTL_genes, D5_sig_eQTL_genes)
D4_D5<-intersect(D4_sig_eQTL_genes, D5_sig_eQTL_genes)
D1_D2_D3<-intersect(intersect(D1_sig_eQTL_genes, D2_sig_eQTL_genes), D3_sig_eQTL_genes)
D1_D2_D4<-intersect(intersect(D1_sig_eQTL_genes, D2_sig_eQTL_genes), D4_sig_eQTL_genes)
D1_D2_D5<-intersect(intersect(D1_sig_eQTL_genes, D2_sig_eQTL_genes), D5_sig_eQTL_genes)
D1_D3_D4<-intersect(intersect(D1_sig_eQTL_genes, D3_sig_eQTL_genes), D4_sig_eQTL_genes)
D1_D3_D5<-intersect(intersect(D1_sig_eQTL_genes, D3_sig_eQTL_genes), D5_sig_eQTL_genes)
D1_D4_D5<-intersect(intersect(D1_sig_eQTL_genes, D4_sig_eQTL_genes), D5_sig_eQTL_genes)
D2_D3_D4<-intersect(intersect(D2_sig_eQTL_genes, D3_sig_eQTL_genes), D4_sig_eQTL_genes)
D2_D3_D5<-intersect(intersect(D2_sig_eQTL_genes, D3_sig_eQTL_genes), D5_sig_eQTL_genes)
D2_D4_D5<-intersect(intersect(D2_sig_eQTL_genes, D4_sig_eQTL_genes), D5_sig_eQTL_genes)
D3_D4_D5<-intersect(intersect(D3_sig_eQTL_genes, D4_sig_eQTL_genes), D5_sig_eQTL_genes)
D1_D2_D3_D4<-intersect(intersect(D1_sig_eQTL_genes, D2_sig_eQTL_genes), 
                       intersect(D3_sig_eQTL_genes, D4_sig_eQTL_genes))
D1_D2_D3_D5<-intersect(intersect(D1_sig_eQTL_genes, D2_sig_eQTL_genes),
                       intersect(D3_sig_eQTL_genes, D5_sig_eQTL_genes))
D1_D2_D4_D5<-intersect(intersect(D1_sig_eQTL_genes, D2_sig_eQTL_genes),
                       intersect(D4_sig_eQTL_genes, D5_sig_eQTL_genes))
D1_D3_D4_D5<-intersect(intersect(D1_sig_eQTL_genes, D3_sig_eQTL_genes),
                       intersect(D4_sig_eQTL_genes, D5_sig_eQTL_genes))
D2_D3_D4_D5<-intersect(intersect(D2_sig_eQTL_genes, D3_sig_eQTL_genes),
                       intersect(D4_sig_eQTL_genes, D5_sig_eQTL_genes))
D1_D2_D3_D4_D5<-intersect(intersect(D1_sig_eQTL_genes, D2_sig_eQTL_genes), 
                          intersect(intersect(D3_sig_eQTL_genes, D4_sig_eQTL_genes), D5_sig_eQTL_genes))




rearrange_eQTL_gene_iterations <- list()
for(j in 1:nIter) {
  if(j%%50==0){
    print(paste("Rearranging iteration #", j, sep=" "))
  }
  temp_name <- paste("Iter", j, sep="_")
  rearrange_eQTL_gene_iterations[[temp_name]] <- list(C1_resampled_eQTL_genes_with_position[[j]],
                                                C2_resampled_eQTL_genes_with_position[[j]],
                                                C3_resampled_eQTL_genes_with_position[[j]],
                                                C4_resampled_eQTL_genes_with_position[[j]],
                                                C5_resampled_eQTL_genes_with_position[[j]],
                                                D1_resampled_eQTL_genes_with_position[[j]],
                                                D2_resampled_eQTL_genes_with_position[[j]],
                                                D3_resampled_eQTL_genes_with_position[[j]],
                                                D4_resampled_eQTL_genes_with_position[[j]],
                                                D5_resampled_eQTL_genes_with_position[[j]])
  
  names(rearrange_eQTL_gene_iterations[[j]]) <- resample_names
  #temp_allC <- 
  #rearrange_eQTL_gene_iterations[[j]][["allC"]] <- union(union(union(rearrange_eQTL_gene_iterations[[j]][[1]], 
  #                                                           rearrange_eQTL_gene_iterations[[j]][[2]]), 
  #                                                     rearrange_eQTL_gene_iterations[[j]][[3]]), 
  #                                               union(rearrange_eQTL_gene_iterations[[j]][[4]], 
  #                                                     rearrange_eQTL_gene_iterations[[j]][[5]]))
  D_names <- c("D1", "D2", "D3", "D4", "D5")
  #   for(i in 1:5){
  #     temp_D <- D_names[i]
  #     #temp_noC_name <- paste(temp_D, "noC", sep="_")
  #     rearrange_eQTL_gene_iterations[[j]][[temp_D]] <- rearrange_eQTL_gene_iterations[[j]][[temp_D]]
  #   }
  combn2 <- combn(D_names, 2)
  combn2a <- apply(combn2, 2, paste, collapse="_")
  for(i in 1:ncol(combn2)){
    temp_Da <- rearrange_eQTL_gene_iterations[[j]][[combn2[1,i]]]
    temp_Db <- rearrange_eQTL_gene_iterations[[j]][[combn2[2,i]]]
    temp_combn_name <- combn2a[i]
    rearrange_eQTL_gene_iterations[[j]][[temp_combn_name]] <- intersect(temp_Da, temp_Db)
  }
  combn3 <- combn(D_names, 3)
  combn3a <- apply(combn3, 2, paste, collapse="_")
  for(i in 1:ncol(combn3)){
    temp_Da <- rearrange_eQTL_gene_iterations[[j]][[combn3[1,i]]]
    temp_Db <- rearrange_eQTL_gene_iterations[[j]][[combn3[2,i]]]
    temp_Dc <- rearrange_eQTL_gene_iterations[[j]][[combn3[3,i]]]
    temp_combn_name <- combn3a[i]
    rearrange_eQTL_gene_iterations[[j]][[temp_combn_name]] <- intersect(temp_Da, intersect(temp_Db, temp_Dc))
  } 
  combn4 <- combn(D_names, 4)
  combn4a <- apply(combn4, 2, paste, collapse="_")
  for(i in 1:ncol(combn4)){
    temp_Da <- rearrange_eQTL_gene_iterations[[j]][[combn4[1,i]]]
    temp_Db <- rearrange_eQTL_gene_iterations[[j]][[combn4[2,i]]]
    temp_Dc <- rearrange_eQTL_gene_iterations[[j]][[combn4[3,i]]]
    temp_Dd <- rearrange_eQTL_gene_iterations[[j]][[combn4[4,i]]]
    temp_combn_name <- combn4a[i]
    rearrange_eQTL_gene_iterations[[j]][[temp_combn_name]] <- intersect(intersect(temp_Da, temp_Db), intersect(temp_Dc, temp_Dd))
  } 
  combn5 <- "D1_D2_D3_D4_D5"
  rearrange_eQTL_gene_iterations[[j]][[combn5]] <- intersect(intersect(rearrange_eQTL_gene_iterations[[j]][["D1"]], 
                                                                 intersect(rearrange_eQTL_gene_iterations[[j]][["D2"]], 
                                                                           rearrange_eQTL_gene_iterations[[j]][["D3"]])), 
                                                       intersect(rearrange_eQTL_gene_iterations[[j]][["D4"]], 
                                                                 rearrange_eQTL_gene_iterations[[j]][["D5"]]))
  C_names <- c("C1", "C2", "C3", "C4", "C5")
  #   for(i in 1:5){
  #     temp_C <- C_names[i]
  #   }
  combn2 <- combn(C_names, 2)
  combn2a <- apply(combn2, 2, paste, collapse="_")
  for(i in 1:ncol(combn2)){
    temp_Ca <- rearrange_eQTL_gene_iterations[[j]][[combn2[1,i]]]
    temp_Cb <- rearrange_eQTL_gene_iterations[[j]][[combn2[2,i]]]
    temp_combn_name <- combn2a[i]
    rearrange_eQTL_gene_iterations[[j]][[temp_combn_name]] <- intersect(temp_Ca, temp_Cb)
  }
  combn3 <- combn(C_names, 3)
  combn3a <- apply(combn3, 2, paste, collapse="_")
  for(i in 1:ncol(combn3)){
    temp_Ca <- rearrange_eQTL_gene_iterations[[j]][[combn3[1,i]]]
    temp_Cb <- rearrange_eQTL_gene_iterations[[j]][[combn3[2,i]]]
    temp_Cc <- rearrange_eQTL_gene_iterations[[j]][[combn3[3,i]]]
    temp_combn_name <- combn3a[i]
    rearrange_eQTL_gene_iterations[[j]][[temp_combn_name]] <- intersect(temp_Ca, intersect(temp_Cb, temp_Cc))
  } 
  combn4 <- combn(C_names, 4)
  combn4a <- apply(combn4, 2, paste, collapse="_")
  for(i in 1:ncol(combn4)){
    temp_Ca <- rearrange_eQTL_gene_iterations[[j]][[combn4[1,i]]]
    temp_Cb <- rearrange_eQTL_gene_iterations[[j]][[combn4[2,i]]]
    temp_Cc <- rearrange_eQTL_gene_iterations[[j]][[combn4[3,i]]]
    temp_Cd <- rearrange_eQTL_gene_iterations[[j]][[combn4[4,i]]]
    temp_combn_name <- combn4a[i]
    rearrange_eQTL_gene_iterations[[j]][[temp_combn_name]] <- intersect(intersect(temp_Ca, temp_Cb), intersect(temp_Cc, temp_Cd))
  } 
  combn5 <- "C1_C2_C3_C4_C5"
  rearrange_eQTL_gene_iterations[[j]][[combn5]] <- intersect(intersect(rearrange_eQTL_gene_iterations[[j]][["C1"]], 
                                                                 intersect(rearrange_eQTL_gene_iterations[[j]][["C2"]], 
                                                                           rearrange_eQTL_gene_iterations[[j]][["C3"]])), 
                                                       intersect(rearrange_eQTL_gene_iterations[[j]][["C4"]], 
                                                                 rearrange_eQTL_gene_iterations[[j]][["C5"]]))
  
  
  
}

saveRDS(rearrange_eQTL_gene_iterations, file="eQTL_hit_lists_including_overlap_from_resampled_SNPs.rds")



pdf("Histograms of simulated eQTL gene overlap, resampling with position.pdf", width=28, height=16)
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
  temp_all_list <- lapply(rearrange_eQTL_gene_iterations, "[[", i+10)
  temp_dist <- sapply(temp_all_list, length)
  
  multiplier <- hist(temp_dist, plot=FALSE)$counts / hist(temp_dist, plot=FALSE)$density
  mydensity <- density(temp_dist)
  mydensity$y <- mydensity$y * multiplier[1] 
  
  hist(temp_dist, main="", ylim=c(0, 380), xlim=c(0, max(temp_dist)*1.5),
       xlab=paste("No. eQTL genes hit in", temp_combn, "overlap", sep=" "))
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




# Do any of the eQTL SNPs hit known candidate genes? (more of a visualisation)
# (visualise on the network diagram made for literature candidates)
# (also visualise on the overlap PPI networks?)

# Are the gene sets hit by eQTL enriched for particular functions/timepoints/body parts?





levels(as.factor(D2_sig_eQTLs$Gene.Symbol))
