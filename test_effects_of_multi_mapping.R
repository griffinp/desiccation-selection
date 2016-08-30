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

setwd("~/Documents/Drosophila Selection Experiment/snp_and_gene_lists/")

Sample_code <- c("C1", "C2", "C3", "C4", "C5",
                 "D1", "D2", "D3", "D4", "D5")

for(i in Sample_code){
  #temp_sig_ex <- get(paste(i, "_sig_ex", sep=""))
  #temp_sig_ex_cols <- temp_sig_ex[,1:2]
  #temp_sig_ex_chrpos <- paste(temp_sig_ex_cols$chr, temp_sig_ex_cols$pos)
  #temp_sig_gene_counts <- get(paste(i, "_sig_gene_counts", sep=""))
  #temp_counts_file_name <- paste(i, "_sig_gene_mapping_counts.txt", sep="")
  temp_single_mapping_SNP_list <- paste(i, "_single_mapping_SNPs", sep="")
  temp_single_mapping_SNP_file_name <- paste(i, "_single_mapping_SNPs.txt", sep="")
  assign(temp_single_mapping_SNP_list, paste(read.table(temp_single_mapping_SNP_file_name,
                                                  sep="\t", header=TRUE)[,1],
                                             read.table(temp_single_mapping_SNP_file_name,
                                                              sep="\t", header=TRUE)[,2]))
  temp_multi_mapping_SNP_list <- paste(i, "_multi_mapping_SNPs", sep="")
  temp_multi_mapping_SNP_file_name <- paste(i, "_multi_mapping_SNPs.txt", sep="")
  assign(temp_multi_mapping_SNP_list, paste(read.table(temp_multi_mapping_SNP_file_name,
                                                        sep="\t", header=TRUE)[,1],
                                             read.table(temp_multi_mapping_SNP_file_name,
                                                        sep="\t", header=TRUE)[,2]))
  temp_single_mapping_gene_list <- paste(i, "_single_mapped_genes", sep="")
  temp_single_mapped_genes_file_name <- paste(i, "_single_mapped_genes.txt", sep="")
  assign(temp_single_mapping_gene_list, read.table(temp_single_mapped_genes_file_name,
                                                  header=FALSE, stringsAsFactors=FALSE)[,1])
  temp_multi_mapping_gene_list <- paste(i, "_multi_mapped_genes", sep="")
  temp_multi_mapped_genes_file_name <- paste(i, "_multi_mapped_genes.txt", sep="")
  assign(temp_multi_mapping_gene_list, read.table(temp_multi_mapped_genes_file_name,
                                                   header=FALSE, stringsAsFactors=FALSE)[,1])
}

### Extract overlaps for single-mapping genes
C_single_mapping_overlap <- list()
for(i in Sample_code[1:5]){
  assign(paste(i, "_sig_genes", sep=""), get(paste(i, "_single_mapped_genes", sep="")))
}


C_single_mapping_overlap$C1_C2<-intersect(C1_sig_genes, C2_sig_genes)
C_single_mapping_overlap$C1_C3<-intersect(C1_sig_genes, C3_sig_genes)
C_single_mapping_overlap$C1_C4<-intersect(C1_sig_genes, C4_sig_genes)
C_single_mapping_overlap$C1_C5<-intersect(C1_sig_genes, C5_sig_genes)
C_single_mapping_overlap$C2_C3<-intersect(C2_sig_genes, C3_sig_genes)
C_single_mapping_overlap$C2_C4<-intersect(C2_sig_genes, C4_sig_genes)
C_single_mapping_overlap$C2_C5<-intersect(C2_sig_genes, C5_sig_genes)
C_single_mapping_overlap$C3_C4<-intersect(C3_sig_genes, C4_sig_genes)
C_single_mapping_overlap$C3_C5<-intersect(C3_sig_genes, C5_sig_genes)
C_single_mapping_overlap$C4_C5<-intersect(C4_sig_genes, C5_sig_genes)
C_single_mapping_overlap$C1_C2_C3<-intersect(intersect(C1_sig_genes, C2_sig_genes), C3_sig_genes)
C_single_mapping_overlap$C1_C2_C4<-intersect(intersect(C1_sig_genes, C2_sig_genes), C4_sig_genes)
C_single_mapping_overlap$C1_C2_C5<-intersect(intersect(C1_sig_genes, C2_sig_genes), C5_sig_genes)
C_single_mapping_overlap$C1_C3_C4<-intersect(intersect(C1_sig_genes, C3_sig_genes), C4_sig_genes)
C_single_mapping_overlap$C1_C3_C5<-intersect(intersect(C1_sig_genes, C3_sig_genes), C5_sig_genes)
C_single_mapping_overlap$C1_C4_C5<-intersect(intersect(C1_sig_genes, C4_sig_genes), C5_sig_genes)
C_single_mapping_overlap$C2_C3_C4<-intersect(intersect(C2_sig_genes, C3_sig_genes), C4_sig_genes)
C_single_mapping_overlap$C2_C3_C5<-intersect(intersect(C2_sig_genes, C3_sig_genes), C5_sig_genes)
C_single_mapping_overlap$C2_C4_C5<-intersect(intersect(C2_sig_genes, C4_sig_genes), C5_sig_genes)
C_single_mapping_overlap$C3_C4_C5<-intersect(intersect(C3_sig_genes, C4_sig_genes), C5_sig_genes)
C_single_mapping_overlap$C1_C2_C3_C4<-intersect(intersect(C1_sig_genes, C2_sig_genes), 
                       intersect(C3_sig_genes, C4_sig_genes))
C_single_mapping_overlap$C1_C2_C3_C5<-intersect(intersect(C1_sig_genes, C2_sig_genes),
                       intersect(C3_sig_genes, C5_sig_genes))
C_single_mapping_overlap$C1_C2_C4_C5<-intersect(intersect(C1_sig_genes, C2_sig_genes),
                       intersect(C4_sig_genes, C5_sig_genes))
C_single_mapping_overlap$C1_C3_C4_C5<-intersect(intersect(C1_sig_genes, C3_sig_genes),
                       intersect(C4_sig_genes, C5_sig_genes))
C_single_mapping_overlap$C2_C3_C4_C5<-intersect(intersect(C2_sig_genes, C3_sig_genes),
                       intersect(C4_sig_genes, C5_sig_genes))
C_single_mapping_overlap$C1_C2_C3_C4_C5<-intersect(intersect(C1_sig_genes, C2_sig_genes), 
                          intersect(intersect(C3_sig_genes, C4_sig_genes), C5_sig_genes))

### Extract overlaps for multi-mapping SNPs
C_multi_mapping_overlap <- list()
for(i in Sample_code[1:5]){
  assign(paste(i, "_sig_genes", sep=""), get(paste(i, "_multi_mapped_genes", sep="")))
}


C_multi_mapping_overlap$C1_C2<-intersect(C1_sig_genes, C2_sig_genes)
C_multi_mapping_overlap$C1_C3<-intersect(C1_sig_genes, C3_sig_genes)
C_multi_mapping_overlap$C1_C4<-intersect(C1_sig_genes, C4_sig_genes)
C_multi_mapping_overlap$C1_C5<-intersect(C1_sig_genes, C5_sig_genes)
C_multi_mapping_overlap$C2_C3<-intersect(C2_sig_genes, C3_sig_genes)
C_multi_mapping_overlap$C2_C4<-intersect(C2_sig_genes, C4_sig_genes)
C_multi_mapping_overlap$C2_C5<-intersect(C2_sig_genes, C5_sig_genes)
C_multi_mapping_overlap$C3_C4<-intersect(C3_sig_genes, C4_sig_genes)
C_multi_mapping_overlap$C3_C5<-intersect(C3_sig_genes, C5_sig_genes)
C_multi_mapping_overlap$C4_C5<-intersect(C4_sig_genes, C5_sig_genes)
C_multi_mapping_overlap$C1_C2_C3<-intersect(intersect(C1_sig_genes, C2_sig_genes), C3_sig_genes)
C_multi_mapping_overlap$C1_C2_C4<-intersect(intersect(C1_sig_genes, C2_sig_genes), C4_sig_genes)
C_multi_mapping_overlap$C1_C2_C5<-intersect(intersect(C1_sig_genes, C2_sig_genes), C5_sig_genes)
C_multi_mapping_overlap$C1_C3_C4<-intersect(intersect(C1_sig_genes, C3_sig_genes), C4_sig_genes)
C_multi_mapping_overlap$C1_C3_C5<-intersect(intersect(C1_sig_genes, C3_sig_genes), C5_sig_genes)
C_multi_mapping_overlap$C1_C4_C5<-intersect(intersect(C1_sig_genes, C4_sig_genes), C5_sig_genes)
C_multi_mapping_overlap$C2_C3_C4<-intersect(intersect(C2_sig_genes, C3_sig_genes), C4_sig_genes)
C_multi_mapping_overlap$C2_C3_C5<-intersect(intersect(C2_sig_genes, C3_sig_genes), C5_sig_genes)
C_multi_mapping_overlap$C2_C4_C5<-intersect(intersect(C2_sig_genes, C4_sig_genes), C5_sig_genes)
C_multi_mapping_overlap$C3_C4_C5<-intersect(intersect(C3_sig_genes, C4_sig_genes), C5_sig_genes)
C_multi_mapping_overlap$C1_C2_C3_C4<-intersect(intersect(C1_sig_genes, C2_sig_genes), 
                                              intersect(C3_sig_genes, C4_sig_genes))
C_multi_mapping_overlap$C1_C2_C3_C5<-intersect(intersect(C1_sig_genes, C2_sig_genes),
                                              intersect(C3_sig_genes, C5_sig_genes))
C_multi_mapping_overlap$C1_C2_C4_C5<-intersect(intersect(C1_sig_genes, C2_sig_genes),
                                              intersect(C4_sig_genes, C5_sig_genes))
C_multi_mapping_overlap$C1_C3_C4_C5<-intersect(intersect(C1_sig_genes, C3_sig_genes),
                                              intersect(C4_sig_genes, C5_sig_genes))
C_multi_mapping_overlap$C2_C3_C4_C5<-intersect(intersect(C2_sig_genes, C3_sig_genes),
                                              intersect(C4_sig_genes, C5_sig_genes))
C_multi_mapping_overlap$C1_C2_C3_C4_C5<-intersect(intersect(C1_sig_genes, C2_sig_genes), 
                                                 intersect(intersect(C3_sig_genes, C4_sig_genes), C5_sig_genes))

### Extract overlaps for single-mapping SNPs in D replicates
D_single_mapping_overlap <- list()
for(i in Sample_code[6:10]){
  assign(paste(i, "_sig_genes", sep=""), get(paste(i, "_single_mapped_genes", sep="")))
}

D_single_mapping_overlap$D1_D2<-intersect(D1_sig_genes, D2_sig_genes)
D_single_mapping_overlap$D1_D3<-intersect(D1_sig_genes, D3_sig_genes)
D_single_mapping_overlap$D1_D4<-intersect(D1_sig_genes, D4_sig_genes)
D_single_mapping_overlap$D1_D5<-intersect(D1_sig_genes, D5_sig_genes)
D_single_mapping_overlap$D2_D3<-intersect(D2_sig_genes, D3_sig_genes)
D_single_mapping_overlap$D2_D4<-intersect(D2_sig_genes, D4_sig_genes)
D_single_mapping_overlap$D2_D5<-intersect(D2_sig_genes, D5_sig_genes)
D_single_mapping_overlap$D3_D4<-intersect(D3_sig_genes, D4_sig_genes)
D_single_mapping_overlap$D3_D5<-intersect(D3_sig_genes, D5_sig_genes)
D_single_mapping_overlap$D4_D5<-intersect(D4_sig_genes, D5_sig_genes)
D_single_mapping_overlap$D1_D2_D3<-intersect(intersect(D1_sig_genes, D2_sig_genes), D3_sig_genes)
D_single_mapping_overlap$D1_D2_D4<-intersect(intersect(D1_sig_genes, D2_sig_genes), D4_sig_genes)
D_single_mapping_overlap$D1_D2_D5<-intersect(intersect(D1_sig_genes, D2_sig_genes), D5_sig_genes)
D_single_mapping_overlap$D1_D3_D4<-intersect(intersect(D1_sig_genes, D3_sig_genes), D4_sig_genes)
D_single_mapping_overlap$D1_D3_D5<-intersect(intersect(D1_sig_genes, D3_sig_genes), D5_sig_genes)
D_single_mapping_overlap$D1_D4_D5<-intersect(intersect(D1_sig_genes, D4_sig_genes), D5_sig_genes)
D_single_mapping_overlap$D2_D3_D4<-intersect(intersect(D2_sig_genes, D3_sig_genes), D4_sig_genes)
D_single_mapping_overlap$D2_D3_D5<-intersect(intersect(D2_sig_genes, D3_sig_genes), D5_sig_genes)
D_single_mapping_overlap$D2_D4_D5<-intersect(intersect(D2_sig_genes, D4_sig_genes), D5_sig_genes)
D_single_mapping_overlap$D3_D4_D5<-intersect(intersect(D3_sig_genes, D4_sig_genes), D5_sig_genes)
D_single_mapping_overlap$D1_D2_D3_D4<-intersect(intersect(D1_sig_genes, D2_sig_genes), 
                       intersect(D3_sig_genes, D4_sig_genes))
D_single_mapping_overlap$D1_D2_D3_D5<-intersect(intersect(D1_sig_genes, D2_sig_genes),
                       intersect(D3_sig_genes, D5_sig_genes))
D_single_mapping_overlap$D1_D2_D4_D5<-intersect(intersect(D1_sig_genes, D2_sig_genes),
                       intersect(D4_sig_genes, D5_sig_genes))
D_single_mapping_overlap$D1_D3_D4_D5<-intersect(intersect(D1_sig_genes, D3_sig_genes),
                       intersect(D4_sig_genes, D5_sig_genes))
D_single_mapping_overlap$D2_D3_D4_D5<-intersect(intersect(D2_sig_genes, D3_sig_genes),
                       intersect(D4_sig_genes, D5_sig_genes))
D_single_mapping_overlap$D1_D2_D3_D4_D5<-intersect(intersect(D1_sig_genes, D2_sig_genes), 
                          intersect(intersect(D3_sig_genes, D4_sig_genes), D5_sig_genes))


### Extract overlaps for multi-mapping SNPs in D replicates
D_multi_mapping_overlap <- list()
for(i in Sample_code[6:10]){
  assign(paste(i, "_sig_genes", sep=""), get(paste(i, "_multi_mapped_genes", sep="")))
}

D_multi_mapping_overlap$D1_D2<-intersect(D1_sig_genes, D2_sig_genes)
D_multi_mapping_overlap$D1_D3<-intersect(D1_sig_genes, D3_sig_genes)
D_multi_mapping_overlap$D1_D4<-intersect(D1_sig_genes, D4_sig_genes)
D_multi_mapping_overlap$D1_D5<-intersect(D1_sig_genes, D5_sig_genes)
D_multi_mapping_overlap$D2_D3<-intersect(D2_sig_genes, D3_sig_genes)
D_multi_mapping_overlap$D2_D4<-intersect(D2_sig_genes, D4_sig_genes)
D_multi_mapping_overlap$D2_D5<-intersect(D2_sig_genes, D5_sig_genes)
D_multi_mapping_overlap$D3_D4<-intersect(D3_sig_genes, D4_sig_genes)
D_multi_mapping_overlap$D3_D5<-intersect(D3_sig_genes, D5_sig_genes)
D_multi_mapping_overlap$D4_D5<-intersect(D4_sig_genes, D5_sig_genes)
D_multi_mapping_overlap$D1_D2_D3<-intersect(intersect(D1_sig_genes, D2_sig_genes), D3_sig_genes)
D_multi_mapping_overlap$D1_D2_D4<-intersect(intersect(D1_sig_genes, D2_sig_genes), D4_sig_genes)
D_multi_mapping_overlap$D1_D2_D5<-intersect(intersect(D1_sig_genes, D2_sig_genes), D5_sig_genes)
D_multi_mapping_overlap$D1_D3_D4<-intersect(intersect(D1_sig_genes, D3_sig_genes), D4_sig_genes)
D_multi_mapping_overlap$D1_D3_D5<-intersect(intersect(D1_sig_genes, D3_sig_genes), D5_sig_genes)
D_multi_mapping_overlap$D1_D4_D5<-intersect(intersect(D1_sig_genes, D4_sig_genes), D5_sig_genes)
D_multi_mapping_overlap$D2_D3_D4<-intersect(intersect(D2_sig_genes, D3_sig_genes), D4_sig_genes)
D_multi_mapping_overlap$D2_D3_D5<-intersect(intersect(D2_sig_genes, D3_sig_genes), D5_sig_genes)
D_multi_mapping_overlap$D2_D4_D5<-intersect(intersect(D2_sig_genes, D4_sig_genes), D5_sig_genes)
D_multi_mapping_overlap$D3_D4_D5<-intersect(intersect(D3_sig_genes, D4_sig_genes), D5_sig_genes)
D_multi_mapping_overlap$D1_D2_D3_D4<-intersect(intersect(D1_sig_genes, D2_sig_genes), 
                                                intersect(D3_sig_genes, D4_sig_genes))
D_multi_mapping_overlap$D1_D2_D3_D5<-intersect(intersect(D1_sig_genes, D2_sig_genes),
                                                intersect(D3_sig_genes, D5_sig_genes))
D_multi_mapping_overlap$D1_D2_D4_D5<-intersect(intersect(D1_sig_genes, D2_sig_genes),
                                                intersect(D4_sig_genes, D5_sig_genes))
D_multi_mapping_overlap$D1_D3_D4_D5<-intersect(intersect(D1_sig_genes, D3_sig_genes),
                                                intersect(D4_sig_genes, D5_sig_genes))
D_multi_mapping_overlap$D2_D3_D4_D5<-intersect(intersect(D2_sig_genes, D3_sig_genes),
                                                intersect(D4_sig_genes, D5_sig_genes))
D_multi_mapping_overlap$D1_D2_D3_D4_D5<-intersect(intersect(D1_sig_genes, D2_sig_genes), 
                                                   intersect(intersect(D3_sig_genes, D4_sig_genes), D5_sig_genes))



C_single_overlap_lengths <- sapply(C_single_mapping_overlap, length)
C_multi_overlap_lengths <- sapply(C_multi_mapping_overlap, length)
D_single_overlap_lengths <- sapply(D_single_mapping_overlap, length)
D_multi_overlap_lengths <- sapply(D_multi_mapping_overlap, length)


pdf(file="Gene-level overlap caused by single- versus multi-mapping SNPs.pdf",
    width=9, height=5)
par(mfrow=c(1,2))
plot(C_multi_overlap_lengths, C_single_overlap_lengths, 
     col=rgb(100, 100, 100, 50, maxColorValue=255), 
     xlim=c(0, 40), ylim=c(0, 10), pch=19, cex=0.7,
     xlab="Gene overlap caused by multi-mapping SNPs",
     ylab="Gene overlap caused by single-mapping SNPs")
lines(x=seq(0, 150, by=1), y=seq(0, 150, by=1), col="grey", lty=2)
text(C_multi_overlap_lengths, C_single_overlap_lengths, 
     labels=names(C_multi_overlap_lengths), cex=0.3)
text(2, 9.3, labels="A", cex=3)

# plot(D_multi_overlap_lengths, D_single_overlap_lengths, col="white",
#      xlim=c(0, 1200), ylim=c(0, 1200))
# lines(x=seq(0, 1500, by=1), y=seq(0, 1500, by=1), col="grey")
# text(D_multi_overlap_lengths, D_single_overlap_lengths, 
#      labels=names(D_multi_overlap_lengths), cex=0.3)

plot(log10(D_multi_overlap_lengths), log10(D_single_overlap_lengths), 
     col=rgb(100, 100, 100, 50, maxColorValue=255), 
     xlim=c(0.5, 3.2), ylim=c(0.5,3.2), pch=19, cex=0.7,
     xlab="log10(gene overlap caused by multi-mapping SNPs)",
     ylab="log10(gene overlap caused by single-mapping SNPs)")
lines(x=seq(0, 100, by=1), y=seq(0, 100, by=1), col="grey", lty=2)
linrel <- lm(log10(D_single_overlap_lengths)~log10(D_multi_overlap_lengths))
linrel_pred <- predict(linrel, 
                       newdata=data.frame(D_multi_overlap_lengths=seq(1, 5000, by=10)))
lines(x=log10(seq(1, 5000, by=10)), y=linrel_pred, col="red")
text(log10(D_multi_overlap_lengths), log10(D_single_overlap_lengths), 
     labels=names(D_multi_overlap_lengths), cex=0.3)
text(0.6, 3, labels="B", cex=3)

dev.off()
