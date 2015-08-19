# Select only significant SNPs and write to new file

alpha_value <- 1e-06

Sample_code <- c("C1", "C2", "C3", "C4", "C5",
                 "D1", "D2", "D3", "D4", "D5")

Comparison <- c("MB_C1.13", "MB_C2.13", "MB_C3.13", "MB_C4.13", "MB_C5.13",
               "MB_D1.13", "MB_D2.13", "MB_D3.13", "MB_D4.13", "MB_D5.13")
#obtained the no. snps by counting lines in the input files
#(e.g. MB_C1_noN_reduced_150818.mpileup.sync)

No_SNPs <- c(2888628, 2814924, 
             2806047, 2797111, 
             2797155, 2815187, 
             2712254, 2889638, 
             2984067, 2803570)
Alpha <- rep(alpha_value, times=10)
Threshold <- Alpha/No_SNPs

bonferroni_thresholds <- data.frame(Sample_code, Comparison, No_SNPs, Alpha, Threshold)

setwd("/Users/pgriffin/Documents/Drosophila Selection Experiment/pileup_and_sync_files")

read_all_scaffolds <- function(sample_code){
  output_name <- paste(sample_code, "_all", sep="")
  temp_input_name <- paste("MB_", sample_code, "_noN_reduced_150818.mpileup.sync", sep="")
  temp_read <- read.table(temp_input_name, header=FALSE)
  assign(temp_input_name, temp_read)
  return(get(temp_input_name))
}

# read_all_scaffolds <- function(sample_code, date){
#   output_name <- paste("MB_", sample_code, "_all", sep="")
# 
#   for(chr in c("2L", "2R", "3L", "3R", "X")){
#     temp_input_name <- paste("MB_", sample_code, "_drift_SNP_test_", date, "_", chr, ".txt", sep="")
#     temp_read <- read.table(temp_input_name, header=TRUE)
#     if(chr=="2L"){
#       output_table <- temp_read
#     } else {
#       output_table <- rbind(output_table, temp_read)
#     }
#   }
#   assign(temp_input_name, output_table)
#   return(get(temp_input_name))
# }

C1_all <- read_all_scaffolds('C1')
C2_all <- read_all_scaffolds('C2')
C3_all <- read_all_scaffolds('C3')
C4_all <- read_all_scaffolds('C4')
C5_all <- read_all_scaffolds('C5')
D1_all <- read_all_scaffolds('D1')
D2_all <- read_all_scaffolds('D2')
D3_all <- read_all_scaffolds('D3')
D4_all <- read_all_scaffolds('D4')
D5_all <- read_all_scaffolds('D5')


C1_all_sig <- subset(C1_all, C1_all$V6<bonferroni_thresholds[bonferroni_thresholds$Sample_code=="C1",'Threshold'])
C2_all_sig <- subset(C2_all, C2_all$V6<bonferroni_thresholds[bonferroni_thresholds$Sample_code=="C2",'Threshold'])
C3_all_sig <- subset(C3_all, C3_all$V6<bonferroni_thresholds[bonferroni_thresholds$Sample_code=="C3",'Threshold'])
C4_all_sig <- subset(C4_all, C4_all$V6<bonferroni_thresholds[bonferroni_thresholds$Sample_code=="C4",'Threshold'])
C5_all_sig <- subset(C5_all, C5_all$V6<bonferroni_thresholds[bonferroni_thresholds$Sample_code=="C5",'Threshold'])
D1_all_sig <- subset(D1_all, D1_all$V6<bonferroni_thresholds[bonferroni_thresholds$Sample_code=="D1",'Threshold'])
D2_all_sig <- subset(D2_all, D2_all$V6<bonferroni_thresholds[bonferroni_thresholds$Sample_code=="D2",'Threshold'])
D3_all_sig <- subset(D3_all, D3_all$V6<bonferroni_thresholds[bonferroni_thresholds$Sample_code=="D3",'Threshold'])
D4_all_sig <- subset(D4_all, D4_all$V6<bonferroni_thresholds[bonferroni_thresholds$Sample_code=="D4",'Threshold'])
D5_all_sig <- subset(D5_all, D5_all$V6<bonferroni_thresholds[bonferroni_thresholds$Sample_code=="D5",'Threshold'])


write.table(C1_all[,1:2], file="/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/C1_all_SNPs.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(C2_all[,1:2], file="/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/C2_all_SNPs.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(C3_all[,1:2], file="/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/C3_all_SNPs.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(C4_all[,1:2], file="/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/C4_all_SNPs.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(C5_all[,1:2], file="/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/C5_all_SNPs.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

write.table(D1_all[,1:2], file="/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/D1_all_SNPs.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(D2_all[,1:2], file="/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/D2_all_SNPs.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(D3_all[,1:2], file="/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/D3_all_SNPs.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(D4_all[,1:2], file="/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/D4_all_SNPs.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(D5_all[,1:2], file="/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/D5_all_SNPs.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

# write.table(C1_all, file="/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/C1_all_SNPs_sync.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
# write.table(C2_all, file="/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/C2_all_SNPs_sync.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
# write.table(C3_all, file="/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/C3_all_SNPs_sync.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
# write.table(C4_all, file="/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/C4_all_SNPs_sync.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
# write.table(C5_all, file="/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/C5_all_SNPs_sync.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
# 
# write.table(D1_all, file="/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/D1_all_SNPs_sync.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
# write.table(D2_all, file="/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/D2_all_SNPs_sync.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
# write.table(D3_all, file="/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/D3_all_SNPs_sync.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
# write.table(D4_all, file="/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/D4_all_SNPs_sync.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
# write.table(D5_all, file="/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/D5_all_SNPs_sync.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
# 


write.table(C1_all_sig[,1:2], file="/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/C1_sig_SNPs.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(C2_all_sig[,1:2], file="/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/C2_sig_SNPs.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(C3_all_sig[,1:2], file="/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/C3_sig_SNPs.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(C4_all_sig[,1:2], file="/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/C4_sig_SNPs.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(C5_all_sig[,1:2], file="/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/C5_sig_SNPs.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

write.table(D1_all_sig[,1:2], file="/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/D1_sig_SNPs.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(D2_all_sig[,1:2], file="/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/D2_sig_SNPs.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(D3_all_sig[,1:2], file="/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/D3_sig_SNPs.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(D4_all_sig[,1:2], file="/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/D4_sig_SNPs.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(D5_all_sig[,1:2], file="/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/D5_sig_SNPs.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

write.table(C1_all_sig, file="/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/C1_sig_SNPs_sync.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(C2_all_sig, file="/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/C2_sig_SNPs_sync.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(C3_all_sig, file="/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/C3_sig_SNPs_sync.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(C4_all_sig, file="/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/C4_sig_SNPs_sync.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(C5_all_sig, file="/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/C5_sig_SNPs_sync.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

write.table(D1_all_sig, file="/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/D1_sig_SNPs_sync.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(D2_all_sig, file="/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/D2_sig_SNPs_sync.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(D3_all_sig, file="/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/D3_sig_SNPs_sync.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(D4_all_sig, file="/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/D4_sig_SNPs_sync.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(D5_all_sig, file="/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/D5_sig_SNPs_sync.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

