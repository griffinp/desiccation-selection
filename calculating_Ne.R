library(ggplot2)

Sample_code <- c("C1", "C2", "C3", "C4", "C5",
                 "D1", "D2", "D3", "D4", "D5")

setwd("~/Documents/Drosophila Selection Experiment/allele_frequency_difference_testing")

query_header <- 'chr\tpos\tref	MB_counts	C1_counts	C1_dec	C1_inc	C1_startfreq	C1_endfreq	chr:1	pos:1	ref:1	MB_counts:1	C2_counts	C2_dec	C2_inc	C2_startfreq	C2_endfreq	chr:2	pos:2	ref:2	MB_counts:2	C3_counts	C3_dec	C3_inc	C3_startfreq	C3_endfreq	chr:3	pos:3	ref:3	MB_counts:3	C4_counts	C4_dec	C4_inc	C4_startfreq	C4_endfreq	chr:4	pos:4	ref:4	MB_counts:4	C5_counts	C5_dec	C5_inc	C5_startfreq	C5_endfreq	chr:5	pos:5	ref:5	MB_counts:5	D1_counts	D1_dec	D1_inc	D1_startfreq	D1_endfreq	chr:6	pos:6	ref:6	MB_counts:6	D2_counts	D2_dec	D2_inc	D2_startfreq	D2_endfreq	chr:7	pos:7	ref:7	MB_counts:7	D3_counts	D3_dec	D3_inc	D3_startfreq	D3_endfreq	chr:8	pos:8	ref:8	MB_counts:8	D4_counts	D4_dec	D4_inc	D4_startfreq	D4_endfreq	chr:9	pos:9	ref:9	MB_counts:9	D5_counts	D5_dec	D5_inc	D5_startfreq	D5_endfreq'
query_head <- unlist(strsplit(query_header, split='\t', fixed=TRUE))

#############
# FUNCTIONS #
#############

useful_subset <- function(database_query_output){
  individual_col_names <- paste(rep(Sample_code, each=4), 
                                rep(c('_dec', '_inc', '_startfreq', '_endfreq'), times=10), 
                                sep="")
  useful_col_names <- c('chr', 'pos', individual_col_names)
  useful_output <- database_query_output[,colnames(database_query_output)%in%useful_col_names]
}

find_and_rearrange_majmin_calls <- function(frequency_table){
  temp_useful <- useful_subset(frequency_table)
  #temp_useful <- temp_useful[temp_useful[,1]!="V1",]
  to_remove <- c()
  all_locus_output <- list()
  for(i in 1:nrow(temp_useful)){
    temp_row <- temp_useful[i,]
    #temp_sample <- j
    if(i%%1000==0){
      print(paste(i, "loci processed: up to", paste(temp_row[1:2], collapse=" ")))
    }
    MB_result <- unlist(c('MB', temp_row[3:5]))
    
    other_result_loc <- c(3, 4, 6)
    other_result_locations <- c(rep(other_result_loc, times=10) + rep(0:9*4, each=3))
    other_results <- temp_row[other_result_locations]
    
    output_table <- matrix(c(MB_result[2:4], other_results), ncol=11, nrow=3)
    colnames(output_table) <- c('MB', Sample_code)
    
    #identify the sample with the largest frequency change (to set the order)
    freqchange_vector <- as.numeric(output_table[3,2:11])-as.numeric(output_table[3,1])
    temp_sample <- colnames(output_table)[which.max(abs(freqchange_vector))+1]
    
    #now identify the columns with alleles in the wrong order and swap
    correct_order <- output_table[1:2,colnames(output_table)==temp_sample]
    temp_output <- output_table
    cols_to_switch <- which(temp_output[2,]==correct_order[[1]])
    output_table[1,cols_to_switch] <- temp_output[2,cols_to_switch]
    output_table[2,cols_to_switch] <- temp_output[1,cols_to_switch]
    output_table[3,cols_to_switch] <- 1-as.numeric(temp_output[3, cols_to_switch])
    
    #test whether there are still issues with maj/min calling
    if(length(which(output_table[1,]==correct_order[[1]][1]))<10){
      to_remove <- c(to_remove, i)
    }
    
    freqchange_vector <- as.numeric(output_table[3,2:11])-as.numeric(output_table[3,1])
    names(freqchange_vector) <- Sample_code
    all_locus_output[[i]] <- list(temp_row[c(1:6)], output_table, freqchange_vector)
  }
  all_locus_output <- all_locus_output[-to_remove]
}

### Calculating Ne ###

# For each SNP set, for each replicate line: 

# k = number of alleles (2)
# p0 = initial allele frequency
# p1 = final allele frequency
# t = number of generations apart (13)
# n0 = number of diploid individuals genotyped in initial generation (200)
# n1 = number of diploid individuals genotyped in final generation (50)

k <- 2
t <- 13
n0 <- 200
n1 <- 50

# Equation 4.18b from Lynch and Walsh Chapter 4 (original reference Nei and Tajima 1981b)

Fhat2 <- function(input_row, p0_location, p1_location){
  p0 <- input_row[p0_location]
  p1 <- input_row[p1_location]
  if((p0 < 0.01 & p1 < 0.04) | (p0 > 0.99 & p1 > 0.96) | is.na(p1)){
    ple <- NA
    #print(input_row[c(p0_location, p1_location)])
  } else {
    ple <- ((p0-p1)^2)/(((p0+p1)/2) - p0*p1)
    #print(input_row[c(p0_location, p1_location)])
#      if(ple < 0.00001){
#        print(input_row[c(p0_location, p1_location)])
#        #print(ple)
#      }
  }
  return(ple)
}

rep_testing <- subset(all_D_se, all_D_se[,1]==all_D_se[,2])

Ne_cross_gen_calc <- function(freq_data_frame, p0_location, p1_location){
  number_loci <- nrow(freq_data_frame)

  expression_output <- unlist(apply(freq_data_frame, MARGIN=1, FUN=Fhat2, 
                                    p0_location=p0_location, p1_location=p1_location))
  
  mean_F_over_all_loci <- mean(expression_output, na.rm=TRUE)
  no_informative_loci <- length(which(is.na(expression_output)==FALSE))
  gen13_allele_freq <- summary(freq_data_frame[which(is.na(expression_output)==FALSE),p1_location])
  print(paste('Average F = ', mean_F_over_all_loci, 'over', no_informative_loci, 'informative loci'))
  print(gen13_allele_freq)
  # Now plug this into equation 4.18a
  Ne <- (t - 2) / ((2*mean_F_over_all_loci) - (1/n0) - (1/n1))

}



########
# MAIN #
########

# import all C replicate SNP frequency tables and merge


for(i in Sample_code[1:5]){
  input_file_name <- paste(i, "_freq_table.txt", sep="")
  object_name <- paste(i, "_freq_table", sep="")
  assign(object_name, read.table(input_file_name, sep="\t", stringsAsFactors=FALSE, header=FALSE))
  if(i == "C1"){
    all_C_data_frame <- get(object_name)
  }
  else{
    all_C_data_frame <- rbind(allC_data_frame, get(object_name))
  }
}
colnames(all_C_data_frame) <- query_head
Cchrpos <- paste(all_C_data_frame$chr, all_C_data_frame$pos)
#remove duplicate rows
all_C_data_frame <- all_C_data_frame[which(duplicated(Cchrpos)==FALSE),]


# import all D replicate SNP frequency tables and merge

for(i in Sample_code[6:10]){
  input_file_name <- paste(i, "_noC_freq_table.txt", sep="")
  object_name <- paste(i, "_freq_table", sep="")
  assign(object_name, read.table(input_file_name, sep="\t", stringsAsFactors=FALSE, header=FALSE))
  if(i == "D1"){
    all_D_data_frame <- get(object_name)
  }
  else{
    all_D_data_frame <- rbind(all_D_data_frame, get(object_name))
  }
}
colnames(all_D_data_frame) <- query_head
Dchrpos <- paste(all_D_data_frame$chr, all_D_data_frame$pos)
# remove duplicate rows
all_D_data_frame <- all_D_data_frame[which(duplicated(Dchrpos)==FALSE),]

# import random SNP table, remove any that overlap with C or D candidates

#random_snps <- read.table("Random_10000_loci_freq.txt", sep="\t", stringsAsFactors=FALSE, header=FALSE)
random_snps <- read.table("Random_100000_loci_freq.txt", sep="\t", stringsAsFactors=FALSE, header=FALSE)
colnames(random_snps) <- query_head

just_endfreq_cols <- random_snps[,paste(Sample_code, "_endfreq", sep="")]
just_endfreq_cols[just_endfreq_cols=="absent"] <- NA
total_na <- rowSums(is.na(just_endfreq_cols))
start_na <- is.na(random_snps$C1_startfreq)
conditions <- cbind(total_na<5, start_na)
#rows_to_keep <- which(total_na<5, start_na==TRUE)
random_snps <- random_snps[which(conditions[,1]&conditions[,2]==FALSE),]

#Find and exclude SNPs that also occur in selected SNP sets
random_chrpos <- paste(random_snps$chr, random_snps$pos)
random_snps <- random_snps[which(random_chrpos%in%Dchrpos==FALSE & random_chrpos%in%Cchrpos==FALSE),]

# for each of the three tables, cut down to useful columns
# and rearrange frequencies so maj/min allele calls are consistent

all_C_rearranged <- find_and_rearrange_majmin_calls(all_C_data_frame)
all_D_rearranged <- find_and_rearrange_majmin_calls(all_D_data_frame)
random_rearranged <- find_and_rearrange_majmin_calls(random_snps)


# Pull out columns containing MB and gen13frequencies for each replicate line

all_C_start_end <- lapply(all_C_rearranged, "[[", 2)
all_C_se <- lapply(all_C_start_end, "[", 3,)
all_C_se <- lapply(all_C_se, unlist)
all_C_se <- lapply(all_C_se, as.numeric)
all_C_se <- as.data.frame(all_C_se)
colnames(all_C_se) <- paste("locus", 1:ncol(all_C_se), sep="")
all_C_se <- t(all_C_se)
colnames(all_C_se) <- c("MB", Sample_code)


all_D_start_end <- lapply(all_D_rearranged, "[[", 2)
all_D_se <- lapply(all_D_start_end, "[", 3,)
all_D_se <- lapply(all_D_se, unlist)
all_D_se <- lapply(all_D_se, as.numeric)
all_D_se <- as.data.frame(all_D_se)
colnames(all_D_se) <- paste("locus", 1:ncol(all_D_se), sep="")
all_D_se <- t(all_D_se)
colnames(all_D_se) <- c("MB", Sample_code)


random_start_end <- lapply(random_rearranged, "[[", 2)
random_se <- lapply(random_start_end, "[", 3,)
random_se <- lapply(random_se, unlist)
random_se <- lapply(random_se, as.numeric)
random_se <- as.data.frame(random_se)
colnames(random_se) <- paste("locus", 1:ncol(random_se), sep="")
random_se <- t(random_se)
colnames(random_se) <- c("MB", Sample_code)


Ne_results <- data.frame(SNP_set=character(), Rep=character(), C_or_D=character(), Ne=numeric(), stringsAsFactors=FALSE)
Ne_results[1,] <- list("lab", "C1", "C", Ne_cross_gen_calc(all_C_se, p0_location=1, p1_location=2))
Ne_results[2,] <- list("lab", "C2", "C", Ne_cross_gen_calc(all_C_se, p0_location=1, p1_location=3))
Ne_results[3,] <- list("lab", "C3", "C", Ne_cross_gen_calc(all_C_se, p0_location=1, p1_location=4))
Ne_results[4,] <- list("lab", "C4", "C", Ne_cross_gen_calc(all_C_se, p0_location=1, p1_location=5))
Ne_results[5,] <- list("lab", "C5", "C", Ne_cross_gen_calc(all_C_se, p0_location=1, p1_location=6))
Ne_results[6,] <- list("lab", "D1", "D", Ne_cross_gen_calc(all_C_se, p0_location=1, p1_location=7))
Ne_results[7,] <- list("lab", "D2", "D", Ne_cross_gen_calc(all_C_se, p0_location=1, p1_location=8))
Ne_results[8,] <- list("lab", "D3", "D", Ne_cross_gen_calc(all_C_se, p0_location=1, p1_location=9))
Ne_results[9,] <- list("lab", "D4", "D", Ne_cross_gen_calc(all_C_se, p0_location=1, p1_location=10))
Ne_results[10,] <- list("lab", "D5", "D", Ne_cross_gen_calc(all_C_se, p0_location=1, p1_location=11))
Ne_results[11,] <- list("des", "C1", "C", Ne_cross_gen_calc(all_D_se, p0_location=1, p1_location=2))
Ne_results[12,] <- list("des", "C2", "C", Ne_cross_gen_calc(all_D_se, p0_location=1, p1_location=3))
Ne_results[13,] <- list("des", "C3", "C", Ne_cross_gen_calc(all_D_se, p0_location=1, p1_location=4))
Ne_results[14,] <- list("des", "C4", "C", Ne_cross_gen_calc(all_D_se, p0_location=1, p1_location=5))
Ne_results[15,] <- list("des", "C5", "C", Ne_cross_gen_calc(all_D_se, p0_location=1, p1_location=6))
Ne_results[16,] <- list("des", "D1", "D", Ne_cross_gen_calc(all_D_se, p0_location=1, p1_location=7))
Ne_results[17,] <- list("des", "D2", "D", Ne_cross_gen_calc(all_D_se, p0_location=1, p1_location=8))
Ne_results[18,] <- list("des", "D3", "D", Ne_cross_gen_calc(all_D_se, p0_location=1, p1_location=9))
Ne_results[19,] <- list("des", "D4", "D", Ne_cross_gen_calc(all_D_se, p0_location=1, p1_location=10))
Ne_results[20,] <- list("des", "D5", "D", Ne_cross_gen_calc(all_D_se, p0_location=1, p1_location=11))
Ne_results[21,] <- list("neutral", "C1", "C", Ne_cross_gen_calc(random_se, p0_location=1, p1_location=2))
Ne_results[22,] <- list("neutral", "C2", "C", Ne_cross_gen_calc(random_se, p0_location=1, p1_location=3))
Ne_results[23,] <- list("neutral", "C3", "C", Ne_cross_gen_calc(random_se, p0_location=1, p1_location=4))
Ne_results[24,] <- list("neutral", "C4", "C", Ne_cross_gen_calc(random_se, p0_location=1, p1_location=5))
Ne_results[25,] <- list("neutral", "C5", "C", Ne_cross_gen_calc(random_se, p0_location=1, p1_location=6))
Ne_results[26,] <- list("neutral", "D1", "D", Ne_cross_gen_calc(random_se, p0_location=1, p1_location=7))
Ne_results[27,] <- list("neutral", "D2", "D", Ne_cross_gen_calc(random_se, p0_location=1, p1_location=8))
Ne_results[28,] <- list("neutral", "D3", "D", Ne_cross_gen_calc(random_se, p0_location=1, p1_location=9))
Ne_results[29,] <- list("neutral", "D4", "D", Ne_cross_gen_calc(random_se, p0_location=1, p1_location=10))
Ne_results[30,] <- list("neutral", "D5", "D", Ne_cross_gen_calc(random_se, p0_location=1, p1_location=11))


ne_plot <- ggplot(data=Ne_results, aes(x=SNP_set, y=Ne)) + 
  geom_point(aes(colour=factor(C_or_D), size=2)) + 
  scale_colour_manual(values=c(rgb(255, 0, 0, 100, maxColorValue=255), rgb(0, 0, 255, 100, maxColorValue=255))) +
  scale_y_log10(breaks=c(10, 100, 1000, 10000)) + 
  theme(legend.position="none")

pdf("Ne_plot_for_neutral_lab_and_desiccation_candidate_SNP_sets.pdf", height=7, width=5)
ne_plot
dev.off()

write.table(Ne_results, file="Ne_calculations.txt", sep="\t", quote=FALSE, row.names=FALSE)

#Look at the Ne values across the gen0—C rep vs gen0—D rep comparisons to test hypothesis
