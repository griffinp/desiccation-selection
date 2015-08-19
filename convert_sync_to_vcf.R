# Convert sync to vcf format

#############
# FUNCTIONS #
#############

split_by_colon <- function(base_call_object){
  object_split <- strsplit(base_call_object, ":")
  object_as_vector <- as.numeric(unlist(object_split))
  names(object_as_vector) <- c("A", "T", "C", "G", "N", "gap")
  return(object_as_vector)
}

find_maj_and_min <- function(base_count_col){
  base_count_col <- base_count_col[1:4]
  #different_condition <- max(base_count_col) != max(base_count_col[-which.max(base_count_col)])
  #if(different_condition){
    #if the overall maj and min have different total read counts
  #actually, this condition is unnecessary as it's sorted out later
  max_location <- base_count_col[which.max(base_count_col)]
  max_name <- names(max_location)
  subtract_max <- base_count_col[-which.max(base_count_col)]
  min_location <- subtract_max[which.max(subtract_max)]
  min_name <- names(min_location)
  #}
  #else {
    #for the case where the overall maj and min have the same total read counts
  #}
  #return(c(max_location, min_location))
  return(c(max_name, min_name))
}

extract_counts_matching_letters <- function(base_call_col){
  # This is to be used on a column from a data frame
  # containing counts for A to gap, 
  # then maj_letter and then min_letter
  # e.g. 
#   
#   [,1] [,2] [,3] [,4]
#   A          "30" "0"  "0"  "0" 
#   T          "5"  "43" "0"  "0" 
#   C          "0"  "22" "30" "93"
#   G          "0"  "1"  "2"  "50"
#   N          "0"  "0"  "0"  "0" 
#   gap        "0"  "0"  "0"  "0" 
#   maj_letter "A"  "T"  "C"  "G" 
#   min_letter "T"  "C"  "G"  "C" 
  #desired_letters <- which(row.names(base_call_col%in%letter_col)
  maj_location <- which(names(base_call_col)==base_call_col['maj_letter'])
  min_location <- which(names(base_call_col)==base_call_col['min_letter'])
  desired_counts <- base_call_col[c(maj_location, min_location)]
  return(as.numeric(desired_counts))
}

find_maj_and_min_for_2_sample_file <- function(base_call_df_1, base_call_df_2){
  split1 <- apply(base_call_df_1, 1, split_by_colon)
  split2 <- apply(base_call_df_2, 1, split_by_colon)
  total_counts <- split1+split2
  overall_maj_min <- apply(total_counts, 2, find_maj_and_min)
  return(list(split1, split2, overall_maj_min))
}

extract_maj_and_min_for_2_sample_file <- function(base_call_df_1, base_call_df_2){
  relevant_info <- find_maj_and_min_for_2_sample_file(base_call_df_1, base_call_df_2)
  rownames(relevant_info[[3]]) <- c("maj_letter", "min_letter")
  base_counts_df_1 <- rbind(relevant_info[[1]], relevant_info[[3]])
  base_counts_df_2 <- rbind(relevant_info[[2]], relevant_info[[3]])  
  allele_calls <- as.data.frame(matrix(NA, ncol=6,
                                       nrow=nrow(base_call_df_1))
                             )
  df1_majmin_counts <- apply(base_counts_df_1, 2, extract_counts_matching_letters)
  df2_majmin_counts <- apply(base_counts_df_2, 2, extract_counts_matching_letters)
  allele_calls[,3:4] <- c(df1_majmin_counts[1,], df1_majmin_counts[2,])
  allele_calls[,5:6] <- c(df2_majmin_counts[1,], df2_majmin_counts[2,])
  allele_calls[,1:2] <- c(relevant_info[[3]][1,], relevant_info[[3]][2,])
  colnames(allele_calls) <- c("overall_maj", "overall_min", "df_1_maj_count",
                              "df_1_min_count", "df_2_maj_count", "df_2_min_count")
  rownames(allele_calls) <- 1:nrow(allele_calls)
  return(allele_calls)
}


rearrange_maj_min <- function(majmin_call_object){
  #This function rearranges the maj and min calls so that the 'maj' allele
  #is the one that decreased in frequency between sample 1 and sample 2
  majmin_call_object <- as.data.frame(majmin_call_object)
  total_df_1_count <- as.numeric(as.character(majmin_call_object$df_1_maj_count))+as.numeric(as.character(majmin_call_object$df_1_min_count))
  df_1_maj_proportion <- as.numeric(as.character(majmin_call_object$df_1_maj_count))/total_df_1_count
  total_df_2_count <- as.numeric(as.character(majmin_call_object$df_2_maj_count))+as.numeric(as.character(majmin_call_object$df_2_min_count))
  df_2_maj_proportion <- as.numeric(as.character(majmin_call_object$df_2_maj_count))/total_df_2_count
  
  maj_has_increased <- df_2_maj_proportion > df_1_maj_proportion
  new_maj <- as.character(majmin_call_object$overall_maj)
  new_min <- as.character(majmin_call_object$overall_min)
  new_maj[which(maj_has_increased)] <- as.character(majmin_call_object$overall_min[which(maj_has_increased)])
  new_min[which(maj_has_increased)] <- as.character(majmin_call_object$overall_maj[which(maj_has_increased)])
  return(data.frame(new_maj, new_min))  
}

#testing

# test_df_1 <- data.frame(base_calls=c('30:5:0:0:0:0', '0:43:22:1:0:0', 
#                                       '0:0:30:2:0:0', '0:0:93:50:0:0'))
# test_df_2 <- data.frame(base_calls=c('23:10:0:0:0:0', '0:40:32:0:0:0', 
#                                      '0:0:52:2:0:0', '0:0:4:48:0:0'))
# relevant_info <- find_maj_and_min_for_2_sample_file(test_df_1, test_df_2)
# test_majmin_calls <- extract_maj_and_min_for_2_sample_file(test_df_1, test_df_2)
# rearrange_maj_min(test_majmin_calls)


########
# MAIN #
########

Sample_code <- c("C1", "C2", "C3", "C4", "C5",
                 "D1", "D2", "D3", "D4", "D5")

chr_vector <- c("dmel_mitochondrion_genome", "2L", "2R", "3L", "3R", "X")

# for (i in Sample_code){
#   temp_filename <- paste("/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/",
#                          i, "_sig_SNPs_sync.txt", sep="")
#   output_filename <- paste("/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/",
#                            i, "_sig_SNPs_150722.vcf", sep="")
#   temp_table <- read.table(temp_filename, header=FALSE, sep="\t",
#                            stringsAsFactors=FALSE)
#   temp_base_call_col_1 <- data.frame(base_calls=temp_table[,4])
#   temp_base_call_col_2 <- data.frame(base_calls=temp_table[,5])
#   majmin_calls <- extract_maj_and_min_for_2_sample_file(temp_base_call_col_1,
#                                                         temp_base_call_col_2)
#   found_alleles <- rearrange_maj_min(majmin_calls)
#   vcf_df <- data.frame(temp_table[,c(1,2)], 
#                        rep('.', times=nrow(temp_table)), 
#                        found_alleles, 
#                        rep(30, times=nrow(temp_table)), 
#                        rep('PASS', times=nrow(temp_table)), 
#                        rep('', times=nrow(temp_table)))
#   colnames(vcf_df) <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
#   write.table(vcf_df, file=output_filename, sep="\t", row.names=FALSE, quote=FALSE)
#   print(paste(nrow(vcf_df), " lines written for sample ", i, sep=""))
# }
# 
# for (i in Sample_code){

#### NB THESE SYNC FILES DO NOT CONTAIN ALL SNPS, ONLY THOSE WITH P < 1E-05

#   temp_filename <- paste("/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/",
#                          i, "_all_SNPs_sync.txt", sep="")
#   output_filename <- paste("/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/",
#                            i, "_all_SNPs_150722.vcf", sep="")
#   temp_table <- read.table(temp_filename, header=FALSE, sep="\t",
#                            stringsAsFactors=FALSE)
#   temp_base_call_col_1 <- data.frame(base_calls=temp_table[,4])
#   temp_base_call_col_2 <- data.frame(base_calls=temp_table[,5])
#   majmin_calls <- extract_maj_and_min_for_2_sample_file(temp_base_call_col_1,
#                                                         temp_base_call_col_2)
#   found_alleles <- rearrange_maj_min(majmin_calls)
#   vcf_df <- data.frame(temp_table[,c(1,2)], 
#                        rep('.', times=nrow(temp_table)), 
#                        found_alleles, 
#                        rep(30, times=nrow(temp_table)), 
#                        rep('PASS', times=nrow(temp_table)), 
#                        rep('', times=nrow(temp_table)))
#   colnames(vcf_df) <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
#   write.table(vcf_df, file=output_filename, sep="\t", row.names=FALSE, quote=FALSE)
#   print(paste(nrow(vcf_df), " lines written for sample ", i, sep=""))
# }

# Version to run on VLSCI for all_SNPs files:

#dir_path <- "/vlsci/VR0267/pgriffin/hsm/sandra/output_wgs/variant_calls/"

dir_path <- "/Users/pgriffin/Documents/Drosophila\ Selection\ Experiment/snp_and_gene_lists/"

for (i in Sample_code){
  #temp_filename <- paste(dir_path, "MB_",
  #                       i, "_noN_reduced.mpileup.sync", sep="")
  temp_filename <- paste(dir_path, i, "_sig_SNPs_sync.txt", sep="")
  print(paste("Reading sync file for", i, sep=" "))
  temp_table_pre <- read.table(temp_filename, header=FALSE, sep="\t",
                               stringsAsFactors=FALSE)
  for(j in chr_vector){
    output_filename <- paste(dir_path,
                             i, "_", j, "_sig_SNPs_150818.vcf", sep="")
    print(paste("Subsetting", i, "for chr", j, sep=" "))
    temp_table <- subset(temp_table_pre, temp_table_pre[,1]==j)
    print(paste("Processing", i, "for chr", j, sep=" "))
    temp_base_call_col_1 <- data.frame(base_calls=temp_table[,4])
    temp_base_call_col_2 <- data.frame(base_calls=temp_table[,5])
    majmin_calls <- extract_maj_and_min_for_2_sample_file(temp_base_call_col_1,
                                                          temp_base_call_col_2)
    found_alleles <- rearrange_maj_min(majmin_calls)
    vcf_df <- data.frame(temp_table[,c(1,2)], 
                         rep('.', times=nrow(temp_table)), 
                         found_alleles, 
                         rep(30, times=nrow(temp_table)), 
                         rep('PASS', times=nrow(temp_table)), 
                         rep('', times=nrow(temp_table)))
    colnames(vcf_df) <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
    write.table(vcf_df, file=output_filename, sep="\t", row.names=FALSE, quote=FALSE)
    print(paste(nrow(vcf_df), " total lines written for sample ", i, " chr ", j, sep=""))
  }
  remove(temp_table_pre)
}


