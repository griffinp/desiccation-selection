setwd("/vlsci/VR0267/pgriffin/hsm/sandra/output_wgs/variant_calls/")
test_file <- read.csv('MB_C1_noN.mpileup.sync', sep = "\t", header=FALSE)
output_file_name <- "MB_C1_all_freq_change_150821.txt"
output_file_prefix <- "MB_C1_all_freq_change_150821_"
# plot1_prefix <- "Minuslog10_p_vs_starting_freq_MB_C1_"
# plot2_prefix <- "End_freq_vs_starting_freq_MB_C1_"

library(polysat)
library("bbmle")


# setwd("/Users/pgriffin/Documents/Drosophila Selection Experiment/pileup_and_sync_files")
# test_file <- read.csv('MB_C1_noN.mpileup.sync', sep = "\t", header=FALSE)
# output_file_name <- "MB_C1_all_freq_change_150821.txt"
# output_file_prefix <- "MB_C1_all_freq_change_150821_"


#############################################
#                                           #
# Functions for reading in data.            #
# Dealing with real data from a two-pop     #
# 'sync' file (created from a pileup file)  #
#                                           #
#############################################

# Considerations
# - have to make sure the maj and min alleles are the same in both pops
# (can use my GWAS pileup file processing code for this)
# - have to make zero- and one-bounded (as Alex did above)

####
####
####

split_by_colon <- function(base_call_object){
  object_split <- strsplit(base_call_object, ":")
  object_as_vector <- as.numeric(unlist(object_split))
  names(object_as_vector) <- c("A", "T", "C", "G", "N", "gap")
  return(object_as_vector)
}

find_maj_and_min <- function(base_count_col){
  base_count_col <- base_count_col[1:4]
  max_location <- base_count_col[which.max(base_count_col)]
  max_name <- names(max_location)
  subtract_max <- base_count_col[-which.max(base_count_col)]
  min_location <- subtract_max[which.max(subtract_max)]
  min_name <- names(min_location)
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


rearrange_maj_min_and_counts <- function(majmin_call_object){
  #This function rearranges the maj and min calls so that the 'maj' allele
  #is the one that decreased in frequency between sample 1 and sample 2
  #It also outputs the counts correctly
  majmin_call_object <- as.data.frame(majmin_call_object)
  total_df_1_count <- as.numeric(as.character(majmin_call_object$df_1_maj_count))+as.numeric(as.character(majmin_call_object$df_1_min_count))
  df_1_maj_proportion <- as.numeric(as.character(majmin_call_object$df_1_maj_count))/total_df_1_count
  total_df_2_count <- as.numeric(as.character(majmin_call_object$df_2_maj_count))+as.numeric(as.character(majmin_call_object$df_2_min_count))
  df_2_maj_proportion <- as.numeric(as.character(majmin_call_object$df_2_maj_count))/total_df_2_count
  
  maj_has_increased <- df_2_maj_proportion > df_1_maj_proportion
  new_df <- majmin_call_object
  #new_maj <- as.character(majmin_call_object$overall_maj)
  #new_min <- as.character(majmin_call_object$overall_min)
  new_df[which(maj_has_increased), 'overall_maj'] <- as.character(majmin_call_object$overall_min[which(maj_has_increased)])
  new_df[which(maj_has_increased), 'overall_min'] <- as.character(majmin_call_object$overall_maj[which(maj_has_increased)])
  new_df[which(maj_has_increased), 'df_1_maj_count'] <- as.character(majmin_call_object$df_1_min_count[which(maj_has_increased)])
  new_df[which(maj_has_increased), 'df_1_min_count'] <- as.character(majmin_call_object$df_1_maj_count[which(maj_has_increased)])
  new_df[which(maj_has_increased), 'df_2_maj_count'] <- as.character(majmin_call_object$df_2_min_count[which(maj_has_increased)])
  new_df[which(maj_has_increased), 'df_2_min_count'] <- as.character(majmin_call_object$df_2_maj_count[which(maj_has_increased)])
  
  
  return(new_df)  
}




#############################################
#                                           #
#     Now read in the real data             #
#                                           #
#############################################

# lowestp <- test_file
# rm(test_file)

chrs <- c("2L", "2R", "3L", "3R", "X", "dmel_mitochondrion_genome")
print("starting loop")

for (tempchr in chrs){
  print(paste("About to subset for", tempchr))
  chr_subset <- subset(test_file, test_file[,1]==tempchr)
  print(paste(nrow(chr_subset), "rows in subset"))
  output_file_name <- paste(output_file_prefix, tempchr, '.txt', sep="")
  print(paste("Processing to produce", output_file_name))
  
  temp_base_call_col_1 <- data.frame(base_calls=chr_subset[,4])
  temp_base_call_col_2 <- data.frame(base_calls=chr_subset[,5])
  majmin_calls <- extract_maj_and_min_for_2_sample_file(temp_base_call_col_1,
                                                        temp_base_call_col_2)
  majmin_counts <- rearrange_maj_min_and_counts(majmin_calls)
  colnames(majmin_counts) <- c("main", "alt", "gen0_maj", "gen0_min", "gen13_maj", "gen13_min")
 # majmin_counts <- sapply(majmin_counts, as.numeric)
  
  startfreq <- as.numeric(majmin_counts$gen0_min)/(as.numeric(majmin_counts$gen0_maj)+as.numeric(majmin_counts$gen0_min))
  endfreq <- as.numeric(majmin_counts$gen13_min)/(as.numeric(majmin_counts$gen13_maj)+as.numeric(majmin_counts$gen13_min))
   
  out_tab <- data.frame(chr_subset, majmin_counts[1:2], startfreq, endfreq)
 
  write.table(out_tab, file=output_file_name, quote=FALSE, row.names=FALSE, sep="\t")
  rm(chr_subset)
}

################################################################

