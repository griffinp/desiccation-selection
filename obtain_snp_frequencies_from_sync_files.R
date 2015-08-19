setwd("~/Documents/Drosophila Selection Experiment/pileup_and_sync_files")

Sample_code <- c("C1", "C2", "C3", "C4", "C5",
                 "D1", "D2", "D3", "D4", "D5")

# Read in table of 'interesting' loci 
# (First dealing with those high-effect SNPs that
# are sig. differentiated in one or more D replicates)

D_high <- read.table('/Users/pgriffin/Documents/Drosophila\ Selection\ Experiment/snpEff_SNP_feature_enrichment/D_replicate_high_SNP_effects.txt',
                     stringsAsFactors=FALSE, sep="\t", header=TRUE)

position_as_search_string <- function(chr, pos){
  output <- paste(chr, pos, sep="\t")
  return(output)
}

all_locus_output <- list()
for(i in 1:nrow(D_high)){
  temp_row <- D_high[i,]
  print(temp_row$gene.name)
  output_table <- matrix(NA, nrow=5, ncol=11)
  colnames(output_table) <- c('MB', Sample_code)
  search_string <- position_as_search_string(temp_row$chr, temp_row$pos)
  grep_command <- paste("-P '", search_string, "' MB_*_noN_reduced.mpileup.sync", sep="")
  grep_output <- system2('grep', args=grep_command, stdout=TRUE)
  split_output <- unlist(strsplit(unlist(strsplit(unlist(strsplit(grep_output, split="\t")), split=":")), split="_"))
  MB_result <- matrix(split_output[c(1,8:11)], ncol=1, nrow=5)
  colnames(MB_result) <- MB_result[1,1]
  output_table[1:4,'MB'] <- MB_result[2:5,]
  
  other_result_loc <- c(2,14:17,20)
  samples_present <- c(1:(length(split_output)%/%20))
  for(j in samples_present){
    if(j < 2){
      other_result_locations <- other_result_loc
    }
    other_result_locations <- c(other_result_locations, other_result_loc+j*20)
  }
  other_results <- split_output[other_result_locations]
  other_results <- matrix(other_results, ncol=max(samples_present), nrow=6)
  colnames(other_results) <- other_results[1,]
  for(k in Sample_code){
    if(k %in% colnames(other_results)){
      output_table[,k] <- other_results[2:6,k]
    }
  }
  all_locus_output[[i]] <- list(temp_row[,c(1:6,11)], output_table)
}


#######################
#######################

# This has worked fine,
# but has revealed a major problem!
# Stupidly set a minimum minor allewle count of 5 reads
# in the FET step, which (I think) has meant a lot of
# loci have been excluded esp. from the lower-coverage reps.
# Will have to repeat EVERYTHING. 

#######################
#######################


