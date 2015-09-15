Sample_code <- c("C1", "C2", "C3", "C4", "C5",
                 "D1", "D2", "D3", "D4", "D5")

setwd("~/Documents/Drosophila Selection Experiment/allele_frequency_difference_testing")

# import all C replicate SNP frequency tables and merge

for(i in Sample_code[1:5]){
  input_file_name <- paste(i, "_freq_table.txt", sep="")
  object_name <- paste(i, "_freq_table", sep="")
  assign(object_name, read.table(input_file_name, sep="\t", stringsAsFactors=FALSE, header=FALSE))
}

# import all D replicate SNP frequency tables and merge

# import random SNP table, remove any that overlap with C or D candidates

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

per_locus_expression <- function(input_row, p0_location, p1_location){
  p0 <- input_row[p0_location]
  p1 <- input_row[p1_location]
  ple <- (p0-p1)^2/((p0+p1)/2 - p0*p1)
  return(ple)
}

number_loci <- nrow(*****)

expression_output <- unlist(apply(****, MARGIN=1, FUN=per_locus_expression, 
                           p0_location=1, p1_location=2))

# expression_output <- vector('NA', length=number_loci)
# for(i in 1:number_loci){
#   expression_output[i] <- per_locus_expression(***p0_column, p1_column***)
# }

sum_over_all_loci <- sum(expression_loci)
Fhat2 <- sum_over_all_loci/k

# Now plug this into equation 4.18a

Ne <- (t - 2) / ((2*Fhat2) - (1/n0) - (1/n1))

#Look at the Ne values across the gen0—C rep vs gen0—D rep comparisons to test hypothesis
