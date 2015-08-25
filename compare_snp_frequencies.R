
library(RSQLite)

#setwd("~/Documents/Drosophila Selection Experiment/snp_and_gene_lists")
setwd("/vlsci/VR0267/pgriffin/hsm/sandra/output_wgs/sql_database")

Sample_code <- c("C1", "C2", "C3", "C4", "C5",
                 "D1", "D2", "D3", "D4", "D5")


position_as_search_string <- function(chr, pos){
  output <- paste(chr, pos, sep="\t")
  return(output)
}

db <- dbConnect(SQLite(), dbname="Freq.sqlite")

# Had to make the SQL database on barcoo as it took up a great deal of memory
# (probably > 100 Gb all up?)

dbWriteTable(conn = db, name = "C1", value = "../variant_calls/C1_all_freq_change_150821.txt",
            sep="\t", row.names = FALSE, header = TRUE)
dbListFields(db, "C1")
dbWriteTable(conn = db, name = "C2", value = "../variant_calls/C2_all_freq_change_150821.txt",
            sep="\t", row.names = FALSE, header = TRUE)
dbWriteTable(conn = db, name = "C3", value = "../variant_calls/C3_all_freq_change_150821.txt",
            sep="\t", row.names = FALSE, header = TRUE)
dbWriteTable(conn = db, name = "C4", value = "../variant_calls/C4_all_freq_change_150821.txt",
            sep="\t", row.names = FALSE, header = TRUE)
dbWriteTable(conn = db, name = "C5", value = "../variant_calls/C5_all_freq_change_150821.txt",
            sep="\t", row.names = FALSE, header = TRUE)
print(paste("imported C tables", Sys.time()))
dbWriteTable(conn = db, name = "D1", value = "../variant_calls/D1_all_freq_change_150821.txt",
            sep="\t", row.names = FALSE, header = TRUE)
dbWriteTable(conn = db, name = "D2", value = "../variant_calls/D2_all_freq_change_150821.txt",
            sep="\t", row.names = FALSE, header = TRUE)
dbWriteTable(conn = db, name = "D3", value = "../variant_calls/D3_all_freq_change_150821.txt",
            sep="\t", row.names = FALSE, header = TRUE)
dbWriteTable(conn = db, name = "D4", value = "../variant_calls/D4_all_freq_change_150821.txt",
            sep="\t", row.names = FALSE, header = TRUE)
dbWriteTable(conn = db, name = "D5", value = "../variant_calls/D5_all_freq_change_150821.txt",
            sep="\t", row.names = FALSE, header = TRUE)
print(paste("imported D tables", Sys.time()))

dbListFields(db, "C1")
#dbReadTable(db, "C1")

dbSendQuery(db, 'CREATE INDEX C1_index ON C1 (`chr`, `pos`)')
dbSendQuery(db, 'CREATE INDEX C2_index ON C2 (`chr`, `pos`)')
dbSendQuery(db, 'CREATE TABLE C1_C2 AS SELECT * FROM C1 LEFT OUTER JOIN C2 ON (C1.chr = C2.chr AND C1.pos=C2.pos)')
dbSendQuery(db, 'DROP TABLE C1')
dbSendQuery(db, 'DROP TABLE C2')
dbSendQuery(db, 'CREATE INDEX C1_C2_index ON C1_C2 (`chr`, `pos`)')
print(paste("merged and indexed C1_C2", Sys.time()))

dbSendQuery(db, 'CREATE INDEX C3_index ON C3 (`chr`, `pos`)')
dbSendQuery(db, 'CREATE INDEX C4_index ON C4 (`chr`, `pos`)')
dbSendQuery(db, 'CREATE TABLE C3_C4 AS SELECT * FROM C3 LEFT OUTER JOIN C4 ON (C3.chr = C4.chr AND C3.pos=C4.pos)')
dbSendQuery(db, 'DROP TABLE C3')
dbSendQuery(db, 'DROP TABLE C4')
dbSendQuery(db, 'CREATE INDEX C3_C4_index ON C3_C4 (`chr`, `pos`)')
print(paste("merged and indexed C3_C4", Sys.time()))


dbSendQuery(db, 'CREATE INDEX C5_index ON C5 (`chr`, `pos`)')
dbSendQuery(db, 'CREATE INDEX D1_index ON D1 (`chr`, `pos`)')
dbSendQuery(db, 'CREATE TABLE C5_D1 AS SELECT * FROM C5 LEFT OUTER JOIN D1 ON (C5.chr = D1.chr AND C5.pos=D1.pos)')
dbSendQuery(db, 'DROP TABLE C5')
dbSendQuery(db, 'DROP TABLE D1')
dbSendQuery(db, 'CREATE INDEX C5_D1_index ON C5_D1 (`chr`, `pos`)')
print(paste("merged and indexed C5_D1", Sys.time()))


dbSendQuery(db, 'CREATE INDEX D2_index ON D2 (`chr`, `pos`)')
dbSendQuery(db, 'CREATE INDEX D3_index ON D3 (`chr`, `pos`)')
dbSendQuery(db, 'CREATE TABLE D2_D3 AS SELECT * FROM D2 LEFT OUTER JOIN D3 ON (D2.chr = D3.chr AND D2.pos=D3.pos)')
dbSendQuery(db, 'DROP TABLE D2')
dbSendQuery(db, 'DROP TABLE D3')
dbSendQuery(db, 'CREATE INDEX D2_D3_index ON D2_D3 (`chr`, `pos`)')
print(paste("merged and indexed D2_D3", Sys.time()))


dbSendQuery(db, 'CREATE INDEX D4_index ON D4 (`chr`, `pos`)')
dbSendQuery(db, 'CREATE INDEX D5_index ON D5 (`chr`, `pos`)')
dbSendQuery(db, 'CREATE TABLE D4_D5 AS SELECT * FROM D4 LEFT OUTER JOIN D5 ON (D4.chr = D5.chr AND D4.pos=D5.pos)')
dbSendQuery(db, 'DROP TABLE D4')
dbSendQuery(db, 'DROP TABLE D5')
dbSendQuery(db, 'CREATE INDEX D4_D5_index ON D4_D5 (`chr`, `pos`)')
print(paste("merged and indexed D4_D5", Sys.time()))


dbSendQuery(db, 'CREATE TABLE C1_C2_C3_C4 AS SELECT * FROM C1_C2 LEFT OUTER JOIN C3_C4 ON (C1_C2.chr = C3_C4.chr AND C1_C2.pos=C3_C4.pos)')
dbSendQuery(db, 'DROP TABLE C1_C2')
dbSendQuery(db, 'DROP TABLE C3_C4')
dbSendQuery(db, 'CREATE INDEX C1_C2_C3_C4_index ON C1_C2_C3_C4 (`chr`, `pos`)')
print(paste("merged and indexed C1 to C4", Sys.time()))

dbSendQuery(db, 'CREATE TABLE C5_D1_D2_D3 AS SELECT * FROM C5_D1 LEFT OUTER JOIN D2_D3 ON (C5_D1.chr = D2_D3.chr AND C5_D1.pos=D2_D3.pos)')
dbSendQuery(db, 'DROP TABLE C5_D1')
dbSendQuery(db, 'DROP TABLE D2_D3')
dbSendQuery(db, 'CREATE INDEX C5_D1_D2_D3_index ON C5_D1_D2_D3 (`chr`, `pos`)')
print(paste("merged and indexed C5 to D3", Sys.time()))

dbSendQuery(db, 'CREATE TABLE C1_to_D3 AS SELECT * FROM C1_C2_C3_C4 LEFT OUTER JOIN C5_D1_D2_D3 ON (C1_C2_C3_C4.chr = C5_D1_D2_D3.chr AND C1_C2_C3_C4.pos=C5_D1_D2_D3.pos)')
dbSendQuery(db, 'DROP TABLE C1_C2_C3_C4')
dbSendQuery(db, 'DROP TABLE C5_D1_D2_D3')
dbSendQuery(db, 'CREATE INDEX C1_to_D3_index ON C1_to_D3 (`chr`, `pos`)')
print(paste("merged and indexed C1 to D3", Sys.time()))


dbSendQuery(db, 'CREATE TABLE all_merged AS SELECT * FROM C1_to_D3 LEFT OUTER JOIN D4_D5 ON (C1_to_D3.chr = D4_D5.chr AND C1_to_D3.pos=D4_D5.pos)')


dbSendQuery(db, 'DROP TABLE C1_to_D3')
dbSendQuery(db, 'DROP TABLE D4_D5')

dbSendQuery(db, 'CREATE INDEX all_merged_index ON all_merged (`chr`, `pos`)')
print(paste("merged and indexed all 10 tables", Sys.time()))


testquery <- dbGetQuery(db, "SELECT * FROM all_merged WHERE chr='2L' AND pos='5010'")
print("Test Query")
print(testquery)

# test_D_high <- D_high[1:5,]
# test_D_high[,1] <- 'C1'
# test_D_high[,2] <- '2L'
# test_D_high[,3] <- c(4864, 4867, 4868, 4875, 4879)
# test_D_high[,5] <- c('T', 'G', 'T', 'G', 'G')
# test_D_high[,6] <- c('A', 'C', 'A', 'T', 'C')

extract_certain_loci <- function(loci_table, chr_col, pos_col){
  output_matrix <- matrix('m', nrow=nrow(loci_table), ncol=90)
  #output_col_names
  for(i in 1:nrow(loci_table)){
    temprow <- loci_table[i,]
    tempchr <- loci_table[i,chr_col]
    temppos <- loci_table[i,pos_col]
    tempquerystring <- paste("SELECT * FROM all_merged WHERE chr='", tempchr, "' AND pos='", temppos, "'", sep="")
    tempqueryresult <- dbGetQuery(db, tempquerystring)
    output_matrix[i,] <- unlist(tempqueryresult)
  }
  colnames(output_matrix) <- unlist(dbListFields(db, 'all_merged'))
  return(output_matrix)
}

extract_random_loci <- function(number_loci){
  output_matrix <- matrix('m', nrow=number_loci, ncol=90)
  querystring <- paste("SELECT * FROM all_merged ORDER BY RANDOM() LIMIT ", number_loci, ";", sep="")
  queryresult <- dbGetQuery(db, querystring)
  output_matrix <- unlist(tempqueryresult)
  colnames(output_matrix) <- unlist(dbListFields(db, 'all_merged'))
  return(output_matrix)
}



D_high <- read.table('/vlsci/VR0267/pgriffin/hsm/sandra/output_wgs/sql_database/D_high_effects.txt',
                     stringsAsFactors=FALSE, sep="\t", header=TRUE)
d_high_table <- extract_certain_loci(D_high, chr_col=2, pos_col=3)
write.table(d_high_table, file="D_high_allele_freq.txt", sep="\t"
            row.names=FALSE, col.names=FALSE)

C_mod_high <- read.table('/vlsci/VR0267/pgriffin/hsm/sandra/output_wgs/sql_database/C_mod_high_effects.txt',
                     stringsAsFactors=FALSE, sep="\t", header=TRUE)
c_mod_high_table <- extract_certain_loci(C_mod_high, chr_col=2, pos_col=3)
write.table(c_mod_high_table, file="C_mod_high_allele_freq.txt", sep="\t"
            row.names=FALSE, col.names=FALSE)

random_100000 <- extract_random_loci(100000)
write.table(random_100000, file="Random_100000_loci_freq.txt", sep="\t",
            row.names=FALSE, col.names=FALSE)

for(i in Sample_code){
  print(paste("Processing sig. locus list for", i))
  temp_filename <- paste(i, "_fet_drift_sig_SNPs_sync.txt", sep="")
  temp_outputname <- paste(i, "_freq_table.txt", sep="")
  temp_table <- read.table(temp_filename, header=FALSE, sep="\t")
  temp_output <- extract_certain_loci(get(temp_table), chr_col=1, pos_col=2)
  write.table(temp_output, file=temp_outputname, sep="\t",
              row.names=FALSE, col.names=FALSE)
}

# other frequency tables to get
# frequency table for each list of sig. SNPs (i.e. 10 in total)
# frequency table for 100,000 randomly-chosen loci
# 




### now try with full files...




#The code below works - for finding each locus in each file using grep -
# but is super slow because the files are so large. Investigating SQL instead.... 

# all_locus_output <- list()
# for(i in 1:nrow(D_high)){
#   temp_row <- D_high[i,]
#   temp_sample <- temp_row[,1]
#   print(paste(temp_row$chr, temp_row$pos))
#   output_table <- matrix(NA, nrow=3, ncol=11)
#   colnames(output_table) <- c('MB', Sample_code)
#   search_string <- position_as_search_string(temp_row$chr, temp_row$pos)
#   grep_command <- paste("-P -m 1 '", search_string, "' *_all_freq_change_150821.txt", sep="")
#   grep_output <- system2('grep', args=grep_command, stdout=TRUE)
#   print(grep_output)
#   split_output <- unlist(strsplit(unlist(strsplit(unlist(strsplit(grep_output, split="\t")), split=":")), split="_"))
#   MB_result <- matrix(c('MB', split_output[21:23]), ncol=4, nrow=1)
#   colnames(MB_result) <- c('sample', 'dec', 'inc', 'freq')
#   output_table[1:3,'MB'] <- MB_result[,2:4]
#   
#   other_result_loc <- c(1,21, 22, 24)
#   samples_present <- c(1:(length(split_output)%/%24)-1)
#   for(j in samples_present){
#     if(j < 2){
#       other_result_locations <- other_result_loc
#     }
#     other_result_locations <- c(other_result_locations, other_result_loc+j*24)
#   }
#   other_results <- split_output[other_result_locations]
#   other_results <- matrix(other_results, byrow=TRUE, nrow=max(samples_present), ncol=4)
#   colnames(other_results) <- colnames(MB_result)
#   
#   for(k in Sample_code){
#     if(k %in% other_results[,1]){
#       output_table[,k] <- other_results[other_results[,1]==k, 2:4]
#     }
#   }
#   
#   #now identify the columns with alleles in the wrong order and swap
#   correct_order <- output_table[1:2,colnames(output_table)==temp_sample]
#   temp_output <- output_table
#   cols_to_switch <- which(temp_output[1,]==correct_order[2])
#   output_table[1,cols_to_switch] <- temp_output[2,cols_to_switch]
#   output_table[2,cols_to_switch] <- temp_output[1,cols_to_switch]
#   output_table[3,cols_to_switch] <- 1-as.numeric(temp_output[3, cols_to_switch])
#   
#   all_locus_output[[i]] <- list(temp_row[,c(1:6)], output_table)
# }


#######################
#######################

