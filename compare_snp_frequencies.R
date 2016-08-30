
#library(RSQLite)
library(likert)

setwd("~/Documents/Drosophila Selection Experiment/snp_and_gene_lists")
#setwd("/vlsci/VR0267/pgriffin/hsm/sandra/output_wgs/sql_database")

Sample_code <- c("C1", "C2", "C3", "C4", "C5",
                 "D1", "D2", "D3", "D4", "D5")


position_as_search_string <- function(chr, pos){
  output <- paste(chr, pos, sep="\t")
  return(output)
}

useful_subset <- function(database_query_output){
  individual_col_names <- paste(rep(Sample_code, each=4), 
                                rep(c('_dec', '_inc', '_startfreq', '_endfreq'), times=10), 
                                sep="")
  useful_col_names <- c('chr', 'pos', individual_col_names)
  useful_output <- database_query_output[,colnames(database_query_output)%in%useful_col_names]
}

query_header <- 'chr pos ref MB_counts C1_counts C1_dec C1_inc C1_startfreq C1_endfreq chr:1 pos:1 ref:1 MB_counts:1 C2_counts C2_dec C2_inc C2_startfreq C2_endfreq chr:2 pos:2 ref:2 MB_counts:2 C3_counts C3_dec C3_inc C3_startfreq C3_endfreq chr:3 pos:3 ref:3 MB_counts:3 C4_counts C4_dec C4_inc C4_startfreq C4_endfreq chr:4 pos:4 ref:4 MB_counts:4 C5_counts C5_dec C5_inc C5_startfreq C5_endfreq chr:5 pos:5 ref:5 MB_counts:5 D1_counts D1_dec D1_inc D1_startfreq D1_endfreq chr:6 pos:6 ref:6 MB_counts:6 D2_counts D2_dec D2_inc D2_startfreq D2_endfreq chr:7 pos:7 ref:7 MB_counts:7 D3_counts D3_dec D3_inc D3_startfreq D3_endfreq chr:8 pos:8 ref:8 MB_counts:8 D4_counts D4_dec D4_inc D4_startfreq D4_endfreq chr:9 pos:9 ref:9 MB_counts:9 D5_counts D5_dec D5_inc D5_startfreq D5_endfreq'
query_head <- unlist(strsplit(query_header, split=' ', fixed=TRUE))


net_stacked <- function(x, j, letter_label_to_add, n_to_add) {
  #Function from here http://statisfactions.com/2012/improved-net-stacked-distribution-graphs-via-ggplot2-trickery/
  ## x: a data.frame or list, where each column is a ordered factor with the same levels
  ## lower levels are presumed to be "negative" responses; middle value presumed to be neutral
  ## returns a ggplot2 object of a net stacked distribution plot
  
  ## Test that all elements of x have the same levels, are ordered, etc.
  all_levels <- levels(x[[1]])
  n <- length(all_levels)
  levelscheck <- all(sapply(x, function(y)
    all(c(is.ordered(y), levels(y) == all_levels))
  ))
  if(!levelscheck)
    stop("All levels of x must be ordered factors with the same levels")
  
  ## Reverse order of columns (to make ggplot2 output look right after coord_flip)
  x <- x[length(x):1]
  
  ## Identify middle and "negative" levels
  if(n %% 2 == 1)
    neutral <- all_levels[ceiling(n/2)]
  else
    neutral <- NULL
  
  negatives <- all_levels[1:floor(n/2)]
  positives <- setdiff(all_levels, c(negatives, neutral))
  
  ## remove neutral, summarize as proportion
  listall <- lapply(names(x), function(y) {
    column <- (na.omit(x[[y]]))
    out <- data.frame(Sample = y, prop.table(table(column)))
    names(out) <- c("Sample", "Frequency_Change", "Freq")
    
    if(!is.null(neutral))
      out <- out[out$Frequency_Change != neutral,]
    
    out
  })
  
  dfall <- do.call(rbind, listall)
  
  ## split by positive/negative
  pos <- dfall[dfall$Frequency_Change %in% positives,]
  neg <- dfall[dfall$Frequency_Change %in% negatives,]
  
  ## Negate the frequencies of negative responses, reverse order
  neg$Freq <- -neg$Freq
  #print(neg$Frequency_Change)
  neg$Frequency_Change <- ordered(neg$Frequency_Change, levels = sort(levels(neg$Frequency_Change), decreasing=TRUE))
  print(summary(pos))
  print(summary(neg))
  
  
  stackedchart <- ggplot() +
    
    geom_bar(data = neg, aes(Sample, Freq, fill = Frequency_Change, order = Frequency_Change), stat = "identity") +
    geom_bar(data = pos, aes(Sample, Freq, fill = Frequency_Change, order = Frequency_Change), stat = "identity") + 
    geom_hline(yintercept=0) +
    # include scale showing % freq change
    scale_y_continuous(name = "",
                       labels = paste0(seq(-100, 100, 20), "%"),
                       limits = c(-1, 1),
                       breaks = seq(-1, 1, .2)) +
    #scale_x_discrete(limits=c(-1, 12)) +
    # include scale showing just 'divergent' and 'convergent'
   # scale_y_continuous(name = "",
  #                     #labels = c("divergent", "convergent"),
  #                     labels= rep("", times=11),
  #                     limits = c(-1, 1),
  #                     breaks = seq(-1, 1, .2)) +
    scale_fill_manual(values=c(colors()[190], colors()[153],
                      colors()[230], colors()[253])) + 
    #ggtitle(paste("sig. SNPs in", j, ", n = ", n_to_add)) + 
    annotate("text", y=-0.95, x=10, label=letter_label_to_add, hjust=0, size=8) +
    annotate("text", y=-1, x=1, label=paste("n =", n_to_add), hjust=0) +
    annotate("text", y=-0.5, x=0, label="divergent") +
    annotate("text", y=0.5, x=0, label="convergent") +
    theme(legend.position="none", panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(), 
          plot.background = element_rect(fill = "transparent",colour = NA)) +
    annotate("rect", xmin=0, xmax=5.5, ymin=-Inf, ymax=Inf, fill="red", alpha=0.1) +
    annotate("rect", xmin=5.5, xmax=11, ymin=-Inf, ymax=Inf, fill = "blue", alpha=0.1) +   
    
    #scale_fill_discrete(limits = c(negatives, positives)) +
    coord_flip()

  stackedchart
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

###############
###############
###############

db <- dbConnect(SQLite(), dbname="Freq.sqlite")

# Had to make the SQL database on barcoo as it took up a great deal of memory
# (probably > 100 Gb all up?)

# dbWriteTable(conn = db, name = "C1", value = "../variant_calls/C1_all_freq_change_150821.txt",
#             sep="\t", row.names = FALSE, header = TRUE)
# dbListFields(db, "C1")
# dbWriteTable(conn = db, name = "C2", value = "../variant_calls/C2_all_freq_change_150821.txt",
#             sep="\t", row.names = FALSE, header = TRUE)
# dbWriteTable(conn = db, name = "C3", value = "../variant_calls/C3_all_freq_change_150821.txt",
#             sep="\t", row.names = FALSE, header = TRUE)
# dbWriteTable(conn = db, name = "C4", value = "../variant_calls/C4_all_freq_change_150821.txt",
#             sep="\t", row.names = FALSE, header = TRUE)
# dbWriteTable(conn = db, name = "C5", value = "../variant_calls/C5_all_freq_change_150821.txt",
#             sep="\t", row.names = FALSE, header = TRUE)
# print(paste("imported C tables", Sys.time()))
# dbWriteTable(conn = db, name = "D1", value = "../variant_calls/D1_all_freq_change_150821.txt",
#             sep="\t", row.names = FALSE, header = TRUE)
# dbWriteTable(conn = db, name = "D2", value = "../variant_calls/D2_all_freq_change_150821.txt",
#             sep="\t", row.names = FALSE, header = TRUE)
# dbWriteTable(conn = db, name = "D3", value = "../variant_calls/D3_all_freq_change_150821.txt",
#             sep="\t", row.names = FALSE, header = TRUE)
# dbWriteTable(conn = db, name = "D4", value = "../variant_calls/D4_all_freq_change_150821.txt",
#             sep="\t", row.names = FALSE, header = TRUE)
# dbWriteTable(conn = db, name = "D5", value = "../variant_calls/D5_all_freq_change_150821.txt",
#             sep="\t", row.names = FALSE, header = TRUE)
# print(paste("imported D tables", Sys.time()))
# 
# dbListFields(db, "C1")
# #dbReadTable(db, "C1")
# 
# dbSendQuery(db, 'CREATE INDEX C1_index ON C1 (`chr`, `pos`)')
# dbSendQuery(db, 'CREATE INDEX C2_index ON C2 (`chr`, `pos`)')
# dbSendQuery(db, 'CREATE TABLE C1_C2 AS SELECT * FROM C1 LEFT OUTER JOIN C2 ON (C1.chr = C2.chr AND C1.pos=C2.pos)')
# dbSendQuery(db, 'DROP TABLE C1')
# dbSendQuery(db, 'DROP TABLE C2')
# dbSendQuery(db, 'CREATE INDEX C1_C2_index ON C1_C2 (`chr`, `pos`)')
# print(paste("merged and indexed C1_C2", Sys.time()))
# 
# dbSendQuery(db, 'CREATE INDEX C3_index ON C3 (`chr`, `pos`)')
# dbSendQuery(db, 'CREATE INDEX C4_index ON C4 (`chr`, `pos`)')
# dbSendQuery(db, 'CREATE TABLE C3_C4 AS SELECT * FROM C3 LEFT OUTER JOIN C4 ON (C3.chr = C4.chr AND C3.pos=C4.pos)')
# dbSendQuery(db, 'DROP TABLE C3')
# dbSendQuery(db, 'DROP TABLE C4')
# dbSendQuery(db, 'CREATE INDEX C3_C4_index ON C3_C4 (`chr`, `pos`)')
# print(paste("merged and indexed C3_C4", Sys.time()))
# 
# 
# dbSendQuery(db, 'CREATE INDEX C5_index ON C5 (`chr`, `pos`)')
# dbSendQuery(db, 'CREATE INDEX D1_index ON D1 (`chr`, `pos`)')
# dbSendQuery(db, 'CREATE TABLE C5_D1 AS SELECT * FROM C5 LEFT OUTER JOIN D1 ON (C5.chr = D1.chr AND C5.pos=D1.pos)')
# dbSendQuery(db, 'DROP TABLE C5')
# dbSendQuery(db, 'DROP TABLE D1')
# dbSendQuery(db, 'CREATE INDEX C5_D1_index ON C5_D1 (`chr`, `pos`)')
# print(paste("merged and indexed C5_D1", Sys.time()))
# 
# 
# dbSendQuery(db, 'CREATE INDEX D2_index ON D2 (`chr`, `pos`)')
# dbSendQuery(db, 'CREATE INDEX D3_index ON D3 (`chr`, `pos`)')
# dbSendQuery(db, 'CREATE TABLE D2_D3 AS SELECT * FROM D2 LEFT OUTER JOIN D3 ON (D2.chr = D3.chr AND D2.pos=D3.pos)')
# dbSendQuery(db, 'DROP TABLE D2')
# dbSendQuery(db, 'DROP TABLE D3')
# dbSendQuery(db, 'CREATE INDEX D2_D3_index ON D2_D3 (`chr`, `pos`)')
# print(paste("merged and indexed D2_D3", Sys.time()))
# 
# 
# dbSendQuery(db, 'CREATE INDEX D4_index ON D4 (`chr`, `pos`)')
# dbSendQuery(db, 'CREATE INDEX D5_index ON D5 (`chr`, `pos`)')
# dbSendQuery(db, 'CREATE TABLE D4_D5 AS SELECT * FROM D4 LEFT OUTER JOIN D5 ON (D4.chr = D5.chr AND D4.pos=D5.pos)')
# dbSendQuery(db, 'DROP TABLE D4')
# dbSendQuery(db, 'DROP TABLE D5')
# dbSendQuery(db, 'CREATE INDEX D4_D5_index ON D4_D5 (`chr`, `pos`)')
# print(paste("merged and indexed D4_D5", Sys.time()))
# 
# 
# dbSendQuery(db, 'CREATE TABLE C1_C2_C3_C4 AS SELECT * FROM C1_C2 LEFT OUTER JOIN C3_C4 ON (C1_C2.chr = C3_C4.chr AND C1_C2.pos=C3_C4.pos)')
# dbSendQuery(db, 'DROP TABLE C1_C2')
# dbSendQuery(db, 'DROP TABLE C3_C4')
# dbSendQuery(db, 'CREATE INDEX C1_C2_C3_C4_index ON C1_C2_C3_C4 (`chr`, `pos`)')
# print(paste("merged and indexed C1 to C4", Sys.time()))
# 
# dbSendQuery(db, 'CREATE TABLE C5_D1_D2_D3 AS SELECT * FROM C5_D1 LEFT OUTER JOIN D2_D3 ON (C5_D1.chr = D2_D3.chr AND C5_D1.pos=D2_D3.pos)')
# dbSendQuery(db, 'DROP TABLE C5_D1')
# dbSendQuery(db, 'DROP TABLE D2_D3')
# dbSendQuery(db, 'CREATE INDEX C5_D1_D2_D3_index ON C5_D1_D2_D3 (`chr`, `pos`)')
# print(paste("merged and indexed C5 to D3", Sys.time()))
# 
# dbSendQuery(db, 'CREATE TABLE C1_to_D3 AS SELECT * FROM C1_C2_C3_C4 LEFT OUTER JOIN C5_D1_D2_D3 ON (C1_C2_C3_C4.chr = C5_D1_D2_D3.chr AND C1_C2_C3_C4.pos=C5_D1_D2_D3.pos)')
# dbSendQuery(db, 'DROP TABLE C1_C2_C3_C4')
# dbSendQuery(db, 'DROP TABLE C5_D1_D2_D3')
# dbSendQuery(db, 'CREATE INDEX C1_to_D3_index ON C1_to_D3 (`chr`, `pos`)')
# print(paste("merged and indexed C1 to D3", Sys.time()))
# 
# 
# dbSendQuery(db, 'CREATE TABLE all_merged AS SELECT * FROM C1_to_D3 LEFT OUTER JOIN D4_D5 ON (C1_to_D3.chr = D4_D5.chr AND C1_to_D3.pos=D4_D5.pos)')
# 
# 
# dbSendQuery(db, 'DROP TABLE C1_to_D3')
# dbSendQuery(db, 'DROP TABLE D4_D5')
# 
# dbSendQuery(db, 'CREATE INDEX all_merged_index ON all_merged (`chr`, `pos`)')
# print(paste("merged and indexed all 10 tables", Sys.time()))


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
  querystring <- paste("SELECT * FROM all_merged ORDER BY RANDOM() LIMIT ", number_loci, ";", sep="")
  queryresult <- dbGetQuery(db, querystring)
  output_matrix <- unlist(tempqueryresult)
  colnames(output_matrix) <- unlist(dbListFields(db, 'all_merged'))
  return(output_matrix)
}



D_high <- read.table('/vlsci/VR0267/pgriffin/hsm/sandra/output_wgs/sql_database/D_high_effects.txt',
                     stringsAsFactors=FALSE, sep="\t", header=TRUE)
d_high_table <- extract_certain_loci(D_high, chr_col=2, pos_col=3)
write.table(d_high_table, file="D_high_allele_freq.txt", sep="\t", 
            row.names=FALSE, col.names=FALSE)

C_mod_high <- read.table('/vlsci/VR0267/pgriffin/hsm/sandra/output_wgs/sql_database/C_mod_high_effects.txt',
                     stringsAsFactors=FALSE, sep="\t", header=TRUE)
c_mod_high_table <- extract_certain_loci(C_mod_high, chr_col=2, pos_col=3)
write.table(c_mod_high_table, file="C_mod_high_allele_freq.txt", sep="\t", 
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

for(i in Sample_code[6:10]){
  print(paste("Processing no-C sig. locus list for", i))
  temp_filename <- paste("noC_", i, "_sig_SNPs.txt", sep="")
  temp_outputname <- paste(i, "_noC_freq_table.txt", sep="")
  temp_table <- read.table(temp_filename, header=FALSE, sep="\t")
  temp_output <- extract_certain_loci(get(temp_table), chr_col=1, pos_col=2)
  write.table(temp_output, file=temp_outputname, sep="\t",
              row.names=FALSE, col.names=FALSE)
}


# other frequency tables to get
# frequency table for each list of sig. SNPs (i.e. 10 in total)
# frequency table for 100,000 randomly-chosen loci
# 

### code to interpret output of database queries ###




D_high <- read.table('/Users/pgriffin/Documents/Drosophila\ Selection\ Experiment/snpeff_SNP_feature_enrichment/D_high_effects.txt',
                     stringsAsFactors=FALSE, sep="\t", header=TRUE)
setwd('/Users/pgriffin/Documents/Drosophila\ Selection\ Experiment/allele_frequency_difference_testing')
d_high <- read.table('D_high_allele_freq.txt', header=FALSE, sep="\t",
                     stringsAsFactors=FALSE)
colnames(d_high) <- query_head



d_useful <- useful_subset(d_high)

d_all_info <- merge(D_high, d_useful, by.x=c('chr', 'pos'), by.y=c('chr', 'pos'))

all_locus_output <- list()
for(i in 1:nrow(d_all_info)){
  temp_row <- d_all_info[i,]
  temp_sample <- temp_row[,'sample']
  if(length(strsplit(temp_sample, ',')[[1]]) > 1){
    temp_sample <- strsplit(temp_sample, ',')[[1]][1]
  }
  print(paste(temp_row$chr, temp_row$pos))

  MB_result <- unlist(c('MB', temp_row[10:12]))
  
  other_result_loc <- c(10, 11, 13)
  other_result_locations <- c(rep(other_result_loc, times=10) + rep(0:9*4, each=3))

  other_results <- temp_row[,other_result_locations]

  output_table <- matrix(c(MB_result[2:4], other_results), ncol=11, nrow=3)
  colnames(output_table) <- c('MB', Sample_code)
  
  #now identify the columns with alleles in the wrong order and swap
  correct_order <- output_table[1:2,colnames(output_table)==temp_sample]
  temp_output <- output_table
  cols_to_switch <- which(temp_output[2,]==correct_order[[1]])
  output_table[1,cols_to_switch] <- temp_output[2,cols_to_switch]
  output_table[2,cols_to_switch] <- temp_output[1,cols_to_switch]
  output_table[3,cols_to_switch] <- 1-as.numeric(temp_output[3, cols_to_switch])
  
  all_locus_output[[i]] <- list(temp_row[,c(1:6)], output_table)
}

pdf(file="Plots of high-effect D loci allele freq changes.pdf", height=36, width=40)
par(mfrow=c(6, 8))
for(i in 1:length(all_locus_output)){
  temp_sample <- as.character(all_locus_output[[i]][[1]][3])
  if(length(strsplit(temp_sample, split=', ')[[1]]) > 1){
    temp_sample <- unlist(strsplit(temp_sample, ', '))
  }
  temp_chr_pos <- all_locus_output[[i]][[1]][1:2]
  temp_freq_vector <- all_locus_output[[i]][[2]][3,]
  start <- as.numeric(temp_freq_vector[1])
  start_end <- data.frame(rep(start, 10), unlist(temp_freq_vector[2:11]), 
                          c(rep('lightblue', times=5), rep('pink', times=5)),
                          rep(1, times=10), rep(2, times=10))
  start_end[,3] <- as.character(start_end[,3])
  start_end[which(rownames(start_end)%in%as.character(temp_sample)), 3] <- 'darkred'
  colnames(start_end) <- c("starty", "endy", "colour", "startx", "endx")
  plot(x=c(1, 2), y=c(0, 1), type="p", col="white",
       main=paste(temp_chr_pos), xaxt="n",
       ylab="Allele frequency", xlab="Change Gen0--Gen13")
  segments(x0=start_end$startx, x1=start_end$endx,
           y0=start_end$starty, y1=start_end$endy,
           col=start_end$colour, lwd=2)
}
dev.off()


#######################
#######################

C_high <- read.table('/Users/pgriffin/Documents/Drosophila\ Selection\ Experiment/snpeff_SNP_feature_enrichment/C_mod_high_effects.txt',
                     stringsAsFactors=FALSE, sep="\t", header=TRUE)
setwd('/Users/pgriffin/Documents/Drosophila\ Selection\ Experiment/allele_frequency_difference_testing')
c_high <- read.table('C_mod_high_allele_freq.txt', header=FALSE, sep="\t",
                     stringsAsFactors=FALSE)
colnames(c_high) <- query_head

useful_subset <- function(database_query_output){
  individual_col_names <- paste(rep(Sample_code, each=4), 
                                rep(c('_dec', '_inc', '_startfreq', '_endfreq'), times=10), 
                                sep="")
  useful_col_names <- c('chr', 'pos', individual_col_names)
  useful_output <- database_query_output[,colnames(database_query_output)%in%useful_col_names]
}

c_useful <- useful_subset(c_high)

c_all_info <- merge(C_high, c_useful, by.x=c('chr', 'pos'), by.y=c('chr', 'pos'))

all_locus_output <- list()
for(i in 1:nrow(c_all_info)){
  temp_row <- c_all_info[i,]
  temp_sample <- temp_row[,'sample']
  if(length(strsplit(temp_sample, ',')[[1]]) > 1){
    temp_sample <- strsplit(temp_sample, ',')[[1]][1]
  }
  print(paste(temp_row$chr, temp_row$pos))
  
  MB_result <- unlist(c('MB', temp_row[10:12]))
  
  other_result_loc <- c(10, 11, 13)
  other_result_locations <- c(rep(other_result_loc, times=10) + rep(0:9*4, each=3))
  
  other_results <- temp_row[,other_result_locations]
  
  output_table <- matrix(c(MB_result[2:4], other_results), ncol=11, nrow=3)
  colnames(output_table) <- c('MB', Sample_code)
  
  #now identify the columns with alleles in the wrong order and swap
  correct_order <- output_table[1:2,colnames(output_table)==temp_sample]
  temp_output <- output_table
  cols_to_switch <- which(temp_output[2,]==correct_order[[1]])
  output_table[1,cols_to_switch] <- temp_output[2,cols_to_switch]
  output_table[2,cols_to_switch] <- temp_output[1,cols_to_switch]
  output_table[3,cols_to_switch] <- 1-as.numeric(temp_output[3, cols_to_switch])
  
  all_locus_output[[i]] <- list(temp_row[,c(1:6)], output_table)
}

pdf(file="Plots of high- and moderate-effect C loci allele freq changes.pdf", height=40, width=40)
par(mfrow=c(7, 8))
for(i in 1:length(all_locus_output)){
  temp_sample <- as.character(all_locus_output[[i]][[1]][3])
  if(length(strsplit(temp_sample, split=', ')[[1]]) > 1){
    temp_sample <- unlist(strsplit(temp_sample, ', '))
  }
  temp_chr_pos <- all_locus_output[[i]][[1]][1:2]
  temp_freq_vector <- all_locus_output[[i]][[2]][3,]
  start <- as.numeric(temp_freq_vector[1])
  start_end <- data.frame(rep(start, 10), unlist(temp_freq_vector[2:11]), 
                          c(rep('lightblue', times=5), rep('pink', times=5)),
                          rep(1, times=10), rep(2, times=10))
  start_end[,3] <- as.character(start_end[,3])
  start_end[which(rownames(start_end)%in%as.character(temp_sample)), 3] <- 'darkblue'
  colnames(start_end) <- c("starty", "endy", "colour", "startx", "endx")
  plot(x=c(1, 2), y=c(0, 1), type="p", col="white",
       main=paste(temp_chr_pos), xaxt="n",
       ylab="Allele frequency", xlab="Change Gen0--Gen13")
  segments(x0=start_end$startx, x1=start_end$endx,
           y0=start_end$starty, y1=start_end$endy,
           col=start_end$colour, lwd=2)
}
dev.off()


##################
#
# Now process each of the per-sample sig loci frequency tables...
#
##################

setwd('/Users/pgriffin/Documents/Drosophila\ Selection\ Experiment/allele_frequency_difference_testing')




stats_table <- data.frame(focal_sample=character(), focal_sample_cat=character(),
                          comp_sample=character(), comp_sample_cat=character(),
                          prop_loci_above_zero=numeric(), stringsAsFactors=FALSE)
t_test_table <- data.frame(focal_sample=character(), groupCmean=numeric(), 
                           groupDmean=numeric(), pval=numeric(), stringsAsFactors=FALSE)
for(j in Sample_code){
  temp_freq <- read.table(paste(j, '_freq_table.txt', sep=""), header=FALSE, sep="\t",
                        stringsAsFactors=FALSE)
  colnames(temp_freq) <- query_head
  
  temp_useful <- useful_subset(temp_freq)

  all_locus_output <- list()
  for(i in 1:nrow(temp_useful)){
    temp_row <- temp_useful[i,]
    temp_sample <- j
    #print(paste(temp_row$chr, temp_row$pos))
    MB_result <- unlist(c('MB', temp_row[3:5]))
    
    other_result_loc <- c(3, 4, 6)
    other_result_locations <- c(rep(other_result_loc, times=10) + rep(0:9*4, each=3))
    other_results <- temp_row[,other_result_locations]
    
    output_table <- matrix(c(MB_result[2:4], other_results), ncol=11, nrow=3)
    colnames(output_table) <- c('MB', Sample_code)
    
    #now identify the columns with alleles in the wrong order and swap
    correct_order <- output_table[1:2,colnames(output_table)==temp_sample]
    temp_output <- output_table
    cols_to_switch <- which(temp_output[2,]==correct_order[[1]])
    output_table[1,cols_to_switch] <- temp_output[2,cols_to_switch]
    output_table[2,cols_to_switch] <- temp_output[1,cols_to_switch]
    output_table[3,cols_to_switch] <- 1-as.numeric(temp_output[3, cols_to_switch])
    
    freqchange_vector <- as.numeric(output_table[3,2:11])-as.numeric(output_table[3,1])
    names(freqchange_vector) <- Sample_code
    all_locus_output[[i]] <- list(temp_row[,c(1:6)], output_table, freqchange_vector)
  }
  temp_freqchange_table <- data.frame(lapply(all_locus_output, "[[", 3))
  info_chr <- lapply(lapply(all_locus_output, "[[", 1), "[[", 1)
  info_pos <- lapply(lapply(all_locus_output, "[[", 1), "[[", 2)
  colnames(temp_freqchange_table) <- paste(info_chr, info_pos, sep="_") 
  temp_freqchange_table <- t(temp_freqchange_table)
  
  #make likert table: convert freq change into pos or neg for the non-focal samples
  br <- c(-1, -0.1, 0, 0.1, 1)
  #prelik <- as.data.frame(temp_freqchange_table)
  prelik <- as.data.frame(apply(temp_freqchange_table, MARGIN=2, 
                  FUN=cut, breaks=br, labels=c("-1 < x < -0.1", "-0.1 < x < 0", "0 < x < 0.1", "0.1 < x < 1"),
                  ordered_result=TRUE))

  for(i in 1:10){
    prelik[,i]<-ordered(prelik[,i],
                        levels=c("-1 < x < -0.1", "-0.1 < x < 0", "0 < x < 0.1", "0.1 < x < 1"))
  }
  row.names(prelik) <- row.names(temp_freqchange_table)
  temp_likert <- net_stacked(prelik, j=j, n_to_add=nrow(prelik), letter_label_to_add=LETTERS[which(Sample_code==j)])
  assign(paste(j, "likert", sep="_"), temp_likert)
  
  #add info to table for statistical test
  if(j%in%Sample_code[1:5]){ focal_category <- "C"} else {focal_category <- "D"}
  temp_cols <- prelik[,-which(colnames(prelik)==j)]
  for(k in 1:9){
    temp_comparison_col <- temp_cols[,k]
    comp_sample <- colnames(temp_cols)[k]
    if(comp_sample%in%Sample_code[1:5]){ focal_comp_category <- "C"} else {focal_comp_category <- "D"}
    pos_proportion <- length(which(temp_comparison_col%in%c("0 < x < 0.1", "0.1 < x < 1")))/length(temp_comparison_col)
    new_row <- list(j, focal_category, comp_sample, focal_comp_category, pos_proportion)
    if(nrow(stats_table) < 1){
      stats_table[1,] <- new_row
    } else {
      stats_table <- rbind(stats_table, new_row)
    }
  }
  stats_for_ttest <- subset(stats_table, stats_table$focal_sample==j)
  ttest <- t.test(prop_loci_above_zero~comp_sample_cat, data=stats_for_ttest)
  groupCmean <- ttest$estimate[1]
  groupDmean <- ttest$estimate[2]
  pval <- ttest$p.value
  new_row2 <- list(j, groupCmean, groupDmean, pval)
  if(nrow(t_test_table) < 1){
    t_test_table[1,] <- new_row2
  } else {
    t_test_table <- rbind(t_test_table, new_row2)
  }  
  assign(paste(j, "sig_loci_freq_list", sep="_"), all_locus_output)
}


# extract freq change for each focal line

freq_change_df <- data.frame(sample=character(), freq_change <- numeric())
for(i in Sample_code){
  temp_freq_list <- get(paste(i, "_sig_loci_freq_list", sep=""))
  temp_focal_freq_change_vector_name <- paste(i, "_focal_freq_change_vector", sep="")
  iloc <- which(Sample_code==i)
  temp_freq_change_vector <- lapply(temp_freq_list, "[[", 3)
  temp_focal_freq_change_vector <- sapply(temp_freq_change_vector, "[[", iloc)
  assign(temp_focal_freq_change_vector_name, temp_focal_freq_change_vector)
  print(i)
  print(summary(temp_focal_freq_change_vector))
  freq_change_df <- rbind(freq_change_df, data.frame(sample=rep(i, times=length(temp_focal_freq_change_vector)),
                                                     freq_change=temp_focal_freq_change_vector))
}

pdf(file="Candidate SNP Frequency Change Boxplots.pdf",
    width=6, height=4)
boxplot(freq_change~sample, data=freq_change_df,
        border=c(rep(rgb(0,0,255,100, maxColorValue=255), times=5), 
                 rep(rgb(255,0,0,100, maxColorValue=255), times=5)),
        ylim=c(0, 1), cex=0.6, pch=19,
        xlab="Replicate", ylab="Frequency change for candidate SNPs")
dev.off()




# write test results

write.table(stats_table, file="Proportion_loci_changing_in_same_dirn_as_focal_sample.txt",
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(t_test_table, file="T_test_results_for_individual_focal_sample_SNP_sets.txt",
            sep="\t", quote=FALSE, row.names=FALSE)


#### Likert plot for all 10 samples

pdf("Likert plots for all sig. SNP sets.pdf", width=11, height=24)
multiplot(C1_likert, C2_likert, C3_likert, C4_likert, C5_likert,
          D1_likert, D2_likert, D3_likert, D4_likert, D5_likert,
          cols=2)
dev.off()

pdf("Likert plots for all sig. SNP sets horizontal.pdf", width=24, height=11)
multiplot(C1_likert, D1_likert, C2_likert, D2_likert, C3_likert,
          D3_likert, C4_likert, D4_likert, C5_likert, D5_likert,
          cols=5)
dev.off()

pdf("Likert plot for just D1.pdf", width=5.5, height=5.5)
D1_likert
dev.off()


# ANOVA on proportion of loci showing positive freq change

aov1 <- aov(prop_loci_above_zero~focal_sample_cat*comp_sample_cat, data=stats_table)

glm1 <- glm(prop_loci_above_zero~focal_sample_cat*comp_sample_cat, data=stats_table)

#################
#
# Now identify the loci with the largest freq change
#
#################

D1_freq_change_vector <- unlist(lapply(lapply(D1_sig_loci_freq_list, "[[", 3), "[", 6))

length(which(D1_freq_change_vector>0.8))
D1_high_freq_change <- unlist(lapply(lapply(D1_sig_loci_freq_list, "[[", 1), "[", c(1,2))[which(D1_freq_change_vector>0.8)])
D1_high_freq_change <- matrix(D1_high_freq_change, ncol=2, byrow=TRUE)
D1_high_freq_change <- apply(D1_high_freq_change, 1, paste, collapse=":")

##################
#
# Now process the randomly-selected loci
#
##################

##### WAIT FOR TABLE TO BE REMADE: ISSUE WITH SELECTING COLUMNS
##### THAT CONTAIN NULL VALUES: DATA TYPE IS SET TO SOMETHING WEIRD
##### BY DEFAULT

setwd("/Users/pgriffin/Documents/Drosophila Selection Experiment/allele_frequency_difference_testing")

temp_freq <- read.table("Random_100000_loci_freq.txt", header=FALSE, sep="\t",
                        stringsAsFactors=FALSE)
#temp_freq <- as.matrix(temp_freq)
#temp_freq <- t(temp_freq)
#temp_freq <- unlist(temp_freq)
#temp_freq <- matrix(temp_freq, ncol=90, byrow=FALSE)
colnames(temp_freq) <- query_head

temp_useful <- useful_subset(temp_freq)
#temp_useful <- temp_useful[temp_useful[,1]!="V1",]

to_remove <- c()
all_locus_output <- list()
for(i in 1:nrow(temp_useful)){
  temp_row <- temp_useful[i,]
  #temp_sample <- j
  if(i%%1000==0){
    print(paste(temp_row[1:2], collapse=" "))
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


temp_freqchange_table <- data.frame(lapply(all_locus_output, "[[", 3))
info_chr <- lapply(lapply(all_locus_output, "[[", 1), "[[", 1)
info_pos <- lapply(lapply(all_locus_output, "[[", 1), "[[", 2)
colnames(temp_freqchange_table) <- paste(info_chr, info_pos, sep="_") 
temp_freqchange_table <- t(temp_freqchange_table)
  

temp_jpg_name <- "Random_loci_freq_change_across_samples.jpg"
jpeg(file=temp_jpg_name, width=45, height=45, units="cm", res=300)
par(mfrow=c(10,10))
for(i in Sample_code){
  temp_base_sample <- i
  for(j in 1:10){
    print(paste("Plotting", i, "versus", Sample_code[j]))
    temp_comparison_sample <- Sample_code[j]
    plot(temp_freqchange_table[,j], 
         temp_freqchange_table[,which(colnames(temp_freqchange_table)==i)], 
         cex=0.4, 
         xlim=c(-1, 1), 
         ylim=c(-1, 1),
         ylab=paste("Frequency change in", temp_base_sample),
         xlab=paste("Frequency change in",  temp_comparison_sample),
         pch=19, col=rgb(100, 100, 100, 50, maxColorValue=255))
    
    lines(y=seq(-1, 1, length=1000), x=seq(-1, 1, length=1000), col="red")
    lines(y=rep(0.3, times=1000), x=seq(0.3, 1, length=1000), col="grey", lty=2)
    lines(x=rep(0.3, times=1000), y=seq(0.3, 1, length=1000), col="grey", lty=2)
  }
}
dev.off()

tes <- subset(temp_freqchange_table, temp_freqchange_table[,7] < -0.8)
head(tes)



temp_jpg_name <- "Random_loci_freq_change_histograms.jpg"
jpeg(file=temp_jpg_name, width=45, height=45, units="cm", res=300)
par(mfrow=c(9, 9))
for(i in Sample_code){
  temp_base_sample <- i
  for(j in (1:10)[-which(Sample_code==i)]){
  hist(temp_freqchange_table[,i], freq=TRUE, breaks=100, 
       col=rgb(100, 100, 100, 100, maxColorValue=255), 
       border=rgb(100, 100, 100, 100, maxColorValue=255),
       main=paste("Frequency change in", temp_base_sample, "sig. loci"),
       ylab="Locus count per bin", xlim=c(-1, 1),
       ylim=c(0, nrow(temp_freqchange_table)/5),
       xlab=paste("vs.", Sample_code[j]))
 hist(temp_freqchange_table[,j], freq=TRUE, breaks=100, 
      col=rgb(100, 0, 0, 100, maxColorValue=255), 
      border=rgb(100, 0, 0, 100, maxColorValue=255),
      add=TRUE)
    #bins <- c(-1, -0.25, 0, 0.25, 1) 
  }
}
dev.off()




