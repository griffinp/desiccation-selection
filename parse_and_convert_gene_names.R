###########################################
# After running SNPEff on each vcf file   #
# and using snpSift to extract gene names #
###########################################

Sample_code <- c("C1", "C2", "C3", "C4", "C5",
                 "D1", "D2", "D3", "D4", "D5")

setwd("~/Documents/Drosophila Selection Experiment/snp_and_gene_lists")

# Now look at the results
colname_vector <- c("chr", "pos", paste(rep("genename", times=88), 1:88, sep=""))

# Object naming is a bit confusing: I was previously using 'all' to refer to 
# 'including SNPs that are in the C-replicates' while also before I was 
# using it to refer to 'both significant and non-significant loci'
# have now tried to change this to 'sig' wherever necessary

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

C1_sig_positions <- paste(C1_sig_ex[,1], C1_sig_ex[,2])
C2_sig_positions <- paste(C2_sig_ex[,1], C2_sig_ex[,2])
C3_sig_positions <- paste(C3_sig_ex[,1], C3_sig_ex[,2])
C4_sig_positions <- paste(C4_sig_ex[,1], C4_sig_ex[,2])
C5_sig_positions <- paste(C5_sig_ex[,1], C5_sig_ex[,2])
D1_sig_positions <- paste(D1_sig_ex[,1], D1_sig_ex[,2])
D2_sig_positions <- paste(D2_sig_ex[,1], D2_sig_ex[,2])
D3_sig_positions <- paste(D3_sig_ex[,1], D3_sig_ex[,2])
D4_sig_positions <- paste(D4_sig_ex[,1], D4_sig_ex[,2])
D5_sig_positions <- paste(D5_sig_ex[,1], D5_sig_ex[,2])

C_sig_positions <-union(union(union(C1_sig_positions, C2_sig_positions), C3_sig_positions), 
                        union(C4_sig_positions, C5_sig_positions))

extract_gene_names_from_table <- function(snpsift_table){
  temp_gene_vector <- c()
  for(i in 3:90){
    temp_gene_vector <- c(temp_gene_vector, snpsift_table[,i])
  }
  gene_names <- temp_gene_vector[temp_gene_vector!=""]
  unique_gene_names <- unique(gene_names)
  #unique_gene_names <- unique_gene_names[unique_gene_names!=""]
  return(unique_gene_names[is.na(unique_gene_names)==FALSE])
}

C1_sig_gene_names <- extract_gene_names_from_table(C1_sig_ex)
C2_sig_gene_names <- extract_gene_names_from_table(C2_sig_ex)
C3_sig_gene_names <- extract_gene_names_from_table(C3_sig_ex)
C4_sig_gene_names <- extract_gene_names_from_table(C4_sig_ex)
C5_sig_gene_names <- extract_gene_names_from_table(C5_sig_ex)
D1_sig_gene_names <- extract_gene_names_from_table(D1_sig_ex)
D2_sig_gene_names <- extract_gene_names_from_table(D2_sig_ex)
D3_sig_gene_names <- extract_gene_names_from_table(D3_sig_ex)
D4_sig_gene_names <- extract_gene_names_from_table(D4_sig_ex)
D5_sig_gene_names <- extract_gene_names_from_table(D5_sig_ex)

### NB SnpEff has annotated (some?) cases where two genes overlap in space as 
### 'Exon_chr_start_end'. I think this is a problem with v6.01 rather than 
### with snpEff. However, for network analysis (etc.) these need to be
### converted to include all genes. 

## I made a 'translation table' by manually adding the gene names from Flybase

exons_for_conv <- read.table("exon_names_gene_list_for_conversion.txt", header=TRUE, 
                             fill=TRUE, stringsAsFactors=FALSE)

# Creating new objects that contain exon-parsed gene name vectors
# for each sample in the format 'C1_sig_gene_names_exons_parsed'
for(i in Sample_code){
  temp_object_name <- paste(i, "_sig_gene_names", sep="")
  temp_object <- get(temp_object_name)
  temp_grep <- grep('Exon', temp_object)
  temp_output_object <- temp_object
  if(length(temp_grep)>0){
    for(j in 1:length(temp_grep)){
      temp_exon_name <- grep('Exon', temp_object, value=TRUE)[j]
      temp_row <- exons_for_conv[exons_for_conv[1]==temp_exon_name]
      temp_output_object <- c(temp_output_object, temp_row[2:4])
    }
    temp_output_object <- temp_output_object[temp_output_object!=""]
    temp_output_object <- temp_output_object[-temp_grep]
  }
  unique_gene_names <- unique(temp_output_object[is.na(temp_output_object)==FALSE])
  assign(paste(i, "_sig_gene_names_exons_parsed", sep=""), unique_gene_names)
}


C_sig_gene_names <- union(union(union(C1_sig_gene_names_exons_parsed, C2_sig_gene_names_exons_parsed), 
                                C3_sig_gene_names_exons_parsed), 
                          union(C4_sig_gene_names_exons_parsed, C5_sig_gene_names_exons_parsed))


# List of D gene names excluding C gene names
D_sig_gene_names <- setdiff(union(union(union(D1_sig_gene_names_exons_parsed, D2_sig_gene_names_exons_parsed), 
                                    D3_sig_gene_names_exons_parsed), 
                              union(D4_sig_gene_names_exons_parsed, D5_sig_gene_names_exons_parsed)), 
                        C_sig_gene_names)


# C_sig_gene_list <- paste(C_sig_gene_names, collapse=" ")
# 
# C_sig_genes_in_D1 <- intersect(D1_sig_gene_names_exons_parsed, C_sig_gene_names)
# C_sig_genes_in_D2 <- intersect(D2_sig_gene_names_exons_parsed, C_sig_gene_names)
# C_sig_genes_in_D3 <- intersect(D3_sig_gene_names_exons_parsed, C_sig_gene_names)
# C_sig_genes_in_D4 <- intersect(D4_sig_gene_names_exons_parsed, C_sig_gene_names)
# C_sig_genes_in_D5 <- intersect(D5_sig_gene_names_exons_parsed, C_sig_gene_names)

D1_noC_sig_gene_names <- setdiff(D1_sig_gene_names_exons_parsed, C_sig_gene_names)
D2_noC_sig_gene_names <- setdiff(D2_sig_gene_names_exons_parsed, C_sig_gene_names)
D3_noC_sig_gene_names <- setdiff(D3_sig_gene_names_exons_parsed, C_sig_gene_names)
D4_noC_sig_gene_names <- setdiff(D4_sig_gene_names_exons_parsed, C_sig_gene_names)
D5_noC_sig_gene_names <- setdiff(D5_sig_gene_names_exons_parsed, C_sig_gene_names)

#############################
# Make conversion table for #
# conversion to FBgn format #
#############################

# Currently just using the 'sig' SNPEff output; will have to repeat
# with 'all' when those files are ready
# NB had to edit these files slightly to remove hash from start of header row
# can't rely on col numbers as these vary among files!

rm(gene_list_for_conversion)
for (i in 1:10){
  temp_sample <- Sample_code[i]
  temp_file_name <- paste("/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/snpEff_genes_", 
                          temp_sample, "_sig.txt", sep="")
  temp_file <- read.table(temp_file_name, header=TRUE, 
                          stringsAsFactors=FALSE, sep="\t",
                          nrows=18000, quote="\"")
  print(paste(nrow(temp_file), "rows in", temp_sample))
  #Fix the 'Exon' rows, which had a column missing
  ## NB apparently no longer necessary once I specified tab-delimited input
#   exon_row_numbers <- grep('Exon', temp_file[,2])
#   fixing_exon_rows <- cbind(temp_file[exon_row_numbers,1:2],
#                             rep(NA, times=length(exon_row_numbers)),
#                             temp_file[exon_row_numbers,3:ncol(temp_file)])
#   temp_file[exon_row_numbers,] <- fixing_exon_rows
  temp_subs <- temp_file[,c('GeneId', 'GeneName', 
                            'Length..DOWNSTREAM.',
                            'Length..GENE.',
                            'Length..UPSTREAM.')]
  if(temp_sample=="C1"){
    gene_list_for_conversion <- temp_subs
  } else {
    gene_list_for_conversion <- rbind(gene_list_for_conversion, temp_subs)
  }
}
colnames(gene_list_for_conversion) <- c("FBgn", "name", "downstream_length", 
                                        "gene_length", "upstream_length")
nrow(gene_list_for_conversion)
gene_list_for_conversion[,3] <- as.integer(gene_list_for_conversion[,3])
gene_list_for_conversion[,4] <- as.integer(gene_list_for_conversion[,4])
gene_list_for_conversion[,5] <- as.integer(gene_list_for_conversion[,5])
gene_list_for_conversion <- unique(gene_list_for_conversion)
nrow(gene_list_for_conversion)
gene_list_for_conversion$total_length <- gene_list_for_conversion$downstream_length + gene_list_for_conversion$gene_length + gene_list_for_conversion$upstream_length

### Have to replace brackets in (some) gene names with underscores
### to match the way these are writen in the snpEff output.
gtest <- gsub("(", "_", gene_list_for_conversion$name, fixed=TRUE)
gtest2 <- gsub(")", "_", gtest, fixed=TRUE)
gene_list_for_conversion$name <- gtest2

write.table(gene_list_for_conversion, file="All_sig_gene_list_for_name_conversion.txt",
            sep="\t", quote=FALSE, row.names=FALSE)


##### NOW REPEATING WITH 'ALL' FILES ######

#rm(gene_list_for_conversion)
for (i in 1:10){
  temp_sample <- Sample_code[i]
  temp_file_name <- paste("/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/snpEff_genes_", 
                          temp_sample, "_all.txt", sep="")
  temp_file <- read.table(temp_file_name, header=TRUE, 
                          stringsAsFactors=FALSE, sep="\t", quote="\"")
  print(paste(nrow(temp_file), "rows in", temp_sample))
  temp_subs <- temp_file[,c('GeneId', 'GeneName', 
                            'Length..DOWNSTREAM.',
                            'Length..GENE.',
                            'Length..UPSTREAM.')]
  if(temp_sample=="C1"){
    gene_list_for_conversion <- temp_subs
  } else {
    gene_list_for_conversion <- rbind(gene_list_for_conversion, temp_subs)
  }
}
colnames(gene_list_for_conversion) <- c("FBgn", "name", "downstream_length", 
                                        "gene_length", "upstream_length")
nrow(gene_list_for_conversion)
gene_list_for_conversion[,3] <- as.integer(gene_list_for_conversion[,3])
gene_list_for_conversion[,4] <- as.integer(gene_list_for_conversion[,4])
gene_list_for_conversion[,5] <- as.integer(gene_list_for_conversion[,5])
gene_list_for_conversion <- unique(gene_list_for_conversion)
nrow(gene_list_for_conversion)
gene_list_for_conversion$total_length <- gene_list_for_conversion$downstream_length + gene_list_for_conversion$gene_length + gene_list_for_conversion$upstream_length

### Have to replace brackets in (some) gene names with underscores
### to match the way these are writen in the snpEff output.
gtest <- gsub("(", "_", gene_list_for_conversion$name, fixed=TRUE)
gtest2 <- gsub(")", "_", gtest, fixed=TRUE)
gene_list_for_conversion$name <- gtest2

write.table(gene_list_for_conversion, file="All_gene_list_for_name_conversion.txt",
            sep="\t", quote=FALSE, row.names=FALSE)




##############################
# Convert gene names into    #
# FBgn format                #
##############################

# To do this you can have already converted the 'Exon_3L_2341122_2353455'
# gene names into their actual genes (see above, lines 76-104)

# In case you haven't, importing again:
exon_for_conv <- read.table("/Users/pgriffin/Documents/Drosophila\ Selection\ Experiment/snp_and_gene_lists/exon_names_gene_list_for_conversion.txt",
                               header=FALSE, stringsAsFactors=FALSE, fill=TRUE)

# Importing the name to FBgn conversion table as well (object created above):

gene_list_for_conversion <- read.table("/Users/pgriffin/Documents/Drosophila\ Selection\ Experiment/snp_and_gene_lists/All_sig_gene_list_for_name_conversion.txt",
                                       header=TRUE, stringsAsFactors=FALSE, sep="\t", quote="\"")

gene_name_to_FBgn <- function(to_translate, conversion_table, exon_translation){
  # A function to translate names from 'gene name' to 'FBgn' format.
  # First checks if there are any untranslated 'Exon' names and deals with them
  # if so
  if(length(grep('Exon', to_translate))>0){
    exons_in_list <- to_translate[grep('Exon', to_translate)]
    relevant_exons <- exon_translation[which(exon_translation[,1]%in%exons_in_list),]
    exon_translation <- unlist(relevant_exons[,2:4])
    non_exons <- to_translate[-(grep('Exon', to_translate))]
    relevant_genes <- conversion_table[which(conversion_table[,2]%in%non_exons),1]
    all_FBgn <- c(relevant_genes, exon_translation)
  }
  else{
    existing_FBgn <- to_translate[grep('FBgn', to_translate)]
    all_FBgn <- unique(c(conversion_table[which(conversion_table[,2]%in%to_translate),1], existing_FBgn))
    #print(paste('started with', length(to_translate), 'gene names',
    #      'and converted to', length(all_FBgn), 'FBgn names', sep=" "))
  }
  return(all_FBgn)
}

for(i in Sample_code[1:5]){
  # Convert C gene lists
  object_name <- paste(i, "sig_gene_names", sep="_")
  temp_list <- get(object_name)
  output_object_name <- paste(i, "sig_gene_names_FBgn", sep="_")
  temp_output <- gene_name_to_FBgn(temp_list,
                                   conversion_table=gene_list_for_conversion,
                                   exon_translation=exon_for_conv)
  assign(output_object_name, temp_output)
}
for(i in Sample_code[6:10]){
  # Convert D ('noC') gene lists
  object_name <- paste(i, "noC_sig_gene_names", sep="_")
  temp_list <- get(object_name)
  output_object_name <- paste(i, "noC_sig_gene_names_FBgn", sep="_")
  temp_output <- gene_name_to_FBgn(temp_list,
                                   conversion_table=gene_list_for_conversion,
                                   exon_translation=exon_for_conv)
  assign(output_object_name, temp_output)
}

allC_gene_list_FBgn_format <- gene_name_to_FBgn(C_sig_gene_names, 
                                                conversion_table=gene_list_for_conversion,
                                                exon_translation=exon_for_conv)

######################
# Saving gene lists  #
######################

write.table(C1_sig_gene_names_exons_parsed, file="C1_sig_gene_list.txt", quote=FALSE,
            row.names=FALSE, col.names=FALSE)
write.table(C2_sig_gene_names_exons_parsed, file="C2_sig_gene_list.txt", quote=FALSE,
            row.names=FALSE, col.names=FALSE)
write.table(C3_sig_gene_names_exons_parsed, file="C3_sig_gene_list.txt", quote=FALSE,
            row.names=FALSE, col.names=FALSE)
write.table(C4_sig_gene_names_exons_parsed, file="C4_sig_gene_list.txt", quote=FALSE,
            row.names=FALSE, col.names=FALSE)
write.table(C5_sig_gene_names_exons_parsed, file="C5_sig_gene_list.txt", quote=FALSE,
            row.names=FALSE, col.names=FALSE)

write.table(C1_sig_gene_names_FBgn, file="C1_sig_gene_names_FBgn_format.txt", quote=FALSE,
            row.names=FALSE, col.names=FALSE)
write.table(C2_sig_gene_names_FBgn, file="C2_sig_gene_names_FBgn_format.txt", quote=FALSE,
            row.names=FALSE, col.names=FALSE)
write.table(C3_sig_gene_names_FBgn, file="C3_sig_gene_names_FBgn_format.txt", quote=FALSE,
            row.names=FALSE, col.names=FALSE)
write.table(C4_sig_gene_names_FBgn, file="C4_sig_gene_names_FBgn_format.txt", quote=FALSE,
            row.names=FALSE, col.names=FALSE)
write.table(C5_sig_gene_names_FBgn, file="C5_sig_gene_names_FBgn_format.txt", quote=FALSE,
            row.names=FALSE, col.names=FALSE)

write.table(C_sig_gene_names, file="Putative_lab_adaptation_gene_list.txt",
            row.names=FALSE, col.names=FALSE, quote=FALSE)

write.table(allC_gene_list_FBgn_format, file="Putative_lab_adaptation_gene_list_FBgn_format.txt",
            row.names=FALSE, col.names=FALSE, quote=FALSE)

write.table(D1_sig_gene_names_exons_parsed, 
            file="D1_sig_withC_gene_list.txt", quote=FALSE, row.names=FALSE,
            col.names=FALSE, sep="\t")
write.table(D2_sig_gene_names_exons_parsed, 
            file="D2_sig_withC_gene_list.txt", quote=FALSE, row.names=FALSE,
            col.names=FALSE, sep="\t")
write.table(D3_sig_gene_names_exons_parsed,  
            file="D3_sig_withC_gene_list.txt", quote=FALSE, row.names=FALSE,
            col.names=FALSE, sep="\t")
write.table(D4_sig_gene_names_exons_parsed, 
            file="D4_sig_withC_gene_list.txt", quote=FALSE, row.names=FALSE,
            col.names=FALSE, sep="\t")
write.table(D5_sig_gene_names_exons_parsed, 
            file="D5_sig_withC_gene_list.txt", quote=FALSE, row.names=FALSE,
            col.names=FALSE, sep="\t")

write.table(D1_noC_sig_gene_names, file="D1_noC_sig_gene_list.txt", quote=FALSE,
            row.names=FALSE, col.names=FALSE)
write.table(D2_noC_sig_gene_names, file="D2_noC_sig_gene_list.txt", quote=FALSE,
            row.names=FALSE, col.names=FALSE)
write.table(D3_noC_sig_gene_names, file="D3_noC_sig_gene_list.txt", quote=FALSE,
            row.names=FALSE, col.names=FALSE)
write.table(D4_noC_sig_gene_names, file="D4_noC_sig_gene_list.txt", quote=FALSE,
            row.names=FALSE, col.names=FALSE)
write.table(D5_noC_sig_gene_names, file="D5_noC_sig_gene_list.txt", quote=FALSE,
            row.names=FALSE, col.names=FALSE)

write.table(D1_noC_sig_gene_names_FBgn, file="D1_noC_sig_gene_names_FBgn_format.txt", quote=FALSE,
            row.names=FALSE, col.names=FALSE)
write.table(D2_noC_sig_gene_names_FBgn, file="D2_noC_sig_gene_names_FBgn_format.txt", quote=FALSE,
            row.names=FALSE, col.names=FALSE)
write.table(D3_noC_sig_gene_names_FBgn, file="D3_noC_sig_gene_names_FBgn_format.txt", quote=FALSE,
            row.names=FALSE, col.names=FALSE)
write.table(D4_noC_sig_gene_names_FBgn, file="D4_noC_sig_gene_names_FBgn_format.txt", quote=FALSE,
            row.names=FALSE, col.names=FALSE)
write.table(D5_noC_sig_gene_names_FBgn, file="D5_noC_sig_gene_names_FBgn_format.txt", quote=FALSE,
            row.names=FALSE, col.names=FALSE)

#saving data to the clipboard
tempdata <- cbind(D5_sig_gene_names_exons_parsed, D5_sig_gene_names_exons_parsed%in%C_sig_gene_names)
clip <- pipe("pbcopy", "w")                       
write.table(tempdata, file=clip, col.names=FALSE, row.names=FALSE)                               
close(clip)
