####################################
# Testing whether observed overlap #
# among gene lists is significant  #
####################################

#############
# FUNCTIONS #
#############

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

########
# MAIN #
########

#Resample simulated "significant" SNPs from the file of all SNPs

Sample_code <- c("C1", "C2", "C3", "C4", "C5", 
               "D1", "D2", "D3", "D4", "D5")

setwd("~/Documents/Drosophila Selection Experiment/snp_and_gene_lists")

# Now look at the results
colname_vector <- c("chr", "pos", paste(rep("genename", times=88), 1:88, sep=""))

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

resample_names <- c("C1", "C2", "C3", "C4", "C5","D1", "D2", "D3", "D4", "D5"
)
nIter=1000

### Resampling SNPs where position
### is taken into account. 

for (i in 1:length(resample_names)) {
  temp_name <- resample_names[i]
  temp_filename <- paste("/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/MB_",
                         temp_name,
                         "_all_genes_extracted.txt", sep="")
  temp_full <- read.table(temp_filename, fill=TRUE, header=FALSE, 
                          colClasses="character", col.names=colname_vector)
  temp_full <- temp_full[temp_full$chr%in%c("2L", "2R", "3L", "3R", "X"),]
  temp_full$chrpos <- paste(temp_full$chr, temp_full$pos, sep="_")
  temp_sig <- get(paste(temp_name, "_sig_ex", sep=""))
  temp_sig$chrpos <- paste(temp_sig$chr, temp_sig$pos, sep="_")
  #temp_length <- nrow(temp_sig)
  find_sig_rows <- which(temp_full$chrpos %in% temp_sig$chrpos)
  resample_snp_table=as.list(rep(NA, length=nIter))
  for (j in 1:nIter) {
    
    temp_resample_index_shift <- sample(1:nrow(temp_full), size=1)
    temp_resample_index <- find_sig_rows + temp_resample_index_shift
    # 'wrap around' indices that would extend off the table
    temp_resample_index_wrapped <- temp_resample_index[temp_resample_index > nrow(temp_full)]
    temp_resample_index_to_use <- c(temp_resample_index_wrapped - nrow(temp_full),
                                    temp_resample_index[temp_resample_index < nrow(temp_full)])
    #print(cbind(find_sig_rows, temp_resample_index, temp_resample_index_to_use))
    temp_resample <- temp_full[temp_resample_index_to_use,]
    temp_gene_list <- extract_gene_names_from_table(temp_resample)
    resample_snp_table[[j]] <- temp_gene_list
    if(j%%50==0){
      print(paste(temp_name, " Resampling iteration #", j, sep=""))
      print(paste(nrow(temp_sig), "sig. SNPs;", 
                  length(find_sig_rows), "sig. SNPs in full data;",
                  temp_resample_index_shift, "is the index shift;",
                  nrow(temp_resample), "rows in resampled table;",
                  "Equates to", length(temp_gene_list), "genes",
                  sep=" "))
    }
  }
  assign(paste(temp_name, "resampled_genes_with_position", sep="_"), resample_snp_table)
  remove(temp_full)
}

setwd("/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/")
saveRDS(resample_snp_table, file="Gene_lists_from_resampled_SNP_table.rds")


#####################

rearrange_iterations <- list()
for(j in 1:nIter) {
  if(j%%50==0){
    print(paste("Rearranging iteration #", j, sep=" "))
  }
  temp_name <- paste("Iter", j, sep="_")
  rearrange_iterations[[temp_name]] <- list(C1_resampled_genes_with_position[[j]],
                                            C2_resampled_genes_with_position[[j]],
                                            C3_resampled_genes_with_position[[j]],
                                            C4_resampled_genes_with_position[[j]],
                                            C5_resampled_genes_with_position[[j]],
                                            D1_resampled_genes_with_position[[j]],
                                            D2_resampled_genes_with_position[[j]],
                                            D3_resampled_genes_with_position[[j]],
                                            D4_resampled_genes_with_position[[j]],
                                            D5_resampled_genes_with_position[[j]])

  names(rearrange_iterations[[j]]) <- resample_names
  temp_allC <- 
    rearrange_iterations[[j]][["allC"]] <- union(union(union(rearrange_iterations[[j]][[1]], 
                                                             rearrange_iterations[[j]][[2]]), 
                                                       rearrange_iterations[[j]][[3]]), 
                                                 union(rearrange_iterations[[j]][[4]], 
                                                       rearrange_iterations[[j]][[5]]))
  D_names <- c("D1", "D2", "D3", "D4", "D5")
  for(i in 1:5){
    temp_D <- D_names[i]
    temp_noC_name <- paste(temp_D, "noC", sep="_")
    rearrange_iterations[[j]][[temp_noC_name]] <- setdiff(rearrange_iterations[[j]][[temp_D]], rearrange_iterations[[j]][["allC"]])
  }
  combn2 <- combn(D_names, 2)
  combn2a <- apply(combn2, 2, paste, collapse="_")
  for(i in 1:ncol(combn2)){
    temp_Da <- rearrange_iterations[[j]][[paste(combn2[1,i], "noC", sep="_")]]
    temp_Db <- rearrange_iterations[[j]][[paste(combn2[2,i], "noC", sep="_")]]
    temp_combn_name <- combn2a[i]
    rearrange_iterations[[j]][[temp_combn_name]] <- intersect(temp_Da, temp_Db)
  }
  combn3 <- combn(D_names, 3)
  combn3a <- apply(combn3, 2, paste, collapse="_")
  for(i in 1:ncol(combn3)){
    temp_Da <- rearrange_iterations[[j]][[paste(combn3[1,i], "noC", sep="_")]]
    temp_Db <- rearrange_iterations[[j]][[paste(combn3[2,i], "noC", sep="_")]]
    temp_Dc <- rearrange_iterations[[j]][[paste(combn3[3,i], "noC", sep="_")]]
    temp_combn_name <- combn3a[i]
    rearrange_iterations[[j]][[temp_combn_name]] <- intersect(temp_Da, intersect(temp_Db, temp_Dc))
  } 
  combn4 <- combn(D_names, 4)
  combn4a <- apply(combn4, 2, paste, collapse="_")
  for(i in 1:ncol(combn4)){
    temp_Da <- rearrange_iterations[[j]][[paste(combn4[1,i], "noC", sep="_")]]
    temp_Db <- rearrange_iterations[[j]][[paste(combn4[2,i], "noC", sep="_")]]
    temp_Dc <- rearrange_iterations[[j]][[paste(combn4[3,i], "noC", sep="_")]]
    temp_Dd <- rearrange_iterations[[j]][[paste(combn4[4,i], "noC", sep="_")]]
    temp_combn_name <- combn4a[i]
    rearrange_iterations[[j]][[temp_combn_name]] <- intersect(intersect(temp_Da, temp_Db), intersect(temp_Dc, temp_Dd))
  } 
  combn5 <- "D1_D2_D3_D4_D5"
  rearrange_iterations[[j]][[combn5]] <- intersect(intersect(rearrange_iterations[[j]][["D1_noC"]], 
                                                             intersect(rearrange_iterations[[j]][["D2_noC"]], 
                                                                       rearrange_iterations[[j]][["D3_noC"]])), 
                                                   intersect(rearrange_iterations[[j]][["D4_noC"]], 
                                                             rearrange_iterations[[j]][["D5_noC"]]))
  C_names <- c("C1", "C2", "C3", "C4", "C5")
  #   for(i in 1:5){
  #     temp_C <- C_names[i]
  #   }
  combn2 <- combn(C_names, 2)
  combn2a <- apply(combn2, 2, paste, collapse="_")
  for(i in 1:ncol(combn2)){
    temp_Ca <- rearrange_iterations[[j]][[combn2[1,i]]]
    temp_Cb <- rearrange_iterations[[j]][[combn2[2,i]]]
    temp_combn_name <- combn2a[i]
    rearrange_iterations[[j]][[temp_combn_name]] <- intersect(temp_Ca, temp_Cb)
  }
  combn3 <- combn(C_names, 3)
  combn3a <- apply(combn3, 2, paste, collapse="_")
  for(i in 1:ncol(combn3)){
    temp_Ca <- rearrange_iterations[[j]][[combn3[1,i]]]
    temp_Cb <- rearrange_iterations[[j]][[combn3[2,i]]]
    temp_Cc <- rearrange_iterations[[j]][[combn3[3,i]]]
    temp_combn_name <- combn3a[i]
    rearrange_iterations[[j]][[temp_combn_name]] <- intersect(temp_Ca, intersect(temp_Cb, temp_Cc))
  } 
  combn4 <- combn(C_names, 4)
  combn4a <- apply(combn4, 2, paste, collapse="_")
  for(i in 1:ncol(combn4)){
    temp_Ca <- rearrange_iterations[[j]][[combn4[1,i]]]
    temp_Cb <- rearrange_iterations[[j]][[combn4[2,i]]]
    temp_Cc <- rearrange_iterations[[j]][[combn4[3,i]]]
    temp_Cd <- rearrange_iterations[[j]][[combn4[4,i]]]
    temp_combn_name <- combn4a[i]
    rearrange_iterations[[j]][[temp_combn_name]] <- intersect(intersect(temp_Ca, temp_Cb), intersect(temp_Cc, temp_Cd))
  } 
  combn5 <- "C1_C2_C3_C4_C5"
  rearrange_iterations[[j]][[combn5]] <- intersect(intersect(rearrange_iterations[[j]][["C1"]], 
                                                             intersect(rearrange_iterations[[j]][["C2"]], 
                                                                       rearrange_iterations[[j]][["C3"]])), 
                                                   intersect(rearrange_iterations[[j]][["C4"]], 
                                                             rearrange_iterations[[j]][["C5"]]))
  
  
  
}

saveRDS(rearrange_iterations, file="Gene_lists_including_overlap_from_resampled_SNPs.rds")

#############################
# Parse 'Exon_2L_start_end' #
# gene names in resampled   #
# gene lists                #
#############################

exons_for_conv <- read.table("/Users/pgriffin/Documents/Drosophila\ Selection\ Experiment/snp_and_gene_lists/exon_names_gene_list_for_conversion.txt", header=TRUE, 
                             fill=TRUE, stringsAsFactors=FALSE)

# Creating new objects that contain exon-parsed gene name vectors


  
parse_exons_in_gene_name_vector <- function(gene_name_vector, exon_conversion_table){
  temp_grep <- grep('Exon', gene_name_vector)
  temp_output_object <- gene_name_vector
  if(length(temp_grep)>0){
    for(j in 1:length(temp_grep)){
      temp_exon_name <- grep('Exon', gene_name_vector, value=TRUE)[j]
      temp_row <- exon_conversion_table[exon_conversion_table[1]==temp_exon_name]
      temp_output_object <- c(temp_output_object, temp_row[2:4])
    }
    temp_output_object <- temp_output_object[temp_output_object!=""]
    temp_output_object <- temp_output_object[-temp_grep]
  }
  unique_gene_names <- unique(temp_output_object[is.na(temp_output_object)==FALSE])
  #assign(paste(i, "_sig_gene_names_exons_parsed", sep=""), unique_gene_names)
  return(unique_gene_names)
}

exons_parsed <- list()

for(i in 1:nIter){
  temp_list <- rearrange_iterations[[i]]
  exons_parsed[[i]] <- lapply(temp_list, parse_exons_in_gene_name_vector, exon_conversion_table=exons_for_conv)
  names(exons_parsed[[i]]) <- names(rearrange_iterations[[i]])
  if(i%%100==0){
    print(paste("Parsing exons in iteration #", i))
  }
}

saveRDS(exons_parsed, file="Gene_lists_exons_parsed_including_overlap_from_resampled_SNPs.rds")


##### Figure S6 #######

setwd("/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/")

# Obtain C genes in each D replicate

C_sig_genes <- read.csv('Putative_lab_adaptation_gene_list.txt', header=FALSE,
                        stringsAsFactors=FALSE)[,1]

for(i in Sample_code[6:10]){
  temp_input_name <- paste(i, '_sig_withC_gene_list.txt', sep="")
  temp_file <- read.csv(temp_input_name, header=FALSE, stringsAsFactors=FALSE)[,1]
  temp_output_name <- paste('C_genes_in_',i, sep="")
  temp_output <- intersect(C_sig_genes, temp_file)
  assign(temp_output_name, temp_output)
}


setwd("/Users/pgriffin/Documents/Drosophila Selection Experiment/gene_list_overlap_testing/")

pdf("Histograms of C gene overlap for each D replicate after resampling with position.pdf", width=8, height=12)
# NB this is (now) using the exon-parsed resimulated gene lists
# which is better than pre-parsing.

par(mfcol=c(3,2))
for(i in 1:5){
  temp_D <- Sample_code[i+5]
  C_genes_in_temp <- paste("C_genes_in_", temp_D, sep="")
  temp_all_list <- lapply(exons_parsed, "[[", i+5)
  temp_all_noC_list <- lapply(exons_parsed, "[[", i+11)
  temp_all_length <- sapply(temp_all_list, length)
  temp_all_noC_length <- sapply(temp_all_noC_list, length)
  temp_dist <- temp_all_length-temp_all_noC_length
  
  multiplier <- hist(temp_dist, plot=FALSE)$counts / hist(temp_dist, plot=FALSE)$density
  mydensity <- density(temp_dist)
  mydensity$y <- mydensity$y * multiplier[1] 
  #offset <- sd(temp_dist)*2.5
  hist(temp_dist, xlim=c(0, 600), ylim=c(0, 320), main="",
       xlab=paste("No.", temp_D, "genes in putative lab-adaptation list", sep=" "))
  lines(mydensity)
  lines(y=c(0, 250), 
        x=c(length(get(C_genes_in_temp)), length(get(C_genes_in_temp))),
        col="red")
  myx <- seq(min(temp_dist), max(temp_dist), length = 100)
  normal <- dnorm(x=myx, mean = mean(temp_dist), sd = sd(temp_dist))
  lines(myx, normal * multiplier[1], col = "blue", lwd = 1)
  normless <- pnorm(q = length(get(C_genes_in_temp)), mean=mean(temp_dist), sd=sd(temp_dist))
  #pval <- (1-normless)*2
  text(x=length(get(C_genes_in_temp)), y=265, labels=signif(normless, 3))
  text(x=20, y=300, labels=LETTERS[i], cex=2)
  
}
dev.off()



###### Check the number of *genes* called from simulated SNPs

#First get real gene lists (again, before removing C genes, and after exon parsing)
for(i in Sample_code[6:10]){
  temp_input_name <- paste('/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/',
                           i, '_sig_withC_gene_list.txt', sep="")
  temp_file <- read.csv(temp_input_name, header=FALSE, stringsAsFactors=FALSE)[,1]
  temp_object_name <- paste(i, '_sig_gene_names', sep="")
  assign(temp_object_name, temp_file)
}



pdf("Histograms of simulated gene number for each replicate after resampling with position.pdf", width=8, height=16)

par(mfcol=c(5,2))
for(i in 1:10){
  temp_rep <- Sample_code[i]
  real_rep_genes <- paste(temp_rep, "_sig_gene_names", sep="")
  temp_all_list <- lapply(exons_parsed, "[[", i)
  temp_all_length <- sapply(temp_all_list, length)
  temp_mean <- mean(temp_all_length)
  temp_7sd <- sd(temp_all_length)*7

  hist(temp_all_length, main="", xlim=c(temp_mean-temp_7sd, temp_mean+temp_7sd),
       xlab=paste(temp_rep, "number of simulated genes", sep=" "))
  abline(v=length(get(real_rep_genes)), col="red")
}
dev.off()

##########################
# Convert to FBgn format
##########################

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

setwd("/Users/pgriffin/Documents/Drosophila Selection Experiment/gene_list_overlap_testing/")

exons_parsed <- readRDS(file="Gene_lists_exons_parsed_including_overlap_from_resampled_SNPs.rds")


iterations_FBgn_format <- list()
for(i in 1:length(exons_parsed)){
  if(i%%100==0){
    print(paste(i, "iterations of renaming done"))
  }
  temp_iteration <- exons_parsed[[i]]
  temp_renamed <- lapply(temp_iteration, FUN=gene_name_to_FBgn, 
                         conversion_table=gene_list_for_conversion, 
                         exon_translation=exon_translation)
  iterations_FBgn_format[[i]] <- temp_renamed
  names(iterations_FBgn_format[[i]]) <- names(temp_iteration)
}

saveRDS(iterations_FBgn_format, file="Resampled_gene_lists_from_SNPs_with_LD_FBgn_format.rds")



