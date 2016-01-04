######################################
# Making lists of D SNPs *excluding* #
# SNPs within 1000 bp of C genes     #
######################################

# Gowinda requires SNPs as input (even when it's working
# in the 'gene' mode). This means until now I've been 
# working with D SNP lists that still include the 
# putative lab-adaptation SNPs. 
# Now want to exclude those SNPs.

library(refGenome)

setwd("~/Documents/Drosophila Selection Experiment/genome_6.01")

dros_genome <- ucscGenome()
dros_gtf <- read.gtf(dros_genome, "dmel-all-r6.01.ensembl.gtf")
gene_positions <- getGenePositions(dros_genome, by="gene_id")
# Not entirely sure this has located *all* genes (17023,
# versus 17262 supposed to be located on the genome according to
# v6.01 release notes... but might be okay)
# Checked a couple of locations with Flybase, they look correct

setwd("~/Documents/Drosophila Selection Experiment/gowinda_gene_category_enrichment")


allC_genes <- unique(read.csv("/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/Putative_lab_adaptation_gene_list_FBgn_format.txt",
                              header=FALSE, stringsAsFactors=FALSE, quote="\"")[,1])

length(which(as.character(gene_positions$gene_id)%in%allC_genes))
length(allC_genes)

# if lengths are not the same... it means I can't access all the genes I need this way.

notfound <- which(as.character(allC_genes%in%gene_positions$gene_id)==FALSE)
allC_genes[notfound]

[1]  "Gene_FBtr0304206"
#this is a micro RNA so not important


# X:17,629,791..17,632,047, 2L:10,388,105..10,393,228
# extra_coords <- data.frame(seqid=c("X", "2L"),
#                            start=c(17629791, 10388105),
#                            end=c(17632047,10393288))
other_coords <- gene_positions[as.character(gene_positions$gene_id)%in%allC_genes,c("seqid", "start", "end")]
# 
# allC_coords <- rbind(other_coords, extra_coords)

allC_coords <- other_coords
allC_coords$up1000 <- allC_coords$start-1000
allC_coords$down1000 <- allC_coords$end+1000

# D1_SNP_file <- read.csv("/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists/D1_sig_SNPs.txt", 
#                         sep="\t",
#                         stringsAsFactors=FALSE, header=FALSE)
# colnames(D1_SNP_file) <- c("chr", "pos")

remove_C_SNPs_from_D_file <- function(input_filename, C_SNP_table){
  temp_D_SNP_table <- read.csv(input_filename, sep="\t",
                               stringsAsFactors=FALSE, header=FALSE)
  colnames(temp_D_SNP_table) <- c("chr", "pos")
  for(i in 1:nrow(C_SNP_table)){
    temp_C_chr <- C_SNP_table[i,1]
    temp_C_up1000 <- C_SNP_table[i,4]
    temp_C_down1000 <- C_SNP_table[i,5]
    cond_vector <- temp_D_SNP_table$chr==temp_C_chr & temp_D_SNP_table$pos>=temp_C_up1000 & temp_D_SNP_table$pos<=temp_C_down1000
    #print(summary(cond_vector))
    temp_D_SNP_table <- temp_D_SNP_table[cond_vector==FALSE,]
  }
  new_file_name <- paste("noC_", input_filename, sep="")
  write.table(temp_D_SNP_table, file=new_file_name, sep="\t", col.names=FALSE,
              row.names=FALSE, quote=FALSE)
}

setwd('/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists')
write.table(allC_coords, file="Putative_lab_adaptation_gene_coords_plus_1000_bp.txt", quote=FALSE, row.names=FALSE)

remove_C_SNPs_from_D_file("D1_sig_SNPs.txt", allC_coords)
remove_C_SNPs_from_D_file("D2_sig_SNPs.txt", allC_coords)
remove_C_SNPs_from_D_file("D3_sig_SNPs.txt", allC_coords)
remove_C_SNPs_from_D_file("D4_sig_SNPs.txt", allC_coords)
remove_C_SNPs_from_D_file("D5_sig_SNPs.txt", allC_coords)
